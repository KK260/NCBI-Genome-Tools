import argparse
import os
import subprocess
import pandas as pd

def process_genbank_files(file_paths, select_group=None):
    """
    Process the GenBank files to extract gene information.
    
    Args:
    - file_paths (list): List of paths to the GenBank files.
    - select_group (str): Optional; Specific group to filter CDS features by.
    
    Returns:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - all_genes (set): Set of all unique gene names found across all files.
    """
    organism_section_gene_map = {}
    all_genes = set()

    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        with open(file_path, 'r') as file:
            organism_name = "Unknown Organism"
            section_name = "Unknown Section"
            genes = {}
            current_gene = None
            in_cds_feature = False
            sequence_lines = []
            
            for line in file:
                if line.startswith("  ORGANISM  "):
                    try:
                        organism_name = " ".join(line.split()[1:3])  # Expect two words for organism name
                    except IndexError:
                        organism_name = "Unknown Organism"

                elif line.startswith("            "):
                    terms = line.split(";")
                    section_name = terms[-2].strip() if len(terms) > 1 else "Unknown Section"
                    if select_group and select_group != section_name:
                        continue  # Skip this file if it doesn't match the selected group

                if line.startswith("     CDS"):
                    in_cds_feature = True
                elif in_cds_feature and "/gene=" in line:
                    current_gene = line.split("=")[1].strip().replace('"', '').lower()
                    genes[current_gene] = None
                    all_genes.add(current_gene)
                elif in_cds_feature and line.startswith("ORIGIN"):
                    in_cds_feature = False
                    sequence_lines = []  # Start capturing sequence lines
                elif sequence_lines is not None:
                    if line.startswith("//"):
                        if current_gene and sequence_lines:
                            genes[current_gene] = ''.join(sequence_lines).replace(' ', '').replace('\n', '')
                        sequence_lines = None  # Stop sequence capture
                    else:
                        sequence_lines.append(''.join(line.strip().split()[1:]))

            organism_section_gene_map[organism_name] = {
                "section": section_name,
                "genes": genes
            }

    return organism_section_gene_map, all_genes

def reorder_columns(df, section_order):
    """
    Reorder the columns of a DataFrame based on a specified section order.
    
    Args:
    - df (pd.DataFrame): DataFrame with gene presence/absence information.
    - section_order (list): List of section names in the desired order.
    
    Returns:
    - pd.DataFrame: DataFrame with reordered columns.
    """
    ordered_columns = ['Gene']
    for section in section_order:
        ordered_columns.extend([col for col in df.columns if col != 'Gene' and section in col])
    return df[ordered_columns]

def write_summary_files(organism_section_gene_map, all_genes, output_base, feature_section_summary, section_order, overwrite):
    """
    Write the summary information to CSV and XLSX files.
    
    Args:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - all_genes (set): Set of all unique gene names found across all files.
    - output_base (str): Base path for the output files (without extension).
    - feature_section_summary (bool): Whether to generate the section-wise summary.
    - section_order (list): List of section names in the desired order.
    - overwrite (bool): Whether to overwrite existing files and directories.
    """
    data = {'Gene': sorted(all_genes)}
    for organism in organism_section_gene_map:
        data[organism] = [gene if gene in organism_section_gene_map[organism]['genes'] else '' for gene in sorted(all_genes)]
    
    df = pd.DataFrame(data)
    
    if section_order:
        df = reorder_columns(df, section_order)
    
    # Write to CSV
    csv_output_file = f"{output_base}.csv"
    if not os.path.exists(csv_output_file) or overwrite:
        df.to_csv(csv_output_file, index=False)
        print(f"Gene summary written to {csv_output_file}")
    
    # Write to XLSX
    xlsx_output_file = f"{output_base}.xlsx"
    if not os.path.exists(xlsx_output_file) or overwrite:
        df.to_excel(xlsx_output_file, index=False)
        print(f"Gene summary written to {xlsx_output_file}")
    
    if feature_section_summary:
        section_gene_map = {}
        for organism in organism_section_gene_map:
            section = organism_section_gene_map[organism]['section']
            if section not in section_gene_map:
                section_gene_map[section] = set()
            section_gene_map[section].update(organism_section_gene_map[organism]['genes'].keys())
        
        section_data = {'Gene': sorted(all_genes)}
        for section in section_gene_map:
            section_data[section] = [gene if gene in section_gene_map[section] else '' for gene in sorted(all_genes)]
        
        section_df = pd.DataFrame(section_data)
        
        if section_order:
            section_df = reorder_columns(section_df, section_order)
        
        # Write section summary to CSV
        section_csv_output_file = f"{output_base}_section_summary.csv"
        if not os.path.exists(section_csv_output_file) or overwrite:
            section_df.to_csv(section_csv_output_file, index=False)
            print(f"Section-wise gene summary written to {section_csv_output_file}")
        
        # Write section summary to XLSX
        section_xlsx_output_file = f"{output_base}_section_summary.xlsx"
        if not os.path.exists(section_xlsx_output_file) or overwrite:
            section_df.to_excel(section_xlsx_output_file, index=False)
            print(f"Section-wise gene summary written to {section_xlsx_output_file}")

def write_combined_gene_sequences(organism_section_gene_map, output_fasta_file, overwrite):
    """
    Write all gene sequences from all organisms into a single combined FASTA file.

    Args:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and gene sequences.
    - output_fasta_file (str): Output path for the combined FASTA file.
    - overwrite (bool): Whether to overwrite existing files.
    """
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_fasta_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if not os.path.exists(output_fasta_file) or overwrite:
        with open(output_fasta_file, 'w') as f:
            for organism, info in organism_section_gene_map.items():
                for gene, sequence in info['genes'].items():
                    if sequence:  # Check if a valid sequence is available
                        # Write the organism and gene name in FASTA format
                        f.write(f">{organism.replace(' ', '_')}|{gene}\n{sequence}\n")
        print(f"Combined gene sequences written to {output_fasta_file}")
    else:
        print(f"Combined FASTA file already exists: {output_fasta_file}. Use --overwrite to regenerate.")

def align_combined_sequences(input_fasta_file, aligned_output_file, overwrite):
    """
    Align combined gene sequences using MAFFT.
    
    Args:
    - input_fasta_file (str): Path to the input combined FASTA file.
    - aligned_output_file (str): Output path for the aligned sequences.
    - overwrite (bool): Whether to overwrite existing alignment files.
    """
    if not os.path.exists(aligned_output_file) or overwrite:
        mafft_cmd = f"mafft --auto {input_fasta_file} > {aligned_output_file}"
        subprocess.run(mafft_cmd, shell=True)
        print(f"Aligned sequences written to {aligned_output_file}")

def run_raxml_ng(aligned_fasta_file, output_dir, overwrite):
    """
    Run RAxML-NG for phylogenetic analysis on the aligned sequences.
    
    Args:
    - aligned_fasta_file (str): Path to the aligned FASTA file.
    - output_dir (str): Directory for RAxML-NG output.
    - overwrite (bool): Whether to overwrite existing RAxML-NG files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, "raxml_output")
    if not os.path.exists(output_path) or overwrite:
        raxml_cmd = f"raxml-ng --all --msa {aligned_fasta_file} --model GTR+G --prefix {output_path}"
        subprocess.run(raxml_cmd, shell=True)
        print(f"RAxML-NG analysis completed for {aligned_fasta_file}")

import os
import subprocess
import urllib.request

def ensure_astral_installed(astral_jar_path):
    """
    Ensure that ASTRAL is installed by checking for the ASTRAL jar file.
    If it's not present, download the ASTRAL jar file.

    Args:
    - astral_jar_path (str): Path to the ASTRAL jar file.
    """
    astral_url = "https://github.com/smirarab/ASTRAL/raw/master/Astral/astral.5.7.8.jar"

    if not os.path.exists(astral_jar_path):
        print("ASTRAL not found, downloading...")
        try:
            urllib.request.urlretrieve(astral_url, astral_jar_path)
            print(f"ASTRAL downloaded and saved at {astral_jar_path}")
        except Exception as e:
            print(f"Error downloading ASTRAL: {e}")
            raise

def check_java_installed():
    """
    Check if Java is installed on the system.
    
    Raises:
    - OSError: If Java is not installed.
    """
    try:
        subprocess.run(["java", "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise OSError("Java is not installed. Please install Java and try again.")

def perform_astral_analysis(input_dir, output_dir, overwrite):
    """
    Perform ASTRAL analysis on RAxML-NG best tree files.

    Args:
    - input_dir (str): Directory containing RAxML-NG best tree files.
    - output_dir (str): Directory to save the ASTRAL output files.
    - overwrite (bool): Whether to overwrite existing files.
    """
    # Full path to the ASTRAL jar file (in the current working directory)
    astral_jar_path = os.path.join(os.getcwd(), "astral.jar")

    # Ensure ASTRAL is installed
    ensure_astral_installed(astral_jar_path)

    # Ensure Java is installed
    try:
        check_java_installed()
    except OSError as e:
        print(e)
        return

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Perform ASTRAL analysis on each RAxML best tree file
    for tree_file in os.listdir(input_dir):
        if tree_file.endswith(".raxml.bestTree"):
            input_path = os.path.join(input_dir, tree_file)
            output_path = os.path.join(output_dir, tree_file.replace(".raxml.bestTree", "_astral.tree"))

            if not os.path.exists(output_path) or overwrite:
                # Command to run ASTRAL
                astral_cmd = f"java -jar {astral_jar_path} -i {input_path} -o {output_path}"

                # Run ASTRAL and check if the command executed successfully
                try:
                    subprocess.run(astral_cmd, shell=True, check=True)
                    print(f"ASTRAL analysis completed for {tree_file}")
                except subprocess.CalledProcessError as e:
                    print(f"Error during ASTRAL analysis: {e}")

def main():
    parser = argparse.ArgumentParser(description="Process GenBank files to extract gene information and perform phylogenetic analyses.")
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Paths to the GenBank files')
    parser.add_argument('-f', '--feature_summary', required=True, help='Base path for feature summary files (CSV and XLSX)')
    parser.add_argument('-g', '--group_feature_summary', action='store_true', help='Whether to generate a section-wise feature summary')
    parser.add_argument('-o', '--group_order', nargs='+', help='Optional: Order of sections in the summary files')
    parser.add_argument('-s', '--generate_gene_sequences', action='store_true', help='Generate gene sequences in FASTA format')
    parser.add_argument('-a', '--align_sequences', action='store_true', help='Align gene sequences using MAFFT')
    parser.add_argument('-r', '--run_raxml', action='store_true', help='Run RAxML-NG for phylogenetic analysis')
    parser.add_argument('-t', '--run_astral', action='store_true', help='Perform ASTRAL analysis on RAxML trees')
    parser.add_argument('--select_group', default=None, help='Limit gene extraction to a specific group')
    parser.add_argument('--output_dir', default='gene_output', help='Output directory for gene sequences and alignment results')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files if they exist')
    
    args = parser.parse_args()

    # Process GenBank files and extract gene information
    organism_section_gene_map, all_genes = process_genbank_files(args.input, select_group=args.select_group)
    
    # Write summary files
    write_summary_files(organism_section_gene_map, all_genes, args.feature_summary, args.group_feature_summary, args.group_order, args.overwrite)
    
    # Generate combined gene sequence FASTA file
    combined_fasta_file = os.path.join(args.output_dir, "combined_gene_sequences.fasta")
    write_combined_gene_sequences(organism_section_gene_map, combined_fasta_file, args.overwrite)

    # Align combined sequences if specified
    aligned_fasta_file = os.path.join(args.output_dir, "aligned_gene_sequences.fasta")
    if args.align_sequences:
        align_combined_sequences(combined_fasta_file, aligned_fasta_file, args.overwrite)
    
    # Run RAxML-NG if specified
    if args.run_raxml:
        run_raxml_ng(aligned_fasta_file, args.output_dir, args.overwrite)
    
    # Run ASTRAL if specified
    if args.run_astral:
        perform_astral_analysis(args.output_dir, args.output_dir, args.overwrite)

if __name__ == "__main__":
    main()
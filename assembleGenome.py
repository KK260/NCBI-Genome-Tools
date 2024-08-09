#!/usr/bin/env python

"""
This script finds all available plastomes, mitogenomes, and nuclear genomes from a given group (e.g., genus, family, or order) in public databases. 

License:
    Copyright 2024 Kevin Karbstein
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import argparse
import os
import pandas as pd

def process_genbank_files(file_paths):
    """
    Process the GenBank files to extract gene information.
    
    Args:
    - file_paths (list): List of paths to the GenBank files.
    
    Returns:
    - organism_section_gene_map (dict): Mapping of organism names to their section names and dictionaries of gene names and sequences.
    - all_genes (set): Set of all unique gene names found across all files.
    """
    organism_section_gene_map = {}
    all_genes = set()

    # Parse each file and extract the necessary information
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            organism_name = None
            section_name = None
            genes = {}
            current_gene = None
            
            for line in lines:
                if line.startswith("  ORGANISM  "):
                    organism_name = line.split()[1] + " " + line.split()[2]
                elif line.startswith("            "):
                    terms = line.split(";")
                    section_name = terms[-2].strip() if len(terms) > 1 else "Unknown"

                if "/gene=" in line and "/CDS" in line:  # consider only CDS features
                    current_gene = line.split("=")[1].strip().replace('"', '').lower()
                    genes[current_gene] = None
                    all_genes.add(current_gene)
                if "ORIGIN" in line:
                    # Extract the DNA sequence
                    dna_sequence = ''.join(line.strip().split()[1:])
                    genes[current_gene] = dna_sequence

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
    
    # Reorder columns based on section order
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
    
    # Generate section-wise summary if requested
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
        
        # Reorder columns based on section order
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

def main():
    """
    Main function to parse command-line arguments and generate the gene summary files.
    """
    parser = argparse.ArgumentParser(description='Process GenBank files and extract gene names and sequences.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input GenBank file paths. Options: FILE1.gb FILE2.gb ...')
    parser.add_argument('--group_order', nargs='+', help='Order of groups for columns in the output. Options: GROUP1 GROUP2 ...')
    parser.add_argument('-o', '--feature_summary', required=True, help='Basic name of the feature summary tables and output files')
    parser.add_argument('-s', '--group_feature_summary', action='store_true', help='Generate group-wise summary')
    parser.add_argument('-g', '--gene_sequences', action='store_true', help='Generate gene sequence FASTA files')
    parser.add_argument('-a', '--align_sequences', action='store_true', help='Align gene sequences using MAFFT')
    parser.add_argument('-r', '--run_raxml', action='store_true', help='Run RAxML-NG analysis on aligned sequences')
    parser.add_argument('-x', '--run_astral', action='store_true', help='Run ASTRAL analysis using RAxML best trees')
    parser.add_argument('--overwrite', type=bool, default=False, help='Overwrite existing files and directories. Options: True/False')

    args = parser.parse_args()
    
    file_paths = args.input
    output_base = args.feature_summary
    feature_section_summary = args.group_feature_summary
    gene_sequences = args.gene_sequences
    align_sequences = args.align_sequences
    run_raxml = args.run_raxml
    run_astral = args.run_astral
    section_order = args.group_order
    overwrite = args.overwrite
    
    # Process the GenBank files to extract gene information
    organism_section_gene_map, all_genes = process_genbank_files(file_paths)
    
    # Write the summary to CSV and XLSX files
    write_summary_files(organism_section_gene_map, all_genes, output_base, feature_section_summary, section_order, overwrite)
    
    # Write gene sequences to FASTA files if requested
    if gene_sequences:
        gene_sequences_dir = os.path.join(output_base + "_gene_sequences")
        write_gene_sequences(organism_section_gene_map, gene_sequences_dir, overwrite)
    
    # Align gene sequences using MAFFT if requested
    if align_sequences:
        aligned_sequences_dir = os.path.join(output_base + "_aligned_sequences")
        align_gene_sequences(gene_sequences_dir, aligned_sequences_dir, overwrite)
    
    # Run RAxML-NG analysis if requested
    if run_raxml:
        raxml_output_dir = os.path.join(output_base + "_raxml_trees")
        run_raxml_ng(aligned_sequences_dir, raxml_output_dir, overwrite)
    
    # Run ASTRAL analysis if requested
    if run_astral:
        astral_output_dir = os.path.join(output_base + "_astral_trees")
        perform_astral_analysis(raxml_output_dir, astral_output_dir, overwrite)

if __name__ == "__main__":
    main()
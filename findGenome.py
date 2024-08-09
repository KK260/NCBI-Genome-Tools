#!/usr/bin/env python

"""
This script finds all available plastomes, mitogenomes, and nuclear genomes from a given group (e.g., genus, family, or order) in public databases. 

License:
    Copyright 2024 Kevin Karbstein and Lara Kösters
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
import os
import logging
import math
from argparse import ArgumentParser
from Bio import Entrez, SeqIO
from tqdm import tqdm
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def validate_inputs(group, org_type):
    """Validate user inputs."""
    valid_org_types = ['chloroplast', 'mitochondrial', 'nuclear_genome']
    if not group:
        raise ValueError("Group must be specified.")
    if org_type not in valid_org_types:
        raise ValueError(f"Invalid organism type. Choose from {valid_org_types}.")
    logging.info("Inputs validated successfully.")

def setup_output_folder(outfolder, overwrite=False):
    """Ensure the output folder exists or is overwritten."""
    if overwrite and os.path.exists(outfolder):
        shutil.rmtree(outfolder)
        logging.info(f"Output folder '{outfolder}' removed and will be recreated.")
    
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
        logging.info(f"Output folder '{outfolder}' created.")
    else:
        logging.info(f"Output folder '{outfolder}' already exists.")

def setup_excluded_folder(outfolder):
    """Ensure the excluded folder exists within the output folder."""
    excluded_folder = os.path.join(outfolder, "excluded")
    if not os.path.exists(excluded_folder):
        os.makedirs(excluded_folder)
        logging.info(f"Excluded folder '{excluded_folder}' created.")
    else:
        logging.info(f"Excluded folder '{excluded_folder}' already exists.")
    return excluded_folder

def search_genomes(group, org_type):
    """Search for genomes in the NCBI database."""
    db = 'nucleotide' if org_type != 'nuclear_genome' else 'genome'
    term = f'("complete genome") AND {group}[Organism]'
    if org_type == 'chloroplast':
        term += ' AND chloroplast'
    else:
        term += ' NOT chloroplast'

    try:
        handle = Entrez.esearch(db=db, term=term, rettype='count')
        retmax = int(Entrez.read(handle)['Count'])
        logging.info(f"Found {retmax} genomes for group '{group}' and type '{org_type}'.")
        return retmax, term
    except Exception as e:
        logging.error(f"Error searching genomes: {e}")
        raise

def estimate_genome_size(term, retmax, db):
    """Estimate the total size of the genomes in gigabytes."""
    try:
        handle = Entrez.esearch(db=db, term=term, idtype="acc", retmax=retmax)
        record = Entrez.read(handle)
        id_list = record['IdList']
        total_size = 0
        for genome_id in id_list:
            handle = Entrez.esummary(db=db, id=genome_id, rettype="gb", retmode="text")
            summary = Entrez.read(handle)[0]
            total_size += int(summary['Length'])
        total_size_gb = total_size / (1024 ** 3)  # Convert from bytes to gigabytes
        return total_size_gb
    except Exception as e:
        logging.error(f"Error estimating genome size: {e}")
        raise

def fetch_genome_ids(term, retmax, batch_size, db):
    """Fetch genome IDs from NCBI in batches."""
    try:
        handle = Entrez.esearch(db=db, term=term, idtype="acc", retmax=retmax)
        record = Entrez.read(handle)
        return record['IdList'], math.ceil(retmax / batch_size)
    except Exception as e:
        logging.error(f"Error fetching genome IDs: {e}")
        raise

def download_genomes(batch_ids, db):
    """Download genomes from NCBI."""
    try:
        handle = Entrez.efetch(db=db, id=batch_ids, rettype="gb", retmode="text")
        return SeqIO.parse(handle, "genbank")
    except Exception as e:
        logging.error(f"Error downloading genomes: {e}")
        raise

def save_genomes(records, outfolder, org_type):
    """Save downloaded genomes to files."""
    for record in records:
        filename = os.path.join(outfolder, f"{record.name}_{org_type}.gb")
        if os.path.exists(filename):
            logging.warning(f"File '{filename}' already exists. Skipping.")
            continue
        with open(filename, "w") as f:
            SeqIO.write(record, f, "genbank")
        logging.info(f"Saved {record.name} to '{filename}'.")

def find_full_organelle(group, outfolder, length_threshold, considered, org_type='chloroplast', batch_size=50, duplicate_removal=False, max_individuals_per_species=None, overwrite=False):
    """Main function to find and download genomes."""
    validate_inputs(group, org_type)
    setup_output_folder(outfolder, overwrite=overwrite)

    retmax, term = search_genomes(group, org_type)
    if retmax == 0:
        logging.warning("No genomes found. Exiting.")
        return

    db = 'nucleotide' if org_type != 'nuclear_genome' else 'genome'
    estimated_size_gb = estimate_genome_size(term, retmax, db)

    # Prompt user to proceed with the download
    user_input = input(f"Do you want to download {retmax} genomes with an estimated size of {estimated_size_gb:.2f} GB? (yes/no): ").strip().lower()
    if user_input != 'yes':
        logging.info("Download aborted by user.")
        return

    id_list, batches = fetch_genome_ids(term, retmax, batch_size, db)

    with tqdm(total=retmax) as pbar:
        for batch_num in range(0, batches):
            i = batch_num * batch_size
            batch_ids = id_list[i:i+batch_size]
            records = download_genomes(batch_ids, db)
            save_genomes(records, outfolder, org_type)
            pbar.update(len(batch_ids))
    
    # First, remove duplicates if specified
    if duplicate_removal:
        remove_duplicates(outfolder)

    # Then, apply the max individuals per species limit
    if max_individuals_per_species:
        apply_max_individuals_per_species(outfolder, max_individuals_per_species)

def remove_duplicates(outfolder):
    """Remove .gb files with duplicate FASTA sequences, keeping the 'NC_*' accession or latest release date."""
    fasta_sequences = {}
    excluded_folder = setup_excluded_folder(outfolder)
    
    for filename in os.listdir(outfolder):
        if filename.endswith(".gb"):
            filepath = os.path.join(outfolder, filename)
            with open(filepath, "r") as file:
                record = SeqIO.read(file, "genbank")
                seq = str(record.seq)
                organism = record.annotations["organism"]

                # Track sequences by organism
                if seq in fasta_sequences:
                    fasta_sequences[seq].append((filepath, organism))
                else:
                    fasta_sequences[seq] = [(filepath, organism)]
    
    for seq, files in fasta_sequences.items():
        if len(files) > 1:
            files.sort(key=lambda f: (
                not os.path.basename(f[0]).startswith("NC_"),  # Prioritize 'NC_' accession
                SeqIO.read(open(f[0]), "genbank").annotations.get("date", "")  # Sort by date if no 'NC_'
            ))
            
            kept_file = files[0]
            for file_to_remove in files[1:]:
                shutil.move(file_to_remove[0], os.path.join(excluded_folder, os.path.basename(file_to_remove[0])))
                logging.info(f"Removed duplicate file '{file_to_remove[0]}' (Organism: {file_to_remove[1]}) and kept the file '{kept_file[0]}' (Organism: {kept_file[1]}).")

def apply_max_individuals_per_species(outfolder, max_individuals_per_species):
    """Remove older files if max_individuals_per_species is specified, keeping only the most recent ones."""
    organisms = {}
    excluded_folder = setup_excluded_folder(outfolder)
    
    for filename in os.listdir(outfolder):
        if filename.endswith(".gb"):
            filepath = os.path.join(outfolder, filename)
            with open(filepath, "r") as file:
                record = SeqIO.read(file, "genbank")
                organism = record.annotations["organism"]

                # Track files by organism
                if organism not in organisms:
                    organisms[organism] = []
                organisms[organism].append(filepath)
    
    for organism, files in organisms.items():
        if len(files) > max_individuals_per_species:
            # Sort by date (PLN field)
            files.sort(key=lambda f: SeqIO.read(open(f), "genbank").annotations.get("date", ""))
            for file_to_remove in files[:-max_individuals_per_species]:
                shutil.move(file_to_remove, os.path.join(excluded_folder, os.path.basename(file_to_remove)))
                logging.info(f"Removed file '{file_to_remove}' due to max individuals per species limit for organism '{organism}'.")

if __name__ == "__main__":
    parser = ArgumentParser(description="Download organelle genomes from NCBI.")
    parser.add_argument("group", type=str, help="Organism group to search for.")
    parser.add_argument("outfolder", type=str, help="Folder to save the downloaded genomes.")
    parser.add_argument("--length_threshold", type=int, default=1000, help="Minimum sequence length to consider.")
    parser.add_argument("--considered", type=int, default=50, help="Number of sequences to consider for analysis.")
    parser.add_argument("--org_type", type=str, choices=['chloroplast', 'mitochondrial', 'nuclear_genome'], default="chloroplast", help="Type of organelle genome to search for.")
    parser.add_argument("--batch_size", type=int, default=50, help="Batch size for downloading.")
    parser.add_argument("--duplicate_removal", action='store_true', help="Whether to remove duplicate sequences. Options: ['TRUE', 'FALSE']. Default is FALSE.")
    parser.add_argument("--max_individuals_per_species", type=int, help="Maximum number of individuals per species to retain.")
    parser.add_argument("--overwrite", action='store_true', help="Whether to overwrite existing output folder. Options: ['TRUE', 'FALSE']. Default is FALSE.")
    args = parser.parse_args()

    try:
        Entrez.email = "your.email@example.com"  # Always provide an email address when using Entrez
        find_full_organelle(args.group, args.outfolder, args.length_threshold, args.considered, args.org_type, args.batch_size, args.duplicate_removal, args.max_individuals_per_species, args.overwrite)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
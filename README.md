# Analyzers-For-Finding-Comparing-Assembing-NCBI-Genomes
Automatically download, compare, and assemble genomes from NCBI


- findGenome.py

This script will automatically download and filter genomes from the NCBI database.

```
python findGenome.py "group" "outfolder" --org_type "chloroplast" --length_threshold INT --batch_size 50 --duplicate_removal --max_individuals_per_species INT

# example
python findGenome.py "Ranunculus" ./output_plastomes --org_type "chloroplast" --length_threshold 170000 --batch_size 50 --duplicate_removal --max_individuals_per_species 2


# usage:
findGenome.py [-h] [--length_threshold LENGTH_THRESHOLD] [--considered CONSIDERED]
                             [--org_type {chloroplast,mitochondrial,nuclear_genome}] [--batch_size BATCH_SIZE]
                             [--duplicate_removal] [--max_individuals_per_species MAX_INDIVIDUALS_PER_SPECIES]
                             group outfolder

Download plastid, mitochondrial, or nuclear genomes from NCBI.

positional arguments:
  group                 The taxonomic group to search for (e.g., genus, family, order).
  outfolder             The output folder where genomes will be saved.

options:
  -h, --help            show this help message and exit
  --length_threshold LENGTH_THRESHOLD
                        Minimum length of the genomes to be considered.
  --considered CONSIDERED
                        Specific genomes to consider.
  --org_type {chloroplast,mitochondrial,nuclear_genome}
                        The type of genome to download.
  --batch_size BATCH_SIZE
                        Batch size for downloading genomes.
  --duplicate_removal   Remove .gb files with duplicate sequences or based on max individuals per species.
  --max_individuals_per_species MAX_INDIVIDUALS_PER_SPECIES
                        Maximum number of individuals per species to retain.
```


### If you use any of the scripts, please cite the following reference until the article is published: 
### Karbstein et al. 2023, BioRxiv (https://doi.org/10.1101/2023.08.08.552429)

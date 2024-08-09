# Find-Compare-Assemble-Genomes-Analyzers
Automatically download, compare, and assemble genomes from NCBI

FindGenome.py

This script will automatically download and filter genomes from the NCBI database.

```
python findGenome.py "group" "outfolder" --org_type "chloroplast" --length_threshold INT --batch_size 50 --duplicate_removal --max_individuals_per_species INT

# example
python findGenome.py "Ranunculus" ./output_plastomes --org_type "chloroplast" --length_threshold 170000 --batch_size 50 --duplicate_removal --max_individuals_per_species 2

```


If you use any of the scripts, please cite the following reference until the article is published: Karbstein et al. 2023 (https://doi.org/10.1101/2023.08.08.552429)

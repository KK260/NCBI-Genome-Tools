# Tools-For-Finding-Comparing-Assembling-NCBI-Genomes
Automatically download, compare, and assemble genomes from NCBI


# - findGenome.py

This script will automatically download and filter genomes from the NCBI database.

```
# basic code:
python findGenome.py "group" "outfolder" --org_type "chloroplast" --length_threshold INT --batch_size 50 --duplicate_removal=FALSE --max_individuals_per_species INT --overwrite=FALSE

# example:
python findGenome_NCBI_v2.py "Ranunculus" ./plastomes_ranunculus --org_type "chloroplast" --batch_size 50 --duplicate_removal=TRUE --max_individuals_per_species 2 --overwrite=TRUE

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
  --max_individuals_per_species MAX_INDIVIDUALS_PER_SPECIES
                        Maximum number of individuals per species to retain.
  --duplicate_removal   Remove .gb files with duplicate sequences or based on max individuals per species.
```

# - assembleGenome.py

This script automatically extracts and compares features, and aligns and assembles annotated CDS regions from .gb files downloaded from NCBI.

```
# basic code:
python assembleGenome.py -i INPUT --group_order GROUP_ORDER -o FEATURE_SUMMARY --s GROUP_FEATURE_SUMMAY -g -a -r -x --overwrite=FALSE

# example:
python assembleGenome.py -i *.gb --group_order Fumarioideae Thalictroideae Delphinieae Ranunculeae Anemoneae -o mitogenome_features -s -g -a -r -x --overwrite=TRUE

# usage:
assembleGenome.py [-h] -i INPUT [FILE1.gb FILE2.gb ...] [--group_order GROUP_ORDER [GROUP1 GROUP2 ...]]
                         -o FEATURE_SUMMARY [-s]
                         [-g] [-a] [-r] [-x]
                         [--overwrite]

Process GenBank files and extract gene names and sequences.

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input GenBank file paths
  --group_order GROUP_ORDER [GROUP_ORDER ...]
                        Order of groups for columns in the output
  -o, --feature_summary FEATURE_SUMMARY
                        Base path for the output files (without extension)
  -s, --group_feature_summary
                        Generate group-wise summary
  -g, --gene_sequences  Generate gene sequence FASTA files
  -a, --align_sequences
                        Align gene sequences using MAFFT
  -r, --run_raxml       Run RAxML-NG analysis on aligned sequences
  -x, --run_astral      Run ASTRAL analysis using RAxML best trees
  --overwrite           Overwrite existing files and directories
```

### If you use any of the scripts, please cite the following reference until the journal article is published: 
Karbstein et al. (2023), BioRxiv (https://doi.org/10.1101/2023.08.08.552429)

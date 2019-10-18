#  Recovering and refining metagenome-assembled genomes (MAGs) from the Svalbard permafrost core samples

![Workflow](./img/Recover_and_Refine_MAGs.jpg)

# Usage
```
usage: Decon_MAG_by_taxa.py [-h] -r R -o O -l L -m M -k K -t T

Extract dominant taxa of subset based on Kaiju annotation results

Mandatory Arguments:
  -h, --help  show this help message and exit
  -r R        Path of your raw fasta folder
  -o O        Path of your output subset folder
  -l L        List of MAG ids you want to process
  -m M        Mapping information for Contig-MAG table
  -k K        Kaiju annotation result
  -t T        Threshold for dominant taxa percentage, default: 0.5

Example:


    python Decon_MAG_by_taxa.py -r ./raw_MAG_fa/ -o ./cleaned_MAG_fa/ -l decon_mag_list.txt -m scaf2bin.txt -k kaiju_anno.txt -t 0.51


```
# Dependencies
The whole workflow 

## input

## output

## Examples: running Decon_MAG_by_taxa.py on sample data.



# Reference

#  Recovering and refining metagenome-assembled genomes (MAGs) from the Svalbard permafrost core samples

![Workflow](./img/Recover_and_Refine_MAGs.jpg)


# Description

The entire workflow includes several bioinformatic steps and tools, including quality control (MOCAT), co-assembly (MegaHit), binning (MaxBin2, MetaBAT2), dereplication and aggregation (DASTool), and quality check (CheckM). In order to recover more MAGs, we suggest you follow the workflow. Detailed information of running above tools will not be introduced here, please check their references. Here we focus on instruction on how to perform further MAG refinement based on taxonomic classification.   

We observed that a large portion of MAGs had a high contamination percentage even after using DASTool. To improve the quality of MAGs, we developed a script, called “Decon_MAG_by_taxa.py”, that will subset each bin to into collections of contigs from the same taxonomic classification. 

By default, Kaiju will return a ‘NA’ if it cannot find a taxonomic classification, which results in loss of hierarchical rank structure. To keep track of the rank, here we considered ‘NA’ in Kaiju annotation as a special taxonomic rank, and sustained the hierarchical structure under the following rules: 1) when ‘NA’ observed in a low taxonomic rank a label is generated via combining higher taxonomic rank information with ‘_NA_’ denotation as a rank identifier (P:Phylum, C:Class, O:Order, F:Family, G:Genus, S:Species) 2) if ‘NA’ appeared at the phylum level a label is generated as ‘P_NA’. For example, if a contig is annotated as: ‘C1; Proteobacteria; Alphaproteobacteria; Rhizobiales; NA;  NA; Unknown species’, then it will be converted to: ‘C1; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiales_NA_F;  Rhizobiales_NA_F_NA_G; Unknown species’. Later, the script calculates the percentage of every taxa label in each rank, and kept labels whose percentage were higher than a user-defined threshold (default=0.5). 


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
  -k K        Kaiju annotation result with taxa labels
  -t T        Threshold for dominant taxa percentage, default: 0.5

Example:


    python Decon_MAG_by_taxa.py -r ./raw_MAG_fa/ -o ./cleaned_MAG_fa/ -l decon_mag_list.txt -m scaf2bin.txt -k kaiju_anno.txt -t 0.51


```
# Dependencies
* Python 3.6 or higher
* Python Package: pandas
* Biopython: https://biopython.org/


# Preparation

* Run Kaiju annotation
```
kaiju -t nodes.dmp -f kaiju_db_nr.fmi -i input_contigs.fa -o kaiju_tax_nr_contigs.out
```
* Add taxa labels
```
addTaxonNames -t nodes.dmp -n names.dmp -i kaiju_tax_nr_contigs.out -r phylum,class,order,family,genus,species -o kaiju_anno.txt 
```
* Generate a tab separated file of Contig_IDs and MAG_IDs (Dastool output file: DASTool_scaffolds2bin.txt)
```
k127_36213	bin_15
k127_42658	bin_15
k127_66274	bin_15
k127_69089	bin_16
k127_167398	bin_16
k127_286487	bin_9
k127_306140	bin_9
```
* Select high-contaminated (in general over 10%) MAG_IDs based on CheckM report.(Of course you could run for all MAGs if you prefer).

# Examples: running Decon_MAG_by_taxa.py on sample data.

For example, you recovered 8 MAGs, you want to run further 'Decon_MAG_by_taxa.pu' for 6 of them ('decon_mag_list.txt').

```
└─[0] <> ls ./raw_MAG_fa
bin_10.fa  bin_15.fa  bin_16.fa  bin_17.fa  bin_21.fa  bin_3.fa  bin_4.fa  bin_9.fa
└─[0] <> cat decon_mag_list.txt
bin_10
bin_15
bin_21
bin_3
bin_4
bin_9
└─[0] <> head kaiju_anno.txt
C	k127_12932	1100719	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Limnohabitans; Limnohabitans sp. Rim11;
C	k127_28115	2182327	Proteobacteria; Betaproteobacteria; Nitrosomonadales; Gallionellaceae; Candidatus Nitrotoga; Candidatus Nitrotoga fabula;
C	k127_44777	83494	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; NA; NA;
C	k127_83189	1736433	Proteobacteria; Betaproteobacteria; Burkholderiales; NA; Rhizobacter; Rhizobacter sp. Root1221;
C	k127_107259	192843	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Rhodoferax; Rhodoferax ferrireducens;
C	k127_108296	28065	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Rhodoferax; NA;
C	k127_131934	281915	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Curvibacter; NA;
C	k127_173258	28065	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Rhodoferax; NA;
C	k127_184610	1797557	Proteobacteria; Betaproteobacteria; Burkholderiales; NA; NA; Burkholderiales bacterium RIFCSPHIGHO2_12_FULL_61_11;
C	k127_202089	192843	Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Rhodoferax; Rhodoferax ferrireducens;
└─[0] <> head scaf2bin.txt
k127_36213	bin_15
k127_42658	bin_15
k127_66274	bin_15
k127_69089	bin_15
k127_167398	bin_15
k127_222707	bin_15
k127_233768	bin_15
k127_247578	bin_15
k127_286487	bin_15
k127_306140	bin_15
```
Also you want to output subset with a higher threshold (0.6):
```
└─[2] <> python Decon_MAG_by_taxa.py -r ./raw_MAG_fa/ -o ./cleaned_MAG_fa/ -l decon_mag_list.txt -m scaf2bin.txt -k kaiju_anno.txt -t 0.6

Parse MAG:  bin_10
['phylum#Proteobacteria#0.98', 'class#_Betaproteobacteria#0.95', 'order#_Burkholderiales#0.87', 'family#_Comamonadaceae#0.77']
Finish!

Parse MAG:  bin_15
['phylum#candidate_division_WPS-2#0.62', 'class#_NA#0.72', 'order#_NA#0.72', 'family#_NA#0.74', 'genus#_NA#0.78', 'species#_candidate_division_WPS-2_bacterium#0.62']
Finish!

Parse MAG:  bin_21
['family#_NA#0.61', 'genus#_NA#0.64']
Finish!

Parse MAG:  bin_3
['phylum#Actinobacteria#0.94', 'class#_Actinobacteria#0.93']
Finish!

Parse MAG:  bin_4
['phylum#Proteobacteria#0.97', 'class#_Betaproteobacteria#0.88', 'order#_Burkholderiales#0.81', 'family#_Oxalobacteraceae#0.65']
Finish!

Parse MAG:  bin_9
['phylum#Acidobacteria#0.92', 'class#_Acidobacteriia#0.88', 'order#_Acidobacteriales#0.88', 'family#_Acidobacteriaceae#0.78']
Finish!

└─[0] <> ls ./cleaned_MAG_fa
bin_10.class___Betaproteobacteria__0.95.cleaned.fa  bin_15.phylum__candidate_division_WPS-2__0.62.cleaned.fa              bin_4.order___Burkholderiales__0.81.cleaned.fa
bin_10.family___Comamonadaceae__0.77.cleaned.fa     bin_15.species___candidate_division_WPS-2_bacterium__0.62.cleaned.fa  bin_4.phylum__Proteobacteria__0.97.cleaned.fa
bin_10.order___Burkholderiales__0.87.cleaned.fa     bin_21.family___NA__0.61.cleaned.fa                                   bin_9.class___Acidobacteriia__0.88.cleaned.fa
bin_10.phylum__Proteobacteria__0.98.cleaned.fa      bin_21.genus___NA__0.64.cleaned.fa                                    bin_9.family___Acidobacteriaceae__0.78.cleaned.fa
bin_15.class___NA__0.72.cleaned.fa                  bin_3.class___Actinobacteria__0.93.cleaned.fa                         bin_9.order___Acidobacteriales__0.88.cleaned.fa
bin_15.family___NA__0.74.cleaned.fa                 bin_3.phylum__Actinobacteria__0.94.cleaned.fa                         bin_9.phylum__Acidobacteria__0.92.cleaned.fa
bin_15.genus___NA__0.78.cleaned.fa                  bin_4.class___Betaproteobacteria__0.88.cleaned.fa
bin_15.order___NA__0.72.cleaned.fa                  bin_4.family___Oxalobacteraceae__0.65.cleaned.fa

```
Our script will print all taxa informtion whose percantage is higher than a user-defined threshod at each rank, and generate subset fasta files. Subset fasta is renamed with following rules: MAG_ID + rank identifier + taxa label + taxa percentage.

As our script provides multiple subsets of fasta corresponding to different ranks for each MAG, you can run CheckM with all of these subsets and evaluate the best tradeoff between completeness and contamination. 


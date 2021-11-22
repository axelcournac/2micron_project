#  µproject

Some codes and scripts to reproduce the analyses and plots. 


### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `hicstuff` / [hicstuff](https://github.com/koszullab/hicstuff)
* `tinyMapper` / [tinyMapper](https://github.com/js2264/tinyMapper)
* `cooler` / [cooler](https://github.com/open2c/cooler)

### Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA can be used to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory
 
```bash
./fastq-dump SRR639031 --split-3 -O /home/axel/data/
```

#### Alignment
To align the reads and generate the contact files in cool format, we hicstuff pipeline: 
```bash
hicstuff pipeline -t 18 -i -D -a bowtie2 --matfmt bg2 --no-cleanup -F -p -o out_Micro-C_WT_log_classic_genome  -g SC288_with_micron SRR7939017.1_1.fastq SRR7939017.1_2.fastq
```
and convert into cool file:
```bash
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 sacCer3.chr_sizes.txt:200 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
```

### Processing of genomic data like Mnase-seq, ChIP-seq or RNA-seq
We used tinyMapper: 
```bash
./tinyMapper.sh -m MNase -s SRR6246290.1 -g SC288_with_micron_SC88 -o results_ATAC-seq

./tinyMapper.sh -m ChIP -s SRR7175393.1 -i SRR7175394.1 -g SC288_with_micron_SC88 -o results_CHIP

./tinyMapper.sh -m RNA -s SRR8503057.1 -g SC288_with_micron_SC88 -o results_RNAseq
```

### Visualisation of contact maps with plasmid signal

#### Alignment 

```bash
hicstuff pipeline -t 18 --read-len=50  -D -a bowtie2 -e DpnII,HinfI --matfmt bg2 --no-cleanup -F -p -o out2_pKan-STB-P -g /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_pKan-STB-P/SC288_with_pKan-STB-P /media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/20210716_FG/FG104_S2_R1_001.fastq.gz /media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/20210716_FG/FG104_S2_R2_001.fastq.gz
```
#### Convertion into cool file:
```bash
cooler cload pairs --zero-based -c1 1 -p1 2 -c2 4 -p2 5 /home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_pKan-STB-P/sacCer3.chr_sizes.txt:2000  out2_pKan-STB-P/tmp/valid_idx_pcrfree.pairs  out2_pKan-STB-P/tmp/valid_idx_pcrfree.pairs.cool 
```

#### Visualisation of contact maps
```bash
python3 /home/axel/Bureau/z_python_scripts_copy/plasmid_micron2_hot_spots_ARG1.py out2_pKan-STB-P/tmp/valid_idx_pcrfree.pairs.cool     pKan-STB-P     pKan-STB-P
```
To have the 1D enrichment plot for contact signal of the plasmid:

```bash
 python3 plasmid_HSC_1D_agglo_norm2.py  valid_idx_pcrfree.pairs.cool plasmid_p2-micron topo2 HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated
 ```
 
To automatically detect the peaks of contact between plasmid and yeast chromosomes, we use plasmid_micron2_hot_spots_ARG1_fig_sup1_norm3.py 



#### Computation of pileup plot according the gene structure: 
creation of bw file with the python code then

```bash

computeMatrix scale-regions -S signal_contact_2micron_plasmid.bw -R /home/axel/Bureau/YEAST/GENES_SC288/long_genes_only_host_chrm.txt2  --beforeRegionStartLength 3000 --regionBodyLength 7000 --afterRegionStartLength 3000 --outFileName signal_2u_genes.gz

plotProfile -m signal_2u_genes.gz -out ExampleProfile1.png --numPlotsPerRow 2  --plotTitle "Test data profile"
```














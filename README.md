#  2µ project

Codes and scripts to reproduce the analyses and plots. 


### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `SRA tools` / [SRA tools](https://github.com/ncbi/sra-tools)
* `Bowtie2` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `hicstuff` / [hicstuff](https://github.com/koszullab/hicstuff)
* `cooler` / [cooler](https://github.com/open2c/cooler)
* `tinyMapper` / [tinyMapper](https://github.com/js2264/tinyMapper)
* `deepTools` / [deepTools](https://deeptools.readthedocs.io/en/develop/)
* `RIdeogram` / [RIdeogram](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html)
* `seqkit`   /   [seqkit](https://bioinf.shenwei.me/seqkit/)
* `dnaglider`   /   [dnaglider](https://github.com/cmdoret/dnaglider)
* `pyBigWig`   /   [pyBigWig](https://github.com/deeptools/pyBigWig)


### Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA toolkit can be used to download and split both mates of a NGS library: 
 
```bash
./fastq-dump SRR639031 --split-3 -O /home/data/
```

#### Alignment of the Hi-C libraries
To align the reads and generate the contact files in cool format, we used hicstuff pipeline: 
```bash
hicstuff pipeline -t 18 --mapping="iterative" -D -a bowtie2 --matfmt bg2 --no-cleanup -F -p -o out_Micro-C_WT_log_classic_genome  -g SC288_with_micron SRR7939017.1_1.fastq SRR7939017.1_2.fastq
```
and convert into cool file:
```bash
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 sacCer3.chr_sizes.txt:200 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
```

### Processing of genomic data like Mnase-seq, RNA-seq or ChIP-seq 
We used tinyMapper: 
```bash
./tinyMapper.sh -m MNase -s SRR5399542.1 -g S288C_with_micron -o results_H3_CC 

./tinyMapper.sh -m RNA -s SRR7692240.1 -g S288C_with_micron -o results_RNAseq

./tinyMapper.sh -m MNase -s SRR11235539.1 -g S288C_with_micron -o results_ATAC-seq

./tinyMapper.sh -m ChIP -s SRR7175393.1 -i SRR7175394.1 -g S288C_with_micron -o results_CHIP_Rpb3

```

### Visualisation of contact maps with plasmid signal

We used home python code to plot contact maps and signal of 2 micron plasmid. 

#### Visualisation of contact maps
```bash
python3 plasmid_micron2_hot_spots_ARG1.py out2_pKan-STB-P/tmp/valid_idx_pcrfree.pairs.cool  pKan-STB-P pKan-STB-P
python3 plasmid_micron_its_Map.py  valid_idx_pcrfree.pairs.cool.200 plasmid_p2-micron micro-C-log
```
To have the 1D enrichment plot for contact signal of the plasmid:

```bash
python3 plasmid_HSC_1D_agglo_norm3.py  valid_idx_pcrfree.pairs.cool plasmid_p2-micron topo2 HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated
 ```
 
To automatically detect the peaks of contact between plasmid and yeast chromosomes, we use plasmid_micron2_hot_spots_ARG1_fig_sup1_norm3.py 

#### Averaged plot of genomic signals around hot spots of contact: 

```bash
python plasmid_micron2_Chip-seq_ARG5.py SRR13736589.bis.fastq.sam.MQ0 SRR13736587.bis.fastq.sam.MQ0 H3_log H3_log

python plasmid_micron2_Chip-seq_ARG5_1Mnase.py SRR7692240.1_1.fastq.gz.sam.MQ0 RNAseq RNAseq # with log 

python plasmid_micron2_Chip-seq_ARG5_1Mnase.py SRR6246290.1.bis_1.fastq.sam.MQ0 ATAC_WT
```

#### Computation of aggregated plots around the 73 loci contacted by 2u plasmid:

```bash
computeMatrix scale-regions  -S  Micro-C_WT_log_interpolated.bw -R HSC_73.bed2  --beforeRegionStartLength 20000  --regionBodyLength 10  --afterRegionStartLength 20000  --sortRegions keep -o heatmap.gz

plotHeatmap -m heatmap.gz -out heat_map_contact_best.pdf --colorMap afmhot_r  --missingDataColor white --sortRegions keep



computeMatrix scale-regions -S SRR14693235.bis^mapped_SC288_with_micron_SC88^ZW95QD.unstranded.CPM.bw -R HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated.bed --beforeRegionStartLength 20000  --regionBodyLength 10  --afterRegionStartLength 20000  --skipZeros -o heatmap.gz

plotHeatmap -m heatmap.gz -out heat_map_pol2_test4.pdf --colorMap Blues --missingDataColor grey  --sortRegions keep  --zMin 0. --zMax 20.0 --interpolationMethod  nearest
```
#### Plot for the gene boxes: 
Use home made python code gene_boxes.py
```bash
pdfjam $(ls -v  window7_gene_*.pdf) --nup 1x73 --landscape --outfile Page12_win_genes7.pdf
```


#### Computation of pileup plot according the gene structure: 
First, we create of bw file for the contact signal with the python code create_big_wig.py then

```bash

computeMatrix scale-regions -S Micro-C_WT_log_redone_plasmid_contact_signal.bw   -R long_genes_only_host_chrm.txt2  --beforeRegionStartLength 7000 --regionBodyLength 7000 --afterRegionStartLength  7000 --outFileName signal_2u_genes.gz 

plotProfile -m signal_2u_genes.gz -out Profile_contact_long_genes.pdf --numPlotsPerRow 2  --plotTitle "Contact signal at long gene (size>7b)"
```

#### 
To plot the agglomerated plots of contact signal around centromeres: 
```bash
python3 plasmid_HSC_1D_agglo_norm3.py  valid_idx_pcrfree.pairs.cool  plasmid_p2-micron log_centros centro1.dat55
```

### Generation and comparison of agglomerated plot of contact signal of 2 micron:
To automate the visualization and comparison of the average contact signal between 2 biological conditions, we use a home made script in bash. 

```bash
file_set_positions="HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated"
cool_file1="out_FG01/tmp/valid_idx_pcrfree.pairs.2000.cool"
cool_file2="out_FG02_4/tmp/valid_idx_pcrfree.pairs.2000.cool"

name1="FG1"
name2="FG2"

python3 plasmid_HSC_1D_agglo_norm2.py $cool_file1 plasmid_p2-micron $name1  $file_set_positions
python3 plasmid_HSC_1D_agglo_norm2.py $cool_file2 plasmid_p2-micron $name2  $file_set_positions

file1='/home/axel/Bureau/compare_agglo_redone_ylim/'$name1'_files/agglomerated_signal_on_HSC__'$name1'_2.0kb_norm2.txt'
file2='/home/axel/Bureau/compare_agglo_redone_ylim/'$name2'_files/agglomerated_signal_on_HSC__'$name2'_2.0kb_norm2.txt'

python3 plot_agglo_HSC2_norm2_arg_n.py  2 $file1 $file2 $name1 $name2 

```

### Generation and visualisation of scatter plot for the ChIP-exo libraries

After alignment with bowtie2 of the 1251 libraries, we filtered the reads with MQ>0 and then process all the files: 

```bash
#!/bin/bash

for i in *.fastq.sam.MQ0 ; 
do j=$(echo $i | grep -E -o  SRR[0-9]+) ; 
if [ ! -f "1D_ENRICHMENT_HSC_"$j"_73Hot_spots_plasmid.pdf" ];
then 
echo $j;
python3 plasmid_micron2_Chip-seq_ARG3.py $i $i; 
fi; done
```

repo: /media/axel/RSG4/diverse_yeast_data/CHip_seq_2018

We use the code density_scatter_plot5.py to generate the plot. 

#### Computation of GC content for regions contacted by the plasmid and the whole genome

```bash
seqkit fx2tab sequences_73HSC.txt -n --length --gc --gc-skew --header-line > gc_content_73HSC.txt

dnaglider-linux -window 10000 -threads 8 -fields "GC,GCSKEW" -fasta SC288_with_micron.fa  -out gc_stats.tsv

```

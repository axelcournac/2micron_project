# project

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
./tinyMapper.sh -m RNA -s SRR8503057.1 -g SC288_with_micron_SC88 -o results_RNAseq

./tinyMapper.sh -m ChIP -s SRR7175393.1 -i SRR7175394.1 -g SC288_with_micron_SC88 -o results_CHIP

./tinyMapper.sh -m MNase -s SRR6246290.1.bis -g SC288_with_micron_SC88 -o results_ATAC-seq
```














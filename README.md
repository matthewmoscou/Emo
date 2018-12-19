# Emo
Genome sequencing and annotation of *Ecdeiocolea monostachya*

## Illumina sequencing of *Ecdeiocolea monostachya*
Sheath tissue from accession E001 was collected, genomic DNA extracted, and DNA sent to Novogene for sequencing. Illumina paired end libraries were generated using inserts of 250bp (DSW66921) and 350bp (DSW66909-V) (150 bp reads). *De novo* assemblies were carried out using several Amazon AWS EC2 instances including (1) m5a.12xlarge (48 CPU, 192 GB RAM), (2) r5.12xlarge (48 CPU, 374 Gb RAM), and (3) XXX.

### Trimmomatic cleaning of Illumina PE reads
Trimmomatic v0.36 was used to clean reads prior to *de novo* assembly. A size limit of 36 bp is set for the retention of paired reads.

```bash
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTGL5CCXY_L6_1.fq.gz EM009_E1_DSW66909-V_HTGL5CCXY_L6_2.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_350_1_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTH7NCCXY_L7_1.fq.gz EM009_E1_DSW66909-V_HTH7NCCXY_L7_2.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_350_2_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTGL5CCXY_L7_1.fq.gz EM009_E1_DSW66909-V_HTGL5CCXY_L7_2.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_350_3_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HWHMWCCXY_L8_1.fq.gz EM009_E1_DSW66909-V_HWHMWCCXY_L8_2.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_350_4_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66921_HTGL5CCXY_L6_1.fq.gz EM009_E1_DSW66921_HTGL5CCXY_L6_2.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_250_1_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66921_HTGL5CCXY_L7_1.fq.gz EM009_E1_DSW66921_HTGL5CCXY_L7_2.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_forward_unpaired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_reverse_paired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > Ecdeiocolea_monostachya_250_2_trimmomatic.run.log 2>&1 &
```

Reads were also cleaned based on the stringent requirement that all 150 bp meet quality requirements prior to *de novo* assembly (needed for `edena`).

```bash
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTGL5CCXY_L6_1.fq.gz EM009_E1_DSW66909-V_HTGL5CCXY_L6_2.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_350_1_edena_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTH7NCCXY_L7_1.fq.gz EM009_E1_DSW66909-V_HTH7NCCXY_L7_2.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_350_2_edena_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HTGL5CCXY_L7_1.fq.gz EM009_E1_DSW66909-V_HTGL5CCXY_L7_2.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_350_3_edena_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66909-V_HWHMWCCXY_L8_1.fq.gz EM009_E1_DSW66909-V_HWHMWCCXY_L8_2.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_350_4_edena_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66921_HTGL5CCXY_L6_1.fq.gz EM009_E1_DSW66921_HTGL5CCXY_L6_2.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_250_1_edena_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 EM009_E1_DSW66921_HTGL5CCXY_L7_1.fq.gz EM009_E1_DSW66921_HTGL5CCXY_L7_2.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_edena_forward_paired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_edena_forward_unpaired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_edena_reverse_paired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_edena_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:150 > Ecdeiocolea_monostachya_250_2_edena_trimmomatic.run.log 2>&1 &
```

## Assessment of *k*-mer distribution
`jellyfish` 1.1.12 was used to assess *k*-mer distribution, determine the presence of heterozygosity, and estimate genome size.

```bash
jellyfish count -t 48 -C -m 24 -s 30G -o emo_jellyfish_24mer Ecdeiocolea_monostachya_250_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_reverse_paired.fq
jellyfish histo -h 3000000 -o emo_jellyfish_24mer.histo emo_jellyfish_24mer_0
```

`R` and `ggplot2` were used to visualize the results.

```R
library(ggplot2)

data = data.frame(data)
data = read.table(file="emo_jellyfish_24mer.histo.ID", header=T)

postscript(file="emo_jellyfish_24mer_distribution.ps", width=6, height=4)
ggplot(data, aes(k, count)) + geom_point() + xlim(c(4,400)) + ylim(c(0,1.5e7)) + xlab("Frequency") + ylab("Total counts")
dev.off()
```

![alt text](figures/emo_jellyfish_24mer_distribution.png "Frequency distribution for Ecdeiocolea monostachya sequence data")

Sequencing of *Ecdeiocolea monostachya* accession E01 exhibits two local maxima, which indicates that this accession is heterozygous and diploid. We used the function `localMaxima` in order to identify the peaks in the distribution.

```
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

localMaxima(data$count[5:400])
```

The results below indicate the peaks are at 77 and 161.

```R
 [1]   1  77 161 247 299 326 329 331 334 389 395
```

Next, we use [findGSE](https://github.com/schneebergerlab/findGSE) from the Schneeberger laboratory to get an estimate of genome size.

```R
library("findGSE")
findGSE(histo="emo_jellyfish_24mer.histo", sizek=24, outdir="emo_jellyfish_24mer_findGSE")
```

```R
Genome size estimate for emo_jellyfish_24mer.histo: 1477061969 bp.
```

The final estimate was 1.47 Gb. This is a slightly lower genome size estimate as compared to flow cytometry of propidium iodide-stained nuclei, of which this accession was found to have 2C = ~2.0 pg = ~2.0 Gb)

`KmerGenie` 1.7048 was used to characterize the *k*-mer distribution based on several *k* and to identify an optimal *k* for genome assembly. Code was downloaded from [KmerGenie](http://kmergenie.bx.psu.edu/).

```bash
~/genome/src/kmergenie-1.7051/kmergenie list_files --diploid -t 46
```

[alt text](figures/histograms.dat.png "Optimal k-mer analysis")

The best predicted *k* = 107 with a predicted genome size of 1.31 Gb. Results can be found in `data/kmerGenie`.

### *De novo* assembly using edena
`edena` 3.131028 was used for *de novo* genome assembly.

```bash
edena -paired Ecdeiocolea_monostachya_250_1_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_edena_reverse_paired.fq -p Emo_edena_v1 -nThreads 48 > Emo.edena.run.log 2>&1 &
edena -e Emo_edena_v1.ovl -m 100 -p Emo_edena_v1_m100
```

This initial run failed (*k*=96), likely due to the need of more memory. Consider running it again at a later date. Also, ensure all files are decompressed.

### *De novo* assembly using ABySS
`ABySS` version 1.3.6 was used for *de novo* genome assembly. Source code was cloned from Github from [ABySS](https://github.com/bcgsc/abyss).

```bash
abyss-pe k=96 name=Emo_abyss_k96 lib='pea peb pec ped pee pef' pea='Ecdeiocolea_monostachya_250_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_reverse_paired.fq' peb='Ecdeiocolea_monostachya_250_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_reverse_paired.fq' pec='Ecdeiocolea_monostachya_350_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_reverse_paired.fq' ped='Ecdeiocolea_monostachya_350_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_reverse_paired.fq' pee='Ecdeiocolea_monostachya_350_3_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_reverse_paired.fq' pef='Ecdeiocolea_monostachya_350_4_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_reverse_paired.fq' np=48 > Emo.abyss.k96.run.log 2>&1 &
```

This initial run failed (*k*=96), likely due to the need of more memory. Consider running it again at a later date. Also, ensure all files are decompressed.


### *De novo* assembly using IDBA-UD
IDBA is the basic iterative de Bruijn graph assembler for second-generation sequencing reads. IDBA-UD, an extension of IDBA, is designed to utilize paired-end reads to assemble low-depth regions and use progressive depth on contigs to reduce errors in high-depth regions is a fuzzy Bruijn graph approach to long noisy reads assembly. Source code was cloned from Github from [IDBA](https://github.com/loneknightpy/idba).

```bash
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_250_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_paired.fa
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_250_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_paired.fa
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_350_1_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_paired.fa
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_350_2_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_paired.fa
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_350_3_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_paired.fa
~/genome/src/idba/bin/fq2fa --filter --merge Ecdeiocolea_monostachya_350_4_gDNA_forward_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_reverse_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_paired.fa

cat Ecdeiocolea_monostachya_250_1_gDNA_paired.fa Ecdeiocolea_monostachya_250_2_gDNA_paired.fa Ecdeiocolea_monostachya_350_1_gDNA_paired.fa Ecdeiocolea_monostachya_350_2_gDNA_paired.fa Ecdeiocolea_monostachya_350_3_gDNA_paired.fa Ecdeiocolea_monostachya_350_4_gDNA_paired.fa > Ecdeiocolea_monostachya_gDNA_paired.fa
~/genome/src/idba/bin/idba_ud -r Ecdeiocolea_monostachya_gDNA_paired.fa -o Emo_idba --num_threads 46
```

### *De novo* assembly using SOAPdenovo2
`SOAPdenovo2` is a de novo de Bruijn graph assembler. Source code was cloned from Github from [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2).

```bash
soapdenovo2-63mer pregraph -s Emo.config -K 63 -R -p 40 -o Emo_soapdenovo2_k63 1>pregraph.log 2>pregraph.err
soapdenovo2-63mer contig -g Emo_soapdenovo2_k63 -R -p 46 1>contig.log 2>contig.err &
soapdenovo2-63mer map -s Emo.config -g Emo_soapdenovo2_k63 -p 46 1>map.log 2>map.err &
soapdenovo2-63mer scaff -g Emo_soapdenovo2_k63 -F -p 46 1>scaff.log 2>scaff.err &
```

### *De novo* assembly using minia
`minia` is a short-read assembler based on a de Bruijn graph with low memory requirements. Source code was cloned from Github from [minia](https://github.com/GATB/minia).

```bash
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 21 -out Emo.minia.k21 > Emo.minia.k21.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 31 -out Emo.minia.k31 > Emo.minia.k31.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 41 -out Emo.minia.k41 > Emo.minia.k41.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 51 -out Emo.minia.k51 > Emo.minia.k51.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 61 -out Emo.minia.k61 > Emo.minia.k61.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 71 -out Emo.minia.k71 > Emo.minia.k71.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 81 -out Emo.minia.k81 > Emo.minia.k81.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 91 -out Emo.minia.k91 > Emo.minia.k91.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 97 -out Emo.minia.k97 > Emo.minia.k97.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 99 -out Emo.minia.k99 > Emo.minia.k99.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 101 -out Emo.minia.k101 > Emo.minia.k101.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 103 -out Emo.minia.k103 > Emo.minia.k103.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 105 -out Emo.minia.k105 > Emo.minia.k105.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 107 -out Emo.minia.k107 > Emo.minia.k107.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 111 -out Emo.minia.k111 > Emo.minia.k111.log 2>&1 &
~/genome/bin/minia -in Ecdeiocolea_monostachya_gDNA_paired.fa -kmer-size 121 -out Emo.minia.k121 > Emo.minia.k121.log 2>&1 &
```

### Assessment of *de novo* genome assemblies using Illumina

```bash
assembly-stats -t Emo.minia.k*.contigs.fa > Emo.minia.stats.txt
```


```bash
./scripts/run_BUSCO.py -i Emo.minia.k21.contigs.fa -o Emo.minia.k21.busco --lineage_path embryophyta_odb9 -m genome -c 16 > Emo.minia.k21.busco.log 2>&1 &
```

```bash
export AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config/
```

```bash
~/emo/src/hisat2-2.1.0/hisat2-build -p 46 Emo.minia.k$1.contigs.fa Emo.minia.k$1.contigs
~/emo/src/hisat2-2.1.0/hisat2 --max-intronlen 20000 -k 1 -p 46 --no-softclip --dta-cufflinks -x Emo.minia.k$1.contigs -1 Ecdeiocolea_monostachya_$2_RNAseq_forward_paired.fq.gz -2 Ecdeiocolea_monostachya_$2_RNAseq_forward_paired.fq.gz -S Emo.minia.k$1.contigs_E01flower_RNAseq.cufflinks.sam
samtools view -F 4 -Shub Emo.minia.k$1.contigs_$2_RNAseq.cufflinks.sam > Emo.minia.k$1.contigs_$2_RNAseq.cufflinks.bam
samtools sort -o Emo.minia.k$1.contigs_$2_RNAseq.cufflinks.sorted.bam Emo.minia.k$1.contigs_$2_RNAseq.cufflinks.bam
cufflinks/cufflinks -p 4 Emo.minia.k$1.contigs_$2_RNAseq.cufflinks.sorted.bam
```

```
cuffmerge -s reconciled_assembly_v4.fa -p 4 gtf_files.txt
gffread merged_asm/merged.gtf -g reconciled_assembly_v4.fa -w transcripts_CGS.fa
TransDecoder.LongOrfs -t transcripts_CGS.fa
./interproscan.sh --output-dir . --input longest_orfs.pep --iprlookup --seqtype p --appl Coils,Gene3D,ProSitePatterns,Pfam,PANTHER,SUPERFAMILY > longest_orfs_interproscan.log 2>&1 &

```

## Nanopore sequencing of *Ecdeiocolea monostachya*

### *De novo* assembly using Wtdbg2 
Wtdbg2 is a fuzzy Bruijn graph approach to long noisy reads assembly. Source code was cloned from Github from [wtdbg2](https://github.com/ruanjue/wtdbg2).

# Emo
Genome sequencing and annotation of *Ecdeiocolea monostachya*

## Illumina sequencing of *Ecdeiocolea monostachya*
Sheath tissue from accession E001 was collected, genomic DNA extracted, and DNA sent to Novogene for sequencing. Illumina paired end libraries were generated using inserts of 250bp (DSW66921) and 350bp (DSW66909-V) (150 bp reads). All *de novo* assemblies were carried out using an m5a.12xlarge (48 CPU, 192 GB RAM) instance on Amazon AWS EC2.

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


### *De novo* assembly using edena
edena 3.131028 was used for *de novo* genome assembly.

```bash
edena -paired Ecdeiocolea_monostachya_250_1_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_250_1_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_250_2_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_1_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_2_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_3_gDNA_edena_reverse_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_edena_forward_paired.fq Ecdeiocolea_monostachya_350_4_gDNA_edena_reverse_paired.fq -p Emo_edena_v1 -nThreads 48 > Emo.edena.run.log 2>&1 &
edena -e Emo_edena_v1.ovl -m 100 -p Emo_edena_v1_m100
```

### *De novo* assembly using ABySS
ABySS version 1.3.6 was used for *de novo* genome assembly.

```bash
abyss-pe k=96 name=Emo_abyss_v1 lib='pea peb pec ped pee pef' pea='Ecdeiocolea_monostachya_250_1_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_250_1_gDNA_reverse_paired.fq.gz' peb='Ecdeiocolea_monostachya_250_2_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_250_2_gDNA_reverse_paired.fq.gz' pec='Ecdeiocolea_monostachya_350_1_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_1_gDNA_reverse_paired.fq.gz' ped='Ecdeiocolea_monostachya_350_2_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_2_gDNA_reverse_paired.fq.gz' pee='Ecdeiocolea_monostachya_350_3_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_3_gDNA_reverse_paired.fq.gz' pef='Ecdeiocolea_monostachya_350_4_gDNA_forward_paired.fq.gz Ecdeiocolea_monostachya_350_4_gDNA_reverse_paired.fq.gz' np=48 > Emo.abyss.run.log 2>&1 &
```

This initial run failed, likely due to the need of more memory. Consider running it again at a later date. Also, ensure all files are decompressed.

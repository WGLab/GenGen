## Introduction

Given a list of genomic regions, the scan_region.pl program can be used to find overlapped regions between two files, or find regions that match specific genomic features (such as being annotated genes, conserved elements, segmental duplications, etc.). The program is highly multi-functional and can do a lot of similar tasks under the same umbrella. Note that for all program arguments, the single dash (-) and double dash (--) has the same effect, and the argument name can be abbreviated to the first a few letters as long as there is no ambiguity (so `--argument`, `-argument`, `--arg` and `-argum` has the same effect).

## Scan two sets of genomic regions

One of the major functionality of the scan_regoin.pl program is used to scan chromosome regions from a query-location-file against a DB-location-file. For example, the query file may contain a set of genomic regions under linkage peaks in one linkage study, yet the DB file contains previously reported linkage regions, and the goal is to identify a set of query linkage regions that have been previously reported. Another example is that the query file may contain a set of CNV calls, whereas the DB file may contain a set of centromeric regions, and the goal is to identify CNV calls in centromeric regions (potential artifacts).

### - File format

A sample query file is shown below (The chr prefix is optional but highly desirable in the beginning of each line):

```
[kai@beta ~/usr/kgenome/trunk/example]$ cat region1 
chr3:100-500 first query region
chr5:2000-3000 second query region
chr8:1000-4000 third query region
```

So the file has one genomic region per line, and each region is shown with a chromosome identifier and a start-end position. The first field in the line is the genomic region, while the following fields can be anything you want to add.

The DB file can be either in the same form as the query file, or in tab-delimited format (with the first three columns as chr, start, end):

```
[kai@beta ~/usr/kgenome/trunk/example]$ cat region2
chr3:400-1000 first DB region
chr8:900-5000 second DB region
```

When scanning two sets of genomic regions, the scan_region.pl program first reads in both files, sort the regions for each file and store them in memory, and then for each DB region, print out query regions that overlap with the DB region. Neither the query file nor the DB file need to be sorted by chromosome or positions. (This is different from scanning UCSC Genome Browser regions, which will be described below in another section.)

### - Basic usage

An example to find the query regions that overlap with DB regions is given below:

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2
chr8:1000-4000 third query region
chr3:100-500 first query region
```

In this example, out of the three regions in the query file, two of them overlap with those in the DB file, so they are printed out to the screen. Note that the strings after the region in each line is also included in the output.

### - the --overlap argument

If the `--overlap` argument is set, then only the overlapped portion of the regions are included in the output:

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -overlap
chr8:1000-4000 third query region
chr3:400-500 first query region
```

In the example above, the chr3:400-500 region in the first query is printed out, because this region is shared by the query (chr3:400-1000) and by the DB (chr3:100-500).

### - the --minquerydbratio and --maxquerydbratio arguments

The user can specify the `--minquerydbratio` and `--maxquerydbratio` argument, to filter the output based on the query/DB length ratio. If the query is too short or too long (with respect to the overlapping DB), it will not be printed out. For example, when scanning one set of CNV calls against another set of calls, setting "-minquerydbratio 0.5 -maxquerydbratio 2" ensure that the size of the query and DB are within 2 fold difference, so that one can claim that any overlap probably represent the same event (CNV calls may have boundary discordances).

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minquerydbratio .5
chr8:1000-4000 third query region
chr3:100-500 first query region

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minquerydbratio .7
chr8:1000-4000 third query region 

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -maxquerydbratio .7
chr3:100-500 first query region
```

The length of the first query region is 401, while the length of the first DB region is 601, with a ratio of ~0.67. So when -minquerydbratio of 0.7 is set, the first query region is no longer printed out. Similarly, when maxquerydbratio is set to 0.7, the third query region is not printed out as it's query/DB ratio is 3001/4101=0.73.

### - the --mindbfrac and --maxdbfrac argument

The user can specify the --mindbfrac and --maxdbfrac argument, so that the fraction of the overlap in DB region must be higher or lower than the criteria for the query to be printed out. These arguments are useful in some scenarios: for example, when scanning CNV calls against known cytogenetic abnormalities (such as a canonical 1.5Mb deletion on 15q13.3), setting the -mindbfrac ensures that a small CNV within this canonical region is not printed out and only relatively large CNV calls that cover most of the bases on 15q13.3 are regarded as the same as the canonical deletion.

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -mindbfrac 0.05
chr8:1000-4000 third query region
chr3:100-500 first query region

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -mindbfrac 0.2
chr8:1000-4000 third query region
```

For the first query region, the overlapped region is 101 bp (chr3:400-500), whereas the DB region is 601 bp, so the fraction of overlap in DB region is 101/601=0.17. When -mindbfrac is set as 0.2, the first query region is not printed out.

### - the --minqueryfrac and --maxqueryfrac argument

The user can specify the --minqueryfrac and --maxqueryfrac argument, so that the fraction of the overlap in query region must be higher or lower than the criteria for the query region to be printed out. These arguments are useful in some scenarios: for example, when scanning CNV calls against "bad regions" that tend to produce artifacts (such as immunoglobulin region, or centromeric regions, or telomeric regions), it is always a good idea to set -minqueryfrac, to ensure that a very large query is not marked as "bad" CNV call just because it by chance overlap with a "bad region".

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minqueryfrac 0.1
chr8:1000-4000 third query region
chr3:100-500 first query region
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minqueryfrac 0.5
chr8:1000-4000 third query region
```

In the above example, the overlapped region is 101 bp, whereas the query region is 401 bp, so the fraction of overlap is 0.25. When -minqueryfrac is set as 0.5, the first query region is not printed out. Note that the program is not smart enough to handle situations where two query regions match the same DB region (one at head and one at tail): both query regions will be considered separately rather than jointly when calculating query fraction.

### - the -minoverlap argument

The -minoverlap argument requires that the overlapped portion must be larger than the specified fraction of either query or DB region. In a certain sense, it combines the -minqueryfrac and -mindbfrac arguments but with an OR operation rather than AND operation.

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minoverlap .1
chr8:1000-4000 third query region
chr3:100-500 first query region

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minoverlap .2
chr8:1000-4000 third query region
chr3:100-500 first query region

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region1 region2 -minoverlap .5
chr8:1000-4000 third query region
```

In the above example, the first query region is 101bp, and the overlap fraction for query and DB is 0.25 and 0.17, respectively. When -minoverlap of 0.2 is set, the criteria is still met for query so the first query is printed out. When -minoverlap of 0.5 is set, the first query does not meet the criteria and does not show up in the output.

## Scan genomic regions against annotated genes

Another functionality of the scan_region.pl program is to scan regions against annotated genes or exons. For example, suppose there is a list of copy number variation (CNV) regions, and a user can use the program to specifically pick out those gene-disrupting CNVs as well as exonic CNVs. Another example is to map SNPs in a given array to either an overlapping gene (if the SNP is located within the gene) or its closest gene (if the SNP is located in intergenic regions).

For this analysis, the user needs to download gene annotations from UCSC Genome Browser. Two types of annotations are widely used: RefGene and UCSC Gene. The former is well annotated but may miss some genuine transcripts, yet the later consists of many computationally predicted genes and is more comprehensive. Note that despite the names, both annotations are "transcript" annotations, rather than real "gene" annotations. (In contrast, the Ensembl does provide gene annotations, in addition to transcript annotations.)

### - Scan region against RefGene annotation

If the user has never downloaded the gene annotation before, below are the commands for downloading Human RefGene annotations (NCBI36 build):

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refLink.txt.gz
gunzip refGene.txt.gz; mv refGene.txt hg18_refGene.txt
gunzip refLink.txt.gz; mv refLink.txt hg18_refLink.txt
```

Note that in the above command, I rename the files with the hg18 prefix, to indicate the version and the organism that the files are for. This is important when dealing with multiple species and multiple genome builds to avoid confusion.

```
[kai@beta ~/usr/kgenome/trunk/example]$ cat region3 
chr3:4000000-4100000 CNV1
chr8:8000000-9000000 CNV2

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region3 hg18_refGene.txt -refgene -reflink hg18_refLink.txt 
chr3:4000000-4100000 CNV1 NOT_FOUND NOT_FOUND
chr8:8000000-9000000 CNV2 CLDN23,MFHAS1,PRAGMIN,THEX1 0
```

In the above example, the input file contains two genomic regions, and we scan these two regions against the RefGene database. The -refgene argument specifies that the DB file is in RefGene format. In the output, the first region does not overlap with any gene, so the "NOT_FOUND" is printed; the second region overlap with four genes (comma-separated). The number "0" indicates the distance between the CNV and the genes, and since they do ovlerap, the distance is zero. (later on we'll see how to expand the query so genes closeby can be scanned.)

The `-reflink hg18_refLink.txt` argument provides a "link" file that annotates the standard gene name. If you do not use the -reflink argument, then something like "NM_001080826,NM_004225,NM_153332,NM_194284" will be printed out as the gene name. However, by adding `--name2` argument to the command, the gene name annotated in the RefGene file will be printed out (they may differ from the gene name in the RefLink file which is more authentic, but usually they are the same).

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region3 ~/lib/ucsc/refgene/hg18_refGene.txt -refgene 
chr3:4000000-4100000 CNV1 NOT_FOUND NOT_FOUND
chr8:8000000-9000000 CNV2 NM_001080826,NM_004225,NM_153332,NM_194284 0

[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region3 ~/lib/ucsc/refgene/hg18_refGene.txt -refgene -name2
chr3:4000000-4100000 CNV1 NOT_FOUND NOT_FOUND
chr8:8000000-9000000 CNV2 CLDN23,MFHAS1,PRAGMIN,THEX1 0
```

Since the first region above does not overlap with any gene, we are interested in knowing what is the closest gene. This requires "expand" the query region to match the closest genes.

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region3 ~/lib/ucsc/refgene/hg18_refGene.txt -refgene -reflink ~/lib/ucsc/refgene/hg18_refLink.txt -expandmax 1m
chr3:4000000-4100000 CNV1 LRRN1 135613
chr8:8000000-9000000 CNV2 CLDN23,MFHAS1,PRAGMIN,THEX1 0
```

The above example tries to expand the query regions to a maximum of 1MB and try to find the closest gene: it found LRRN1 which is 135kb away from the first region.

The command can be very useful in finding the closest genes for each SNP in a given SNP array: simply prepare an input file with two fields per line (region and SNP identifier), then scan against the RefGene database with -expandmax argument and then with some simple reformatting you'll have a file mapping SNP to genes as well as their distances.

Sometimes it may be useful to use the -expanddb argument as well: this basically expand the DB region to be scanned: for example, a DB region of chr3:400-1000 will be treated as if it were chr3:300-1100 when "-expanddb 100" is used. A SNP at chr3:350-350 may be outside of an annotated RefGene transcript region but can be treated as having distance of 0 with the gene if "-expanddb 100" is used.

### - Scan region against the UCSC Gene annotation

To use the UCSC Gene annotation rather than RefGene annotation, the command is slightly different (note that --knowngene flag is used, and the --kgxref is used to specify the "link" file):

```
[kai@beta ~/usr/kgenome/trunk/example]$ scan_region.pl region3 ~/lib/ucsc/kg/hg18_knownGene.txt -knowngene -kgxref ~/lib/ucsc/kg/hg18_kgXref.txt -expandmax 1m
chr3:4000000-4100000 CNV1 UNQ3037 0
chr8:8000000-9000000 CNV2 BC035792,CLDN23,CR593780,FLJ00269,MFHAS1,THEX1 0
```

As we can see here, the UCSC Gene annotation thinks that the first region actually overlaps with a gene, and the second region overlaps with six genes (rather than four genes).

### - Scan RefGene coding sequence

In some cases, one might be particularly interested in finding overlap between a query and the coding sequence of a gene (from the first exon to the last exon), rather than the entire transcript. In this case, a refGene file can be provided, but the `--refcds` argument (as opposed to `--refgene`) should be specified.

### - Scan RefGene exons

In many cases, especially when analyzing copy number variations, one might be particularly interested in finding overlap between a query and the exon of a gene, excluding all introns. In this case, a refGene file can be provided, but the `--refexon` argument (as opposed to `--refgene`) should be specified.

## Scan genomic features annotated in UCSC genome browser databases

### - General overview of UCSC genome browser tables

When processing DB files as UCSC tables, this program works by first reading all information from the query file and store them in memory, then scan each line of the DB file to find overlaps. THEREFORE, THE DB FILE CANNOT CONTAIN TWO LOCATIONS THAT SHARE OVERLAPS. IN ADDITION, THE TEMPLATE FILE MUST BE SORTED FIRST BY CHROMOSOME THEN BY START LOCATION. The chromosome can be sorted alphabetically (since there will usually be X, Y and MT chromosomes), but start location should be always sorted numerically. Also note that by convention, positions annotated in UCSC databases are 0-based rather than 1-based; this has been taken care of by this program.

These DB files can be downloaded from the UCSC Genome Browser Annotation Database. Generally speaking, I recommended saving the DB file with the same file name as the original database file, but prefixing with the version and organism of the corresponding genome. For example, for human genome, the DB-file name should start with hg17_ or hg18_, to eliminate confusions when dealing with different genome builds. For example, after downloading the file evofold.txt.gz for human May 2004 assembly from the UCSC genome browser database, you can use the sort command to rename the file as hg17_evofold.sorted (see details below).

### - MCE format

The `--mce_flag` argument specifies that template file is in MCE (most conserved element) format from the UCSC genome browser. You can download a MCE template file using a simple command:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/phastConsElements17way.txt.gz
gunzip phastConsElements17way.txt.gz
```

This will generate the phastConsElements17way.txt file. Next perform sorting:

```
sort -k 2,2 -k 3,3n phastConsElements17way.txt > hg18_phastConsElements17way.sorted
```

A sample sorted MCE file is shown below:

```
585 chr1 1865 1943 lod=31 307
585 chr1 2039 2096 lod=104 444
585 chr1 2473 2566 lod=179 505
585 chr1 2873 2917 lod=107 447
585 chr1 3081 3135 lod=58 378
```

The `--score_threshold` argument operates on the fifth field in the line (such as "lod=31"), and the `--normscore_threshold` argumet operates on the sixth field (such as 307).

To scan the queryfile against a MCE file:

```
scan_region.pl queryfile hg18_phastConsElements17way.sorted -mce
```

### - EvoFold format

The EvoFold file from UCSC genome browser contains genomic segments that are predicted to retain stable RNA structural folds in different species. You can download a EvoFold DB file for human NCBI36 genome assembly using

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/evofold.txt.gz
gunzip evofold.txt.gz
sort -k 2,2 -k 3,3n evofold.txt > hg18_evofold.sorted
```

The first a few lines of a sorted EvoFold table is shown below:

```
591 chr1 886773 886798 608_0_-_96 96 - 25 (((((.....(((...))).))))) 0.97,0.98,0.99,0.99,0.9,0.78,0.87,0.99,0.98,0.98,0.21
591 chr1 888417 888435 617_0_+_156 156 + 18 ((((((......)))))) 0.97,0.99,0.99,1.0,0.99,0.92,1.0,1.0,0.99,1.0,1.0
```

The `--score_threshold` and `--normscore_threshold` arguments both operate on the sixth field in the line (such as "96").

### - TFBS format

The TFBS file from UCSC genome browser contains predicted transcription factor binding sites (TFBS), based on positional weight matrices (PWM) from the TRANFAC database. You can download a TFBS DB-file for human NCBI36 genome assembly using:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/tfbsConsSites.txt.gz
gunzip tfbsConsSites.txt.gz
sort -k 2,2 -k 3,3n tfbsConsSites.txt > hg18_tfbsConsSites.sorted
```

The first a few lines of a sorted TFBS table is shown below:

```
585 chr1 1832 1848 V$ARP1_01 818 - 2.01
585 chr1 3211 3232 V$NRSF_01 749 - 2.19
585 chr1 3603 3623 V$YY1_02 790 - 1.93
585 chr1 3949 3956 V$NKX25_01 1000 - 2.48
585 chr1 3996 4014 V$CART1_01 798 - 1.75
```

The `--score_threshold` argument operates on the eighth field in the line (such as "2.01"), and the `--normscore_threshold` argumet operates on the sixth field (such as 818).

To decode the strange values like V$ARP1_01 and V$NRSF_01, you can additionally download the TF annotation file at http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/tfbsConsFactors.txt.gz and examine the file for detailed information on the transcription factors binding to the sites.

### - wgRna format

The wgRna table from UCSC genome browser contains microRNA and small nucleolar RNA information. You can download a wgRna DB-file for human NCBI36 genome assembly using:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/wgRna.txt.gz
gunzip wgRna.txt.gz
sort -k 2,2 -k 3,3n wgRna.txt > hg18_wgRna.sorted
```

The first a few lines of a sorted wgRna table is shown below:

```
593 chr1 1092346 1092441 hsa-mir-200b 960 + 1092402 1092425 miRna
593 chr1 1093105 1093195 hsa-mir-200a 960 + 1093120 1093142 miRna
593 chr1 1093105 1093195 hsa-mir-200a 960 + 1093158 1093180 miRna
593 chr1 1094247 1094330 hsa-mir-429 960 + 1094297 1094319 miRna
611 chr1 3467118 3467214 hsa-mir-551a 480 - 3467133 3467154 miRna
```

The `--score_threshold` and `--normscore_threshold` arguments both operate on the sixth field in the line (such as "960").

### - SegDup format

The segDup table from UCSC genome browser contains chromosome regions with segmental duplications, as well as the "target" regions that match these duplications. You can download a segdup DB-file for human NCBI36 genome assembly using:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/genomicSuperDups.txt.gz
gunzip genomicSuperDups.txt.gz
sort -k 2,2 -k 3,3n genomicSuperDups.txt.gz > hg18_genomicSuperDups.sorted
```

The first a few lines of a sorted segdup table is shown below:

```
585 chr1 465 30596 No.1139,chr2:114046528 570211 - chr2 114046528 114076216 243018229 1139 1000 N/A Filtered N/A N/A build35/align_both/0012//both064929 30176 43 531 29645 29282 363 128 235 0.987755 0.986324 0.012346 0.0123574
585 chr1 486 30596 No.2251,chr9:844 582086 + chr9 844 30515 138429268 2251 1000 N/A Filtered N/A N/A build35/align_both/0013//both065185 30137 14 491 29646 29474 172 64 108 0.994198 0.993729 0.00582435 0.00582657
```

The `--score_threshold` and `--normscore_threshold` arguments both operate on the sixth field in the line (such as "570211").

### - Generic format

The program can also take a generic UCSC genome browser format, where the first four columns are index, chromosome, start and end position, with zero-based (half-open) position counting schemes. This is indeed how the MCE, EvoFold and other tracks are organized in the database for the newest genome build (but not necessarily the case for older genome build, and not necessarily true for non-human species). The --anno_flag argument should be used to flag that a generic format is expected for the DB file.

### - phastCons format

Despite its name, do not confuse it with the MCE (most conserved element) track of the genome browser. The phastcons file merely contains conservation scores for each genome position that can be aligned to other species (in comparison, the MCE file contains the conservation scores for the top 5% most conserved genomic regions in the genome). For human genome, much more than 5% of the genome (probably ~50%) can be actually aligned, although only 5% of the genome can be considered as conserved during evolution.

In other word, the MCE annotates conserved sequence segments, but the phastcons annotates individual base pairs that can be aligned.

Since the file contains a score for each genomic position, the file size is extremely large (typically >10GB), and the scanning takes a lot of time (up to a whole day). In most circumstances, I recommended using the MCE file for sequence conservation analysis, since the phastcons score is less meaningful than the MCE score, which use a hidden Markov model to identify and summarize a genomic region that is conserved during evolution.


## Introduction

A few people asked me the best way to identify SNPs in a particular array (for example, Illumina array with 550K markers), that are best proxy for a list of candidate SNPs (for example, dozens of SNPs reported in a meta-analysis paper using imputation data). The `find_ld_snp.pl` program is designed to address this issue. It finds the best proxy markers for a given list of candidate markers based on linkage disequilibrium from a GWAS data set (in PLINK format).

## Instruction

The program requires that PLINK be installed in the system first, since it calls PLINK for LD calculation. Typing the command without arguments will print out a simple help message.

```
[kaiwang@biocluster ~/]$ find_ld_snp.pl 
Usage:
     find_ld_snp.pl [arguments] <query-SNP-list-file> <candidate-SNP-list-file> <PLINK-binary-prefix>

     Optional arguments:
            -v, --verbose                   use verbose output
            -h, --help                      print help message
            -m, --man                       print complete documentation
        
     Function: find SNPs in candidate list that are best proxy for SNPs in query list
 
     Example: find_ld_snp.pl querylist humanhap550.snplist hapmap_ceu_r23a
```

It takes three input files: a query SNP list file which has one SNP per line, a larger  candidate SNP list file file containing all candidate SNPs to be selected, and a PLINK file prefix. For example, if you have `hapmap_CEU_r23a.bim`, `hapmap_CEU_r23a.fam`, `hapmap_CEU_r23a.bed` file, you just need to specify `hapmap_CEU_r23a`.

Now issue the command:

```
[kaiwang@cc ~/]$ find_ld_snp.pl speliotes.snplist snplist.illumina ~/lib/hapmap/plink1/hapmap_CEU_r23a
```

it will tell what are the best SNPs to use in the Illumina array, their r2 measures, and their distance to each other. Use the "-m" argument to read the manual for more information.
## GenGen main package

The stable version of GenGen main package can be downloaded below. The main package DOES NOT include the necessary library files.

[gengen.tar.gz](https://github.com/WangGenomicsLab/GenGen/archive/v1.0.0.tar.gz)

Additional update 2010Jun24: A few people asked me the best way to identify SNPs in a particular array (for example, Illumina array with 550K markers), that are best proxy for a list of candidate SNPs (for example, dozens of SNPs reported in a meta-analysis paper using imputation data). The [find_ld_snp.pl](http://www.openbioinformatics.org/gengen/download/find_ld_snp.pl) program is designed to do this: just run "find_ld_snp.pl snplist.candidate snplist.550karray hapmap_CEU_r23a" and it will tell what are the best SNPs to use in the Illumina array, their r2 measures, and their distance to each other. Use the "-m" argument to read the manual first.

## Supplementary library files

SNP-gene mapping files for various Illumina and Affymetrix SNP arrays, as well as SNPTable file for commonly used arrays, can be downloaded below.

[gengenlib.tar.gz](http://www.openbioinformatics.org/gengen/download/gengenlib.tar.gz)

Among the SNP-gene mapping files, hh550* are for HumanHap550 array, hhall\* are for other Illumina Infinium arrays, affy500k* are for Affymetrix 500K arrays, gw6* are for Affymetrix genome-wide 6.0 arrays. You can contact me directly to get SNP-gene mapping annotation for other arrays. Or, you can follow the instructions to generate a custom annotation file using the scan_region.pl program included in the package. Among the SNPTable files, hh550_610.snptable is for Illumina 550K or 610K arrays (all versions), and hh1mv13.snptable is for Illumina 1M arrays (version 1 or 3).

## SNPTable file for Illumina arrays and Affymetrix arrays

- [CNV370 version 1 array](http://www.openbioinformatics.org/gengen/download/cnv370v1_snptable.txt.gz)
- [CNV370 version 3 array](http://www.openbioinformatics.org/gengen/download/cnv370v3_snptable.txt.gz)
- [Cyto12 version 1 array](http://www.openbioinformatics.org/gengen/download/cyto12v1_snptable.txt.gz)
- [HumanHap1M version 1 array](http://www.openbioinformatics.org/gengen/download/hh1mv1_snptable.txt.gz)
- [HumanHap1M version 3 array](http://www.openbioinformatics.org/gengen/download/hh1mv3_snptable.txt.gz)
- [HumanHap550 version 3 array](http://www.openbioinformatics.org/gengen/download/hh550v3_snptable.txt.gz)
- [Human610 version 1 array](http://www.openbioinformatics.org/gengen/download/hh610v1_snptable.txt.gz)
- [Human660 array](http://www.openbioinformatics.org/gengen/download/hh660_snptable.txt.gz)
- [Omni version 1 array](http://www.openbioinformatics.org/gengen/download/ho1v1_snptable.txt.gz)

Contact me with issues/bugs in these files as they were not tested by me.



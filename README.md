![Logo](http://genomics.usc.edu/images/lab/logo/gengen.png)

## GenGen

GenGen is a suite of free software tools to facilitate the analysis of high-throughput genomics data sets. The package is currently a work-in-progress and infrequently updated.

Some of the most useful and more stable programs in this package are:

- calculate_association.pl: perform genome-wide association analysis for cohorts with case-control design and trio design, with ability to handle PLINK-formatted TPED/TFAM files and handle phenotype permutations.
- calculate_gsea.pl, combine_gsea.pl: use pathway-based approaches to test the significane of genes within a pathway for genome-wide association studies (see reference below).
- scan_region.pl: scan two sets of genomic regions to find overlap (for example, scan most conserved genomic elements, scan transcription factor binding sites, scan segmental duplication regions, scan predicted microRNA sites, scan pseudogenes, scan EvoFold regions, etc.
- convert_mach.pl, convert_mach2snptest.pl, combine_snptest.pl: convert MACH-imputed genotype data sets into standard PED/MAP format or SNPTEST format for subsequent association test.
- convert_bim_allele.pl: convert SNP allele coding between Illumina A/B alleles, Illumina 1/2 alleles, Illumina TOP strand alleles and dbSNP forward strand alleles.
- retrieve_ensembl.pl, retrieve_variation.pl, retrieve_compara.pl: query the Ensembl databases to get the desired information on genes, transcripts, variations, cross-references, domains, predicted regulatory elements and so on.

If you like this repository, please click on the "Star" button on top of this page, to show appreciation to the repository maintainer. If you want to receive notifications on changes to this repository, please click the "Watch" button on top of this page.

## Reference

- Wang K, Li M, Bucan M. [Pathway-based approaches for analysis of genome-wide association studies](http://www.cell.com/ajhg/abstract/S0002-9297(07)63775-6). _**American Journal of Human Genetics**_, 81:1278-1283, 2007
- Wang K, Li M, Hakonarson H. [Analysing biological pathways from genome-wide association studies](http://www.nature.com/nrg/journal/v11/n12/full/nrg2884.html). _**Nature Reviews Genetics**_, 11:843-854, 2010

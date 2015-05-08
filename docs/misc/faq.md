1. **What is the difference between the various snpgenemap file with or without the "exp20k" suffix?**

    A regular snpgenemap file gives the SNP, it's closest RefSeq transcript and their distances in three fields each line. The \*.exp20k file defines a gene as the corresponding transcript plus the upstream and downstream 20kb, so the distance in this file will be 20kb less than those in the regular snpgenemap file. That's the only difference. For additional details see here.

2. **The alleles are coded as 1 and 2 in my genotype data, so how can I convert them to ACGT alleles?**

    The file may be coded as Illumina's A/B alleles, or be coded by PLINK using `--recode12` argument (see additional explanation [here](../tutorial/coding.md)). In the former case, you can use the `convert_bim_allele.pl` program for conversion to ACGT; in the latter case, you can do noting with it and it is best to ask whoever generated the data for original genotypes.

3. **I have a new SNP array, so how can I find the closest genes for each of the SNPs?**

    The `scan_region.pl` program can solve the problem. Detailed commands can be found [here](../tutorial/scan.md).

4. **Can I convert MACH-imputed data to some kind of dosage format for association analysis?**

    Unlike the PED format, there is no standard format for "allele dosage". It is probably best to convert and combine all the batches of imputed data to SNPTEST format, and then write a simple script to sum the dosage for each SNP in cases and controls (by processing each line in the SNPTEST file). I do have such a script but it is not mature enough to be included in the package and if you really want to use it, you can email me to get it.



## Introduction

MACH is one of the most popular genotype imputation software. Typically, it takes phased haplotype data from a reference population, then take genotype data for a group of subjects genotyped on a particular SNP array, then try to impute or predict the genotypes on markers that were not on the SNp arrays. The MACH software produces several output files, including mlinfo, mlgeno, mldose, mlqc, mlprob files. These files cannot be directly used in subsequent association analysis, unless converted to appropriate file formats, such as the standard PED/MAP format, or the SNPTEST format (which takes into account of genotype imputation uncertainty when performing association test).

The command line arguments for the `convert_mach.pl` program is listed below:

```
[kai@beta ~]$ convert_mach.pl
Usage:
convert_mach.pl [arguments] <geno-file> <info-file>

Optional arguments:
-v, --verbose use verbose output
-h, --help print help message
-m, --man print complete documentation
--legendfile <file> HapMap legend file for the chromosome
--pedinfo <file> ped information file or PLINK fam file (first 6 columns will be used)
--prefix <string> output prefix (default:<geno-file>
--rsq <float> rsq threshold to define high-qulaity SNP and write to extract file (default: 0.3)
--qc_threshold <float> genotype posterior probability threshold to output (default: 0.9)

Function: convert MACH imputation results to standard PED and MAP files for association analysis.

Example: convert_mach.pl chr22_step2.mlgeno chr22_step2.mlinfo chr22_step2.mlqc -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt -ped chr22.ped -prefix chr22_step2
```

The command line arguments for the `convert_mach2snptest.pl` program is listed below:

```
[kai@beta ~]$ convert_mach2snptest.pl
Usage:
convert_mach2snptest.pl [arguments] <prob-file> <info-file>

Optional arguments:
-v, --verbose use verbose output
-h, --help print help message
-m, --man print complete documentation
--prefix <string> output prefix (default:<geno-file>)
--rsq <float> rsq threshold to write SNPs to exclude file (default: 0.3)
--legendfile <file> HapMap legend file used in MACH imputation
--mapfile <file> MAP file (only useful when imputation results contains markers not annotated in legend file)
--keepfile <file> a two-column file specifying subjects to output
--blocksize <int> size of each reading block (default=100,000,000, reduce if out-of-memory)

Function: convert MACH imputation results to SNPTEST input files (a gen file, a sample file, and an exclude file).

Example: convert_mach2snptest.pl chr22.mlprob chr22_step2.mlinfo -prefix output -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt
convert_mach2snptest.pl chr22.mlprob chr22_step2.mlinfo -prefix output -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt -rsq 0.6 -keep caseid
```

The command line argument for the `combine_snptest.pl` program is listed below:

```
[kai@beta ~]$ combine_snptest.pl
Usage:
combine_snptest.pl [arguments] <prefix1 | ...>

Optional arguments:
-v, --verbose use verbose output
-h, --help print help message
-m, --man print complete documentation
--prefix <string> output prefix (default:<combine_prefix1>)
--keepfile <file> a two-column file specifying subjects to output
--gen_suffix <string> gen file suffix (default: gen)
--sample_suffix <string> sample file suffix (default: sample)

Function: combine multiple SNPTEST files (*.gen and *.sample) into one single file

Example: combine_snptest.pl file1 file2 -keep caseid.keep -prefix combine
#the above command combines file1.gen, file1,sample, file2.gen, file2.sample, and output combine.gen and combine.sample based on subject in caseid.keep file

combine_snptest.pl file1 -keep caseid.keep -prefix caseonly
#the above command retrieve a subset of subjects (based on caseid.keep file) from file1.gen and file2.sample, then write to new file caseonly.gen and caseonly.sample

``` 

## Convert MACH-imputed data to PED/MAP format

The program uses three files from MACH imputation runs: mlgeno, mlinfo and mlqc. In addition, it needs to read the orginal haplotype file used in imputation analysis (to figure out marker chromosome and location information), as well as the original PED file for genotype markers (to figure out subject information such as sex and disease status). The `--prefix` argument specifies the output file prefix.

An example file is shown below:

```
[kai@beta ~/project/]$ convert_mach.pl chr22_step2.mlgeno chr22_step2.mlinfo chr22_step2.mlqc -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt -ped chr22.ped -prefix imp
```

After running the command, a imp.ped and imp.map file will be generated and can be used in association test. The genotype calls are based on a QC threshold of 0.9, that is, a call is generated if the posterior probability is higher than 0.9 (otherwise non-call will be generated for the corresponding genotype). Furthermore, any marker that has Rsq measure less than 0.3 will be written to an imp.exclude file (which can be used by `--exclude` argument in PLINK during association analysis). The Rsq measure estimates the squared correlation between imputed and true genotypes, with smaller values indicating more difficult SNPs to impute due to lack of correlation information with surrounding markers. The QC threshold and Rsq threshold can be controlled by the `--qc_threshold` and `--rsq` argument, respectively.

Note that in the above command, the --ped argument is used (chr22.ped is the file used in imputation), which serves the only purpose of providing the first six columns to the final imp.ped file. Therefore, you can also use PLINK-formatted fam file for this argument as well, as only the first six columns are used.

In a typical imputation run, an input file is split into multiple parts (each with maybe ~500 subjects for one chromosome), and then imputation is performed on each part separately. Once all imputation is done, the PED/MAP file can be generated for each part, and then one can simply concatenate all PED file for a given chromosome together into one single large PED file (by the "cat" command in Linux/Unix system), since all the MAP files for the same chromosome are identical. Then one can also cat all the \*.exclude file together, take the unique ones and generate a new exclude file to be used in association analysis. As this step of combining data can be easily done by a simple system command, there is no need to write additional programs to handle it; however, as we will see later, this step does not work for the SNPTEST files, so we'll need the combine_snptest.pl program for the purpose of combining imputed genotypes for many subjects.

## Convert MACH-imputed data to SNPTEST format

The command line arguments for the convert_mach2snptest.pl program are almost identical to that of `convert_mach.pl` program without the `--pedinfo` argument, but the user can additionally use the `--keep` argument, to select only a subset of samples (such as cases or controls) to be written to output file. The reason is that SNPTEST requires a separate genotype file for cases and for controls, so essentially, one has to run the same program twice but with differetn `--keep` files to generate genotypes for cases and controls, respectively.

```
[kai@beta ~/project/]$ convert_mach2snptest.pl chr22.mlprob chr22_step2.mlinfo -prefix imp -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt
```

The imp.gen and imp.sample file will be generated as a result of the above command. Note that this program may generate many temporary files, in the name of `prefix.block1.gen`, `prefix.block1.sample`, etc. Each of these temporary files contains information for one "block" of genotypes, with the default value of 100Mb for one block. If there is an out-of-memory error, the user can try to reduce the blocksize to 50000000 or even smaller. Typically, the default size works well for a computer with 2GB memory. All the temporary files will be removed after the proram finishes.

The output files have prefix of \*.gen and \*.sample, representing genotype data (with confidence scores) and sample data, respectively. Note that it is possible to use the gzip system command in Linux/Unix system to compress the large \*.gen files, as SNPTEST can accept gzipped input files.

When performing association analysis by SNPTEST, it is important to turn on the -proper argument so that genotype imputation uncertainty is taken into account. Otherwise SNPTEST will simply used the most likely genotype calls for the association test.

## Merge multiple SNPTEST files together

As mentioned above, most imputation runs will split the input genotype files into multiple parts, run imputation on each parts, then combine the results together, due to memory issues and also due to speed issues. SNPTEST files are not in "one-sample-per-line" format, so they cannot be directly concatenated into one single file and must require specific software for combining them. The combine_snptest.pl program can handle this easily. For example,

```
[kai@beta ~]$ combine_snptest.pl file1 file2 file3 file4 -prefix combine
```

The above command will attempt to combine genotype information from `file1.gen`, `file1.sample`, `file2.gen`, `file2.sample`, `file3.gen`, `file3.sample`, `file4.gen` and `file4.sample` together, and generate the combine.gen file and combine.sample file. The `--keep` argument can be also used here to extract a specific group of samples (such as cases or controls).

Furthermore, if only file1 is supplied to the program together with the `--keep` argument, then the program will effectively extract specific samples from the `file1.gen` and `file1.sample` file and generate a new file.


# Example: the QTL-MAS 2012 data
Refer to https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-8-S5-S1 for how the data set was simulated.

# Construct a GRM
Input files: genotype files (geno.bin, geno.mrk, and geno.indi) and SNP info file (all.snp_info.csv)
```
./bfmap --compute_grm 1 --binary_genotype_file geno --snp_info_file all.snp_info.csv --output grm1 --num_threads 10
```
Output files: GRM (grm1.grm.bin and grm1.grm.indi)

# Estimate heritability
Input files: phenotype file (phen.csv) and GRM files (grm1.\*)
```
./bfmap --varcomp --phenotype phen.csv --trait milk --binary_grm_file grm1 --output milk --num_threads 10
```
Output file: milk.varcomp.csv

Expected: heritability = (proportion,0.307879)

# GWAS
Input files: phenotype file, SNP info file, genotype files, and GRM files
```
./bfmap --assoc --phenotype phen.csv --trait milk --snp_info_file all.snp_info.csv --snp_set single --binary_genotype_file geno --binary_grm grm1 --heritability 0.307879 --output milk.single.assoc --num_threads 10
```
Output file: milk.single.assoc.csv

Expected: top SNP = (1006499,1,40.1925,-0.498007,40.6905,4.6119e-23,0.840652|)

# Finemapping
Additional input file: SNP info file for a QTL region of interest (e.g. topQTL.snp_info.csv)
```
./bfmap --sss --phenotype phen.csv --trait milk --snp_info_file topQTL.snp_info.csv --snp_weight weight --binary_genotype_file geno --binary_grm grm1 --heritability 0.307879 --output milk.topQTL --num_threads 10
```
Output files: milk.topQTL.pip.csv and milk.topQTL.model.csv

Expected: (Posterior number of causal variants is 2.09235)

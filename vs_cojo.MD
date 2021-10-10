### download the GigaScience pig data to ~/giga_pig/
```sh
plink --bfile ~/giga_pig/3000_gwas_ok --from-mb 20 --to-mb 25 --chr 1 --maf 0.01 --make-bed --out chr1-20-25

plink --bfile chr1-20-25 --fill-missing-a2 --make-bed --out filled

plink --bfile filled --recode A --out add
```

### R code for simulating phenotypes
Randomly selecting SNPs that are in high positive LD, e.g., 1:22201760 (idx=447) and 1:22247692 (idx=464) in Chr1:20-25Mb

```R
# Copied from https://cnsgenomics.com/software/gcta/#MakingaGRM
# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  closeAllConnections()
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

# read GRM
bin = ReadGRMBin( "pig_grm" )
np = length(bin$diag)
G = matrix(0, nrow=np, ncol=np)
G[upper.tri(G)] = bin$off
G = G + t(G)
 
# simulate phenotypes with h2=0.5
raw = read.table("~/finemap/add.raw",head=T)
raw = raw[,-(1:6)]
diag(G) = bin$diag + 1
 
pheno = t(chol(G)) %*% rnorm(np)
pheno = pheno + raw[,447] * 0.3 + raw[,464] * 0.3
 
dat = cbind(bin$id, pheno)
colnames(dat) = c("IID","FID","QT")
write.table(dat, file="~/finemap/gcta.pheno.txt", row.names=F, col.names=F, quote=F)
write.csv(dat, file="~/finemap/bfmap.pheno.csv", row.names=F, quote=F)

```

### BFMAP analysis
```
./bfmap --compute_grm 1 --binary_genotype ~/giga_pig/3000_gwas_ok --snp_info gwas_ok.snp_info.csv --output bfmap_grm

./bfmap --varcomp --phenotype bfmap.pheno.csv --trait QT --binary_grm bfmap_grm --output hsq --num_threads 10

./bfmap --phenotype bfmap.pheno.csv --trait QT --snp_info chr1-20-25.snp_info.csv --binary_genotype 3000_gwas_ok --output test --num_threads 10 --binary_grm bfmap_grm --heritability 0.51
```

### GCTA LMM analysis
```
gcta64 --mlma --bfile filled --grm ~/giga_pig/pig_grm --pheno gcta.pheno.txt --out lmm --thread-num 10
```

### R code for producing .ma file
```R
assoc = read.table("~/finemap/lmm.mlma",head=T)
summary(assoc$p)
write.table( cbind(assoc[,-c(1,3)], 2797), file="lmm.ma", row.names=F, quote=F)
assoc[which(assoc$p == min(assoc$p)),]
assoc[c(447,464),]
```
### GCTA-COJO analysis
```sh
gcta64 --cojo-slct --bfile filled --cojo-file lmm.ma --cojo-p 1e-6 --out cojo
```

### Z-scores-based methods
```R
assoc = read.table("~/finemap/lmm.mlma",head=T)

# Computing Bayes factors using simple correlations
# when 447 and 464 are included
pve = 0.01
R = matrix(c(1,0.8247017,0.8247017,1),nrow=2)
z0 = assoc$b[c(447,464)] / assoc$se[c(447,464)]
z = z0
0.5 * (sum(z * solve(R + diag(2)/np/pve, z)) - log(det(R*np*pve + diag(2))) ) 
 
# when only 465 is included
R = matrix(c(1),nrow=1)
z = z0[2]
0.5 * (sum(z * solve(R + diag(1)/np/pve, z)) - log(det(R*np*pve + diag(1))) ) 

# check the difference between the two logBFs above

# Computing Bayes factors using relatedness-adjusted correlations
# when 447 and 464 are included
L = t(chol(G))
X = as.matrix(raw[,c(447,464)])
X = scale(X, center=T, scale=F)
Xstar = solve(L, X)
Rstar = t(Xstar) %*% Xstar
D = sqrt( diag(Rstar) )
z = z0
0.5 * (sum(z*D * solve(Rstar + diag(2)/pve, D*z)) - log(det(Rstar*pve + diag(2))) ) 

# when only 465 is included
z = z0[2]
Rstar = Rstar[2,2]
D = D[2]
0.5 * (sum(z*D * solve(Rstar + diag(1)/pve, D*z)) - log(det(Rstar*pve + diag(1))) ) 

# check the difference between the two logBFs above

```

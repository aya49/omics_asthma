---
title: "asthma"
author: "aya43@sfu.ca"
output:
  rmarkdown::github_document:
  number_sections: true
  # pandoc_args: [
  #   "--output=README.md"
  #   ]
---

```
# Cell Specific eQTL Analysis without Sorting Cells
# y ~ g*c # y=rnaseq, g=dna, c=cells
# NOTE: Yi = β0 + Xiβ1 + Ziβ2 + Wiβ3 + XiZiβ4 + XiWiβ5 + ZiWiβ6 + XiZiWiβ7 + ei
# Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
# Y ~ X * Z * W
# Y ~ (X + Z + W)^3
```

a few extra notes about the data

## feature types

### rnaseq

#### input raw fastq > output sample x base/isoform/gene matrix of counts/abundence

1. install packages used in order: **fastqc** > **star** > **rsem**
  - one of the easiest ways to do this is via [python](https://docs.python.org/3/using/unix.html#getting-and-installing-the-latest-version-of-python) and then [conda](https://conda.io/docs/user-guide/install/index.html#installing-conda-on-a-system-that-has-results/enrichr-python-installations-or-packages) after which you can create an environment to install your packages in. click on the links to see how to do this in your respective operating system. linux example below.

```{bash}
#install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

#create an environment to install the packages in
conda env create --name myenv_py3 --file environment.yaml

#install the packages via conda; searchable via https://anaconda.org/
conda install -c bioconda fastqc
conda install -c bioconda star
conda install -c bioconda rsem

#conda install -c r r #r packages for putting together the final matrix; optional
#conda install -c r r-dplyr 

#load the environment
source activate myenv_py3
```

2. given a folder of raw fastq files, run fastQC (if fastQC not done already), STAR (map transcripts to the organism's genome e.g. [GRC](https://www.ncbi.nlm.nih.gov/grc)), and RSEM (get abundence). to do so, first **write the raw fastq file paths into [config.yaml](./config.yaml)** (we write all the file names to keep a record, you can also specify a folder) and then run [Snakefile](./Snakefile) (adapted from the (Sakemake pipeline)[https://snakemake.readthedocs.io/en/stable/]). note: if there is more than one run of a sample (i.e. more than one fastq file per sample), put them together! linux example below.

```{bash}
#write raw fastq file paths into config.yaml
#run snakemake; its configurations and parameters are in config.yaml
snakemake -np --cores 6 #dry run
snakemake -p --cores 6 #actual run; change the number of cores as needed

```

3. given the <sample>.<genes/isoforms>.results rsem output files for each fastq file, ([nice link](https://ycl6.gitbooks.io/rna-seq-data-analysis/quantification_using_rsem1.html)) consolidate them into a count matrix. linux example below.

```{bash}
#run count1.sh
chmod u+x count1.sh
./count1.sh
```


### dna (affymetrix axiom)
- converting probe id to SNP id: [genechip array annotation files](https://www.thermofisher.com/ca/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-annotation-files.html)


#### imputation (incomplete)
we use **ricopili** to infer continuous probabilistic values for missing SNPs using reference HapMap data; original tutorial [here](https://sites.google.com/a/broadinstitute.org/ricopili/). there are results/enrichr options such as splink2, whichever one is suitable.

ensure that you have already installed software and reference data listed [here](https://sites.google.com/a/broadinstitute.org/ricopili/installation/external-software#TOC-External-Software-Packages). In linux:

```{bash}
# Liftover
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver

# METAL
wget http://www.sph.umich.edu/csg/abecasis/Metal/download/Linux-metal.tar.gz
tar -xvzf Linux-metal.tar.gz

# SHAPEIT
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz
tar -xvzf shapeit.v2.r904.glibcv2.17.linux.tar.gz

# IMPUTE2
wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
tar -xvzf impute_v2.3.2_x86_64_static.tgz

# PLINK (v2.0)
wget http://s3.amazonaws.com/plink2-assets/alpha1/plink2_linux_avx2.zip
unzip plink2_linux_avx2.zip

# EIGENSOFT
wget https://github.com/DReichLab/EIG/archive/v7.2.1.tar.gz
tar -xvzf v7.2.1.tar.gz
cd EIG-7.2.1
cd src/
# may need to edit the Makefile to specify where the correct site for LAPACK and BLAS are ***
# e.g. "LAPACK = -llapack" to "LAPACK = -L/com/extra/lapack/3.5.0/lib -llapack -lblas"
# where -L specifies the path to the LAPACK library.
# LAPACK can be downloaded and installed via:
apt install liblapack-dev liblapack3 libopenblas-base libopenblas-dev
# Once LAPACK has been installed, run the following commands:
make clobber
make install

# EAGLE
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.tar.gz
tar -xvzf Eagle_v2.4.tar.gz
```

OR with miniconda

```{bash}
conda install -c bioconda ucsc-liftover 

# no conda for metal
wget http://www.sph.umich.edu/csg/abecasis/Metal/download/Linux-metal.tar.gz
tar -xvzf Linux-metal.tar.gz

conda install -c bioconda shapeit
conda install -c bioconda impute2
conda install -c bioconda plink2
conda install -c bioconda eigensoft
conda install -c bioconda eagle 
conda install -c soil eagle-phase 
```

- install recopili into the  (replace *<version>* with latest file from the [downloads page](https://sites.google.com/a/broadinstitute.org/ricopili/download)). 

```{bash}
#download and install ricopili
#wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.<version>.tar.gz
wget https://sites.google.com/a/broadinstitute.org/ricopili/download/rp_bin.2018_Jun_6.001.tar.gz
tar -xvzf rp_bin.2018_Jun_6.001.tar.gz

tar -xvzf rp_bin.<version>.tar.gz
#tar -xvzf rp_bin.2018_May_28.001.tar.gz
```

- download [hapmap reference data](http://zzz.bwh.harvard.edu/plink/res.shtml) (linked 20180531) for whichever race / all races you prefer and unzip into the **hapmap** folder.


### rnaelements, rnapancancer, metab, cell, cellseqgenes

adapted from amrits' paper; rna counts in log10 scale, cells are presented as original counts




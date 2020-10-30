

# asthma ER vs DR

[midterm poster](./20180808_poster.pdf), [midterm presentation](./20180723_presentation.pdf), [abstract](201806_prelim-abstract.doc), follow up paper DOI: 10.1183/23120541.00107-2019

**Purpose**: human asthma patients experience an asthma attack or constriction in their airways upon contact with allergens. a subgroup of these patients will recover after their initial attack. however, some will suffer from inflammation some hours later.
- hence this project aims to find biomarkers that differentiates between patents who
  - will (dual response DR) and
  - will not suffer from this second wave of inflammation (early response ER)
- from blood sampled
  - before allergen contact and
  - a few hours after the initial asthma attack (theoretically some time before a possible second wave of inflammation).
  
Note some terms:
- class: phenotype
- flippers: people who have been tested to be ER and DR at some point and therefore have uncertain phenotypes
- features: columns of feature matrices in folder ```feat```
- feature types: what machines each sample was run through e.g. dna = microarray genotype, rnaseq = RNAseq
- scripts / src / code: a .R/.Rmd file in this project
- paper / amrits' paper: see ```src/paper```

Folder structure: src, data (download: https://vault.sfu.ca/index.php/s/x1jV2ky6WKxvIBR)

```{bash}
├─── src (source code)
├──+ data (raw data)
|  ├─+ dna
|  | ├── plink (plink software for imputation)
|  | ├── hapmap (hapmap ref for imputation)
|  | └── reference_SNP
|  ├─+ metab
|  ├─+ other
|  ├─+ rnaelements
|  ├─+ rnapancancer
|  └─+ rnaseq
|    └── fastq (folder containing .fastq files)
├──+ result
|  ├── meta
|  ├── feat
|  ├─+ stat
|  | ├── pca
|  | ├── gwas
|  | └── eqtl
|  └─+ supervised
|    └── blocksplsda
```

Using the scripts:
- all code files are labelled with a number indicating the order in which it should be ran
- set variable ```root``` on all scripts to this directory


## Prepare data: import, normalize, format

Scripts starting with ```00``` imports raw data and normalize/reshape them into feature matrices (feat); each feature has two meta matrices that describe its row (file / sample / subject) and columns (col / feature)

**input**: ```data```

**output**: ```result/meta/file```, ```result/meta/col-<feature type>```, ```result/feat/<feature type>.<time>```


### `meta/file`

This section documents definitions in the `meta/file`.

```result/meta/file``` (subject x subject meta data); columns include:
- ```id``` = row names of feature matrices: indicates subject name and is unique; all duplicate samples are removed, those removed are samples with the least amount of feature types made on it and its response does not conform with a majority of the samples made on that patient
- ```response```: a sample is **ER** if its associated minimum FEV (forced expiratory volume) over times 180-420 min after asthmatic allergen challenge is over -15, results/enrichrwise its a **DR** i.e. min(c("F180L","F240L","F300L","F360L","F420L")) > -15 is ER, <= -15 is DR. if there is no FEV information, response_mac is used (phenotyping done on clinical side)
- ```flipper_calc```: TRUE if a subject (a unique 'id') has multiple samples taken and those samples have inconsistent 'response'; based on calculation above
- ```cohort```: **Discovery** = used to derive panel; **Validation** = used to test panel on
- ```centre```, ```date```: centre (site) and date sample was collected at
- ```race```, ```sex```, ```weight```, ```height```, ```age```, ```bmi```, ...: other useful things
- various columns: indicate what feature types are available for each subject
- ```filename_dna```: microarray well id
- ```filename_rnaseq.<time>```: .bam filename

- ```result/meta/file.raw``` (sample x subject + sample meta data); same as above except columns:
  - ```time```: sample collected **Pre** (before) or **Post** (after) allergen challenge
  - ```id```: indicates subject id; each subject has one or two unique samples (column ```time``` = **Pre** & **Post**)
  
```meta/file``` **indices** used to run analysis on a subset of subjects:
- ```goodppl```: 36 subjects with certain phenotypes used in the discovery cohort of amrits' paper
- ```flipperdr```: the ER's from ```goodppl``` + the **flippers** from the discovery cohort of amrits' paper; comparable with ```goodppl```
- ```all```: everyone!


### `meta/col(s)`

This section documents definitions in the `meta/col(s)`.

```result/meta/col-<feature>``` (feature x feature meta data); feature meta data
- ```id``` = column names of feature matrices
- for **rna** features, important columns include: gene ```symbol```, ```start``` and ```end``` positions, ```chr```omosome
- for **dna** features, important columns include: ```dbSNP``` (snp rs numbers), ```pos_phys``` position, ```chromosome```

```meta/col``` **indices** used to run analysis on a subset of features: these are further split into feature types; for now, only applies to ```dna```
- ```asthma```: asthma related snps or snps close to an asthma related gene derived from (most) databases
- ```asthma-st```: asthma related snps discussed with casey and scott
- ```asthma-rod```: asthma related snps as found by rod
- ```asthma-gwrns```: asthma related snps derived from gwrns database
- ```asthma-ebiclinvaromim```: asthma related snps derived from ebi database
- ```all```: all features!

### feat(ures)

This section documents the features generated for each data set.

- ```result/feat/<feature>.<pre/post/diff>``` (subject x feature): feature matrices with row/colnames coinciding with the 'id' column of meta/file and meta/col; split into whether samples were collected **pre**/**post** allergen test; **diff** is the difference between the two times; 
  - features are split into ```.<pre/post>``` indicating whether the sample was taken before or after the asthma allergen challenge based on availability

below lists the feature types
- ```dna/<01/12>``` (86 subjects x 287624 SNPs): microarray genetics; 01 = dominant / 12 = recessive model; SNPs that has the same genotype across all subjects, or has less than 3 subjectsof an alternate genotype is deleted
- ```rnaseq<genes/isoforms><annotated using ucsc/ens(embl); mapped using star to reference genome / trinity>``` (36 subjects x 9858/18104 (STAR mapping, ensembl annotation), 15209/31688 (STAR mapping, ucsc annotation), /147816 (trinity assembly) genes/isoforms): RNAseq transcriptomics done for the **discovery** cohort
- ```rnaseq<genes/isoforms>``` (36 subjects x 9454/13359/76277 filtered genes/isoforms/isopct): same as ```rnaseq<genes/isoforms>ensstar```, but filtered on TTM normalized data rather than total sum normalization
- ```rnaelements<genes>``` (74 subjects x 75 genes): nanostring elements panel, genes derived from analysis of rnaseq data; (done for the **discovery** cohort (to further filter derived genes) and the **validation**)
- ```rnapc<genes>```  (35 subjects x 600 genes): nanostring pancancer panel, genes related to immune system
- ```metab```olites (35 subjects x 163 metabolites)
- ```cell``` (71 subjects x 8 cell types): white blood cell counts
- ```cellseqgenes``` (32 subjects x 82 cell types): cell counts derived from rnaseq data, see paper for more details (related genes are removed if rnaseq and cellseqgenes are analyzed together)


## Previous pilot analysis

The following folder contains results from a previous pilot project on the same research topic.

**output**: ```result/<script name>.html```, ```result/data```

scripts starting with ```01``` are mostly analyses done by amrit for each data type and on the meta/file from his paper; see html files for more details


## Calculate statistics and plot

**input**: ```result/meta```, ```result/feat```

**output**: ```result/stat```, ```result/supervised```

### PCA

Dimensionality reduction on the features are done using PCA.

```result/stat/pca/<feature type>-<meta/file index>X<meta/col index>_pca-iso```: PCA plots for some demographic features for each subject to see if any confounding factors should be removed, overall subjects are pretty homogenous


### Association studies

```result/stat/gwas/<feature type>-<meta/file index>X<meta/col index>_class-response_DR<# of DRs>vER<# of ER>_test-<chi2/lmbayes type of test used>_<cauc/none -- caucasians only or all races>```: each feature is compared to our class to see if there are any associations (chi2 test done for categorical features, linear monel done for continuous features)
- ```..._id.csv```: subjects used
- ```..._de.png```: differential expression
- ```.../logfc.vavg```: log fold change versus average count
- ```.../qq```: qq plot for unadjusted p values
- ```.../manhattan```: manhattan plots
- ```.../reg```: plots features vs class for the features with the lowest p values
- ```.../dose-gene```: enrichment analysis on known diseases related to the genes / snps
  
the main ```.csv``` files' columns include:
- ```<stat test>_p_<p value adjustment method>```: p values calculated using what test and adjusted using what method
- ```log10de```: log10(differential expression) for continuous features
- ...other meta data


### Quantitative trait loci

Note: only most significant p values are shown

```result/stat/eqtl/<feature type 1>-<feature type 2>-<meta/file index>X<meta/col index>_class-response<could also be sex etc, indicates additional interactions>_DR<# of DRs>vER<# of ER>_<model used>_cisdist-<within how many bases is counted as 'local'>``` (optional: significant features in common between different times ```_<pre/post/diff><# of significant features of the time>v<pre/post/diff><# of significant features of the time>```): just look at the ```modelLINEAR_CROSS```
- ```..._id.csv```: subjects used
- ```.../qq```: qq plot for unadjusted p values

The main ```.csv``` files' columns include:
- ```<stat test>_p_<p value adjustment method>```: p values calculated using what test and adjusted using what method
- ```<feature type 1>```: id for feature type 1
- ```<feature type 2>```: id for feature type 2
- ```pvalue```: unadjusted p value
- ```FDR```: fdr adjusted p value
- ```cis_trans```: whether two features are local (cis) or not (trans); features without genome positions will all be (cis)
- ```gwas.p_<feature type 1>```: association study unadjusted p value for feature type 1
- ```gwas.p_<feature type 2>```: association study unadjusted p value for feature type 2
- ...other meta data; will append ```<feature type 1>``` to indicate which feature type meta data is for
- note: if p value is missing, it wasn't significant > 0.01




### Blocksplsda (diablo): sparse partial least squares differential analysis using regression

```result/supervised/blocksplsda/<feature type 1...n>-<meta/file index>X<meta/col index>_class-response>_DR<# of DRs>vER<# of ER>_<model used>_pthres-<p value threshold for dna, only significant snps are included>_<tune - is number of factors tuned?>_constrains-<from 0 - 1 how constrained should each feature type be with each other>_<bind/none; binded means that every feature of the same type (e.g. rnaseq, rna pancancer) are merged into the same block>``` :
- ```..._id.csv```: subjects used
- ```....Rdata```: result, includes the following (see mixOmics for more details in section block.plsda())
  - ```.../loading_<feat type>.csv```: features used towards model
  - ```.../comp_<feat type>.csv```: component / variate values for each subject
  - ```.../ev_<feat type>.csv```: Percentage of explained variance for each component and each block
  - ```.../weights.csv```: correlation between the components of each feature type and outcome; used to weigh outcome


- ```.../predictscore-<test type: loo=leave one out, Mfold = 10 fold cross validation>```: performance of a feature set in the model
- ```.../tune```: number of factors used vs error rate (tries all number of factors to tune out best number of factors to use); only exists if model was tuned

Other plots: blue = ERs, orange = DRs (see sections in [mixOmics](https://cran.r-project.org/web/packages/mixOmics/mixOmics.pdf) for details)
- ```.../arrows``` plotArrow(): each sample is plotted as an arrow, start/end of arrow is its position in factors of feature type 1/2 (if there's more than 2 feature Types: start/end of arrow is its position as calculated by the median of factors of all feature types / factors of each block); short arrows mean strong agreement between feature types :)
- ```.../circos``` circosPlot(): indicates whether or not features between feature types are correlated
- ```.../heatmap``` cim(): rows = class, columns = feature type; features can be found in the loading csv
- ```.../loadings``` plotLoadings): indicates contribution / coefficient of each feature towards model
- ```.../network``` network(): similarity between features is obtained by calculating the sum of the correlations between the original features and each of the latent components of the model
- ```.../plots``` plotDiablo(): plots correlation between components of different feature types
- ```.../var``` plotVar(): plots the correlations between each feature and components with concentric circles of radius one et radius given by rad.in



# OTHER


```{r}
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





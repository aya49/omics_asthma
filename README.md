# asthma ER vs DR
human asthma patients experience constriction in their airways upon contact with allergen. Amongst them, some also experience a second wave of inflammation after their first reaction.

This project tries to find biomarkers that differentiate these two groups of patients using data from blood sampled
- before allergen contact
- a few hours after initial allergen contact and theoretically some time before a possible second wave of inflammation


Data types used include
- RNAseq
- Genotyping (Affymetrix Axiom)





## Data pre-processing



### RNAseq
input raw fastq; output sample x base/isoform/gene matrix of counts/abundence
- install packages required: fastqc, star, rsem
  - one of the easiest ways to do this is via [python](https://docs.python.org/3/using/unix.html#getting-and-installing-the-latest-version-of-python) and then [conda](https://conda.io/docs/user-guide/install/index.html#installing-conda-on-a-system-that-has-other-python-installations-or-packages) after which you can create an environment to install your packages in. click on the links to see how to do this in your respective operating system. see below for an example in linux.

```bash

#install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

#create an environment to install the packages in
conda env create --name myenv_py3 --file environment.yaml

#install the packages via conda; searchable via https://anaconda.org/
conda install -c bioconda fastqc
conda install -c bioconda star
conda install -c bioconda rsem

#load the environment
source activate myenv_py3

```


- given a folder of raw fastq files, run fastQC (optional if fastQC is already done), STAR (map transcripts to the organism's genome e.g. [GRC](https://www.ncbi.nlm.nih.gov/grc)), and RSEM (get abundence). to do so, first write the raw fastq file paths into [config.yaml](./config.yaml) and then run [Snakefile](./Snakefile) (adapted from the (Sakemake pipeline)[https://snakemake.readthedocs.io/en/stable/]). see below for an example in linux. note: if there is more than one run of a sample (i.e. more than one fastq file per sample), put them together!

```bash
#run snakemake; its configurations and parameters are in config.yaml
snakemake -np --cores 6 #dry run
snakemake -p --cores 6 #actual run; change the number of cores as needed

```














# by Yolanda Yang & Casey P Shannon
# 
# Edit the "config.yaml" file to change the parameter values
#    snakemake -p --cores 8      # -np instead of -p for dry run

configfile: "config.yaml"

# Generate the final outputs: RSEM outputs and fastQC dignostics
rule all:
    input:
        expand("fastQC/{sample}/", sample = config["samples"]),
        expand("rsem/{sample}.genes.results", sample = config["samples"]),
        expand("rsem/{sample}.isoforms.results", sample = config["samples"])
    shell:
        "conda env export > environment.yaml"

rule downloadGTF:
    output:
        expand("{downloadDir}gtf/gencode.v27.annotation.gtf", downloadDir = config["downloadDir"])
    params:
        url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
        downloadDir = config["downloadDir"]
    shell:
        "wget -O - {params.url} | gunzip -c > {output}"

rule downloadFASTA:
    output:
        expand("{downloadDir}fasta/GRCh38.p10.genome.fa", downloadDir = config["downloadDir"])
    params:
        url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz",
        downloadDir = config["downloadDir"]
    shell:
        "wget -O - {params.url} | gunzip -c > {output}"

rule generateRSEMIndex:
    input:
        gtf = expand("{downloadDir}gtf/gencode.v27.annotation.gtf", downloadDir = config["downloadDir"]),
        fasta = expand("{downloadDir}fasta/GRCh38.p10.genome.fa", downloadDir = config["downloadDir"])
    output:
        expand("{downloadDir}index/rsem/", downloadDir = config["downloadDir"])
    shell:
        "rsem-prepare-reference \
        -p 32 \
        --star \
        --gtf {input.gtf} \
        {input.fasta} \
        {output}GRCh38"

rule fastQC:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "fastQC/{sample}/"
    benchmark:
        "benchmarks/fastQC/{sample}.benchmark.txt"
    threads: 32
    params:
        outDir = "fastQC/{sample}/"
    shell:
        "zcat -c {input} | fastqc -t {threads} --outdir {params.outDir} stdin"

rule RSEM:
    input:
        fastqs = lambda wildcards: config["samples"][wildcards.sample],
        ref = expand("{downloadDir}index/rsem/", downloadDir = config["downloadDir"])
    output:
        "rsem/{sample}.genes.results",
        "rsem/{sample}.isoforms.results"
    benchmark:
        "benchmarks/rsem/{sample}.benchmark.txt"
    threads: 32
    params:
        readsType = config["readsType"],
        transcriptBam = config["transcriptBam"],
        genomeBam = config["genomeBam"],
        starGenomeBam = config["starGenomeBam"],
        sortBam = config["sortBam"],
        id = "{sample}"
    shell:
        "rsem-calculate-expression {params.readsType} \
                                   --star \
                                   --star-gzipped-read-file \
                                   -p {threads} \
                                   {params.transcriptBam} {params.genomeBam} {params.starGenomeBam}\
                                   {params.sortBam} \
                                   {input.fastqs} \
                                   {input.ref}GRCh38 \
                                   {params.id} && \
        mv -v *.genes.results *.isoforms.results *.stat/ -t rsem/" #add a *.bam* between results and *.stat if you want the bam file outputs from STAR
        mv -r fastQC ../data/RNAseq/fastQC
        mv -r rsem ../data/RNAseq/rsem
        mv -r benchmarks ../data/RNAseq/benchmarks




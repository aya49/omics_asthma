## input: RSEM output .isoforms/genes.results files
## output: RSEM expected count matrices and meta
## aya43@sfu.ca
## created 20180509
## last modified 20180522

# assumes: this file is in the directory where snakefile etc were ran; genes/transcripts in each file are in same order!
# usage: chmod u+x count1.sh; ./count1.sh

# columns in .results file: gene_id   transcript_id(s)   length   effective_length   expected_count   TPM   FPKM (IsoPct; for isoforms)



paste rsem/*.genes.results > genes.rsem.all
paste rsem/*.isoforms.results > isoforms.rsem.all

awk -v col=5 -v totalcol=7 -f col.awk < genes.rsem.all > genes.rsem.counts
awk -v col=5 -v totalcol=8 -f col.awk < isoforms.rsem.all > isoforms.rsem.counts
awk -v col=8 -v totalcol=8 -f col.awk < isoforms.rsem.all > isoforms.rsem.counts

rm genes.rsem.all
rm isoforms.rsem.all

cut -f 1,2 rsem/HI.1157.003.Index_25.48.genes.results > genes.rsem.genetrans
cut -f 2,1 rsem/HI.1157.003.Index_25.48.isoforms.results > isoforms.rsem.genetrans

cd rsem
printf '%s\t' *.genes.results > ../genes.rsem.files
printf '%s\t' *.isoforms.results > ../isoforms.rsem.files
cd ..

# move everything made to data folder
mkdir rsem_tables
mv *.rsem.* rsem_tables/
mv -tr ../data/RNAseq .snakemake benchmarks fastQC rsem rsem_tables 



# sed -i "1s/.*/$rowname/" genes.rsem.counts


# PASTE joins .results files side-by-side
# CUT choose columns with expected_count info and put them into output

# paste RNASEQ_data/rsem_*/rsem.genes.results | tail -n+2 | \
# cut -f1,5,12,19,26 > RNASEQ_data/edgeR.genes.rsem.txt

# paste RNASEQ_data/rsem_*/rsem.isoforms.results | tail -n+2 | \
# cut -f1,5,13,21,29 > RNASEQ_data/edgeR.isoforms.rsem.txt

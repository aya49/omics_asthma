## define input/output file paths
## aya43@sfu.ca
## created 20180509
## last modified 20180816

## HOW TO USE: in the main script, run the following to set home directory
# root = "~/projects/asthma"
# source(paste0(root, "/code/_dirs.R"))


## root directory ---------------------------------------
setwd(root); cat("\nset directories: ", root)


## input directory ------------------------------------
data_dir = paste0(root,"/data")

# dna
dna_dir = paste0(data_dir, "/dna")
gt1_call_dir = paste0(dna_dir, "/AxiomGT1.calls.txt") # data
meta_snp_temp_dir = paste0(dna_dir, "/Axiom_PMRA.na35.annot.csv") # feature snp annotations
meta_snp_idrod_dir = paste0(dna_dir, "/rod/rod.csv")

# rnaseq
rnaseq_dir = paste0(data_dir,"/rnaseq")
rsem_dir = paste0(rnaseq_dir,"/rsem")
meta_col_temp_dir = paste0(rnaseq_dir,"/HuGene-2_1-st-v1.na36.hg19.transcript.csv")
meta_col_tr_temp_dir = paste0(rnaseq_dir,"/HuGene-2_1-st-v1.na36.hg19.probeset.csv")

trinity_map_dir = paste0(rnaseq_dir, "/asthma.trinity.blastx.outfmt6.txt")
rnaseqa_data_dir = paste0(rnaseq_dir, "/allRnaseqDatasets_rawData.RDATA")
rnaseqa_datanew_dir = paste0(rnaseq_dir, "/allRnaseqDatasets.RDATA")
rnaseqa_datanorm_dir = paste0(rnaseq_dir, "/allRnaseqDatasets_normalized.RDATA")
rnaseqa_trinbl_dir = paste0(rnaseq_dir,"/asthma.trinity.blastx.outfmt6.txt")
rnaseq_ucsc_dir = paste0(rnaseq_dir, "/geneCounts")
rnaseq_iso_dir = paste0(rnaseq_dir, "/geneIsoCounts")
rnaseq_trin_dir = paste0(rnaseq_dir, "/trinity")
rnaseq_control_dir = paste0(rnaseq_dir, "/ERCC_Controls_Analysis.txt")
rnaseq_grch38_dir = paste0(rnaseq_dir, "/ensembl/asthma.GRCh38_ERCC.gencode21.trimmedCounts.txt")

rnapc_preeperfcov_dir = paste0(rnaseq_dir, "/data/prePerfNoCov.rds")
rnapc_panelperfresult_dir = paste0(rnaseq_dir, "/results/rnaseq_biomarkerPanelPerformanceResults.rds")

# rnaelements nanostring rnaelements_data.Rmd
rnaelements_dir = paste0(data_dir,"/rnaelements")
rnaseq_to_ns_dir = paste0(rnaelements_dir, "/results/biomarkerCandidateSelection/rnaseq_to_nanoString_biomarkers.RDATA") # biomarker gene candidates
hba_attenuation_discovery_dirs =
  paste0(rnaelements_dir, "/data/validation/HBA2_attenuation_rediscovery/",
         c("20160114_test0 max fovjan14-2016 _RCC",
           "20160116_set2 jan 16-2016 max fov_RCC",
           "20160116_set 3 3-1 to 3-12 jan16 2016 max fov_RCC",
           "20160116_set 4 4-1 to 4-12 maxfov jan 16-2016_RCC"))
hba_attenuation_discovery_dir =
  paste0(rnaelements_dir, "/data/HBA2Attenutation_asthmaBiomarkersDiscovery.csv")
hba_attenuation_discovery_demo_dir =
  paste0(rnaelements_dir, "/data/validationStudy/reCalibration/attenuation/AttenuationnanoElementsDemo_labSheet.csv")
hba_attenuation_discovery_data_dir =
  paste0(rnaelements_dir, "/data/validation/HBA2_attenuation_rediscovery/HBA2Attenutation_asthmaBiomarkersDiscovery.csv")
hba_attenuation_confirmatory_data_dirs =
  paste0(rnaelements_dir, "/data/validation/HBA2_attenuation_replication/",
         c("20160210_asthma biomarkersvalidation set1 v1-v12 maxFOV Feb-10-2016_RCC",
           "20160210_asthma biomarkersvalidation set2 v13-v24 maxFOV Feb-10-2016_RCC",
           "20160210_asthma biomarkersvalidation set3 v25-v36 maxFOV Feb-10-2016_RCC",
           "20160211_asthma biomarkersvalidation set4 v37-v48 maxFOV Feb-11-2016_RCC",
           "20160211_asthma biomarkersvalidation set5 v49-v60 maxFOV Feb-11-2016_RCC",
           "20160211_asthma biomarkersvalidation set6 v61-v72 maxFOV Feb-11-2016_RCC"))
hba_attenuation_confirmatory_data_dir =
  paste0(rnaelements_dir, "/data/HBA2Attenutation_asthmaBiomarkersValidation_reAnaysis.csv")
hba_attenuation_confirmatory_demo_dir =
  paste0(rnaelements_dir, "/data/validation/HBA2_attenuation_replication/validationCohort_mappngFile.csv")

asthma_templt_dir = paste0(rnaelements_dir, "/data/asthmaTestplate.Jan16_2015.csv")

# rnapancancer
rnapc_dir = paste0(data_dir, "/rnapancancer")
rnapc_data_dir <- paste0(rnapc_dir, "/data/rawdata/nanoString/",
                         c("20150703_asingh - pancancer MAXset A_RCC",
                           "20150709_setB max fov - ykw n amrit_RCC",
                           "20150709_Set3-C-AMRIT-KYW-MAX_RCC"))
rnapc_raw_dir = paste0(rnapc_dir, "/data/preprocessedData/cleanedRawData_panCancerDiscovery.csv")
rnapc_datanew_dir = paste0(rnapc_dir, "/data/preprocessedData/panCancerDatasets.RDATA")
rnapc_panel_dir = paste0(rnapc_dir, "/results/latePhaseBiomarkerPanel.rds")
rnapc_panelperf_dir = paste0(rnapc_dir, "/results/panelPerformance.rds")
rnapc_panelloocv_dir = paste0(rnapc_dir, "/results/panelLOOCV.rds")
rnapc_geo_dir = paste0(rnapc_dir, "/data/preprocessedData/geneExp_geo_plosone.RDATA")
dingo_dir = paste0(rnapc_dir, "/results/dingoResult.rds")

# cell genes
cell_pc_immune_dir = paste0(rnapc_dir, "/Cells_nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.csv")
cell_gene_dir = paste0(data_dir,"/cellgene_dmap.iris.lm22.pc.Rdata")
cell_gene2_dir = paste0(data_dir, "/LM22.txt")

# metab
metab_dir = paste0(data_dir,"/metab")
metab_names_dir = paste0(metab_dir, "/data/metNames.csv")
asthma_discov_dir = paste0(metab_dir, "/data/2015-02-10_Conc_asthmaDiscoveryPlateFeb10.2015.csv")

# demographics meta files
meta_file_temp1_dir = paste0(rnaseq_dir,"/RNASeq.asthma.clinical_sequencing_merged.csv") # cell_data.R
meta_file_temp2_dir = paste0(rnaseq_dir,"/rnaseq_demo.Rdata")
meta_file_extra_dir = paste0(rnaseq_dir, "/asthmaDemo_allsite.xlsx")
gt1_meta_dir = paste0(dna_dir, "/additional_sample_data.txt")
meta_file_temp_dir = paste0(dna_dir, "/meta_file_temp.csv")
meta_file_rnapc_dir = paste0(rnaelements_dir,"/data/allsitesDemo_clean.rds")
meta_file_temp3_dir = paste0(data_dir, "/allsitesDemo.csv")
meta_file_rnaseqa_dir = paste0(data_dir, "/allsitesDemo.rds")
meta_file_data_dir = paste0(data_dir,"/asthmaDemo_allsite_withSampleInfo_DH_v5.csv")

meta_file_rnae_dir = paste0(rnaelements_dir, "/data/demo_eset_nanoString.RDATA")
meta_hiseq_to_id_dir = paste0(rnaseq_dir, "/demo/rnaseqHiSeqInfo.csv")

# other
other_dir = paste0(data_dir, "/other")

other_difbc_dir = paste0(other_dir, "/diffexp/BioCarta_2016_table.txt")
other_difkegg_dir = paste0(other_dir, "/diffexp/KEGG_2016_table.txt")

other_erhga_dir = paste0(other_dir, "/enrichr/biomarkerpanel/cdMolecules_Human_Gene_Atlas_table.txt")
other_erwiki_dir = paste0(other_dir, "/enrichr/biomarkerpanel/innImmuneGenes_WikiPathways_2016_table.txt")
other_erreact_dir = paste0(other_dir, "/enrichr/asthmaExacerbationGenes/Reactome_2016_table.txt")
other_erexach_dir = paste0(other_dir, "/enrichr/asthmaExacerbationGenes/LINCS_L1000_Chem_Pert_down_table.txt")
other_erexali_dir = paste0(other_dir, "/enrichr/asthmaExacerbationGenes/LINCS_L1000_Ligand_Perturbations_up_table.txt")


## output directory ------------------------------------

result_dir = paste0(root, "/result"); dir.create(result_dir, showWarnings=F)
meta_dir = paste0(result_dir,"/meta"); dir.create(meta_dir, showWarnings=F)
feat_dir = paste0(result_dir,"/feat"); dir.create(feat_dir, showWarnings=F)

# meta
meta_file_dir = paste0(meta_dir,"/file")
meta_fileall_dir = paste0(meta_dir,"/file.raw")
meta_filepcdiff_dir = paste0(meta_file_dir,".rnapancancer.setdiff") #one biological sample we don't account for, we use another sample from the same patient L_ST_005 that got other non-rnapancancer transcriptomics done on it as well

grch38_dir = paste0(meta_dir,"/grch38") #full grch genes annotation

meta_col_dir = paste0(meta_dir,"/col-")
meta_col_gt_dir = paste0(meta_col_dir,"dna")
meta_col_rnapc_dir = paste0(meta_col_dir,"rnapcgenes")
meta_col_rnaseq_dir = paste0(meta_col_dir,"rnaseq")
meta_col_rnaelements_dir = paste0(meta_col_dir,"rnaelements")
meta_col_pretty_dir = paste0(meta_dir,"/meta_col_asthma_pretty.csv") # pretty-fied asthma related snps
meta_col_rnacells_dir = paste0(meta_col_dir,"rnaseqstarensgenes_id_cell.Rdata") #genes indicitive of cell populations... derived from amrit's files, should be good

# features
feat_cell_dir = paste0(feat_dir,"/cell")
feat_cellseqgenes_dir = paste0(feat_dir,"/cellseqgenes")
feat_dna_dir = paste0(feat_dir,"/dna")
feat_rnaelements_dir = paste0(feat_dir,"/rnaelements")
feat_rnapc_dir = paste0(feat_dir,"/rnapcgenes")
feat_rnaseq_dir = paste0(feat_dir,"/rnaseq")
feat_metab_dir = paste0(feat_dir,"/metab")

# results
stat_dir = paste0(result_dir,"/stat"); dir.create(stat_dir, showWarnings=F)
other_res_dir = paste0(result_dir,"/other"); dir.create(other_res_dir, showWarnings=F)
colpeek_dir = paste0(other_res_dir,"/col"); dir.create(colpeek_dir, showWarnings=F)
preprocess_dir = paste0(other_res_dir,"/preprocessing"); dir.create(preprocess_dir, showWarnings=F)
pca_dir = paste0(stat_dir,"/pca"); dir.create(pca_dir, showWarnings=F)
eqtl_dir = paste0(stat_dir,"/eqtl"); dir.create(eqtl_dir, showWarnings=F)
gwas_dir = paste0(stat_dir,"/gwas"); dir.create(gwas_dir,showWarnings=F)
# gwasc_dir = paste0(stat_dir,"/gwas_compare"); dir.create(gwasc_dir,showWarnings=F)

super_dir = paste0(result_dir,"/supervised"); dir.create(super_dir, showWarnings=F)
blocksplsda_dir = paste0(super_dir,"/blocksplsda"); dir.create(blocksplsda_dir,showWarnings=F)

# old plots
result_data_dir = paste0(result_dir,"/data"); dir.create(result_data_dir, showWarnings=F)

rnaseq_res_dir = paste0(result_data_dir,"/rnaseq"); dir.create(rnaseq_res_dir, showWarnings=F)
rnaseq_cell_count_dir = paste0(rnaseq_res_dir, "/cell_count.pdf")
rnaseq_demo_dir = paste0(rnaseq_res_dir, "/demographics.pdf")
rnaseq_2b_dir = paste0(rnaseq_res_dir, "/Figure2b.pdf")
rnaseq_2c_dir = paste0(rnaseq_res_dir, "/Figure2c.pdf")
rnaseq_2d_dir = paste0(rnaseq_res_dir, "/Figure2d.pdf")
rnaseq_3a_dir = paste0(rnaseq_res_dir, "/Figure3a.pdf")
rnaseq_3c_dir = paste0(rnaseq_res_dir, "/Figure3c.pdf")
rnaseq_3d_dir = paste0(rnaseq_res_dir, "/Figure3d.pdf")
rnaseq_sig_dir = paste0(rnaseq_res_dir, "/significantgenes_DEanalysis.csv")

rnaelements_res_dir = paste0(result_data_dir,"/rnaelements"); dir.create(rnaelements_res_dir, showWarnings=F)
cell_type_fig_dir = paste0(rnaelements_res_dir, "/cellTypes_")
rnae_5a_dir = paste0(rnaelements_res_dir, "/Figure5a.pdf")
rnae_5b_dir = paste0(rnaelements_res_dir, "/Figure5b.pdf")
rnae_5c_dir = paste0(rnaelements_res_dir, "/Figure5c.pdf")
rnae_5d_dir = paste0(rnaelements_res_dir, "/Figure5d.pdf")
meta_fev_dir = paste0(rnaelements_res_dir, "/Figure4a_fev.pdf")
meta_rev_dir = paste0(rnaelements_res_dir, "/Figure4b_revised.pdf")
meta_5a_dir = paste0(rnaelements_res_dir, "/Figure5a.pdf")
meta_5b_dir = paste0(rnaelements_res_dir, "/Figure5b.pdf")
meta_5c_dir = paste0(rnaelements_res_dir, "/Figure5c.pdf")
meta_5d_dir = paste0(rnaelements_res_dir, "/Figure5d.pdf")

rnapc_res_dir = paste0(result_data_dir,"/rnapancancer"); dir.create(rnapc_res_dir, showWarnings=F)
rnapc_pairs_dir = paste0(rnapc_res_dir, "/FigureS1.pdf")

metab_res_dir = paste0(result_data_dir,"/metab"); dir.create(metab_res_dir, showWarnings=F)
metab_hist_dir = paste0(metab_res_dir, "/hist_metab.pdf")
meta_metab_1a_dir = paste0(metab_res_dir, "/Figure1A.pdf")
meta_metab_1b_dir = paste0(metab_res_dir, "/Figure1B.jpeg")
meta_metab_1c_dir = paste0(metab_res_dir, "/Figure1C.pdf")
meta_metab_1d_dir = paste0(metab_res_dir, "/Figure1D.pdf")
meta_metab_1da_dir = paste0(metab_res_dir, "/Figure1Da.pdf")
meta_metab_1e_dir = paste0(metab_res_dir, "/Figure1E.pdf")

## input: a whole bunch of raw demographics / met afiles
## output: well formatted meta files (subject x subject meta data)
## aya43@sfu.ca
## created 20180509
## last modified 20180509


## logistics
root = "~/projects/asthma" # root directory, used for _dirs.R
source(paste0(root, "/src/_dirs.R"))
source(paste0(root, "/src/_func.R"))
source(paste0(root, "/src/_func-asthma.R"))
libr(append(pkgs(),c("pracma")))


## options
writecsv = T # write results as csv on top of Rdata?

responsetype_col = "_calc" # which response column to use? "_gen" for genetics reponse, "_mac" for mcmaster response label (most complete)?


## start
start = Sys.time()


## meta 1: dna ----------------------------------

# load dna filenames
gt1_calls_colnames = names(read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1, nrows=1))
wells = gsub(".CEL","",sapply(str_split(gt1_calls_colnames,"_"), function(x) x[length(x)]))

# load and bind dna affymetrix & lab meta files (same rows)
gt1_meta = fread(gt1_meta_dir, data.table=F)
gt1lab_meta = fread(meta_file_temp_dir, data.table=F)
meta_file0 =
  cbind(gt1lab_meta,
        gt1_meta[match(gt1lab_meta[,"Sample ID"],gt1_meta[,"Name"]), !colnames(gt1_meta)%in%c("Name")])

# display column names and how many unique elements in each, of meta_file; delete cols with only one unique element
ucol = col_probe(meta_file0)
meta_file0 = meta_file0[,-ucol$u1]

# remove control samples
meta_file0 = meta_file0[meta_file0[, "Sample Type"]=="Normal",]

# adjust values
meta_file0[meta_file0[,"Kit"]=="Mini kit","Batch"] = 0
meta_file0 = gsub(" ","",as.matrix(meta_file0))
# # get rid of batch 0; there's a copy of it in results/enrichr batches
# row_ind = meta_file0[,"Batch"]!="0" & !grepl("Normal",meta_file0[,"Sample ID"],ignore.case=T) # make patient name the unique id
# meta_file = meta_file0[row_ind,]


# delete columns
meta_file0 =
  meta_file0[, !colnames(meta_file0) %in%
               c("Sample Type", "dishQC",
                 "Cluster CR (Axiom inlier cluster call rate)",
                 "Array", "Replicate Info", "Number", "Random order",
                 "Inferred Gender","Random order", "Number",
                 "1000 Genome Concordance", "Reproducibility", "Replicate Info",
                 "Cluster CR (Axiom inlier cluster call rate)", "dishQC") &
               !grepl("[/]|Plate|Het|Rate|Concordance|Reproducibility|Filename",colnames(meta_file0))]

# rename columns (temporarily, to match latter meta files -- Gen_response is class!)
colnames(meta_file0) =
  c("filename_dna","SEX","SITE","NAME",
    "Gen_Response","kit_dna", "batch_dna","RACE")
# meta_file0 = meta_file0[,-c("kit")] #batch = 0 is a minikit; else it's a DNeasy

# order rows according to dna matrix
meta_file0 = as.data.frame(meta_file0)
roworder = match(wells,meta_file0[,paste0("filename_dna")])
meta_file01 = meta_file0[roworder[!is.na(roworder)],]

# check for flippers
response_table = rowSums(table(meta_file01$NAME, meta_file01$Gen_Response) != 0)
flippers_gen <- names(response_table[response_table >1])
flippers_non_gen <- names(response_table[response_table ==1])


# isolate one sample only for each subject
meta_file010 = meta_file01
u_id = sapply(unique(meta_file01$NAME), function(x) {
  return(ifelse(sum(meta_file01$NAME==x)==1, which(meta_file01$NAME==x), which(meta_file01$NAME==x & meta_file01$kit_dna=="Minikit")))
})
meta_file01 = meta_file01[u_id,]


## meta 2: general ------------------------------
# meta_file02a = read.csv(meta_file_data_dir, check.names=T); colnames(meta_file02a)[8] = "NAME"
# meta_file02a = meta_file02a[!grepl("Normal",meta_file02a$NAME),]
meta_file02a = read.xls(meta_file_extra_dir, check.names=T)
meta_file02a = meta_file02a[!grepl("Normal",meta_file02a$NAME),]

# adjust id
meta_file02a$NAME[grepl("WRF",meta_file02a$NAME)] = "WRF"

#check if any rows missing for the pancancer data set, if so, merge
meta_file02_id = paste(meta_file02a$NAME, meta_file02a$AIC_YMD, meta_file02a$Time)


## meta 3: all (pan-cancer too) -------------------------
meta_file03 = readRDS(meta_file_rnapc_dir)
meta_file03 = meta_file03[!grepl("Normal",meta_file03$NAME),] #delete controls
meta_file03$NAME[grepl("WRF",meta_file03$NAME)] = "WRF"
meta_file03_id = paste(meta_file03$NAME, meta_file03$AIC_YMD, meta_file03$Time)

flippers_pc = unique(meta_file03$NAME)[sapply(unique(meta_file03$NAME), function(x) any(meta_file03$ranseq_flipper[meta_file03$NAME==x]=="Y"))]
flippers_non_pc = meta_file03$NAME[!meta_file03$NAME%in%flippers_pc]

#check if there are subjects in pan-cancer meta not in general meta
pc_id = paste0(meta_file03$NAME, "_", meta_file03$AIC_YMD)
extra_id = paste0(meta_file02a$NAME, "_", meta_file02a$AIC_YMD)
meta_file02 = meta_file02a
if (sum(!meta_file03_id%in%meta_file02_id)>0)
  meta_file02 = rbind(meta_file02a,meta_file03[
    !meta_file03_id%in%meta_file02_id,
    match(tolower(colnames(meta_file02a)), tolower(colnames(meta_file03)))])




## meta 4: rnaseq ---------------------------------
meta_file040 = get(load(meta_file_temp2_dir)) # contains rnaseq ids

meta_file040 = data.frame(lapply(meta_file040, as.character), stringsAsFactors=FALSE)
meta_file040[grepl("WRF",meta_file040[,"NAME"]),"NAME"] = "WRF"


## meta 4.5: cell ---------------------------------
# meta_file1_temp0 = fread(meta_file_temp1_dir, data.table=F)
# meta_file1_temp0[grepl("WRF",meta_file1_temp0[,"Name"]),"Name"] = "WRF"
# meta_file1_temp0[,"Filename"] = gsub(".bam","",meta_file1_temp0[,"Filename"])
# meta_file04.5 = meta_file1_temp0
# meta_file04.5$Sample.ID[is.na(meta_file04.5$Sample.ID)] = meta_file04.5$Subject[is.na(meta_file04.5$Sample.ID)]

load(rnaseqa_data_dir)
meta_file04 = demoRNASeq
meta_file04$NAME = as.character(meta_file04$NAME)
meta_file04[grepl("WRF",meta_file04[,"NAME"]),"NAME"] = "WRF"

# all(meta_file040$Filename.Prefix%in%meta_file04.5$Filename) #T!
# all(meta_file04.5$Filename%in%meta_file040$Filename.Prefix) #F!
# meta_file040 = meta_file040[match(meta_file04.5$Filename, meta_file040$Filename.Prefix),]

# # combine 04, 04.5
# intercol = intersect(colnames(meta_file04.5), colnames(meta_file040))
# meta_file04.5[,intercol] = sapply(intercol, function(x) {
#   vals = meta_file04.5[,x]
#   misval = is.na(vals) | vals==""
#   if(sum(misval)>0) vals[misval] = meta_file040[misval,x]
#   vals
# })

# meta_file04 = cbind(meta_file040[,!colnames(meta_file040)%in%append(colnames(meta_file04.5),c("Filename.Prefix"))],meta_file04.5)
# meta_file04$MST_LST_numbers[is.na(meta_file04$MST_LST_numbers)] = meta_file04$Sample.ID[is.na(meta_file04$MST_LST_numbers)]
# meta_file04$MST_LST_numbers[grepl("L_ST_025",meta_file04$MST_LST_numbers)] = "L_ST_025"

#check for missing responses, L_ST_025 has no response
meta_file04$rp04 = meta_file04$Response_usedInRNAseqBiomarkerAnalysis
meta_file04$rp04[is.na(meta_file04$rp04)] = meta_file04$Mac_Response[is.na(meta_file04$rp04)]
meta_file04$rp04[is.na(meta_file04$rp04)] = meta_file04$CorrectResponse[is.na(meta_file04$rp04)]
meta_file04$rp04[is.na(meta_file04$rp04)] = meta_file04$calculated_Response[is.na(meta_file04$rp04)]



## adjust values ------------------------------------

pc_ind = paste0(meta_file03$NAME, "_", meta_file03$AIC_YMD)
extraa_ind = paste0(meta_file02a$NAME, "_", meta_file02a$AIC_YMD)
if (sum(!pc_ind %in% meta_file02)>0)
meta_file02 = rbind(meta_file02a,
                        meta_file03[!pc_ind %in% meta_file02,
                          match(tolower(colnames(meta_file02a)), tolower(colnames(meta_file03)))])

meta_file02$UniqueID[grepl("L_ST_025",meta_file02$UniqueID)] = "L_ST_025"
meta_file02$UniqueID = gsub(" ","_",meta_file02$UniqueID)


#bmi
meta_file02$bmi = as.numeric(meta_file02[,"Wt..Kg."])/(meta_file02[,"HT.cm."]^2)

# id_namedate = subject names _ date
meta_file02$id_namedate = paste("uniqueID",
                                    as.numeric(factor(as.character(meta_file02$NAME))),
                                    as.numeric(factor(as.character(meta_file02$AIC_YMD))), sep="_")

# id_name = subject names
meta_file02$id_name = factor(paste("ID", as.numeric(factor(as.character(meta_file02$NAME))), sep="_"))

# scale breathing test at diff times columns based on blfev
dropFEV1 = as.data.frame(t(100*scale(t(
  meta_file02[, c("BLFEV","F10L","F20L","F30L","F45L","F60L","F90L","F120L",
                      "F180L","F240L","F300L","F360L","F420L")]),
  center=meta_file02$BLFEV, scale=meta_file02$BLFEV)))

# calculate response
ear = suppressWarnings( apply(dropFEV1, 1, function(x) {min(x[1:8], na.rm=T)}) )
ear[ear == "Inf"] = NA
lar = suppressWarnings( apply(dropFEV1, 1, function(x) {min(x[9:13], na.rm=T)}) )
lar[lar == "Inf"] = NA
meta_file02$EAR = round(ear, 1)
meta_file02$LAR = round(lar, 1)

time = gsub("L", "", substring(colnames(dropFEV1), 2))
# time = vapply(str_split(colnames(dropFEV1), "_"),'[',3, FUN.VALUE=character(1)) %>% as.numeric
time[1] = 0;
time = as.numeric(time)/60

# AUR trapz: area of a function for each patient, all and LAS time vs drops
meta_file02$AUC = apply(dropFEV1, 1, function(y) trapz(x=time, y=y) )
meta_file02$AUC_LAR = apply(dropFEV1, 1, function(y) trapz(x=time[9:13], y=y[9:13]) )

#AIS: PC20 of pre/post
ais = meta_file02[meta_file02$Time == "Pre", "PC20"] / meta_file02[meta_file02$Time == "Post", "PC20"]
names(ais) = meta_file02$id_namedate[meta_file02$Time == "Pre"]
meta_file02$AIS = ais

# calculate a reponse
meta_file02$response_calc = rep(NA, nrow(meta_file02))
meta_file02$response_calc[meta_file02$LAR > -15] = "ER"
meta_file02$response_calc[meta_file02$LAR <= -15] = "DR"
meta_file02$Mac_Response[!as.character(meta_file02$Mac_Response) %in% c("ER","DR")] = NA
meta_file02$CorrectResponse[!as.character(meta_file02$CorrectResponse) %in% c("ER","DR")] = NA

# check for flippers
response_table = rowSums(table(meta_file02$NAME, meta_file02$response_calc) != 0)
flippers_calc <- names(response_table[response_table >1])
flippers_non_calc <- names(response_table[response_table ==1])

response_table = rowSums(table(meta_file02$NAME, meta_file02$Mac_Response) != 0)
flippers_mac <- names(response_table[response_table >1])
flippers_non_mac <- names(response_table[response_table ==1])

# delete people without a response
meta_file02 = meta_file02[apply(meta_file02[,grep("esponse", colnames(meta_file02))], 1, function(x) !all(is.na(x)|x=="")),]

response_col = paste0("response",responsetype_col)

meta_file02$response = meta_file02[,response_col]
meta_file02$response[is.na(meta_file02$response)] = meta_file02$Mac_Response[is.na(meta_file02$response)]

meta_file021 = meta_file02

# reconcile duplicate values
meta_file021$SITE[meta_file021$SITE=="MAC"] = "McMaster"
meta_file021$SITE = tolower(meta_file021$SITE)
# meta_file021$centre = tolower(meta_file021$SITE)
# meta_file021$SITE[is.na(meta_file021$SITE)] =
#   meta_file021$centre[is.na(meta_file021$SITE)]


# meta_file021$race[meta_file021$race==""] = "Unknown"
#
# meta_file021$sex[is.na(meta_file021$sex)] = meta_file021$SEX[is.na(meta_file021$sex)]
# meta_file021$sex[meta_file021$sex=="Female"] = "F"
# meta_file021$sex[meta_file021$sex=="Male"] = "M"
# meta_file021$sex[meta_file021$sex==""] = NA


## merge meta dna and meta general with democlean ---------------------------------
meta_file022 = meta_file021
for (meta_file_tmp in list(meta_file01,meta_file03,meta_file04)) {
  print(sum(!meta_file_tmp$NAME%in%meta_file022$NAME)) #0
  
  if (sum(duplicated(meta_file_tmp$NAME))==0 & !is.null(meta_file_tmp$NAME)) { #meta_file01
    meta_file_tmpt = meta_file_tmp[match(meta_file022$NAME,meta_file_tmp$NAME),]
    
  } else if (sum(rownames(meta_file_tmp)%in%paste0(meta_file022$id_namedate, ".", meta_file022$Time))>0) { #meta_file03
    meta_file_tmpt = meta_file_tmp[match(paste0(meta_file022$id_namedate, ".", meta_file022$Time),rownames(meta_file_tmp)),]
    
  } else if (!is.null(meta_file_tmp$Filename.Prefix)) { #meta_file04; 2 samples not in 022; no worrites they don't have responses!
    # meta_file_tmpt = meta_file_tmp[match(meta_file022$rnaseqID_quebec,meta_file_tmp$rnaseqID_quebec),] # SHOULD BE THIS ONE, but if i do this... i need to put back in samples with less of results/enrichr data types done on, so, ill stick to this for now :)
    tmpid = paste0(meta_file_tmp$NAME,"_",meta_file_tmp$Time,"_",meta_file_tmp$AIC_YMD)
    m022id = paste0(meta_file022$NAME,"_",meta_file022$Time,"_",meta_file022$AIC_YMD)
    meta_file_tmpt = meta_file_tmp[match(m022id,tmpid),]
  }
  meta_file022 = cbind(meta_file022, meta_file_tmpt[,!colnames(meta_file_tmpt)%in%colnames(meta_file022)])
  # meta_file022 = merge(meta_file_tmpt[,!colnames(meta_file_tmpt)%in%colnames(meta_file022) | colnames(meta_file01)%in%"NAME"], meta_file022, by="NAME", all=T, suffixes=c("",""))
  for (x in intersect(colnames(meta_file_tmpt), colnames(meta_file022))) {
    narow = meta_file022[,x]=="" | is.na(meta_file022[,x])
    meta_file022[narow,x] = meta_file_tmpt[narow,x]
  }
}
# # add rows not merged from meta_file04
# meta_file022 = rbind(meta_file022, 
#                      t(t(as.matrix(meta_file04))[match(colnames(meta_file022),colnames(meta_file04)),
#                                       !meta_file04$MST_LST_numbers%in%meta_file022$MST_LST_numbers]))

# colnames(meta_file03)[colnames(meta_file03)%in%"calculated_Response"] = "response_calc"
# meta_file03$ID = meta_file03$NAME
#
# match_col = match(tolower(colnames(meta_file022)), tolower(colnames(meta_file03)))
#
# meta_file022 = rbind(meta_file022,meta_file03[!meta_file03$ID%in%meta_file022$id,match_col])
# colnames(meta_file022)[colnames(meta_file022)%in%"response_calc"] = "response"

## adjust values again ------------------------------
meta_file022$RACE[grepl("Unknown",meta_file022$RACE) | meta_file022$RACE=="" | is.na(meta_file022$RACE)] =
  meta_file022$RACE[grepl("Unknown",meta_file022$RACE) | meta_file022$RACE=="" | is.na(meta_file022$RACE)]
meta_file022$RACE[meta_file022$RACE=="White"] = "Caucasian"
meta_file022$RACE = tolower(meta_file022$RACE)
namenorace = meta_file022$NAME[meta_file022$RACE]

meta_file022$SEX[meta_file022$SEX=="Female"] = "F"
meta_file022$SEX[meta_file022$SEX=="Male"] = "M"
meta_file022$SEX[!meta_file022$SEX%in%c("F","M")] = NA

meta_file022 = meta_file022[!is.na(meta_file022[,response_col]),]

meta_file022_temp = meta_file022

# isolate one sample only for each subject's pre/post; based on response_calc
# backup_col = "response_mac"
u_ind = sapply(unique(meta_file022$NAME), function(x) {
  n_ind = which(meta_file022$NAME==x)
  if (max(table(meta_file022[n_ind,"Time"]))<2) return(n_ind)
  y = meta_file022[n_ind,"response"]
  n_ind1 = n_ind[y==mode_set(y)]
  yl = lapply(unique(meta_file022[n_ind,"AIC_YMD"]), #dates
              function(z) meta_file022[meta_file022$NAME==x & meta_file022$AIC_YMD==z,])
  names(yl) = unique(meta_file022[n_ind,"AIC_YMD"])
  n_date = names(yl)[which.max(sapply(yl, function(yli) ifelse(is.na(yli$Cohort), sum(yli!="" & !is.na(yli)), Inf) ))]
  n_ind[meta_file022[n_ind,"AIC_YMD"]==n_date]
})
meta_file022 = meta_file022[as.vector(u_ind),]
meta_file022 = meta_file022[,apply(meta_file022,2,function(x) !all(is.na(x) | x==""))]

## extract cells feature -------------------------------------------
meta_file1_temp0 = meta_file022
cells = meta_file022[,c(29:34,129:142,149:177)]
crn = meta_file022$NAME
# cells = cells[,order(colnames(cells))]
tentimes = as.numeric(gsub("[.]","e",gsub("^x|[.]$","",str_extract(colnames(cells),"x10[.e0-9]+"))) )

# use counts, not percentages
cells$NEU.x10e9.L.[is.na(cells$NEU.x10e9.L.)] = 
  (cells$NEU.WBC. * cells$WBC.x10e9.L.) [is.na(cells$NEU.x10e9.L.)]
cells$LYM.x10e9.L.[is.na(cells$LYM.x10e9.L.)] = 
  (cells$LYM.WBC. * cells$WBC.x10e9.L.) [is.na(cells$LYM.x10e9.L.)]
cells$Mono.x10e9.L.[is.na(cells$Mono.x10e9.L.)] = 
  (cells$Mono.WBC. * cells$WBC.x10e9.L.) [is.na(cells$Mono.x10e9.L.)]
cells$Eos.x10e9.L.[is.na(cells$Eos.x10e9.L.)] = 
  (cells$Eos.WBC. * cells$WBC.x10e9.L.) [is.na(cells$Eos.x10e9.L.)]
cells$Baso.x10e9.L.[is.na(cells$Baso.x10e9.L.)] = 
  (cells$Baso.WBC. * cells$WBC.x10e9.L.) [is.na(cells$Baso.x10e9.L.)]
cells$Baso.x10e9.L.[is.na(cells$Baso.x10e9.L.)] = 
  (cells$Baso.WBC. * cells$WBC.x10e9.L.) [is.na(cells$Baso.x10e9.L.)]
ccn = colnames(cells)
cells = sapply (1:ncol(cells), function(x) { #doesn't word for column "TotalCC..x10.6."
  if (is.na(tentimes[x])) return(as.numeric(cells[,x]))
  return(as.numeric(cells[,x]) * as.numeric(tentimes[x]))
})
colnames(cells) = ccn

# split data into pre/post allergic attack samples for each subject
cells_pre = cells[meta_file1_temp0$Time=="Pre",]
rownames(cells_pre) = crn[meta_file1_temp0$Time=="Pre"]
cells_pre = cells_pre[apply(cells_pre,1,function(x)!all(is.na(x))),
                      apply(cells_pre,2,function(x)!all(is.na(x)))]
cells_pre = cells_pre[,grepl(".L.",colnames(cells_pre))]

cells_post = cells[meta_file1_temp0$Time=="Post",]
rownames(cells_post) = crn[meta_file1_temp0$Time=="Post"]
cells_post = cells_post[apply(cells_post,1,function(x)!all(is.na(x))),
                        apply(cells_post,2,function(x)!all(is.na(x)))]
cells_post = cells_post[,grepl(".L.",colnames(cells_post))]

# save
save(cells, file=paste0(feat_cell_dir,".raw.Rdata"))
if (writecsv) write.csv(cells, file=paste0(feat_cell_dir,".raw.csv"))
save(cells_pre, file=paste0(feat_cell_dir,".pre.Rdata"))
if (writecsv) write.csv(cells_pre, file=paste0(feat_cell_dir,".pre.csv"))
save(cells_post, file=paste0(feat_cell_dir,".post.Rdata"))
if (writecsv) write.csv(cells_post, file=paste0(feat_cell_dir,".post.csv"))



## delete columns from merged meta_file & final response annotation ---------------------------------------
# (save raw version with both response=pre/post, anresults/enrichr one with just pre)
meta_file2 = meta_file022[,c(1:14,16:ncol(meta_file022))] # data availability
colnames(meta_file2)[c(1:3,5:6,8:12,15,35:36,40,42,46,72,107)] =
  c("sponsor","centre","drug","id","date","race","sex",
    "weight","height","age","allergen","time","response_mac",
    "cohort","id_mstlst","flipper_rnaseq","response_dna","filename_rnaseq")
# Read.Set.Id = filename_rnaseq
meta_file2 = meta_file2[,apply(meta_file2,2,function(x) any(!is.na(x)))]
# convert data to right format: char, numeric, factor
for (ci in 1:ncol(meta_file2)) {
  x = meta_file2[,ci]
  if (sum(!is.na(as.numeric(as.character(x[!is.na(x)]))))/sum(!is.na(x))>.7 & sum(!is.na(x))>0) {
    meta_file2[,ci] = as.numeric(x)
  } else { meta_file2[,ci] = as.character(x) }
}


# annotate response based on majority; does everything twice... but not computationally expensive so alright
meta_file2$response_dna = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file010$NAME, NA, mode_set(meta_file010$Gen_Response)) )
meta_file2$response_mac = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file02$NAME, NA, mode_set(meta_file02$Mac_Response)) )
meta_file2$response_calc = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file02$NAME, NA, mode_set(meta_file02$response_calc)) )

# add flipper info
meta_file2$flipper_calc = meta_file2$flipper_gen = meta_file2$flipper_mac = meta_file2$flipper_rnaseq = NA
meta_file2$flipper_calc[meta_file2$id%in%flippers_calc] =
  meta_file2$flipper_gen[meta_file2$id%in%flippers_gen] =
  meta_file2$flipper_rnaseq[meta_file2$id%in%flippers_pc] =
  meta_file2$flipper_mac[meta_file2$id%in%flippers_mac] = TRUE
meta_file2$flipper_calc[meta_file2$id%in%flippers_non_calc] =
  meta_file2$flipper_gen[meta_file2$id%in%flippers_non_gen] =
  meta_file2$flipper_rnaseq[meta_file2$id%in%flippers_non_pc] =
  meta_file2$flipper_mac[meta_file2$id%in%flippers_non_mac] = FALSE


# isolate meta table with id as subject name
meta_file = meta_file2[meta_file2$time=="Pre",]
colnames(meta_file)[colnames(meta_file)%in%"filename_rnaseq"] = "filename_rnaseq.pre"
meta_file$filename_rnaseq.post = meta_file2$filename_rnaseq[grep("Post",meta_file2$time)[
  match(meta_file$id,meta_file2$id[meta_file2$time=="Post"])]]


## final response annotation (take majority) ------------------------------
# not doing anymore, just use calculated response + Mac response


## save --------------------------------------------
save(meta_file2, file=paste0(meta_fileall_dir,".Rdata"))
if (writecsv) write.csv(meta_file2, file=paste0(meta_file_dir,".raw.csv"))
save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))


# isolate non-flippers: let's just use the ones from the previous paper! there are 6 flippers by calculation here though, and there are at least two flippers based on all the columns
demo = get(load(rnaseqa_datanorm_dir))
goodppl = unique(as.character(demo$NAME))
goodppl[grepl("WRF",goodppl)] = "WRF"
# goodppl = meta_file$id[!meta_file[,paste0("flipper",responsetype_col)]] # based on calculated data
save(goodppl, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))

flipperdr = meta_file$id[meta_file$id%in%goodppl & meta_file$response=="ER"]
flipperdr = append(flipperdr, meta_file$id[!meta_file$id%in%goodppl & meta_file$cohort%in%"Discovery" & meta_file$response=="DR"])
# flipperdr = meta_file$id[meta_file[,flipper_col] | meta_file$response=="ER"] # based on calculated data
save(flipperdr, file=paste0(meta_file_dir,"_id_flipperdr.Rdata"))


## make fev plots ----------------------

fltime_colind = colnames(meta_file)[grepl("F[0-9]+L",colnames(meta_file))]
# fltimes = as.numeric(gsub("[A-Z]","",colnames(meta_file2)[fltime_colind]))
# fldf = Reduce('rbind', lapply(1:length(fltime_colind), function(i)
#   data.frame(id=meta_file2$id, response=meta_file2$response,
#              time=rep(fltimes[i],nrow(meta_file2)),
#              fev=meta_file2[,fltime_colind[i]]) )) #the manual way of melt()
fldf = melt(meta_file[append(c("id","response"),fltime_colind)], measure.vars=fltime_colind)
fldf$variable = as.numeric(gsub("[A-Z]","",fldf$variable))


#### this part doesn't work, it was for testing new pretty plots, don't worry
fldfa = ddply(fldf, .(response, variable), function(x)
  c(mean=mean(x$value), sd = sd(x$value),
    lower=as.numeric(quantile(x$value, .25)), upper=as.numeric(quantile(x$value, .75))) )

pd <- position_dodge(width = 0.2) # move them .2 to the left and right
gbase  = ggplot(fldfa, aes(y=mean, colour=response)) +
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=pd) +
  # geom_point(position=pd) +
  geom_line() # + facet_grid(variable ~ sex)
# gline = gline %+% fldfa

ggsave(filename=paste0(stat_dir,"/fev.png"),
       gbase + aes(x=variable) +
         xlab("time (min after allergen exposure)") +
         ylab("mean FEV (forced expiratory volume)"),
       scale = 1.5, width = 5, height = 3, units = c("in"))

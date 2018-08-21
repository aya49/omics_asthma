## input: create meta file for object / subject
## output: well formatted meta files
## aya43@sfu.ca
## created 20180509
## last modified 20180509


## logistics
root = "~/projects/asthma"; commandArgs <- function(...) root  # root directory, used for _dirs.R
source(paste0(root, "/code/_dirs.R"))
source(paste0(root, "/code/_func.R"))
source(paste0(root, "/code/_func-asthma.R"))
libr(pkgs())


## options
writecsv = T # write results as csv on top of Rdata?

responsetype_col = "_calc" # which response column to use? "_gen" for genetics reponse, "_mac" for mcmaster response label (most complete)?


## start
start = Sys.time()


## meta 1: genotype ----------------------------------

# load genotype filenames
gt1_calls_colnames = names(read.table(gt1_call_dir, sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1, nrows=1))
wells = gsub(".CEL","",sapply(str_split(gt1_calls_colnames,"_"), function(x) x[length(x)]))

# load and bind genotype affymetrix & lab meta files (same rows)
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
# # get rid of batch 0; there's a copy of it in other batches
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
  c("filename_genotype","SEX","SITE","NAME",
    "Gen_Response","kit_genotype", "batch_genotype","RACE")
# meta_file0 = meta_file0[,-c("kit")] #batch = 0 is a minikit; else it's a DNeasy

# order rows according to genotype matrix
meta_file0 = as.data.frame(meta_file0)
roworder = match(wells,meta_file0[,paste0("filename_genotype")])
meta_file01 = meta_file0[roworder[!is.na(roworder)],]

# check for flippers
response_table = rowSums(table(meta_file01$NAME, meta_file01$Gen_Response) != 0)
flippers_gen <- names(response_table[response_table >1])
flippers_non_gen <- names(response_table[response_table ==1])


# isolate one sample only for each subject
meta_file010 = meta_file01
u_id = sapply(unique(meta_file01$NAME), function(x) {
  return(ifelse(sum(meta_file01$NAME==x)==1, which(meta_file01$NAME==x), which(meta_file01$NAME==x & meta_file01$kit_genotype=="Minikit")))
})
meta_file01 = meta_file01[u_id,]


## meta 2: general ------------------------------
meta_file02a = read.xls(meta_file_extra_dir, check.names=T)
meta_file02a = meta_file02a[!grepl("Normal",meta_file02a$NAME),]

# adjust id
meta_file02a$NAME[grepl("WRF",meta_file02a$NAME)] = "WRF"

#check if any rows missing for the pancancer data set, if so, merge
meta_file02_id = paste(meta_file02a$NAME, meta_file02a$AIC_YMD, meta_file02a$Time)


## meta 3: pan-cancer -------------------------
meta_file03 = readRDS(meta_file_rnapc_dir)
meta_file03 = meta_file03[!grepl("Normal",meta_file03$NAME),] #delete controls
meta_file03$NAME[grepl("WRF",meta_file03$NAME)] = "WRF"
meta_file03_id = paste(meta_file03$NAME, meta_file03$AIC_YMD, meta_file03$Time)

#check if there are subjects in pan-cancer meta not in general meta
pc_id = paste0(meta_file03$NAME, "_", meta_file03$AIC_YMD)
extra_id = paste0(meta_file02a$NAME, "_", meta_file02a$AIC_YMD)
meta_file02 = meta_file02a
if (sum(!meta_file03_id%in%meta_file02_id)>0)
  meta_file02 = rbind(meta_file02a,meta_file03[
    !meta_file03_id%in%meta_file02_id,
    match(tolower(colnames(meta_file02a)), tolower(colnames(meta_file03)))])




## meta 4: rnaseq ---------------------------------
meta_file04 = get(load(meta_file_temp2_dir)) # contains rnaseq ids

meta_file04 = data.frame(lapply(meta_file04, as.character), stringsAsFactors=FALSE)
meta_file04[grepl("WRF",meta_file04[,"NAME"]),"NAME"] = "WRF"




## adjust values ------------------------------------
pc_ind = paste0(meta_file03$NAME, "_", meta_file03$AIC_YMD)
extraa_ind = paste0(meta_file02a$NAME, "_", meta_file02a$AIC_YMD)
if (sum(!pc_ind %in% meta_file02)>0)
meta_file02 = rbind(meta_file02a,
                        meta_file03[!pc_ind %in% meta_file02,
                          match(tolower(colnames(meta_file02a)), tolower(colnames(meta_file03)))])

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


# isolate one sample only for each subject's pre/post; based on response_calc
backup_col = "response_mac"
meta_file02_id = paste0(meta_file02$NAME,"_",meta_file02$Time)
u_ind = sapply(unique(meta_file02_id), function(x) {
  n_ind = which(meta_file02_id==x)
  y = meta_file02[n_ind,"response"]
  n_ind1 = n_ind[y==mode_set(y)]
  n_ind1[which.max(apply(meta_file02[n_ind1,],1, function(y) sum(y!="" & !is.na(y))))]
})
meta_file021 = meta_file02[u_ind,]


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


## merge meta genotype and meta general with democlean ---------------------------------
meta_file022 = meta_file021
for (meta_file_tmp in list(meta_file01,meta_file03,meta_file04)) {
  print(sum(!meta_file_tmp$NAME%in%meta_file022$NAME)) #0
  if ("AIC_YMD"%in%colnames(meta_file_tmp)) {
    meta_file_tmpt = meta_file_tmp[match(paste0(meta_file022$NAME,meta_file022$AIC_YMD),
                        paste0(meta_file_tmp$NAME,meta_file_tmp$AIC_YMD)),]
  } else {
    meta_file_tmpt = meta_file_tmp[match(meta_file022$NAME,meta_file_tmp$NAME),]
  }
  meta_file022 = cbind(meta_file022, meta_file_tmpt[,!colnames(meta_file_tmpt)%in%colnames(meta_file022)])
  # meta_file022 = merge(meta_file_tmpt[,!colnames(meta_file_tmpt)%in%colnames(meta_file022) | colnames(meta_file01)%in%"NAME"], meta_file022, by="NAME", all=T, suffixes=c("",""))
  for (x in intersect(colnames(meta_file_tmpt), colnames(meta_file022))) {
    narow = meta_file022[,x]=="" | is.na(meta_file022[,x])
    meta_file022[narow,x] = meta_file_tmpt[narow,x]
  }
}

# colnames(meta_file03)[colnames(meta_file03)%in%"calculated_Response"] = "response_calc"
# meta_file03$ID = meta_file03$NAME
#
# match_col = match(tolower(colnames(meta_file022)), tolower(colnames(meta_file03)))
#
# meta_file022 = rbind(meta_file022,meta_file03[!meta_file03$ID%in%meta_file022$id,match_col])
# colnames(meta_file022)[colnames(meta_file022)%in%"response_calc"] = "response"


meta_file022$RACE[grepl("Unknown",meta_file022$RACE) | meta_file022$RACE=="" | is.na(meta_file022$RACE)] =
  meta_file022$RACE[grepl("Unknown",meta_file022$RACE) | meta_file022$RACE=="" | is.na(meta_file022$RACE)]
meta_file022$RACE[meta_file022$RACE=="White"] = "Caucasian"
meta_file022$RACE = tolower(meta_file022$RACE)
namenorace = meta_file022$NAME[meta_file022$RACE]

meta_file022$SEX[meta_file022$SEX=="Female"] = "F"
meta_file022$SEX[meta_file022$SEX=="Male"] = "M"
meta_file022$SEX[!meta_file022$SEX%in%c("F","M")] = NA

meta_file022 = meta_file022[!is.na(meta_file022[,response_col]),]


## delete columns from merged meta_file & final response annotation ---------------------------------------
# (save raw version with both response=pre/post, another one with just pre)
meta_file2 = meta_file022[,c(1:3,5:14,16:116)] # data availability
colnames(meta_file2)[c(1:5,7:11,14,34:35,39,41,45,71)] =
  c("sponsor","centre","drug","id","date","race","sex",
    "weight","height","age","allergen","time","response_mac",
    "cohort","id_mstlst","flipper_rnaseq","response_gen")
# Read.Set.Id = filename_rnaseq
meta_file2 = meta_file2[,apply(meta_file2,2,function(x) any(!is.na(x)))]
for (ci in 1:ncol(meta_file2)) {
  x = meta_file2[,ci]
  if( sum(!is.na(as.numeric(x[!is.na(x)])))/sum(!is.na(x))>.7 & sum(!is.na(x))>0)
    meta_file2[,ci] = as.numeric(x)
}


# annotate response based on majority; does everything twice... but not computationally expensive so alright
meta_file2$response_gen = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file010$NAME, NA, mode_set(meta_file010$Gen_Response)) )
meta_file2$response_mac = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file02$NAME, NA, mode_set(meta_file02$Mac_Response)) )
meta_file2$response_calc = sapply(meta_file2$id, function(x)
  ifelse (!x%in%meta_file02$NAME, NA, mode_set(meta_file02$response_calc)) )

# add flipper info
meta_file2$flipper_calc = meta_file2$flipper_gen = meta_file2$flipper_mac = NA
meta_file2$flipper_calc[meta_file2$id%in%flippers_calc] =
  meta_file2$flipper_gen[meta_file2$id%in%flippers_gen] =
  meta_file2$flipper_mac[meta_file2$id%in%flippers_mac] = T
meta_file2$flipper_calc[meta_file2$id%in%flippers_non_calc] =
  meta_file2$flipper_gen[meta_file2$id%in%flippers_non_gen] =
  meta_file2$flipper_mac[meta_file2$id%in%flippers_non_mac] = F


# isolate meta table with id as subject name
meta_file = meta_file2[meta_file2$time=="Pre",]
colnames(meta_file)[colnames(meta_file)%in%"filename_rnaseq"] = "filename_rnaseq.pre"
meta_file$filename_rnaseq.post = meta_file2$filename_rnaseq[grep("Post",meta_file2$time)[
  match(meta_file$id,meta_file2$id[meta_file2$time=="Post"])]]

## final response annotation (take majority) ------------------------------


## save --------------------------------------------
save(meta_file2, file=paste0(meta_fileall_dir,".Rdata"))
if (writecsv) write.csv(meta_file2, file=paste0(meta_file_dir,".raw.csv"))
save(meta_file, file=paste0(meta_file_dir,".Rdata"))
if (writecsv) write.csv(meta_file, file=paste0(meta_file_dir,".csv"))

goodppl = meta_file$id[!meta_file[,paste0("flipper",responsetype_col)]]
save(goodppl, file=paste0(meta_file_dir,"_id_goodppl.Rdata"))



## make fev plots ----------------------

fltime_colind = colnames(meta_file)[grepl("F[0-9]+L",colnames(meta_file))]
# fltimes = as.numeric(gsub("[A-Z]","",colnames(meta_file2)[fltime_colind]))
# fldf = Reduce('rbind', lapply(1:length(fltime_colind), function(i)
#   data.frame(id=meta_file2$id, response=meta_file2$response,
#              time=rep(fltimes[i],nrow(meta_file2)),
#              fev=meta_file2[,fltime_colind[i]]) )) #the manual way of melt()
fldf = melt(meta_file[append(c("id","response"),fltime_colind)], measure.vars=fltime_colind)
fldf$variable = as.numeric(gsub("[A-Z]","",fldf$variable))

#### need to debug...
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







me = Matrix_eQTL_main(
  snps = f1_sd,
  gene = f2_sd,
  cvrt = cvrt, # SlicedData object with additional covariates. Can be an empty SlicedData object in case of no covariates. The constant is always included in the model and would cause an error if included in cvrt. The order of columns must match those in snps and gene.
  
  output_file_name = eqtl_tra_dir, # significant associations (all significant associations if pvOutputThreshold=0 or only distant if pvOutputThreshold>0). If the file with this name exists, it is overwritten.
  output_file_name.cis = eqtl_cis_dir, #output_file_name_cis=tempfile(); significant local associations
  
  pvOutputThreshold = pvalthres_tra, 
  pvOutputThreshold.cis = pvalthres_cis,
  
  useModel = useModel, 
  # errorCovariance = errorCovariance, # numeric. The error covariance matrix. Use numeric() for homoskedastic independent errors.
  verbose = F, 
  
  snpspos = f1_pos,
  genepos = f2_pos, 
  cisDist = cisDist,
  pvalue.hist = T, # logical, numerical, or "qqplot" (faster if false); To record information for a histogram set pvalue.hist to the desired number of bins of equal size. Finally, pvalue.hist can also be set to a custom set of bin edges.
  min.pv.by.genesnp = T, # record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis (faster if false)
  noFDRsaveMemory = F # save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTLs are not recorded
) 

rt_analysis = function(ctrl_dir, sample_dir, annotation_path, ctrl, sample, GOI, libType, gen_info, batch, norm='1', viewSamples=c('ARS2', 'CPSF73', 'INTS11', 'RRP40'), ctrl_sample='EGFP', statistic='median', usExt=5000, plot=TRUE, verbose=TRUE, verbose2=TRUE, pdf=FALSE){
  ##########################################################
  #### (1) read in gene specific data
  ########################################################## 
  
  # Output the index of the gene of interest (GOI) and its name if verbose or verbose2 is TRUE
  if (verbose | verbose2){
    cat(which(rownames(gen_info)==GOI), GOI, sep='\t')
  }
  
  # Output a tab or newline character based on the value of verbose
  cat(ifelse(verbose, '\t', '\n'))
  
  # Create an empty list to store data for each replicate
  rt_list = list()
  
  ###################### KAOS-KODE #########################
  
  strand_options = c("+"="plus", "-"="minus")
  
  
  # Annotation importing
  annot_gr = import(annotation_path) # update filename
  
  
  # Subsetting annotation on gene name
  gene_annot = subset(annot_gr, gene_name == GOI)
  chrom_no = as.character(seqnames(gene_annot)@values)
  strand_sign = as.character(strand(gene_annot)@values)
  
  # Coordinates
  start_coord = min(start(gene_annot))
  end_coord = max(end(gene_annot))
  
  # Empty DF with correct length
  indexes <- seq(start_coord, end_coord)
  GOI_data <- data.frame(index = indexes)
  
  # Loop for loading replicates and saving to DF
  for (fname in ctrl) {
    ctrl_fname = paste0(ctrl_dir, fname, "_", strand_options[strand_sign], ".bw")
    ctrl_bw = import(ctrl_fname, 
                     which = GenomicRanges::GRanges(seqnames = chrom_no, 
                                                    ranges = IRanges::IRanges(start = start_coord, end = end_coord)), 
                     as = "NumericList")[[1]]
    GOI_data[[fname]] = ctrl_bw
  }
  
  # Loop for loading replicates and saving to DF
  for (fname in sample) {
    sample_fname = paste0(sample_dir, fname, "_", strand_options[strand_sign], ".bw")
    sample_bw = import(sample_fname, 
                       which = GenomicRanges::GRanges(seqnames = chrom_no, 
                                                      ranges = IRanges::IRanges(start = start_coord, end = end_coord)), 
                       as = "NumericList")[[1]]
    GOI_data[[fname]] = sample_bw
  }
  
  # Setting the DF indices to the positions of the reads and deleting index column
  rownames(GOI_data) = GOI_data$index
  GOI_data$index = NULL
  
  # Find coord of next gene
  next_gene_annot = subset(annot_gr, seqname = chrom_no, stand = strand_sign)
  if (strand_sign == "+") {
    ss = sort(start(next_gene_annot))
    next_gene_coord = which(ss > end_coord)[1]
  } else {
    ss = sort(start(next_gene_annot))
    next_gene_coord = which(ss < end_coord)[1]
  }
  
  ###################### KAOS-KODE #########################
  
  
  
  ##########################################################  
  #### (2) log2-transform data (+ pseudocount of 1)
  ##########################################################  
  
  # Perform a log2 transformation on the gene of interest data, adding a pseudocount of 1 to avoid taking the log of zero
  log2_GOI_data = log2(GOI_data + 1)
  
  # Remove the original gene of interest data from memory to save space
  rm(GOI_data)
  
  ###############################################################
  #### (3) remove batch effects in log2-transformed L_sample data
  ###############################################################
  
  # Remove batch effects from the log2-transformed gene of interest data using the removeBatchEffect function from the limma package
  log2_GOI_data.nobatch = limma::removeBatchEffect(log2_GOI_data, batch=batch)
  
  # Remove the original log2-transformed gene of interest data from memory to save space
  rm(log2_GOI_data)
  
  ###############################################################
  #### (4) Normalize to gene body signal
  ###############################################################
  
  # Determine the upstream transcript start site (uTSS) and downstream transcript end site (dTES) for the GOI
  uTSS = ifelse(strand=='+', gen_info[GOI, 'start'], gen_info[GOI, 'end'])
  uTSS_row = which(rownames(log2_GOI_data.nobatch)==as.character(uTSS))
  dTES = ifelse(strand=='+', gen_info[GOI, 'end'], gen_info[GOI, 'start'])
  dTES_row = which(rownames(log2_GOI_data.nobatch)==as.character(dTES))
  
  # Extract transcript start sites (TSSs) and determine the main TSS
  TSSs = extract_TSSs(gen_info, GOI) #@ extract_TSSs is a home-made function, see separate annotation
  if (strand == '+'){
    TSS = ifelse(length(TSSs[TSSs >= uTSS & TSSs < dTES]) > 0, TSSs[TSSs >= uTSS & TSSs < dTES][1], uTSS)
  }else{
    TSS = ifelse(length(TSSs[TSSs <= uTSS & TSSs > dTES]) > 0, TSSs[TSSs <= uTSS & TSSs > dTES][1], uTSS)
  }
  TSS_row = which(rownames(log2_GOI_data.nobatch)==as.character(TSS))
  
  # Extract transcript end sites (TESs) and determine the main TES
  TESs = extract_TESs(gen_info, GOI) #@ extract_TESs is a home-made function, see separate annotation
  TESs = unique(TESs)
  if (strand == '+'){
    TESs = TESs[TESs > TSS & TESs <= dTES]
  }else{
    TESs = TESs[TESs < TSS & TESs >= dTES]
  }
  if (length(TESs) == 0){
    TESs = dTES
  }
  TES_rows = unlist(sapply(TESs, function(x) {which(rownames(log2_GOI_data.nobatch)==as.character(x))}))
  
  # Calculate standard deviations for multiple TESs and select the one with the minimum standard deviation
  if (length(TESs) > 1){
    stdevs = rep(NA, length(TESs))
    for (i in 1:length(TES_rows)){
      TES_row = TES_rows[i]
      if (TES_row>TSS_row){
        body.norm.factors = apply(log2_GOI_data.nobatch[TSS_row:TES_row,], 2, median)
        log2_GOI_data.nobatch.bodynorm = t(t(log2_GOI_data.nobatch) - body.norm.factors)
        stdevs[i] = sd(log2_GOI_data.nobatch.bodynorm[TSS_row:TES_row, ])
      }
    }
    TES = TESs[which(stdevs == min(stdevs, na.rm=TRUE))]
  }else{
    TES = TESs
  }
  TES_row = which(rownames(log2_GOI_data.nobatch)==as.character(TES))
  
  # Calculate body normalization factors and adjust data accordingly
  body.norm.factors = apply(log2_GOI_data.nobatch[TSS_row:TES_row,], 2, median)
  for (curr.sample in paste0(libType, '_', viewSamples)){
    ctrl_rep1 = paste0(libType, '_', ctrl_sample, '_', 'rep1')
    sample_rep1 = paste0(curr.sample, '_', 'rep1')
    sample_rep1_diff = as.numeric(body.norm.factors[sample_rep1] - body.norm.factors[ctrl_rep1])
    ctrl_rep2 = paste0(libType, '_', ctrl_sample, '_', 'rep2')
    sample_rep2 = paste0(curr.sample, '_', 'rep2')
    sample_rep2_diff = as.numeric(body.norm.factors[sample_rep2] - body.norm.factors[ctrl_rep2])
    body_diff = mean(c(sample_rep1_diff, sample_rep2_diff))
    rt_list[[curr.sample]] = list('TSS'=TSS_row-usExt, 'TES'=TES_row-usExt, 'rep1_diff'=sample_rep1_diff, 'rep2_diff'=sample_rep2_diff, 'body_diff'=body_diff)
  }
  
  # Provide information about the normalization process if specified
  if (verbose){
    cat(paste(norm_types[[norm]], 'normalization'), '\n')
  }
  
  # Adjust normalization factors based on the specified method
  if (norm=='0'){
    body.norm.factors = rep(0, length(body.norm.factors))
  }
  if (norm=='2'){   
    upTES_row = ifelse(TES_row-TSS_row+1 >= 500, TES_row-500+1, TSS_row) ## up to 500 bp upstream of TES
    body.norm.factors = apply(log2_GOI_data.nobatch[upTES_row:TES_row,], 2, median)
  }
  
  # Adjust data based on the calculated body normalization factors
  log2_GOI_data.nobatch.bodynorm = t(t(log2_GOI_data.nobatch) - body.norm.factors)
  
  # Calculate mean coverage for each sample
  log2_GOI_data.nobatch.bodynorm.means = matrix(0, nrow=nrow(log2_GOI_data.nobatch.bodynorm), ncol=5)
  colnames(log2_GOI_data.nobatch.bodynorm.means) = paste0(libType, '_', c('EGFP', 'ARS2', 'CPSF73', 'INTS11', 'RRP40'))
  rownames(log2_GOI_data.nobatch.bodynorm.means) = rownames(log2_GOI_data.nobatch)
  for (sample in paste0(libType, '_', c('EGFP', 'ARS2', 'CPSF73', 'INTS11', 'RRP40'))){
    rep_1 = paste0(sample, '_rep1')
    rep_2 = paste0(sample, '_rep2')
    mean.cov = rowMeans(log2_GOI_data.nobatch.bodynorm[, c(rep_1, rep_2)])
    log2_GOI_data.nobatch.bodynorm.means[, sample] = mean.cov
  }
  
  # Remove the original log2-transformed gene of interest data from memory to save space
  rm(log2_GOI_data.nobatch)
  
  ###############################################################
  #### (5) Check upstream signal
  ###############################################################
  
  # Determine if the major transcript start site (TSS) is the upstream transcription start site (uTSS) or within the larger gene
  if (TSS_row > uTSS_row){
    within.gene = TRUE
    us.TSS_rows = uTSS_row:(TSS_row-1)
  }else{
    within.gene = FALSE
  }
  
  # Create a vector of positions representing upstream regions
  us_rows = 1:usExt
  
  # Calculate the mean signal within the gene body and upstream regions for each sample
  for (curr.sample in paste0(libType, '_', c(ctrl_sample, viewSamples))){
    body_mean  = mean(log2_GOI_data.nobatch.bodynorm.means[TSS_row:TES_row, curr.sample])
    if (within.gene){
      us.TSS_mean = mean(log2_GOI_data.nobatch.bodynorm.means[us.TSS_rows, curr.sample])
      us.TSS_diff = body_mean - us.TSS_mean
    }else{
      us.TSS_diff = NA
    }
    us_mean = mean(log2_GOI_data.nobatch.bodynorm.means[us_rows, curr.sample])
    us_diff = body_mean - us_mean
    
    # Store the differences between gene body and upstream signal for each sample in the result list
    if (curr.sample != paste0(libType, '_', ctrl_sample)){
      rt_list[[curr.sample]][['us_diff']] = us_diff           
      rt_list[[curr.sample]][['us.TSS_diff']] = us.TSS_diff   
    }
  }
  
  ###############################################################
  #### (6) Subtract control signal
  ###############################################################
  
  # Define the names and columns of the control sample and its replicates
  curr.ctrl = paste0(libType, '_', ctrl_sample)
  curr.ctrl_col = which(colnames(log2_GOI_data.nobatch.bodynorm.means) == curr.ctrl)
  curr.ctrl_rep1 = paste0(curr.ctrl, '_rep1')
  curr.ctrl_rep1_col = which(colnames(log2_GOI_data.nobatch.bodynorm) == curr.ctrl_rep1)
  curr.ctrl_rep2 = paste0(curr.ctrl, '_rep2')
  curr.ctrl_rep2_col = which(colnames(log2_GOI_data.nobatch.bodynorm) == curr.ctrl_rep2)
  
  # Subtract control signal from the gene body and remove control columns
  log2_GOI_data.nobatch.bodynorm.means.ctrlsubtract = log2_GOI_data.nobatch.bodynorm.means - log2_GOI_data.nobatch.bodynorm.means[, curr.ctrl_col]
  log2_GOI_data.nobatch.bodynorm.means.ctrlsubtract = log2_GOI_data.nobatch.bodynorm.means.ctrlsubtract[, -curr.ctrl_col]
  
  # Subtract control signal from the replicate columns
  log2_GOI_data.nobatch.bodynorm.ctrlsubtract = log2_GOI_data.nobatch.bodynorm
  log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, c(1, 3, 5, 7, 9)] = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, c(1, 3, 5, 7, 9)] - log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, curr.ctrl_rep1_col]
  log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, c(2, 4, 6, 8, 10)] = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, c(2, 4, 6, 8, 10)] - log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, curr.ctrl_rep2_col]
  log2_GOI_data.nobatch.bodynorm.ctrlsubtract = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, -c(curr.ctrl_rep1_col, curr.ctrl_rep2_col)]
  
  # Note: All dataframes/matrices up until this point contain usExt (max 5000) positions upstream of locus of interest, all positions from locus of interest, and from 'allowed' downstream region.
  # uTSS starts at row usExt + 1 (default 5001)

  ###############################################################
  #### (7) Find states using HMM and fit double-sigmoidal to data
  ###############################################################
  # This code segment iterates through samples, detects putative readthrough events using Hidden Markov Models (HMM), and fits sigmoidal or double sigmoidal functions to the readthrough data to estimate its length and intensity. It also handles cases where no readthrough is detected.
  
  for (curr.sample in paste0(libType, '_', viewSamples)){
    rt_list[[curr.sample]][['rt']] = FALSE              ###@@@ 2) was readthrough detected by one or other method (HMM, sigmoidal or doublesigmoidal fitting)
    rt_list[[curr.sample]][['max_rt_length']] = 0      ###@@@ 3) length of the 'allowed' readthrough region (i.e. data supported region)
    rt_list[[curr.sample]][['rt_TES']] = 0            ###@@@ 4a) TES from which rt starts (real coordinates) determined by HMM
    rt_list[[curr.sample]][['rt_start']] = 0           ###@@@ 4b) readthrough start (relative to uTSS) determined by HMM - also used for sigmoidal or doublesigmoidal fitting
    rt_list[[curr.sample]][['rt_end']] = 0             ###@@@ 5) readthrough end (relative to uTSS) determined by HMM
    rt_list[[curr.sample]][['rt_int']] = 0             ###@@@ 6) readthrough intensity (mean or median signal) in the HMM determined readthrough
    rt_list[[curr.sample]][['rt_sum']] = 0             ###@@@ 7) readthrough intensity integrated in the HMM determined readthrough
    rt_list[[curr.sample]][['rt_max']] = NA             ###@@@ 8) readthrough intensity (mean or median signal) in the HMM determined state with highest signal
    rt_list[[curr.sample]][['rt_max_iv']] = NA             ###@@@ 9) interval with highest readthrough intensity
    rt_list[[curr.sample]][['dsfit']] = FALSE           ###@@@ 10) double sigmoidal fit (TRUE/FALSE)
    rt_list[[curr.sample]][['sfit']] = FALSE            ###@@@ 11) sigmoidal fit (TRUE/FALSE)
    rt_list[[curr.sample]][['rt_end_fitted']] = 0      ###@@@ 12) readthrough end (relative to uTSS) determined by sigmoidal or doublesigmoidal fitting
    rt_list[[curr.sample]][['extrapolated']] = FALSE    ###@@@ 13) is rt_end_fitted determined by extrapolation beyond the data region (TRUE/FALSE)
    rt_list[[curr.sample]][['best_fit']] = NA           ###@@@ 14) data for best fit (used for plotting at the end of function - will be removed from output)
    rt_list[[curr.sample]][['best_fit_asymp']] = NA     ###@@@ 15) asymptotic value for the fitted double sigmoidal curve (used for calculation of rt_end_fitted) (only relevant for doublesigmoidal fitting)
    rt_list[[curr.sample]][['best_fit_Rsq']] = NA       ###@@@ 16) R-squared value for the best fitted curve (only relevant for sigmoidal or doublesigmoidal fitting)
    rt_list[[curr.sample]][['endDeclinePoint_x']] = NA  ###@@@ 17) 'end' of fitted doublesigmoidal curve (used to cut down the stored data if possible) (only relevant for doublesigmoidal fitting)
    rt_list[[curr.sample]][['rt_int_fitted']] = 0      ###@@@ 18) readthrough intensity (mean or median signal) in the fitted determined readthrough
    rt_list[[curr.sample]][['rt_sum_fitted']] = 0      ###@@@ 19) readthrough intensity integrated in the fitted determined readthrough
    
    samples_reps = c(paste0(curr.sample, '_rep1'), paste0(curr.sample, '_rep2'))
    curr.sample_data = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[, samples_reps]  ##@@ same nrow as above
    hmm_result = sample_subtract_hmm(curr.data.list = list('GOI' = curr.sample_data[TSS_row:nrow(curr.sample_data), ]), TSS_row, TES_row, TES_rows, col = cols_trans[[curr.sample]], usExt = usExt, statistic = statistic, cutoff = 2, plot = FALSE, ymax = 0.5, max_cv = 200)
    
    ## sample_subtract_hmm analyses a dataframe starting at the chosen TSS (TSS_row) otherwise containing the rest of the curr.sample_data dataframe: 
    ## rownumbers (rownum) in this dataframe relates to the original dataframe as TSS_row = 1 -> orig_rownum = rownum + TSS_row - 1
    ## the output readthrough interval ('rt_iv' ]TES;end_rt]) is given in rownumbers relating to original dataframes: conversion relative to uTSS: pos = orig_rownum - usExt
    
    if (length(hmm_result) > 0){
      rt_start.row = hmm_result[['rt_iv']][1] + 1
      rt_end.row = hmm_result[['rt_iv']][2]
      reg_end.row = nrow(curr.sample_data)
      if (reg_end.row > rt_start.row){
        if (verbose){
          cat(paste(paste0(curr.sample, ':'), 'putative readthrough detected'), '\t')
        }
        rt_region_means = rowMeans(curr.sample_data[rt_start.row:reg_end.row, ])
        
        ## this dataframe starts at the chosen rt_start (rt_start.row) otherwise containing the rest of the curr.sample_data dataframe: 
        ## rownumbers (rownum) in this dataframe relates to the original dataframe as rt_start.row =  1 -> orig_rownum = rownum + rt_start.row - 1
        ## the output readthrough interval ('rt_iv' ]TES;end_rt]) is given in rownumbers relating to original dataframes: conversion relative to uTSS: pos = orig_rownum - usExt
        
        max_rt_length = length(rt_region_means)
        rt_list[[curr.sample]][['rt']] = TRUE
        rt_list[[curr.sample]][['max_rt_length']] = max_rt_length
        rt_list[[curr.sample]][['rt_TES']] = rownames(curr.sample_data)[hmm_result[['rt_iv']][1]]   ### real coordinates
        rt_list[[curr.sample]][['rt_start']] = rt_start.row - usExt  ### uTSS relative coordinates
        rt_list[[curr.sample]][['rt_end']] = rt_end.row - usExt     ### uTSS relative coordinates
        rt_list[[curr.sample]][['rt_int']] = hmm_result[['rt_int']]
        rt_list[[curr.sample]][['rt_max']] = hmm_result[['max_int']]
        rt_list[[curr.sample]][['rt_max_iv']] = hmm_result[['max_int_iv']] - usExt
        rt_list[[curr.sample]][['rt_sum']] = sum(rt_region_means[1:(rt_end.row - rt_start.row + 1)])   ## orig_rownum = rt_length + rt_start.row - 1 -> rt_length = orig_rownum - rt_start.row + 1
        message = paste(paste0('estimated rt-length based on HMM fitting', ':'), rt_end.row - rt_start.row + 1)
        ## fit a (double) sigmoidal function to the data
        binning = FALSE
        if (length(rt_region_means) > 10000){  ## convert data to 100 bp bins
          bin.size = 100
          no.bins = as.integer(length(rt_region_means) / bin.size)
          binning = TRUE
        } else {
          if (length(rt_region_means) > 100){  ## convert data to 100 bins !!
            no.bins = 100
            bin.size = as.integer(length(rt_region_means) / no.bins)
            binning = TRUE
          }
        }
        if (binning){
          bin.means = rep(0, no.bins)
          bin.xvals = rep(0, no.bins)
          for (bin.no in 1:no.bins){
            iv = ((bin.no - 1) * bin.size + 1):(bin.no * bin.size)
            bin.xvals[bin.no] = rt_start.row + as.integer(mean(iv)) - 1
            bin.means[bin.no] = ifelse(statistic == 'median', median(rt_region_means[iv]), mean(rt_region_means[iv]))
          }
        } else {
          bin.xvals = rt_start.row:reg_end.row
          bin.means = rt_region_means
        }
        tss = sum((rt_region_means - mean(rt_region_means))^2)  ## total sum of squares used for calculating Rsq
        sig = tryCatch({
          suppressWarnings(fitAndCategorize(data.frame('time' = bin.xvals, 'intensity' = bin.means)))
        }, error = function(e) {
          return(list('summaryVector' = list('decision' = 'error')))
        })
        if (sig$summaryVector$decision == 'double_sigmoidal' | sig$summaryVector$decision == 'ambiguous'){
          parameterVector = sig$doubleSigmoidalModel
          xvals_fit = rt_start.row + 1:parameterVector$endDeclinePoint_x - 1
          fit_int_standard = doublesigmoidalFitFormula(xvals_fit, finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate, maximum = parameterVector$maximum_Estimate,
                                                       slope1Param = parameterVector$slope1Param_Estimate,
                                                       midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                       slope2Param = parameterVector$slope2Param_Estimate,
                                                       midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate)
          xvals_Rsq_length = min(length(rt_region_means), length(xvals_fit))
          Rsq_for_fits = rep(NA, 8)
          for (i in 1:8){
            fit_int = adj_function(i, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
            rss = sum((rt_region_means[1:xvals_Rsq_length] - fit_int[1:xvals_Rsq_length])^2)  ## residual sum of squared errors
            Rsq_for_fits[i] = 1 - rss / tss  ## R squared
          }
          best_fit_index = which(Rsq_for_fits == max(Rsq_for_fits, na.rm = TRUE))[1]
          best_fit_Rsq = Rsq_for_fits[best_fit_index]
          best_fit = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
          best_fit_asymp = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = TRUE, fit_asymptote_standard = parameterVector$finalAsymptoteIntensity)
          max_x = which(best_fit == max(best_fit))
          if (length(max_x) > 1){
            max_x = as.integer(median(max_x))
          }
          best_fit_zero = best_fit - best_fit_asymp
          integrated_best_fit_zero = sum(best_fit_zero[max_x:length(best_fit_zero)])
          rel_cumsum_best_fit_zero = cumsum(best_fit_zero[max_x:length(best_fit_zero)]) / integrated_best_fit_zero
          intersect = which(rel_cumsum_best_fit_zero >= 0.95)[1]
          if (!is.na(intersect)){
            new_row_intersect = which(rel_cumsum_best_fit_zero >= 0.95)[1] + max_x - 1
          } else {
            new_row_intersect = max_rt_length
          }
          rt_end_fitted = new_row_intersect + rt_start.row - usExt - 1   ### uTSS relative coordinates ## orig_rownum = rownum + rt_start.row - 1 
          rt_list[[curr.sample]][['dsfit']] = TRUE
          rt_list[[curr.sample]][['sfit']] = FALSE
          rt_list[[curr.sample]][['rt_end_fitted']] = rt_end_fitted
          rt_list[[curr.sample]][['extrapolated']] = ifelse(new_row_intersect + rt_start.row - 1 > reg_end.row, TRUE, FALSE)
          rt_list[[curr.sample]][['best_fit']] = best_fit
          rt_list[[curr.sample]][['best_fit_asymp']] = best_fit_asymp
          rt_list[[curr.sample]][['best_fit_Rsq']] = best_fit_Rsq
          rt_list[[curr.sample]][['endDeclinePoint_x']] = as.integer(round(parameterVector$endDeclinePoint_x + rt_list[[curr.sample]][['rt_start']]))
          rt_fitted_end_for_sum = ifelse(new_row_intersect <= max_rt_length, new_row_intersect, max_rt_length)
          rt_list[[curr.sample]][['rt_int_fitted']] = ifelse(statistic == 'median', median(rt_region_means[1:rt_fitted_end_for_sum]), mean(rt_region_means[1:rt_fitted_end_for_sum]))       
          rt_list[[curr.sample]][['rt_sum_fitted']] = sum(rt_region_means[1:rt_fitted_end_for_sum])
          message = paste(paste0('estimated rt-length based on double sigmoidal fitting', ':'), new_row_intersect)
        }
        if (sig$summaryVector$decision == 'sigmoidal' | sig$summaryVector$decision == 'ambiguous'){
          sigmoidal = TRUE
          parameterVector = sig$sigmoidalModel
          xvals_fit = rt_start.row:reg_end.row
          fit_int_standard = sigmoidalFitFormula(xvals_fit, maximum = parameterVector$maximum_Estimate, slopeParam = parameterVector$slopeParam_Estimate, midPoint = parameterVector$midPoint_Estimate)
          Rsq_for_fits = rep(NA, 8)
          for (i in 1:8){
            fit_int = adj_function(i, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
            rss = sum((rt_region_means - fit_int)^2)  ## residual sum of squared errors
            Rsq_for_fits[i] = 1 - rss / tss  ## R squared
          }
          best_fit_index = which(Rsq_for_fits == max(Rsq_for_fits, na.rm = TRUE))[1]
          best_fit_Rsq = Rsq_for_fits[best_fit_index]
          best_fit = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
          if (!is.na(rt_list[[curr.sample]][['best_fit_Rsq']])){
            if (best_fit_Rsq >= rt_list[[curr.sample]][['best_fit_Rsq']]){
              sigmoidal = TRUE
            } else {
              sigmoidal = FALSE
            }
          }
          if (sigmoidal){  
            rt_list[[curr.sample]][['dsfit']] = FALSE
            rt_list[[curr.sample]][['sfit']] = TRUE
            rt_list[[curr.sample]][['rt_end_fitted']] = reg_end.row - usExt  ### uTSS relative coordinates 
            rt_list[[curr.sample]][['extrapolated']] = FALSE
            rt_list[[curr.sample]][['best_fit']] = best_fit
            rt_list[[curr.sample]][['best_fit_asymp']] = NA
            rt_list[[curr.sample]][['best_fit_Rsq']] = best_fit_Rsq
            rt_list[[curr.sample]][['endDeclinePoint_x']] = reg_end.row - usExt  ### uTSS relative coordinates 
            rt_list[[curr.sample]][['rt_int_fitted']] = ifelse(statistic == 'median', median(rt_region_means[1:max_rt_length]), mean(rt_region_means[1:max_rt_length]))       
            rt_list[[curr.sample]][['rt_sum_fitted']] = sum(rt_region_means[1:max_rt_length])
            message = paste0('estimated rt-length based on sigmoidal fitting', ': >', max_rt_length)
          }
        }
        if (verbose){
          cat(message, '\n')
        }
      } else {
        if (verbose){
          cat(paste(paste0(curr.sample, ':'), 'no readthrough detected'), '\n')
        }
      }
    } else {
      if (verbose){
        cat(paste(paste0(curr.sample, ':'), 'no readthrough detected'), '\n')
      }
    }
  }  
  
  
  ###############################################################
  #### (8) preparing for plotting
  ###############################################################
  # This code segment prepares the data for plotting by binning it and organizing it into appropriate data structures, such as vectors and lists, for easier handling and visualization.

  xvals = 1:nrow(log2_GOI_data.nobatch.bodynorm.means)-usExt
  xmin = min(xvals)
  xmax = max(xvals)
  legend.text = c()    # Initialize legend text vector
  legend.cols = c()    # Initialize legend color vector
  xmaxs = c()          # Initialize vector for storing maximum x-values
  
  # Loop through each sample in rt_list
  for (curr.sample in names(rt_list)){
    pre_tab = 16  # Set the prefix tab size for legend text
    
    # Check conditions for generating legend text based on readthrough data
    if (rt_list[[curr.sample]][['rt']] & 
        rt_list[[curr.sample]][['rt_max']] >= 1 & 
        rt_list[[curr.sample]][['rt_int']] >= 0.25){
      
      # Check if the fitted readthrough data meets certain criteria
      if ((rt_list[[curr.sample]][['dsfit']] | rt_list[[curr.sample]][['sfit']]) & 
          (rt_list[[curr.sample]][['rt_sum_fitted']] >= rt_list[[curr.sample]][['rt_sum']])) {
        
        if (rt_list[[curr.sample]][['sfit']]){
          leg.text = paste0('>', rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']]+1)
          tab = pre_tab - nchar(leg.text)
          legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
        } else {
          if (rt_list[[curr.sample]][['extrapolated']]){
            leg.text = paste0('*', rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']]+1)
            tab = pre_tab - nchar(leg.text)
            legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
          } else {
            leg.text = rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']] + 1
            tab = pre_tab - nchar(leg.text)
            legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
          }
        }
        
        # Calculate extra length and store maximum x-values
        xtra_length = rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1
        xmaxs = c(xmaxs, max(rt_list[[curr.sample]][['endDeclinePoint_x']], rt_list[[curr.sample]][['rt_end']]+xtra_length))
        
      } else {
        # Calculate extra length and store maximum x-values if criteria not met
        xtra_length = rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1
        xmaxs = c(xmaxs, rt_list[[curr.sample]][['rt_end']]+xtra_length)
        leg.text = paste0('#', rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1)
        tab = pre_tab - nchar(leg.text)
        legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
      }
    } else {
      # Set legend text to 'NA' if conditions not met
      leg.text = 'NA'
      tab = pre_tab - nchar(leg.text)
      legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
    }
    legend.cols = c(legend.cols, cols[[curr.sample]])  # Store legend color for the current sample
  }
  
  # Adjust xmax based on maximum x-values and filter xvals accordingly
  if (length(xmaxs) > 0){
    xmax = min(xmax, max(xmaxs))
    xvals = xvals[xvals<=xmax]
  }
  
  # Calculate indices for uTSS and dTES
  uTSS_plot = which(rownames(log2_GOI_data.nobatch.bodynorm.means)==as.character(uTSS)) - usExt
  dTES_plot = which(rownames(log2_GOI_data.nobatch.bodynorm.means)==as.character(dTES)) - usExt
  
  # Calculate the length of the readthrough region and plot width
  rt_region_length = xmax - dTES_plot
  plot_width = usExt + gen_info[GOI, 'length'] + rt_region_length
  
  # Determine if data needs to be binned based on plot width
  binning = FALSE
  if (plot_width > 10000){
    bin.size = 100
    no.bins = as.integer(plot_width/bin.size)
    binning = TRUE
  } else {
    if (length(plot_width) > 100){
      no.bins = 100
      bin.size = as.integer(length(plot_width)/no.bins)
      binning = TRUE
    }
  }
  
  # Bin data if necessary
  if (binning){
    binned.xvals = binning_function(0, no.bins, bin.size, bin.xvals=TRUE, usExt=usExt, statistic=statistic)
  } else {
    binned.xvals = xvals
  }
  
  # Adjust tick marks and labels on x-axis based on readthrough region length
  if (as.integer(rt_region_length/1E6) >= 4){
    stepsize = as.integer(rt_region_length/1E6/4)
    at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E6, dTES_plot+2*stepsize*1E6, dTES_plot+3*stepsize*1E6, dTES_plot+4*stepsize*1E6)
    label_vector = c('', '', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'mb')) 
  } else {
    if (as.integer(rt_region_length/1E5) >= 4){
      stepsize = as.integer(rt_region_length/1E5/4)
      at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E5, dTES_plot+2*stepsize*1E5, dTES_plot+3*stepsize*1E5, dTES_plot+4*stepsize*1E5)
      label_vector = c('', '', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), '00kb')) 
    } else {
      if (as.integer(rt_region_length/1E4) >= 4){
        stepsize = as.integer(rt_region_length/1E4/4)
        at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E4, dTES_plot+2*stepsize*1E4, dTES_plot+3*stepsize*1E4, dTES_plot+4*stepsize*1E4)
        label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), '0kb')) 
      } else {
        if (as.integer(rt_region_length/1E3) >= 4){
          stepsize = as.integer(rt_region_length/1E3/4)
          at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E3, dTES_plot+2*stepsize*1E3, dTES_plot+3*stepsize*1E3, dTES_plot+4*stepsize*1E3)
          label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'kb')) 
        } else {
          if (as.integer(rt_region_length/1E2) >= 4){
            stepsize = as.integer(rt_region_length/1E2/4)
            at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E2, dTES_plot+2*stepsize*1E2, dTES_plot+3*stepsize*1E2, dTES_plot+4*stepsize*1E2)
            label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+0.', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'kb')) 
          } else {
            at_vector = c(-usExt, 1, dTES_plot)
            label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES')
          }
        }
      }
    }
  }
  
  # Define control sample names
  ctrl_name = paste0(libType, '_', ctrl_sample)
  ctrl_rep_1 = paste0(ctrl_name, '_rep1')
  ctrl_rep_2 = paste0(ctrl_name, '_rep2')
  
  # Bin control sample data if necessary
  if (binning){
    ctrl_log2_binned.means = binning_function(log2_GOI_data.nobatch.bodynorm.means[,ctrl_name], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.rep1 = binning_function(log2_GOI_data.nobatch.bodynorm[,ctrl_rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.rep2 = binning_function(log2_GOI_data.nobatch.bodynorm[,ctrl_rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
  } else {
    ctrl_log2_binned.means = log2_GOI_data.nobatch.bodynorm.means[,ctrl_name]
    ctrl_log2_binned.rep1 = log2_GOI_data.nobatch.bodynorm[,ctrl_rep_1]
    ctrl_log2_binned.rep2 = log2_GOI_data.nobatch.bodynorm[,ctrl_rep_2]
  }
  
  # Loop through each sample in rt_list for further processing
  for (curr.sample in names(rt_list)){
    if (!is.na(rt_list[[curr.sample]][['best_fit']][1])){
      rt_start = rt_list[[curr.sample]][['rt_start']]
      rt_list[[curr.sample]][['best_fit']] = rt_list[[curr.sample]][['best_fit']][binned.xvals[binned.xvals >= rt_start]-rt_start+1]
    }
    rep_1 = paste0(curr.sample, '_rep1')
    rep_2 = paste0(curr.sample, '_rep2')
    
    # Bin sample data if necessary
    if (binning){
      log2_binned.means = binning_function(log2_GOI_data.nobatch.bodynorm.means[,curr.sample], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      log2_binned.rep1 = binning_function(log2_GOI_data.nobatch.bodynorm[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      log2_binned.rep2 = binning_function(log2_GOI_data.nobatch.bodynorm[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      subtr_log2_binned.means = binning_function(log2_GOI_data.nobatch.bodynorm.means.ctrlsubtract[,curr.sample], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      subtr_log2_binned.rep1 = binning_function(log2_GOI_data.nobatch.bodynorm.ctrlsubtract[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      subtr_log2_binned.rep2 = binning_function(log2_GOI_data.nobatch.bodynorm.ctrlsubtract[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    } else {
      log2_binned.means = log2_GOI_data.nobatch.bodynorm.means[,curr.sample]
      log2_binned.rep1 = log2_GOI_data.nobatch.bodynorm[,rep_1]
      log2_binned.rep2 = log2_GOI_data.nobatch.bodynorm[,rep_2]
      subtr_log2_binned.means = log2_GOI_data.nobatch.bodynorm.means.ctrlsubtract[,curr.sample]
      subtr_log2_binned.rep1 = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[,rep_1]
      subtr_log2_binned.rep2 = log2_GOI_data.nobatch.bodynorm.ctrlsubtract[,rep_2]
    }
    
    # Store processed data in rt_list for the current sample
    rt_list[[curr.sample]][['binned.xvals']] = binned.xvals
    rt_list[[curr.sample]][['at_vector']] = at_vector
    rt_list[[curr.sample]][['label_vector']] = label_vector
    rt_list[[curr.sample]][['log2_binned.means']] = log2_binned.means
    rt_list[[curr.sample]][['log2_binned.rep1']] = log2_binned.rep1
    rt_list[[curr.sample]][['log2_binned.rep2']] = log2_binned.rep2
    rt_list[[curr.sample]][['ctrl_log2_binned.means']] = ctrl_log2_binned.means
    rt_list[[curr.sample]][['ctrl_log2_binned.rep1']] = ctrl_log2_binned.rep1
    rt_list[[curr.sample]][['ctrl_log2_binned.rep2']] = ctrl_log2_binned.rep2
    rt_list[[curr.sample]][['subtr_log2_binned.means']] = subtr_log2_binned.means
    rt_list[[curr.sample]][['subtr_log2_binned.rep1']] = subtr_log2_binned.rep1
    rt_list[[curr.sample]][['subtr_log2_binned.rep2']] = subtr_log2_binned.rep2
  }
  
  #######################################################
  #### (9) plotting
  #######################################################
  if (plot){  # Check if plotting is enabled
    if (pdf){  # Check if PDF output is requested
      pdf_name = paste0(libType, '_', GOI, '.pdf')  # Define the PDF file name
      pdf(file=pdf_name, width=8, height=12)  # Open the PDF file for plotting
    }else{  # If PDF output is not requested
      dev.new()  # Open a new plotting device
    }
    par(mfrow=c(2,1), family = '')  # Set up a layout with 2 rows and 1 column for the plots
    if (norm=='1'){  # Check the normalization method
      plot.title1 = paste(GOI, '(normalized genebody)')  # Define plot title for the first panel
      plot.title2 = paste(GOI, '(normalized genebody - ctrl subtracted)')  # Define plot title for the second panel
    }else{
      if (norm=='2'){  # Check the normalization method
        plot.title1 = paste(GOI, '(TESnormalized genebody)')  # Define plot title for the first panel
        plot.title2 = paste(GOI, '(TESnormalized genebody - ctrl subtracted)')  # Define plot title for the second panel
      }else{
        plot.title1 = paste(GOI, '(unnormalized genebody)')  # Define plot title for the first panel
        plot.title2 = paste(GOI, '(unnormalized genebody - ctrl subtracted)')  # Define plot title for the second panel
      }
    }
    # Plot 1
    max.val = round(max(log2_GOI_data.nobatch.bodynorm), 1) + 1  # Calculate maximum y-axis value
    min.val = round(min(log2_GOI_data.nobatch.bodynorm), 1) - 1  # Calculate minimum y-axis value
    plot(0,0,type='n', xlim=c(xmin, xmax), ylim=c(min.val, max.val), main=plot.title1, xlab='position', ylab='norm. log2(cov)', las=1, xaxt='n')  # Set up the plot
    axis(side=1, at=at_vector, labels=label_vector)  # Add x-axis labels
    # Plot gene body coverage for each sample
    for (sample in paste0(libType, '_', c(ctrl_sample, viewSamples))){
      rep_1 = paste0(sample, '_rep1')  # Construct sample names for replicates
      rep_2 = paste0(sample, '_rep2')
      if (binning){  # Check if binning is enabled
        binned.means = binning_function(log2_GOI_data.nobatch.bodynorm.means[,sample], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)  # Bin gene body coverage data
        binned.rep1 = binning_function(log2_GOI_data.nobatch.bodynorm[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
        binned.rep2 = binning_function(log2_GOI_data.nobatch.bodynorm[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
      }else{
        binned.means = log2_GOI_data.nobatch.bodynorm.means[,sample]  # Extract gene body coverage data
        binned.rep1 = log2_GOI_data.nobatch.bodynorm[,rep_1]
        binned.rep2 = log2_GOI_data.nobatch.bodynorm[,rep_2]
      }
      # Plot gene body coverage
      lines(binned.xvals, binned.means[1:length(binned.xvals)], lwd=2, col=cols[[sample]])  # Plot mean coverage
      polygon(c(binned.xvals, rev(binned.xvals)), c(binned.rep1[1:length(binned.xvals)], rev(binned.rep2[1:length(binned.xvals)])), col=cols_trans[[sample]], border = NA)  # Plot replicates
    }
    # Add vertical lines for TSS and TES
    abline(v=TSS_row-usExt, lty='dotted', col='black')
    abline(v=TES_row-usExt, lty='dotted', col='black')
    ###
    n = 0  # Initialize counter for additional regions of interest
    for (curr.sample in names(rt_list)){  # Iterate over each sample in the list of regions of interest
      # Check if the sample has valid region of interest data
      if (rt_list[[curr.sample]][['rt']] & rt_list[[curr.sample]][['rt_max']] >= 1 & rt_list[[curr.sample]][['rt_int']] >= 0.25){
        rt_start = rt_list[[curr.sample]][['rt_start']]  # Get start position of the region of interest
        rt_end = rt_list[[curr.sample]][['rt_end']]  # Get end position of the region of interest
        abline(v=rt_start, lty='dashed', col=cols_trans[[curr.sample]])  # Add dashed line for start position
        lines(c(rt_start, rt_end), rep(min.val+n*0.25, 2), col=cols_trans[[curr.sample]], lwd=5, lend=1)  # Plot region of interest
        # Check if the region of interest has fitted data
        if ((rt_list[[curr.sample]][['dsfit']] | rt_list[[curr.sample]][['sfit']]) & (rt_list[[curr.sample]][['rt_sum_fitted']] >= rt_list[[curr.sample]][['rt_sum']])){
          rt_end_fitted = rt_list[[curr.sample]][['rt_end_fitted']]  # Get end position of the fitted region
          # Plot fitted region if it falls within the plot limits
          if (rt_end_fitted <= xmax){
            abline(v=rt_end_fitted, lty='dashed', col=cols_trans[[curr.sample]])  # Add dashed line for end position of fitted region
            if (rt_end_fitted > rt_end){  # Check if end of fitted region exceeds end of observed region
              lines(c(rt_end, rt_end_fitted), rep(min.val+n*0.25, 2), col=cols_trans[[curr.sample]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
            }
          }else{
            if (rt_end < xmax & rt_end_fitted > rt_end){  # Check if observed region ends before plot limit and fitted region exceeds it
              lines(c(rt_end, xmax), rep(min.val+n*0.25, 2), col=cols_trans[[curr.sample]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
            }
          }
        }
        n = n + 1  # Increment counter for additional regions of interest
      }
    }
    # Add legend to the plot
    par(family = 'mono')
    legend('topright', fill=legend.cols, legend=legend.text, cex=0.75)  # Add legend
    par(family = '')
    # Add horizontal line at y=0
    abline(h=0, col='gray')
  }
  if (pdf){  # Check if PDF output was requested
    dev.off()  # Close the PDF file
  }
  
  #######################################################
  #### (10) done
  #######################################################
  return(rt_list)
}

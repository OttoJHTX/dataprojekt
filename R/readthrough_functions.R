sample_subtract_hmm = function(curr.data.list, TSS_row, TES_row, TES_rows, col='gray', usExt=5000, statistic='median', cutoff=2, plot=FALSE, ymax=-2, max_cv=200){
  ## Model initialization
  hmm = tryCatch({suppressWarnings(initHMM(curr.data.list, nStates=2, "IndependentGaussian", sharedCov=TRUE))}, error = function(e) {return(NA)})
  ## Model fitting
  if (class(hmm)[1]=='HMM'){
    hmm_fitted = tryCatch({suppressWarnings(fitHMM(curr.data.list, hmm, maxIters=50))}, error = function(e) {return(NA)})
  }else{
    hmm_fitted = NA
  }
  if (class(hmm_fitted)[1]=='HMM'){
    ## Calculate state path
    viterbi = getViterbi(hmm_fitted, curr.data.list)
    states = as.integer(viterbi[['GOI']])
    data.means = rowMeans(curr.data.list[['GOI']])
    if (plot){  ###@@@
      plot(1:length(data.means), data.means, type='l', col=col, main=GOI)
    }
    state1 = ifelse(statistic=='median', median(data.means[states==1]), mean(data.means[states==1]))
    state2 = ifelse(statistic=='median', median(data.means[states==2]), mean(data.means[states==2]))
    on.state = ifelse(state1 > state2, 1, 2)
    off.state = ifelse(state1 > state2, 2, 1)
    off.state.int = ifelse(statistic=='median', median(data.means[states==off.state]), mean(data.means[states==off.state]))
    states_iv = rle(states)
    #determine one unbroken on-state
    states.list = list()
    states.list[['state']] = rep(NA, length(states_iv$values))
    states.list[['mean']] = rep(NA, length(states_iv$values))
    states.list[['median']] = rep(NA, length(states_iv$values))
    states.list[['max']] = rep(NA, length(states_iv$values))
    states.list[['sd']] = rep(NA, length(states_iv$values))
    states.list[['window']] = list()
    state.start = 1
    for (i in 1:length(states_iv$values)){
      state.length = states_iv$lengths[i]
      state.end = state.start-1+state.length
      state.iv = state.start:state.end
      states.list[['state']][i] = ifelse(states_iv$values[i]==on.state, 1, 0)
      states.list[['mean']][i] = mean(data.means[state.iv])
      states.list[['median']][i] = median(data.means[state.iv])
      states.list[['max']][i] = max(data.means[state.iv])
      states.list[['sd']][i] = sd(data.means[state.iv])
      states.list[['window']][[i]] = c(state.start, state.end)
      if (plot){
        if (states_iv$values[i]==on.state){
          lines(c(state.start, state.end-1), rep(ymax, 2), col=col, lwd=5, lend=1)
        }else{
          lines(c(state.start, state.end-1), rep(ymax-0.2, 2), col='black', lwd=5, lend=1)
        }
      }
      state.start = state.end + 1
    }
    states.on = which(states.list[['state']]==1)
    eligible.on.states = c()
    if (length(states.on) > 0){
      for (state.on in states.on){
        state.end = states.list[['window']][[state.on]][2]
        if (state.end > min(TES_rows)-TSS_row+1){         ###@@@ 190415 change (state.end > TES_row-TSS_row+1)  the on-state end needs to go beyond the most proximal TES
          eligible.on.states = c(eligible.on.states, state.on)
        }
      }
    }
    if (length(eligible.on.states) > 0){
      max.state = which(states.list[[statistic]]==max(states.list[[statistic]][eligible.on.states]))
      if (length(max.state) > 1){
        max.states = c()
        for (max.state.no in max.state){
          max.state.iv = states.list[['window']][[max.state.no]]
          if (max.state.iv[2] > TES_row-TSS_row+1){
            max.states = c(max.states, max.state.no)
          }
        }
        if (length(max.states) > 0){
          max.state = max.states[1]
        }else{
          max.state = rev(max.state)[1]
        }
      }
      max.state.iv = states.list[['window']][[max.state]]
      if (plot){
        lines(max.state.iv, rep(ymax+0.2, 2), col='red', lwd=5, lend=1)
      }
      ext.state.conf = rep(NA, length(states_iv$values))
      ext.state.cvs = rep(NA, length(states_iv$values))
      if (length(states.on[states.on <= max.state]) > 0){
        ext.state.end = max.state.iv[2]
        for (state.on in states.on[states.on <= max.state]){
          ext.state.start = states.list[['window']][[state.on]][1]
          ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)   ### coefficient of variation
          ext.state.conf[state.on] = ext.state.mean - cutoff*ext.state.sd
          ext.state.cvs[state.on] = ext.state.cv
        }
      }
      if (length(states.on[states.on > max.state]) > 0){
        ext.state.start = max.state.iv[1]
        for (state.on in states.on[states.on > max.state]){
          ext.state.end = states.list[['window']][[state.on]][2]
          ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)
          ext.state.conf[state.on] = ext.state.mean - cutoff*ext.state.sd
          ext.state.cvs[state.on] = ext.state.cv
        }
      }
      combined.on.states = which(ext.state.conf > 0 | (ext.state.cvs<cutoff*ext.state.cvs[max.state] & ext.state.cvs>0 & ext.state.cvs<max_cv))
      if (length(combined.on.states) > 0){
        min_int_ds_rt = ifelse(length(states_iv$lengths) > rev(combined.on.states)[1], min(states.list[[statistic]][rev(combined.on.states)[1]:length(states_iv$lengths)]), NA)
        combined.on.state.end = states.list[['window']][[max(combined.on.states)]][2]    ## relative to current data matrix
        if (combined.on.state.end > min(TES_rows)-TSS_row+1){        ###@@@ the end of the combined read-through should be beyond the end of most proximal TES
          combined.on.state.start = states.list[['window']][[min(combined.on.states)]][1]  ## relative to current data matrix
          combined.on.state.start = ifelse(combined.on.state.start > min(TES_rows)-TSS_row, combined.on.state.start, min(TES_rows)-TSS_row+1) ###@@@ 
          combined.region.mean = mean(curr.data.list[['GOI']][combined.on.state.start:combined.on.state.end, ])
          combined.region.sd = sd(curr.data.list[['GOI']][combined.on.state.start:combined.on.state.end, ])
          combined.region.cv = round(100*combined.region.sd/combined.region.mean, 2)
          rt = TRUE
          if (plot){
            lines(c(combined.on.state.start, combined.on.state.end), rep(ymax-0.5, 2), col='darkgreen', lwd=5, lend=1)
          }
          ## determine which TES this on-state most likely relates to
          new.TES_rows = TES_rows-TSS_row+1
          new.TES_row = TES_row-TSS_row+1
          new.TES_rows = sort(new.TES_rows[new.TES_rows <= combined.on.state.start & new.TES_rows > 0])  ###@@@ sort(new.TES_rows[new.TES_rows <= combined.on.state.start & new.TES_rows > 0 & new.TES_rows >= new.TES_row])
          new.TES_rows = setdiff(new.TES_rows-1, new.TES_rows) + 1  ###@@@ remove TESs that are only shifted by 1 position
          if (plot){
            abline(v=new.TES_rows, lty='dotted', col='green')
          }
          TES.ivs = list()
          n = 0
          if (length(new.TES_rows) > 1){
            for (i in 1:(length(new.TES_rows)-1)){
              TES.region.start = new.TES_rows[i]
              TES.region.end = new.TES_rows[i+1] - 1
              TES.region.mean = mean(curr.data.list[['GOI']][TES.region.start:TES.region.end, ])
              TES.region.sd = sd(curr.data.list[['GOI']][TES.region.start:TES.region.end, ])
              TES.region.cv = round(100*TES.region.sd/TES.region.mean, 2)
              if (TES.region.mean > 0){  
                n = n + 1
                TES.ivs[[n]] = c(TES.region.start, TES.region.end)
              }
            }
          }
          ext.state.end = combined.on.state.end
          if (length(TES.ivs) > 0){
            ext.TES.states.cv = rep(NA, length(TES.ivs))
            ext.TES.states.conf = rep(NA, length(TES.ivs))
            for (i in length(TES.ivs):1){
              ext.state.start = TES.ivs[[i]][1]
              ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
              ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
              ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)
              ext.TES.states.cv[i] = ext.state.cv 
              ext.TES.states.conf[i] = ext.state.mean - cutoff*ext.state.sd
              #cat(paste(paste0(ext.state.start, '-', ext.state.end, ':'), ext.state.mean, ext.state.sd, ext.state.cv), '\n')
            }
            combined.ext.TES.states = which(ext.TES.states.conf > 0 | (ext.TES.states.cv<cutoff*combined.region.cv & ext.TES.states.cv>0  & ext.TES.states.cv<max_cv)) # combined.ext.TES.states = which(ext.TES.states > 0)
            if (length(combined.ext.TES.states) != 0){
              combined.ext.TES.states.start = TES.ivs[[min(combined.ext.TES.states)]][1]  ## relative to current data matrix
            }else{
              if (length(new.TES_rows) > 0){
                TESs.dists = combined.on.state.start - new.TES_rows
                us.TESs = new.TES_rows[which(TESs.dists >= 0)]
                combined.ext.TES.states.start = us.TESs[which((combined.on.state.start - us.TESs) == min(combined.on.state.start - us.TESs))]
              }else{
                combined.ext.TES.states.start = new.TES_row
              }
            }
          }else{
            if (length(new.TES_rows) > 0){
              TESs.dists = combined.on.state.start - new.TES_rows
              us.TESs = new.TES_rows[which(TESs.dists >= 0)]
              combined.ext.TES.states.start = us.TESs[which((combined.on.state.start - us.TESs) == min(combined.on.state.start - us.TESs))]
            }else{
              combined.ext.TES.states.start = new.TES_row
            }
          }
          rt_int = ifelse(statistic=='median', median(data.means[(combined.ext.TES.states.start+1):combined.on.state.end]), mean(data.means[(combined.ext.TES.states.start+1):combined.on.state.end]))  ## rt-region starts one nucleotide downstream of chosen TES
          if (rt_int <= 0){
            rt = FALSE
          }
          if (plot){
            lines(c(combined.ext.TES.states.start, combined.on.state.end), rep(ymax+0.5, 2), col='purple', lwd=5, lend=1)
            abline(v=combined.ext.TES.states.start, lty='dashed', col='gray')
            abline(v=combined.on.state.end, lty='dashed', col='gray')
          }
        }else{
          rt = FALSE
        }
      }else{
        rt = FALSE
      }
    }else{
      rt = FALSE
    }
  }else{
    rt = FALSE
  }
  if (rt){
    output = list('rt_iv'=c(combined.ext.TES.states.start, combined.on.state.end)+TSS_row-1, 'rt_int'=rt_int, 'min_int'=min(states.list[[statistic]]), 'max_int'=states.list[[statistic]][max.state], 'max_int_iv'=states.list[['window']][[max.state]]+TSS_row-1, 'off.state_int'=off.state.int, 'min_int_ds_rt'=min_int_ds_rt)  ## the rownumber interval for which readthrough is detected
  }else{
    output = list()
  }
  return(output)
}  

adj_function = function(n, fit_int_standard, rt_region_means, sample_states, asymptote=FALSE, fit_asymptote_standard=0){
  fit_int_zero = fit_int_standard-min(fit_int_standard)
  fit_asymptote_zero = fit_asymptote_standard-min(fit_int_standard)
  if (n==1){
    # 1 (the unadjusted fitted double sigmoidal) 
    fit_int = fit_int_standard
    fit_asymptote = fit_asymptote_standard
  }
  if (n==2){
    ## 2 (adjust using max and min values of data)
    max.int = max(rt_region_means)
    min.int = min(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==3){
    ## 3 adjust using max and 'min_int_ds_rt' (else 'min_int') values of data
    max.int = max(rt_region_means)
    min.int = ifelse(is.na(sample_states[['min_int_ds_rt']]), sample_states[['min_int']], sample_states[['min_int_ds_rt']])    
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==4){
    ## 4 adjust using max_int and 'min_int_ds_rt' (else 'min_int') values of data
    max.int = sample_states[['max_int']]
    min.int = ifelse(is.na(sample_states[['min_int_ds_rt']]), sample_states[['min_int']], sample_states[['min_int_ds_rt']])
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==5){
    ## 5 adjust using max and median values of data
    max.int = max(rt_region_means)
    min.int = median(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==6){
    ## 6 adjust using max_int and median values of data
    max.int = sample_states[['max_int']]
    min.int = median(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==7){
    ## 7 adjust using median value of data
    median.int = median(rt_region_means)
    fit_int = median.int*fit_int_standard/median(fit_int_standard)
    fit_asymptote = median.int*fit_asymptote_standard/median(fit_int_standard)
  }
  if (n==8){
    ## 8 adjust using mean value of data
    mean.int = mean(rt_region_means)
    fit_int = mean.int*fit_int_standard/mean(fit_int_standard)
    fit_asymptote = mean.int*fit_asymptote_standard/mean(fit_int_standard)
  }
  if (asymptote){
    return(fit_asymptote)
  }else{
    return(fit_int)
  }
}


binning_function =  function(data.vector, no.bins, bin.size, bin.xvals=FALSE, usExt=5000, statistic='median'){
  binned.data = rep(0, no.bins)
  for (bin.no in 1:no.bins){
    iv = ((bin.no-1)*bin.size+1):(bin.no*bin.size)
    if (bin.xvals){
      binned.data[bin.no] = -usExt + as.integer(mean(iv))
    }else{
      binned.data[bin.no] = ifelse(statistic=='median', median(data.vector[iv]), mean(data.vector[iv]))
    }
  }
  return(binned.data)
}

norm_types = list('0'='without', '1'='genebody', '2'='TES')

rt_analysis = function(ctrl_dir, sample_dir, annot_gr, ctrl, sample, GOI, estimate_TES=FALSE, batch=NULL, norm='1',  statistic='median', usExt=5000, plot_data=TRUE, verbose=TRUE, verbose2=TRUE, pdf=FALSE){
  ##########################################################
  #### (1) read in gene specific data
  ########################################################## 
  
  # Output the index of the gene of interest (GOI) and its name if verbose or verbose2 is TRUE
  if (verbose | verbose2){
    cat( GOI, sep='\n')
  }
  
  # Output a tab or newline character based on the value of verbose
  if (verbose){
    cat(ifelse(verbose, '\t', '\n'))
  }
  
  # Create an empty list to store data for each replicate
  rt_list = list()
  
  ###################### KAOS-KODE #########################
  
  strand_options = c("+"="plus", "-"="minus")
  
  # Subsetting annotation on gene name
  gene_annot = subset(annot_gr, gene_name == GOI)
  chrom_no = as.character(seqnames(gene_annot)@values)
  strand_sign = as.character(strand(gene_annot)@values)
  chrom_no_loop = paste0("chr", chrom_no) # TEMP FIX
  
  # Coordinates
  start_coord = min(start(gene_annot))
  end_coord = max(end(gene_annot))
  
  # Find coord of next gene
  next_gene_annot = sort(subset(annot_gr, seqnames %in% c(chrom_no, chrom_no_loop) & strand == strand_sign))
  if (strand_sign == "+") {
    s = start(next_gene_annot)
    e = end(next_gene_annot)
    next_gene_coord = s[which(s > end_coord)[1]]
    prev_gene_coord = e[rev(which(e < start_coord))[1]]
  } else {
    s = end(next_gene_annot)
    e = start(next_gene_annot)
    next_gene_coord = s[rev(which(s < start_coord))[1]]
    prev_gene_coord = e[which(e > end_coord)[1]]
  }
  
  start_coord_ext =  ifelse(strand_sign == "+", min(start(gene_annot)) - usExt, next_gene_coord + 1)
  end_coord_ext = ifelse(strand_sign == "+", next_gene_coord - 1, max(end(gene_annot)) + usExt)
  
  
  # Empty DF with correct length
  indexes <- seq(start_coord_ext, end_coord_ext)
  GOI_data <- data.frame(index = indexes)
  
  # Loop for loading replicates and saving to DF
  for (fname in ctrl) {
    ctrl_fname = paste0(ctrl_dir, fname, "_", strand_options[strand_sign], ".bw")
    ctrl_bw = import(ctrl_fname, 
                     which = GenomicRanges::GRanges(seqnames = chrom_no_loop, 
                                                    ranges = IRanges::IRanges(start = start_coord_ext, end = end_coord_ext)), 
                     as = "NumericList")[[1]]
    GOI_data[[fname]] = ctrl_bw
  }
  
  # Loop for loading replicates and saving to DF
  for (fname in sample) {
    sample_fname = paste0(sample_dir, fname, "_", strand_options[strand_sign], ".bw")
    sample_bw = import(sample_fname, 
                       which = GenomicRanges::GRanges(seqnames = chrom_no_loop, 
                                                      ranges = IRanges::IRanges(start = start_coord_ext, end = end_coord_ext)), 
                       as = "NumericList")[[1]]
    GOI_data[[fname]] = sample_bw
  }
  
  # Setting the DF indices to the positions of the reads and deleting index column
  if (strand_sign == "-"){
    GOI_data = GOI_data[order(GOI_data$index, decreasing=TRUE), ]
  }
  rownames(GOI_data) = GOI_data$index
  GOI_data$index = NULL

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
  if (!is.null(batch)) {
    log2_GOI_data = limma::removeBatchEffect(log2_GOI_data, batch=batch)
  } 
  # Remove the original log2-transformed gene of interest data from memory to save space
  
  
  ###############################################################
  #### (4) Normalize to gene body signal
  ###############################################################
  
  #prev_gene_coord
  
  transcript_annot = subset(gene_annot, type == "transcript")
  if (estimate_TES){
      if (length(transcript_annot) > 1) {
      outer_left = max(min(next_gene_coord, prev_gene_coord), min(as.integer(rownames(log2_GOI_data))))
      outer_right = min(max(next_gene_coord, prev_gene_coord), max(as.integer(rownames(log2_GOI_data)))) # outer_right = max(as.integer(rownames(log2_GOI_data))
      lefts = sort(unique(start(transcript_annot)))
      rights = sort(unique(end(transcript_annot)))
      left_diffs = rep(NA, length(lefts))
      right_diffs = rep(NA, length(rights))
      for (i in 1:length(lefts)){
        if ( any(lefts > lefts[i] + 50) ){
          iv_down = (lefts[i]+1):lefts[which(lefts > lefts[i] + 50)[1]]
        }else{
          iv_down = (lefts[i]+1):rights[which(rights > lefts[i] + 50)[1]]
        }
        downleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
        
        left_diff_vector = c()
        if ( any(lefts < lefts[i] - 50) ){
          for (idx in which(lefts < lefts[i] - 50)){
            iv_up = lefts[idx]:lefts[i]
            upleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
            left_diff_vector = c(left_diff_vector, downleft_ctrl_log2_GOI_mean - upleft_ctrl_log2_GOI_mean)
          }
        }else{
          iv_up = outer_left:lefts[i]
          upleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
          left_diff_vector = c(left_diff_vector, downleft_ctrl_log2_GOI_mean - upleft_ctrl_log2_GOI_mean)
        }
        left_diffs[i] = min(left_diff_vector)
      }
      for (j in 1:length(rights)){
        if ( any(rights < rights[j] - 50) ){
          iv_up = rights[rev(which(rights < rights[j] - 50))[1]]:(rights[j]-1)
        }else{
          iv_up = lefts[rev(which(lefts < rights[j] - 50))[1]]:(rights[j]-1)
        }
        upright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
        
        right_diff_vector = c()
        if ( any(rights > rights[j] + 50) ){
          for (idx in which(rights > rights[j] + 50)){
            iv_down = rights[j]:rights[idx]
            downright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
            right_diff_vector = c(right_diff_vector, upright_ctrl_log2_GOI_mean - downright_ctrl_log2_GOI_mean)
          }
        }else{
          iv_down = rights[j]:outer_right
          downright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
          right_diff_vector = c(right_diff_vector, upright_ctrl_log2_GOI_mean - downright_ctrl_log2_GOI_mean)
        }
        right_diffs[j] = min(right_diff_vector)
      }
  
      max_diffs_left = which(left_diffs == max(left_diffs))
      max_diffs_right = which(right_diffs == max(right_diffs))
    }
  }
  final_range = range(transcript_annot)
  if (strand_sign == "+") {
    uTSS = min(start(transcript_annot))
    dTES = max(end(transcript_annot))
    if (estimate_TES){
      TESs = rights[max_diffs_right]
      TSSs = lefts[max_diffs_left]
      start(final_range) = min(TSSs)
      end(final_range) = max(TESs)
    }
  } else {
    uTSS = max(end(transcript_annot))
    dTES = min(start(transcript_annot))
    if (estimate_TES){
      TESs = lefts[max_diffs_left]
      TSSs = rights[max_diffs_right]
      start(final_range) = min(TESs)
      end(final_range) = max(TSSs)
    }
  }
  
  if (strand_sign == "+") {
    TSS = start(final_range)
    TES = end(final_range)
  }else{
    TES = start(final_range)
    TSS = end(final_range)
  }
  uTSS_row = which(rownames(log2_GOI_data)==as.character(uTSS))
  dTES_row = which(rownames(log2_GOI_data)==as.character(dTES))
  TSS_row = which(rownames(log2_GOI_data)==as.character(TSS))
  TES_row = which(rownames(log2_GOI_data)==as.character(TES))
  
  # Calculate body normalization factors and adjust data accordingly
  body.norm.factors = apply(log2_GOI_data[TSS_row:TES_row,], 2, median)
  
  sample_rep_diffs = rep(NA, length(ctrl))
  for (i in 1:length(ctrl)) {
    sample_rep_diffs[i] = as.numeric(body.norm.factors[sample[i]] - body.norm.factors[ctrl[i]])
  }
  body_diff = mean(sample_rep_diffs)
  rt_list[["sample"]] = list('TSS'=TSS_row-usExt, 'TES'=TES_row-usExt, 'rep_diffs'=sample_rep_diffs, 'body_diff'=body_diff)
  
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
    body.norm.factors = apply(log2_GOI_data[upTES_row:TES_row,], 2, median)
  }
  
  # Adjust data based on the calculated body normalization factors
  log2_GOI_data.bodynorm = t(t(log2_GOI_data) - body.norm.factors)
  
  # Calculate mean coverage for each sample
  log2_GOI_data.bodynorm.means = matrix(0, nrow=nrow(log2_GOI_data.bodynorm), ncol=2)
  colnames(log2_GOI_data.bodynorm.means) = c("ctrl", "sample")
  rownames(log2_GOI_data.bodynorm.means) = rownames(log2_GOI_data)
  log2_GOI_data.bodynorm.sd = log2_GOI_data.bodynorm.means
  log2_GOI_data.bodynorm.means[, "ctrl"] = rowMeans(log2_GOI_data.bodynorm[, ctrl])
  log2_GOI_data.bodynorm.means[, "sample"] = rowMeans(log2_GOI_data.bodynorm[, sample])
  log2_GOI_data.bodynorm.sd[, "ctrl"] = apply(log2_GOI_data.bodynorm[, ctrl], 1, sd)
  log2_GOI_data.bodynorm.sd[, "sample"] = apply(log2_GOI_data.bodynorm[, sample], 1, sd)
  
  # Remove the original log2-transformed gene of interest data from memory to save space
  rm(log2_GOI_data)
  
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
  
  for (curr.sample in c("ctrl", "sample")){
    body_mean  = mean(log2_GOI_data.bodynorm.means[TSS_row:TES_row, curr.sample])
    if (within.gene){
      us.TSS_mean = mean(log2_GOI_data.bodynorm.means[us.TSS_rows, curr.sample])
      us.TSS_diff = body_mean - us.TSS_mean
    }else{
      us.TSS_diff = NA
    }
    us_mean = mean(log2_GOI_data.bodynorm.means[us_rows, curr.sample])
    us_diff = body_mean - us_mean
    
    # Store the differences between gene body and upstream signal for each sample in the result list
    if (curr.sample == "sample"){
      rt_list[[curr.sample]][['us_diff']] = us_diff           
      rt_list[[curr.sample]][['us.TSS_diff']] = us.TSS_diff   
    }
  }
  
  ###############################################################
  #### (6) Subtract control signal
  ###############################################################
  
  # Subtract control signal from the gene body and remove control columns
  log2_GOI_data.bodynorm.means.ctrlsubtract = data.frame("sample" = log2_GOI_data.bodynorm.means[, "sample"] - log2_GOI_data.bodynorm.means[, "ctrl"])
  
  # Subtract control signal from the replicate columns
  log2_GOI_data.bodynorm.ctrlsubtract = log2_GOI_data.bodynorm
  
  for (i in 1:length(ctrl)) {
    log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] = log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] - log2_GOI_data.bodynorm.ctrlsubtract[, i] 
    log2_GOI_data.bodynorm.ctrlsubtract[, i] = log2_GOI_data.bodynorm.ctrlsubtract[, i] - log2_GOI_data.bodynorm.ctrlsubtract[, i] 
  }
  log2_GOI_data.bodynorm.sd.ctrlsubtract = data.frame("sample" = apply(log2_GOI_data.bodynorm.ctrlsubtract[,sample], 1, sd))
  
  # Note: All dataframes/matrices up until this point contain usExt (max 5000) positions upstream of locus of interest, all positions from locus of interest, and from 'allowed' downstream region.
  # uTSS starts at row usExt + 1 (default 5001)
  
  ###############################################################
  #### (7) Find states using HMM and fit double-sigmoidal to data
  ###############################################################
  # This code segment iterates through samples, detects putative readthrough events using Hidden Markov Models (HMM), and fits sigmoidal or double sigmoidal functions to the readthrough data to estimate its length and intensity. It also handles cases where no readthrough is detected.
  
  for (curr.sample in c("sample")){
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
    
    samples_reps = sample
    curr.sample_data = log2_GOI_data.bodynorm.ctrlsubtract[, samples_reps]  ##@@ same nrow as above
    hmm_result = sample_subtract_hmm(curr.data.list = list('GOI' = curr.sample_data[TSS_row:nrow(curr.sample_data), ]), TSS_row, TES_row, TES_rows=TES_row, col = cols_trans[[curr.sample]], usExt = usExt, statistic = statistic, cutoff = 2, plot = FALSE, ymax = 0.5, max_cv = 200)
    
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
  
  xvals = 1:nrow(log2_GOI_data.bodynorm.means)-usExt
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
  uTSS_plot = which(rownames(log2_GOI_data.bodynorm.means)==as.character(uTSS)) - usExt
  dTES_plot = which(rownames(log2_GOI_data.bodynorm.means)==as.character(dTES)) - usExt
  
  # Calculate the length of the readthrough region and plot width
  rt_region_length = xmax - dTES_plot
  #@plot_width = usExt + (TES_row - TSS_row + 1) + rt_region_length
  plot_width = usExt + (dTES_row - uTSS_row + 1) + rt_region_length
  
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
  
  ###@@@
  # Bin control sample data if necessary
  if (binning){
    ctrl_log2_binned.means = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.max = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"] + log2_GOI_data.bodynorm.sd[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.min = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"] - log2_GOI_data.bodynorm.sd[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    
    #ctrl_log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm[,ctrl[1]], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #ctrl_log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm[,ctrl[2]], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
  } else {
    ctrl_log2_binned.means = log2_GOI_data.bodynorm.means[,"ctrl"]
    ctrl_log2_binned.max = log2_GOI_data.bodynorm.means[,"ctrl"] + log2_GOI_data.bodynorm.sd[,"ctrl"]
    ctrl_log2_binned.min = log2_GOI_data.bodynorm.means[,"ctrl"] - log2_GOI_data.bodynorm.sd[,"ctrl"]
    
    #ctrl_log2_binned.rep1 = log2_GOI_data.bodynorm[,ctrl[1]]
    #ctrl_log2_binned.rep2 = log2_GOI_data.bodynorm[,ctrl[2]]
  }
  
  # ... "sample" in rt_list for further processing
  if (!is.na(rt_list[["sample"]][['best_fit']][1])){
    rt_start = rt_list[["sample"]][['rt_start']]
    rt_list[["sample"]][['best_fit']] = rt_list[["sample"]][['best_fit']][binned.xvals[binned.xvals >= rt_start]-rt_start+1]
  }
  # rep_1 = sample[1]
  # rep_2 = sample[2]
  
  # Bin sample data if necessary
  if (binning){
    log2_binned.means = binning_function(log2_GOI_data.bodynorm.means[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    log2_binned.max = binning_function(log2_GOI_data.bodynorm.means[,"sample"] + log2_GOI_data.bodynorm.sd[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    log2_binned.min = binning_function(log2_GOI_data.bodynorm.means[,"sample"] - log2_GOI_data.bodynorm.sd[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    
    
    #log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.means = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.max = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] + log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.min = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] - log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #subtr_log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm.ctrlsubtract[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #subtr_log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm.ctrlsubtract[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
  } else {
    log2_binned.means = log2_GOI_data.bodynorm.means[,"sample"]
    log2_binned.max = log2_GOI_data.bodynorm.means[,"sample"] + log2_GOI_data.bodynorm.sd[,"sample"]
    log2_binned.min = log2_GOI_data.bodynorm.means[,"sample"] - log2_GOI_data.bodynorm.sd[,"sample"]
    
    #log2_binned.rep1 = log2_GOI_data.bodynorm[,rep_1]
    #log2_binned.rep2 = log2_GOI_data.bodynorm[,rep_2]
    subtr_log2_binned.means = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"]
    subtr_log2_binned.max = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] + log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"]
    subtr_log2_binned.min = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] - log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"]
    
    #subtr_log2_binned.rep1 = log2_GOI_data.bodynorm.ctrlsubtract[,rep_1]
    #subtr_log2_binned.rep2 = log2_GOI_data.bodynorm.ctrlsubtract[,rep_2]
  }
  
  # Store processed data in rt_list for the current sample
  rt_list[["sample"]][['binned.xvals']] = binned.xvals
  rt_list[["sample"]][['at_vector']] = at_vector
  rt_list[["sample"]][['label_vector']] = label_vector
  rt_list[["sample"]][['log2_binned.means']] = log2_binned.means
  rt_list[["sample"]][['log2_binned.max']] = log2_binned.max
  rt_list[["sample"]][['log2_binned.min']] = log2_binned.min
  rt_list[["sample"]][['ctrl_log2_binned.means']] = ctrl_log2_binned.means
  rt_list[["sample"]][['ctrl_log2_binned.max']] = ctrl_log2_binned.max
  rt_list[["sample"]][['ctrl_log2_binned.min']] = ctrl_log2_binned.min
  rt_list[["sample"]][['subtr_log2_binned.means']] = subtr_log2_binned.means
  rt_list[["sample"]][['subtr_log2_binned.max']] = subtr_log2_binned.max
  rt_list[["sample"]][['subtr_log2_binned.min']] = subtr_log2_binned.min
  
  #######################################################
  #### (9) plotting
  #######################################################
  if (plot_data){  # Check if plotting is enabled
    if (pdf){  # Check if PDF output is requested
      pdf_name = paste0(wd, GOI, '.pdf')  # Define the PDF file name
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
    max.val = round(max(c(ctrl_log2_binned.max, log2_binned.max)), 1) + 1  # Calculate maximum y-axis value
    min.val = round(min(c(ctrl_log2_binned.min, log2_binned.min)), 1) - 1  # Calculate minimum y-axis value
    plot(0,0,type='n', xlim=c(xmin, xmax), ylim=c(min.val, max.val), main=plot.title1, xlab='position', ylab='norm. log2(cov)', las=1, xaxt='n')  # Set up the plot
    axis(side=1, at=at_vector, labels=label_vector)  # Add x-axis labels
    # Plot gene body coverage for "each sample "ctrl"
    lines(binned.xvals, ctrl_log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["ctrl"]])  # Plot mean coverage
    polygon(c(binned.xvals, rev(binned.xvals)), c(ctrl_log2_binned.max[1:length(binned.xvals)], rev(ctrl_log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["ctrl"]], border = NA)  # Plot replicates
    # Plot gene body coverage for "each sample "sample"
    lines(binned.xvals, log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["sample"]])  # Plot mean coverage
    polygon(c(binned.xvals, rev(binned.xvals)), c(log2_binned.max[1:length(binned.xvals)], rev(log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["sample"]], border = NA)  # Plot replicates
    # Add vertical lines for TSS and TES
    abline(v=TSS_row-usExt, lty='dotted', col='black')
    abline(v=TES_row-usExt, lty='dotted', col='black')
    ###
    # Check if the sample has valid region of interest data
    if (rt_list[["sample"]][['rt']] & rt_list[["sample"]][['rt_max']] >= 1 & rt_list[["sample"]][['rt_int']] >= 0.25){
      rt_start = rt_list[["sample"]][['rt_start']]  # Get start position of the region of interest
      rt_end = rt_list[["sample"]][['rt_end']]  # Get end position of the region of interest
      abline(v=rt_start, lty='dashed', col=cols_trans[["sample"]])  # Add dashed line for start position
      lines(c(rt_start, rt_end), rep(min.val, 2), col=cols_trans[["sample"]], lwd=5, lend=1)  # Plot region of interest
      # Check if the region of interest has fitted data
      if ((rt_list[["sample"]][['dsfit']] | rt_list[["sample"]][['sfit']]) & (rt_list[["sample"]][['rt_sum_fitted']] >= rt_list[["sample"]][['rt_sum']])){
        rt_end_fitted = rt_list[["sample"]][['rt_end_fitted']]  # Get end position of the fitted region
        # Plot fitted region if it falls within the plot limits
        if (rt_end_fitted <= xmax){
          abline(v=rt_end_fitted, lty='dashed', col=cols_trans[["sample"]])  # Add dashed line for end position of fitted region
          if (rt_end_fitted > rt_end){  # Check if end of fitted region exceeds end of observed region
            lines(c(rt_end, rt_end_fitted), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
          }
        }else{
          if (rt_end < xmax & rt_end_fitted > rt_end){  # Check if observed region ends before plot limit and fitted region exceeds it
            lines(c(rt_end, xmax), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
          }
        }
      }
    }
    # Add legend to the plot
    par(family = 'mono')
    legend('topright', fill=as.character(unlist(cols)), legend=c("ctrl", legend.text), cex=0.75)  # Add legend
    par(family = '')
    # Add horizontal line at y=0
    abline(h=0, col='gray')
    
    # plot 2
    max.val = round(max(subtr_log2_binned.max), 1) + 1
    min.val = round(min(subtr_log2_binned.min), 1) - 1
    plot(0,0,type='n', xlim=c(xmin, xmax), ylim=c(min.val, max.val), main=paste(GOI, '(normalized genebody - ctrl subtracted)'), xlab='position', ylab='norm. log2(cov)', las=1, xaxt='n')
    axis(side=1, at=at_vector, labels=label_vector)
    lines(binned.xvals, subtr_log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["sample"]])
    polygon(c(binned.xvals, rev(binned.xvals)), c(subtr_log2_binned.max[1:length(binned.xvals)], rev(subtr_log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["sample"]], border = NA)
    abline(v=TSS_row-usExt, lty='dotted', col='black')
    abline(v=TES_row-usExt, lty='dotted', col='black')
    ###
    if (rt_list[["sample"]][['rt']] & rt_list[["sample"]][['rt_max']] >= 1 & rt_list[["sample"]][['rt_int']] >= 0.5){    ## require that the max hmm-state has a value >= 1 (i.e. 2-fold over control)
      rt_start = rt_list[["sample"]][['rt_start']]
      rt_end = rt_list[["sample"]][['rt_end']]
      abline(v=rt_start, lty='dashed', col=cols_trans[["sample"]])
      lines(c(rt_start, rt_end), rep(min.val, 2), col=cols_trans[["sample"]], lwd=5, lend=1)
      if ((rt_list[["sample"]][['dsfit']] | rt_list[["sample"]][['sfit']]) & (rt_list[["sample"]][['rt_sum_fitted']] >= rt_list[["sample"]][['rt_sum']])){
        lines(binned.xvals[binned.xvals >= rt_start], rt_list[["sample"]][['best_fit']][1:length(which(binned.xvals >= rt_start))], col=cols_trans[["sample"]], lwd=4)
        lines(binned.xvals[binned.xvals >= rt_start], rt_list[["sample"]][['best_fit']][1:length(which(binned.xvals >= rt_start))], col='black', lwd=1, lty='dotted')
        #lines(xvals[xvals >= rt_start], rt_list[["sample"]][['best_fit']][xvals[xvals >= rt_start]-rt_start+1], col=cols_trans[["sample"]], lwd=4)
        #lines(xvals[xvals >= rt_start], rt_list[["sample"]][['best_fit']][xvals[xvals >= rt_start]-rt_start+1], col='black', lwd=1, lty='dotted')
        rt_end_fitted = rt_list[["sample"]][['rt_end_fitted']]
        if (rt_end_fitted <= xmax){
          abline(v=rt_end_fitted, lty='dashed', col=cols_trans[["sample"]])
          if (rt_end_fitted > rt_end){
            lines(c(rt_end, rt_end_fitted), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')
          }
        }else{
          if (rt_end < xmax & rt_end_fitted > rt_end){
            lines(c(rt_end, xmax), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')
          }
        }
      }
    }
    par(family = 'mono')
    legend('topright', fill=legend.cols, legend=legend.text, cex=0.75) #, bg=NA
    par(family = '')
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


#!/usr/bin/env Rscript

################# (1) libraries                   #################

if(!require("BiocManager")) {install.packages("BiocManager")}
if(!require("rtracklayer")) {BiocManager::install("rtracklayer")}
if(!require("STAN")) {BiocManager::install("STAN")}
if(!require("limma")) {BiocManager::install("limma")}
if(!require("sicegar")) {BiocManager::install("sicegar")}
if(!require("doParallel")) {install.packages("doParallel")}

library(BiocManager)
library(rtracklayer)
library(STAN)
library(limma)
library(sicegar)
library(doParallel, quietly = TRUE, warn.conflicts = FALSE)

source(func_path)
source(rt_analysis_path)


################# (2) command line arguments      #################

args = commandArgs(trailingOnly=TRUE)

ncores = as.integer(args[1])  ## one less than booked!
ctrl_dir = args[2]         ## default EGFP
sample_dir = args[3]
file_cat_batch = read.csv(args[4])

ctrl = file_cat_batch$file_name[file_cat_batch$cat == "ctrl"]
sample = file_cat_batch$file_name[file_cat_batch$cat == "sample"]
ctrl_batch = file_cat_batch$batch[file_cat_batch$cat == "ctrl"]
sample_batch = file_cat_batch$batch[file_cat_batch$cat == "sample"]
batch = c(ctrl_batch, sample_batch)

GOIs = readLines(args[5])
normalization = args[6]       ## 0, 1, 2
usExt = args[7]
rt_analysis_path = args[8]
func_path = args[9]
output_path = args[10]
annot_path = args[11]

annot = import(annot_path)

if (FALSE){
  ## TESTINGS
  ncores = 1  ## one less than booked!
  ctrl_dir = wd         ## default EGFP
  sample_dir = wd
  ctrl = 'EGFP'
  sample = 'CPS'
  GOIs = readLines("C:/Users/sejeo/Documents/Uni/4. semester/Dataprojekt/GOIs.txt")
  normalization = 1       ## 0, 1, 2
  usExt = 5000
  rt_analysis_path = wd
  func_path = "C:/Users/sejeo/Downloads/50a_readthru_analysis_2024_revisit.R"
  output_path = "C:/Users/sejeo/Downloads/bob.Rdata"
  annot_path = "C:/Users/sejeo/Documents/Uni/4. semester/Dataprojekt/annotation_subset.gtf"
}

################# (4) color schemes               #################
cols = list()
cols[["ctrl"]] = rgb(8/255, 76/255, 97/255, 1)
cols[["sample"]] = rgb(219/255, 58/255, 52/255, 1)

cols_trans = list()
cols_trans[['ctrl']] = rgb(8/255, 76/255, 97/255, 1)
cols_trans[['sample']] = rgb(219/255, 58/255, 52/255, 1)



################# (5) general info table          #################


################# (6) run the script              #################
#use 36 cores and 384 GB 
#ncores = 35

if (ncores > 1){
  registerDoParallel(cores = ncores)
  readthrough.list = foreach(GOI = GOIs) %dopar% tryCatch({rt_analysis(ctrl_dir, sample_dir, annot, ctrl, sample, GOI, batch, normalization, statistic='median', usExt, plot=FALSE, verbose=FALSE, verbose2=TRUE)}, error = function(e) {return('error')})
}else{
  readthrough.list = lapply(GOIs, function(GOI) tryCatch({rt_analysis(ctrl_dir, sample_dir, annot, ctrl, sample, GOI, batch, normalization, statistic='median', usExt, plot=FALSE, verbose=FALSE, verbose2=TRUE)}, error = function(e) {return('error')}) )
}
names(readthrough.list) = GOIs

# save list of outputs
save(readthrough.list, file = output_path)



save_vector_to_csv = function(file_path, new_vector) {
  
  if (file.exists(file_path)) {
    existing_data = read.csv(file_path)
    new_data = rbind(existing_data, new_vector)
    write.csv(new_data, file = file_path, row.names = FALSE)
    
  } else {
    write.csv(new_vector, file = file_path, row.names = FALSE)
  }
}
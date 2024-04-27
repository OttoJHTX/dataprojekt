GOI = 'GAPDH'

# Directory specification
ctrl_dir = "/Users/alexanderengell-hansen/Desktop/git/dataprojekt/data/"
sample_dir = "/Users/alexanderengell-hansen/Desktop/git/dataprojekt/data/"
ctrl = c("control")
sample = c("cps")

strand_options = c("+"="plus", "-"="minus")


# Whatever, temp
annotation = annotation_path

annot_gr = import(annotation_path)


# Subsetting annotation on gene name
gene_annot = subset(annot_gr, gene_name == GOI)
chrom_no = as.character(seqnames(gene_annot)@values)
strand_sign = as.character(strand(gene_annot)@values)

# Coordinates
start_coord = min(start(gene_annot))
end_coord = max(end(gene_annot))

# Empty DF with correct length
indexes <- seq(start_coord, end_coord)
ctrl_df <- data.frame(index = indexes)

# Loop for loading replicates and saving to DF
for (fname in ctrl) {
  ctrl_fname = paste0(ctrl_dir, fname, "_", strand_options[strand_sign], ".bw")
  ctrl_bw = import(ctrl_fname, 
                  which = GenomicRanges::GRanges(seqnames = chrom_no, 
                                                 ranges = IRanges::IRanges(start = start_coord, end = end_coord)), 
                  as = "NumericList")[[1]]
  ctrl_df[[fname]] = ctrl_bw
}

# Setting the DF indices to the positions of the reads and deleting index column
rownames(ctrl_df) = ctrl_df$index
ctrl_df$index = NULL

next_gene_annot = subset(annot_gr, seqname = chrom_no, stand = strand_sign)
if (strand_sign == "+") {
  ss = sort(start(next_gene_annot))
  next_gene_coord = which(ss > end_coord)[1]
} else {
  ss = sort(start(next_gene_annot))
  next_gene_coord = which(ss < end_coord)[1]
}




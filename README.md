# Dataprojekt

-.. . - / . .-. / ..-. . -.. -

... -.- -.-- -.. / -- .. --.


- 




To do:

-- Fremlæggelse 19/02/2024
- Biologisk process
- Data
- Model

Spørgsmål til Søren:
 - Hvad bruges bigwig filen til?; Når vi har ranges og kromosomer i den anden fil, skal vi så kun bruge bigwig'en til scoren?
 - Hvad bruges Scoren til?

<br>
GenomeDK paths: <br>
/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw <br>
/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf <br>

Når R er startet: <br>
annot_gr = import('/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf') <br>
bw = import('/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw') <br>

annot_gr <- subset(annot_gr, gene_name == 'GAPDH')

ts_annot_gr <- subset(annot_gr, type == "transcript") <br>
bw_12 <- subset(bw, seqnames == "chr12") <br>

start_coords <- start(ts_annot_gr) <br>
end_coords <- end(ts_annot_gr) <br>
new_start <- min(start_coords) - 500 <br>
new_end <- max(end_coords) + 100000 <br> 
 <br> 
new_ranges <- GRanges(
  seqnames = "chr12",
  ranges = IRanges(start = new_start, end = new_end)
) <br>
<br>
overlaps <- findOverlaps(new_ranges, bw_12) <br>


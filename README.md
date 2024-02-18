# dataprojekt

20240209 - møde:
- Kasper arbejder torsdag fast. 
- Otto arbejder mandag og torsdag fast.
- Alexander arbejder torsdag fast.
- Som udgangspunkt, mødes vi (med Søren) efter databaser TØ, kl. 11. 
- Vi kan dog også mødes om tirsdagen efter kl. 14. 



To do:

-- Få adgang til Genome cluster:
- Otto - done
- Alex - done
- Kasper - done
-- Connected til Genome cluster:
- Otto - done
- Alex - done
- Kasper - done
-- Setup af conda env og R på cluster:
- Otto - done
- Alex - done
- Kasper - WIP

-- Fremlæggelse 19/02/2024
- Biologisk process
- Data
- Model
<br>
GenomeDK paths: <br>
/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw <br>
/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf <br>

Når R er startet: <br>
annot_gr = import('/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf') <br>
bw = import('/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw') <br>
start_coords <- start(ts_annot_gr) <br>
end_coords <- end(ts_annot_gr) <br>
new_start <- min(start_coords) - 500 <br>
new_end <- max(end_coords) + 100000 <br> 
 <br> 
new_ranges <- GRanges( <br> 
  seqnames = seqnames(bw_12), <br> 
  ranges = IRanges(start = new_start, end = new_end), <br>
  strand = "*" <br>
) <br>
<br>
overlaps <- findOverlaps(new_ranges, bw_12) <br>


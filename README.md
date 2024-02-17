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

GenomeDK paths:
/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw
/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf

Når R er startet:
annot_gr = import('/home/engellalex28/KAOs_project/annotation/Homo_sapiens.GRCh38.108.gtf')
bw = import('/home/engellalex28/KAOs_project/data/L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw')
start_coords <- start(ts_annot_gr)
end_coords <- end(ts_annot_gr)
new_start <- min(start_coords) - 500
new_end <- max(end_coords) + 100000

new_ranges <- GRanges(
  seqnames = seqnames(ts_annot_gr),
  ranges = IRanges(start = new_start, end = new_end),
  strand = strand(ts_annot_gr)
)

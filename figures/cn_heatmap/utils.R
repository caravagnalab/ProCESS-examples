library(dplyr)
library(tidyr)
library(IRanges)
library(purrr)

segment_fixed_windows <- function(x, window_size = 1e6) {
  
  # list of chromosomes present
  chrs <- sort(unique(x$chr))
  
  # process per chromosome to avoid cross-chr mixing
  out_chr <- map_df(chrs, function(chr_i) {
    
    # chromosome-specific segments
    segs_chr <- x %>% filter(chr == chr_i)
    
    # if no segments for this chr, skip
    if (nrow(segs_chr) == 0) return(tibble())
    
    # build windows for this chromosome using max end observed
    chr_max <- max(segs_chr$to)
    wins_from <- seq(1, chr_max, by = window_size)
    wins_to   <- pmin(wins_from + window_size - 1, chr_max)
    wins <- tibble(chr = chr_i, from = wins_from, to = wins_to)
    
    # prebuild IRanges for windows and segments (per-chr)
    wins_ir <- IRanges(start = wins$from, end = wins$to)
    
    # iterate samples for this chromosome
    samples <- unique(segs_chr$sample_id)
    map_df(samples, function(sid) {
      segs_s <- segs_chr %>% filter(sample_id == sid)
      
      # if sample has no segments on this chr, return rows with NAs
      if (nrow(segs_s) == 0) {
        tibble(
          chr = wins$chr,
          from = wins$from,
          to = wins$to,
          sample_id = sid,
          Major = NA_character_,
          ratio = NA_real_,
          segment_id = NA_character_
        )
      } else {
        seg_ir <- IRanges(start = segs_s$from, end = segs_s$to)
        
        # find overlaps between each window and sample segments
        hits <- findOverlaps(wins_ir, seg_ir)
        
        # prepare an empty result and fill it
        res <- tibble(
          chr = wins$chr,
          from = wins$from,
          to = wins$to,
          sample_id = sid,
          Major = NA,
          minor = NA,
          ratio = NA,
          segment_id = NA
        )
        
        if (length(hits) > 0) {
          # For each window, determine overlapping segments and choose best by max overlap
          # hits are pairs (query=window index, subject=segment index)
          hits_df <- as_tibble(as.data.frame(hits)) %>%
            dplyr::rename(win_idx = queryHits, seg_idx = subjectHits)
          
          # compute overlap widths for each hit
          hits_df <- hits_df %>%
            mutate(
              win_start = wins$from[win_idx],
              win_end   = wins$to[win_idx],
              seg_start = segs_s$from[seg_idx],
              seg_end   = segs_s$to[seg_idx],
              ov_width  = pmax(0, pmin(win_end, seg_end) - pmax(win_start, seg_start) + 1)
            )
          
          # pick, per window, the segment with largest overlap
          best_hits <- hits_df %>%
            group_by(win_idx) %>%
            slice_max(order_by = ov_width, n = 1, with_ties = FALSE) %>%
            ungroup()
          
          # fill results for windows that had hits
          res[best_hits$win_idx, "Major"]      <- segs_s$Major[best_hits$seg_idx]
          res[best_hits$win_idx, "minor"]      <- segs_s$minor[best_hits$seg_idx]
          res[best_hits$win_idx, "ratio"]      <- segs_s$ratio[best_hits$seg_idx]
          res[best_hits$win_idx, "segment_id"] <- segs_s$segment_id[best_hits$seg_idx]
        }
        
        res
      }
    })
  })
  
  out_chr
}

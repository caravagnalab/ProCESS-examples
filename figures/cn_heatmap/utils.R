library(dplyr)
library(tidyr)
library(IRanges)
library(purrr)

library(dplyr)

add_missing_chromosomes <- function(
    sequenza,
    chr_reference,
    value_cols,
    caller_name = "sequenza",
    segment_type = "clonal"
) {
  
  # Identify missing chromosomes
  missing_chr <- chr_reference %>%
    anti_join(
      sequenza %>% distinct(chromosome),
      by = "chromosome"
    )
  
  if (nrow(missing_chr) == 0) {
    return(sequenza)
  }
  
  # Create rows for missing chromosomes
  new_rows <- missing_chr %>%
    transmute(
      chromosome = chromosome,
      from = start,
      to = end,
      caller = caller_name,
      is_match = "match",
      type = segment_type,
      !!!setNames(
        replicate(length(value_cols), 1, simplify = FALSE),
        value_cols
      )
    )
  
  # Add any remaining columns as NA
  missing_cols <- setdiff(names(sequenza), names(new_rows))
  new_rows[missing_cols] <- NA
  
  # Bind and return
  bind_rows(sequenza, new_rows) %>%
    arrange(chromosome, from)
}



segment_fixed_windows <- function(x, window_size = 1e6) {
  
  # list of chromosomes present
  chrs <- sort(unique(x$chromosome))
  
  # process per chromosome to avoid cross-chr mixing
  out_chr <- map_df(chrs, function(chr_i) {
    # chromosome-specific segments
    segs_chr <- x %>% filter(chromosome == chr_i)
    
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
          caller = NA,
          is_match = NA,
          type = NA,
          TRUE_ratio = NA,
          TRUE_Major1 = NA,
          TRUE_minor1 = NA,
          INFERRED_Major1= NA,
          INFERRED_minor1 = NA,
          TRUE_Major2= NA,
          TRUE_minor2= NA,
          INFERRED_Major2 = NA,
          INFERRED_minor2 = NA,
          INFERRED_ratio = NA
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
          caller = NA,
          is_match = NA,
          type = NA,
          TRUE_ratio = NA,
          TRUE_Major1 = NA,
          TRUE_minor1 = NA,
          INFERRED_Major1= NA,
          INFERRED_minor1 = NA,
          TRUE_Major2= NA,
          TRUE_minor2= NA,
          INFERRED_Major2 = NA,
          INFERRED_minor2 = NA,
          INFERRED_ratio = NA
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
          res[best_hits$win_idx, "caller"] <- segs_s$caller[best_hits$seg_idx]
          res[best_hits$win_idx, "is_match"] <- segs_s$is_match[best_hits$seg_idx]
          res[best_hits$win_idx, "type"] <- segs_s$type[best_hits$seg_idx]
          res[best_hits$win_idx, "TRUE_ratio"] <- segs_s$TRUE_ratio[best_hits$seg_idx]
          res[best_hits$win_idx, "TRUE_Major1"] <- segs_s$TRUE_Major1[best_hits$seg_idx]
          res[best_hits$win_idx, "TRUE_minor1"] <- segs_s$TRUE_minor1[best_hits$seg_idx]
          res[best_hits$win_idx, "INFERRED_Major1"] <- segs_s$INFERRED_Major1[best_hits$seg_idx]
          res[best_hits$win_idx, "INFERRED_minor1"] <- segs_s$INFERRED_minor1[best_hits$seg_idx]
          res[best_hits$win_idx, "TRUE_Major2"] <- segs_s$TRUE_Major2[best_hits$seg_idx]
          res[best_hits$win_idx, "TRUE_minor2"] <- segs_s$TRUE_minor2[best_hits$seg_idx]
          res[best_hits$win_idx, "INFERRED_Major2"] <- segs_s$INFERRED_Major2[best_hits$seg_idx]
          res[best_hits$win_idx, "INFERRED_minor2"] <- segs_s$INFERRED_minor2[best_hits$seg_idx]
          res[best_hits$win_idx, "INFERRED_ratio"] <- segs_s$INFERRED_ratio[best_hits$seg_idx]
        }
        
        res
      }
    })
  })
  
  out_chr
}


segment_fixed_windows_all_split <- function(x, window_size = 1e6) {
  chrs <- sort(unique(x$chr))
  
  map_df(chrs, function(chr_i) {
    segs_chr <- x %>% filter(chr == chr_i)
    if (nrow(segs_chr) == 0) return(tibble())
    
    chr_max <- max(segs_chr$to)
    wins_from <- seq(1, chr_max, by = window_size)
    wins_to   <- pmin(wins_from + window_size - 1, chr_max)
    wins <- tibble(chr = chr_i, from = wins_from, to = wins_to)
    
    wins_ir <- IRanges(start = wins$from, end = wins$to)
    samples <- unique(segs_chr$sample_id)
    
    map_df(samples, function(sid) {
      segs_s <- segs_chr %>% filter(sample_id == sid)
      
      # no segmentation rows for this sample on this chr -> return NA rows for all windows
      if (nrow(segs_s) == 0) {
        tibble(
          chr = wins$chr,
          from = wins$from,
          to = wins$to,
          sample_id = sid,
          Major = NA_integer_,
          minor = NA_integer_,
          ratio = NA_real_,
          segment_id = NA_character_
        )
      } else {
        seg_ir <- IRanges(start = segs_s$from, end = segs_s$to)
        
        hits <- findOverlaps(wins_ir, seg_ir)
        # If no overlaps for ANY window, return NA rows for each fixed window
        if (length(hits) == 0) {
          tibble(
            chr = wins$chr,
            from = wins$from,
            to = wins$to,
            sample_id = sid,
            Major = NA_integer_,
            minor = NA_integer_,
            ratio = NA_real_,
            segment_id = NA_character_
          )
        } else {
          # build tibble of all window<->segment overlaps
          hits_df <- as_tibble(as.data.frame(hits)) %>%
            dplyr::rename(win_idx = queryHits, seg_idx = subjectHits) %>%
            mutate(
              # compute the precise intersection range between the window and the segment
              int_range = purrr::map2(win_idx, seg_idx, ~ IRanges::pintersect(wins_ir[.x], seg_ir[.y])),
              int_from = purrr::map_int(int_range, IRanges::start),
              int_to   = purrr::map_int(int_range, IRanges::end),
              chr = chr_i,
              sample_id = sid,
              Major = segs_s$Major[seg_idx],
              minor = segs_s$minor[seg_idx],
              ratio = segs_s$ratio[seg_idx],
              segment_id = segs_s$segment_id[seg_idx]
            ) %>%
            select(chr, from = int_from, to = int_to, sample_id, Major, minor, ratio, segment_id) %>%
            arrange(from)
          
          # There may be fixed windows that had no overlap at all (they won't appear in hits_df).
          # We need to add rows for those windows with NA CN fields.
          overlapped_win_idx <- unique(hits_df %>% mutate(win_from = from, win_to = to) %>% 
                                         mutate(win_label = paste0(from, ":", to)) %>% 
                                         pull(win_label))
          
          # create labels for all fixed windows
          all_win_labels <- paste0(wins$from, ":", wins$to)
          non_overlapping_windows <- which(!(all_win_labels %in% overlapped_win_idx))
          
          if (length(non_overlapping_windows) > 0) {
            na_rows <- tibble(
              chr = wins$chr[non_overlapping_windows],
              from = wins$from[non_overlapping_windows],
              to   = wins$to[non_overlapping_windows],
              sample_id = sid,
              Major = NA_integer_,
              minor = NA_integer_,
              ratio = NA_real_,
              segment_id = NA_character_
            )
            hits_df <- bind_rows(hits_df, na_rows) %>% arrange(from)
          }
          
          # Return the per-sample per-chr non-overlapping subwindows with CN annotations
          hits_df
        }
      }
    })
  })
}


segment_fixed_windows_all <- function(x, window_size = 1e6) {
  
  chrs <- sort(unique(x$chr))
  
  map_df(chrs, function(chr_i) {
    segs_chr <- x %>% filter(chr == chr_i)
    if (nrow(segs_chr) == 0) return(tibble())
    
    chr_max <- max(segs_chr$to)
    wins_from <- seq(1, chr_max, by = window_size)
    wins_to   <- pmin(wins_from + window_size - 1, chr_max)
    wins <- tibble(chr = chr_i, from = wins_from, to = wins_to)
    
    wins_ir <- IRanges(start = wins$from, end = wins$to)
    
    samples <- unique(segs_chr$sample_id)
    
    map_df(samples, function(sid) {
      segs_s <- segs_chr %>% filter(sample_id == sid)
      
      if (nrow(segs_s) == 0) {
        tibble(
          chr = wins$chr,
          from = wins$from,
          to = wins$to,
          sample_id = sid,
          Major = NA_integer_,
          minor = NA_integer_,
          ratio = NA_real_,
          segment_id = NA_character_
        )
      } else {
        seg_ir <- IRanges(start = segs_s$from, end = segs_s$to)
        
        hits <- findOverlaps(wins_ir, seg_ir)
        if (length(hits) == 0) {
          tibble(
            chr = wins$chr,
            from = wins$from,
            to = wins$to,
            sample_id = sid,
            Major = NA_integer_,
            minor = NA_integer_,
            ratio = NA_real_,
            segment_id = NA_character_
          )
        } else {
          # For **all overlapping segments**, expand rows
          hits_df <- as_tibble(as.data.frame(hits)) %>%
            dplyr::rename(win_idx = queryHits, seg_idx = subjectHits) %>%
            mutate(
              win_start = wins$from[win_idx],
              win_end   = wins$to[win_idx],
              seg_start = segs_s$from[seg_idx],
              seg_end   = segs_s$to[seg_idx]
            ) %>%
            mutate(
              chr = chr_i,
              from = win_start,
              to   = win_end,
              sample_id = sid,
              Major = segs_s$Major[seg_idx],
              minor = segs_s$minor[seg_idx],
              ratio = segs_s$ratio[seg_idx],
              segment_id = segs_s$segment_id[seg_idx]
            ) %>%
            select(chr, from, to, sample_id, Major, minor, ratio, segment_id)
          
          hits_df
        }
      }
    })
  })
}



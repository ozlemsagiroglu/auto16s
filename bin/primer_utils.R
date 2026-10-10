# Primer matching shared by primer_detect.R and primer_trim.R
suppressPackageStartupMessages({ library(Biostrings); library(ShortRead) })

MAX_OFFSET <- 12L   # primer may start up to 12 bases into the read (heterogeneity spacers)

# allowed mismatches: 1 for primers < 20 nt, 2 for >= 20 nt
mm_for <- function(p) max(1L, floor(nchar(p) / 10))

# End position of the first primer hit that starts within the first MAX_OFFSET + 1 bases; NA if none.
# IUPAC codes in the primer are honoured (fixed = FALSE).
primer_end <- function(reads, primer, mm = mm_for(primer)) {
  reads <- as(reads, "DNAStringSet")
  out <- rep(NA_integer_, length(reads))
  if (!length(reads)) return(out)
  w   <- width(reads)
  pre <- subseq(reads, 1, pmin(w, MAX_OFFSET + nchar(primer)))
  m   <- vmatchPattern(DNAString(primer), pre, max.mismatch = mm, fixed = FALSE)
  st  <- startIndex(m); en <- endIndex(m)
  ix  <- rep(seq_along(st), lengths(st)); s <- unlist(st); e <- unlist(en)
  if (!length(ix)) return(out)
  ok  <- s >= 1 & e <= width(pre)[ix]       # hit fully inside the read start (mismatch hits can overhang)
  ix <- ix[ok]; s <- s[ok]; e <- e[ok]
  o  <- order(ix, s); first <- !duplicated(ix[o])
  out[ix[o][first]] <- e[o][first]
  out
}

# Trim a chunk of read pairs given the primer pair. Pairs in reverse orientation (R1 starts with the
# reverse primer) are swapped so every output pair is forward/reverse. Pairs without both primers are dropped.
trim_chunk <- function(a, b, fp, rp) {
  a1 <- primer_end(sread(a), fp); b1 <- primer_end(sread(b), rp)
  norm <- !is.na(a1) & !is.na(b1)
  a2 <- primer_end(sread(a), rp); b2 <- primer_end(sread(b), fp)
  sw <- !norm & !is.na(a2) & !is.na(b2)
  F <- append(narrow(a[norm], start = a1[norm] + 1L), narrow(b[sw], start = b2[sw] + 1L))
  R <- append(narrow(b[norm], start = b1[norm] + 1L), narrow(a[sw], start = a2[sw] + 1L))
  ok <- width(F) > 0 & width(R) > 0
  offs <- c(a1[norm] - nchar(fp), b2[sw] - nchar(fp))       # bases before the primer (spacer length, approx.)
  list(F = F[ok], R = R[ok], n_norm = sum(norm), n_swap = sum(sw), offsets = offs)
}

# ---------------------------------------------------------------- 3' end: read-through
# When the amplicon is shorter than the read, the read runs past the end of the amplicon into the
# reverse complement of the opposite primer and then the sequencing adapter. Everything from that point
# on is cut off. The opposite primer works for any library type; the adapter sequences catch reads where
# the primer is damaged or absent (status "none").
ADAPTERS <- c(truseq = "AGATCGGAAGAGC", nextera = "CTGTCTCTTATACACATCT")
MIN_PARTIAL <- 10L   # a partial match at the very end of the read must cover at least 10 bases
MIN_KEEP    <- 20L   # hits closer than this to the read start are ignored (not a read-through)

rc_iupac <- function(p) as.character(reverseComplement(DNAString(p)))

# first position (1-based) where `pattern` starts in each read: full match anywhere (start >= MIN_KEEP),
# or a prefix of the pattern (>= MIN_PARTIAL bases) at the very end of the read. NA if none.
first_hit <- function(reads, pattern) {
  out <- rep(NA_integer_, length(reads))
  if (!length(reads)) return(out)
  w <- width(reads); P <- nchar(pattern)
  m <- vmatchPattern(DNAString(pattern), reads, max.mismatch = mm_for(pattern), fixed = FALSE)
  st <- startIndex(m); ix <- rep(seq_along(st), lengths(st)); s <- unlist(st)
  if (length(ix)) {
    ok <- s >= MIN_KEEP & s + P - 1 <= w[ix]
    ix <- ix[ok]; s <- s[ok]
    if (length(ix)) { o <- order(ix, s); f <- !duplicated(ix[o]); out[ix[o][f]] <- s[o][f] }
  }
  for (L in seq(min(P - 1L, 30L), MIN_PARTIAL)) {          # longest partial first
    todo <- which(is.na(out) & w >= L + MIN_KEEP)
    if (!length(todo)) break
    tail_seq <- subseq(reads[todo], start = w[todo] - L + 1L)
    hit <- vcountPattern(DNAString(substr(pattern, 1, L)), tail_seq, max.mismatch = floor(L / 10), fixed = FALSE) > 0
    out[todo[hit]] <- w[todo[hit]] - L + 1L
  }
  out
}

# Cut reads (ShortReadQ) before the opposite primer (rc) or an adapter, whichever comes first.
# Returns the trimmed reads and which reads were cut.
trim_readthrough <- function(x, opposite_primer = NULL) {
  r <- sread(x)
  pats <- c(if (!is.null(opposite_primer)) rc_iupac(opposite_primer), ADAPTERS)
  hits <- vapply(pats, function(p) first_hit(r, p), integer(length(r)))
  if (is.null(dim(hits))) hits <- matrix(hits, nrow = length(r))
  cut <- suppressWarnings(apply(hits, 1, function(h) if (all(is.na(h))) NA_integer_ else min(h, na.rm = TRUE)))
  cut_rows <- which(!is.na(cut))
  if (length(cut_rows)) x[cut_rows] <- narrow(x[cut_rows], start = 1L, end = cut[cut_rows] - 1L)
  list(reads = x, cut = !is.na(cut))
}

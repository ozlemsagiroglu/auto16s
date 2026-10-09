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

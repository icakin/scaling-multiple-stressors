# ======================================================================
# fill_runs_from_delivery.R
#
# Fill the FASTQ file names and MD5 checksums into runs_paired_fastq.tsv
# from a Novogene delivery folder, and report anything that does not match.
#
#   Rscript submission/ena/fill_runs_from_delivery.R <delivery_dir> [PRJEBxxxxx] [instrument]
#
# The delivery folder is searched recursively for *.fq.gz / *.fastq.gz and
# for an MD5 file (MD5.txt, md5.txt, *.md5). Each sample is matched by
# either its SampleID (e.g. M1.20.C1) or the sequencing-provider label
# (the numeric "Community" column, e.g. 1) appearing in the file name.
# Nothing is uploaded; this only writes the sheet.
# ======================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: fill_runs_from_delivery.R <delivery_dir> [PRJEB] [instrument_model]")
dir  <- args[1]; prj <- if (length(args) >= 2) args[2] else "PRJEB_TO_FILL"
inst <- if (length(args) >= 3) args[3] else "INSTRUMENT_TO_FILL"

sheet <- "submission/ena/runs_paired_fastq.tsv"
hdr   <- readLines(sheet, n = 1)
runs  <- read.delim(sheet, skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
meta  <- read.csv("data/Sequencing_metadata2.csv", stringsAsFactors = FALSE)
label <- setNames(meta$Community, meta$SampleID)

fq <- list.files(dir, pattern = "\\.(fq|fastq)\\.gz$", recursive = TRUE, full.names = TRUE)
cat("FASTQ files found:", length(fq), "\n")
md5f <- list.files(dir, pattern = "(?i)md5", recursive = TRUE, full.names = TRUE)
md5 <- character(0)
for (f in md5f) {
  l <- readLines(f, warn = FALSE); l <- l[nzchar(l)]
  p <- strsplit(trimws(l), "\\s+")
  ok <- vapply(p, function(x) length(x) >= 2 && grepl("^[0-9a-f]{32}$", x[1]), logical(1))
  md5 <- c(md5, setNames(vapply(p[ok], `[`, "", 1), basename(vapply(p[ok], function(x) x[length(x)], ""))))
}
cat("MD5 entries found:", length(md5), "\n")
if (!length(md5)) message("No MD5 file found; computing checksums (slow for large files) ...")

find_pair <- function(sid) {
  keys <- unique(c(sid, gsub("\\.", "_", sid), gsub("\\.", "-", sid), label[[sid]]))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  hits <- fq[basename(fq) %in% unlist(lapply(keys, function(k)
    basename(fq)[grepl(paste0("(^|[^A-Za-z0-9])", gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", k), "([^A-Za-z0-9]|$)"), basename(fq))]))]
  f1 <- hits[grepl("_1\\.|_R1[_.]|\\.R1\\.|_1\\.f", basename(hits))]
  f2 <- hits[grepl("_2\\.|_R2[_.]|\\.R2\\.|_2\\.f", basename(hits))]
  if (length(f1) == 1 && length(f2) == 1) c(f1, f2) else NULL
}
sum_md5 <- function(f) { b <- basename(f); if (b %in% names(md5)) md5[[b]] else as.character(tools::md5sum(f)) }

unmatched <- character(0)
for (i in seq_len(nrow(runs))) {
  sid <- runs$sample[i]; pr <- find_pair(sid)
  if (is.null(pr)) { unmatched <- c(unmatched, sid); next }
  runs$forward_file_name[i] <- basename(pr[1]); runs$forward_file_md5[i] <- sum_md5(pr[1])
  runs$reverse_file_name[i] <- basename(pr[2]); runs$reverse_file_md5[i] <- sum_md5(pr[2])
}
runs$study <- prj; runs$instrument_model <- inst
writeLines(hdr, sheet)
suppressWarnings(write.table(runs, sheet, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))

cat("\nmatched:", nrow(runs) - length(unmatched), "of", nrow(runs), "\n")
if (length(unmatched)) {
  cat("UNMATCHED (fix by hand or check the delivery folder):\n"); print(unmatched)
}
paper_files <- unique(c(runs$forward_file_name, runs$reverse_file_name)); paper_files <- paper_files[nzchar(paper_files)]
writeLines(fq[basename(fq) %in% paper_files], "submission/ena/files_to_upload.txt")
cat("wrote submission/ena/files_to_upload.txt (", length(paper_files), "files for Webin upload )\n")

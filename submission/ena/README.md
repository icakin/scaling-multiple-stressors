# Raw-read deposit (ENA) and archive (Zenodo)

Everything here is prepared; the two accounts and the two upload steps are yours.

## A. ENA (raw 16S amplicon reads)

What goes up: the 172 community samples used in the paper plus the 3 extraction
blanks. The 88 evolved-line samples and the 13 monoculture samples on the same
sequencing run are NOT included (different study; confirm with Gab).

0. The Novogene delivery is archived by Hebe at https://doi.org/10.5281/zenodo.8289513
   (open access). Download two files from that page into submission/ena/raw/
   (git-ignored): 00.RawData.zip (7.7 GB) and report_X204SC23022054-Z01-F001.zip
   (26 MB, has the instrument model and the MD5 list). Unzip both there.

1. Point the filler at the unzipped raw data. It matches each sample to its
   pair of FASTQ files by SampleID or by the numeric tube label, pulls the
   MD5s from Novogene's MD5 file, and lists the files to upload:

       Rscript submission/ena/fill_runs_from_delivery.R submission/ena/raw

   It prints any sample it could not match. Fix those by hand in
   runs_paired_fastq.tsv before going on. Samples with "E" in the name are the
   evolved-I20 communities (Hebe's note) and are not in the sheet.

2. Log in at https://www.ebi.ac.uk/ena/submit/webin (create a Webin account
   if you do not have one; use the group/lab account if Gab has one).

3. Register the study: "Register study" -> fill from study.txt. Note the
   PRJEB accession it gives you.

4. Register samples: "Register samples" -> "Download spreadsheet template" is
   not needed; choose "Upload filled spreadsheet" and give it
   samples_ERC000011.tsv. It reports errors per row; the usual one is the
   collection date, which is "not provided" here. If you know the month the
   communities were harvested, put it in as YYYY-MM for every row first.

5. Upload the FASTQ files to your Webin file area. Simplest is the Webin
   File Uploader (Java) from the same portal; FTP to webin2.ebi.ac.uk with
   your Webin credentials also works. Upload only the files listed in
   files_to_upload.txt (the filler writes it).

6. Submit reads: "Submit reads" -> "Upload filled spreadsheet" ->
   runs_paired_fastq.tsv, after re-running step 1 with the study accession
   and instrument model:

       Rscript submission/ena/fill_runs_from_delivery.R /path/to/raw_data PRJEBxxxxx "Illumina NovaSeq 6000"

   The instrument model is in Novogene's delivery report. Webin checks every
   MD5 against the uploaded file, so a mismatch means a bad upload, not a bad
   sheet.

7. Set the study release date to hold until publication. It appears as
   PRJEBxxxxx immediately and the reads stay private until release.

Put PRJEBxxxxx in the Data availability placeholder in manuscript.qmd.

## B. Zenodo (code, processed data, model fits: one DOI)

The repo carries .zenodo.json (authors, funding NE/Y000889/1, licence, keywords,
link to the Carmichael et al. deposit) and CITATION.cff, so the archived record
is filled in automatically.

1. Gab enables the GitHub integration on https://zenodo.org/account/settings/github/
   for GabYvonDurocher/scaling-multiple-stressors (it must be the repo owner).
2. On GitHub, create a release, e.g. tag v1.0.0, title "Submission to Nature
   Communications". Zenodo archives the release and mints a DOI within minutes.
3. Paste the DOI into the Code availability placeholder in manuscript.qmd, and
   add it to .zenodo.json related_identifiers as "isVersionOf" for later
   releases if you want.

Every later release (revision, acceptance) gets its own version DOI under the
same concept DOI, so cite the concept DOI in the paper.

If Gab prefers a manual upload instead: zip the repo at the release commit
(without .git) and upload at https://zenodo.org/uploads/new; the metadata from
.zenodo.json has to be typed in by hand in that case.

# Raw-read deposit at the European Nucleotide Archive

This folder records the deposit of the raw 16S rRNA gene amplicon reads that underlie the amplicon sequence variant tables used in the paper.

## Deposit

- **ENA study accession:** PRJEB127267.
- **Samples deposited:** the 172 endpoint community samples analysed in the paper, plus 3 extraction blanks (175 samples in total).
- **Data:** raw paired-end FASTQ reads (Illumina NovaSeq 6000, PE250), one forward and one reverse file per sample (350 files), each registered with its MD5 checksum. Reads are as delivered by the sequencing provider, with primers not removed.
- **Excluded:** the 88 evolved-line samples and 13 monoculture samples sequenced on the same run belong to a separate study and were not deposited here.
- **Release:** the study is held private and will be released on publication.
- **Provenance:** the reads come from the Novogene sequencing delivery archived at <https://doi.org/10.5281/zenodo.8289513>.

## Files

| File | Contents |
|---|---|
| `study.txt` | Study (project) registration: title, description, study type and centre. |
| `samples_ERC000011.tsv` | Sample registration sheet under the ENA default checklist ERC000011: one row per sample, with community identity, richness, stress regime, temperature, pH, salinity and DNA concentration. |
| `runs_paired_fastq.tsv` | Run registration sheet: one row per sample, with library details and the forward and reverse FASTQ file names and MD5 checksums. |
| `files_to_upload.txt` | Manifest of the 350 FASTQ files uploaded to ENA. |
| `fill_runs_from_delivery.py` | Script that matched each sample to its FASTQ pair in the sequencing delivery, took the MD5 checksums from the delivery's checksum file, and wrote the run sheet and the manifest. |

The raw delivery itself is not part of this repository; it is available from the Zenodo record above.

## Code and processed data

This repository, including the processed data, analysis code and model fits, is archived on Zenodo; the DOI is given in the paper's Code availability statement.

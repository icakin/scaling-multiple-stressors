#!/usr/bin/env python3
"""
Fill FASTQ file names and MD5 checksums into runs_paired_fastq.tsv from the
Novogene delivery folder, and list the files to upload.

    python3 submission/ena/fill_runs_from_delivery.py <raw_data_dir> [PRJEBxxxxx] [instrument_model]

Novogene's 16S delivery gives each sample a folder holding five files:
<sid>.raw_1/raw_2.fastq.gz (demultiplexed reads as sequenced, primers on),
<sid>_1/_2.fastq.gz (primer- and adapter-trimmed) and <sid>.extendedFrags
(FLASH-merged). ENA gets the raw pair. MD5s come from Rawdata_MD5.txt.
Nothing is uploaded; this only writes the sheet and files_to_upload.txt.
"""
import csv, os, re, sys, hashlib

if len(sys.argv) < 2:
    sys.exit(__doc__)
raw = sys.argv[1]
prj = sys.argv[2] if len(sys.argv) > 2 else "PRJEB_TO_FILL"
inst = sys.argv[3] if len(sys.argv) > 3 else "INSTRUMENT_TO_FILL"
sheet = "submission/ena/runs_paired_fastq.tsv"

# index every fastq under raw, and every md5 entry
fq = {}
for d, _, files in os.walk(raw):
    for f in files:
        if f.endswith((".fastq.gz", ".fq.gz")):
            fq[f] = os.path.join(d, f)
md5 = {}
for d, _, files in os.walk(raw):
    for f in files:
        if "md5" in f.lower():
            for line in open(os.path.join(d, f)):
                p = line.split()
                if len(p) >= 2 and re.fullmatch(r"[0-9a-f]{32}", p[0]):
                    md5[os.path.basename(p[-1])] = p[0]
print(f"FASTQ files found: {len(fq)}   MD5 entries: {len(md5)}")

def checksum(path):
    b = os.path.basename(path)
    if b in md5:
        return md5[b]
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

meta = {r["SampleID"]: r for r in csv.DictReader(open("data/Sequencing_metadata2.csv"))}

def pair_for(sid):
    cands = [sid, sid.replace(".", "_"), sid.replace(".", "-"), meta.get(sid, {}).get("Community", "")]
    for k in [c for c in cands if c]:
        for f1, f2 in ((f"{k}.raw_1.fastq.gz", f"{k}.raw_2.fastq.gz"),
                       (f"{k}.raw_1.fq.gz",    f"{k}.raw_2.fq.gz"),
                       (f"{k}_1.fastq.gz",     f"{k}_2.fastq.gz"),
                       (f"{k}_1.fq.gz",        f"{k}_2.fq.gz")):
            if f1 in fq and f2 in fq:
                return fq[f1], fq[f2]
    return None

lines = open(sheet).read().splitlines()
head, cols, rows = lines[0], lines[1].split("\t"), [l.split("\t") for l in lines[2:] if l.strip()]
ix = {c: i for i, c in enumerate(cols)}
unmatched, upload = [], []
for r in rows:
    while len(r) < len(cols): r.append("")
    pr = pair_for(r[ix["sample"]])
    if not pr:
        unmatched.append(r[ix["sample"]]); continue
    r[ix["forward_file_name"]] = os.path.basename(pr[0]); r[ix["forward_file_md5"]] = checksum(pr[0])
    r[ix["reverse_file_name"]] = os.path.basename(pr[1]); r[ix["reverse_file_md5"]] = checksum(pr[1])
    r[ix["study"]] = prj; r[ix["instrument_model"]] = inst
    upload += list(pr)

with open(sheet, "w") as f:
    f.write(head + "\n" + "\t".join(cols) + "\n")
    for r in rows: f.write("\t".join(r) + "\n")
open("submission/ena/files_to_upload.txt", "w").write("\n".join(upload) + "\n")

print(f"matched: {len(rows) - len(unmatched)} of {len(rows)}")
if unmatched:
    print("UNMATCHED (fix by hand or check the delivery folder):"); print("\n".join("  " + u for u in unmatched))
print(f"wrote submission/ena/files_to_upload.txt ({len(upload)} files for Webin upload)")

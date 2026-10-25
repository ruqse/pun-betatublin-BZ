#!/usr/bin/env python3
"""
Compare a pipeline run with the published intermediate files.

Usage:
    python3 pipeline/tests/compare_to_published.py <pipeline_outdir>

Checks, per sample or gene, against the published files in data/ (prints PASS/FAIL,
exits 1 on any FAIL):
  1. haplotypes   headers and sequences == data/haplotypes/Pun-bt-{gene}_{farm}_haplotypes.fasta
  2. matrices     identity_matrices/{gene}_matrix.txt == data/haplotypes/{gene}_matrix.txt
                  (same names, identities within 1e-6)
  3. cds          CDS == data/haplotype_cds/ (nucleotides), and translated CDS == ungapped
                  protein in data/cds_alignments/Pun-bt-{gene}-align.fas (curated alignment)
  4. vcf          CHROM POS REF ALT == data/variants/vcf/
  5. chimera_qc   number of reads flagged chimeric by vsearch
Standard library only.
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
GENES = "ABCDEFG"
FARMS = ("Resistant", "Susceptible")
CODON = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
         for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}
failures = []


def report(check, item, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}\t{check}\t{item}\t{detail}")
    if not ok:
        failures.append((check, item))


def read_fasta(path):
    records, name = {}, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                name = line[1:]
                records[name] = []
            elif name is not None:
                records[name].append(line.strip())
    return {k: "".join(v) for k, v in records.items()}


def cluster_key(header):
    return re.search(r"cluster-(\d+)_ReadCount-(\d+)", header.replace(" ", "_")).groups()


def check_haplotypes(outdir):
    for gene in GENES:
        for farm in FARMS:
            path = outdir / "haplotypes" / f"Pun-bt-{gene}_{farm}.consensus_oriented_sequences.fasta"
            if not path.exists():
                continue
            got = read_fasta(path)
            want = read_fasta(DATA / "haplotypes" / f"Pun-bt-{gene}_{farm}_haplotypes.fasta")
            ok = {k: v.upper() for k, v in got.items()} == {k: v.upper() for k, v in want.items()}
            report("haplotypes", f"Pun-bt-{gene}_{farm}", ok, f"{len(got)} vs {len(want)} published")


def read_matrix(path):
    rows = [l.split() for l in open(path).read().splitlines()[1:] if l.strip()]
    return {r[0]: dict(zip([x[0] for x in rows], map(float, r[1:]))) for r in rows}


def check_matrices(outdir):
    for gene in GENES:
        path = outdir / "identity_matrices" / f"{gene}_matrix.txt"
        if not path.exists():
            continue
        got, want = read_matrix(path), read_matrix(DATA / "haplotypes" / f"{gene}_matrix.txt")
        ok = set(got) == set(want) and all(abs(got[a][b] - want[a][b]) < 1e-6 for a in want for b in want)
        report("matrices", f"Pun-bt-{gene}", ok, f"{len(got)} x {len(got)}")


def check_cds(outdir):
    for gene in GENES:
        curated = {("Resistant" in k,) + cluster_key(k): v.replace("-", "").upper()
                   for k, v in read_fasta(DATA / "cds_alignments" / f"Pun-bt-{gene}-align.fas").items()
                   if "cluster-" in k}
        for farm in FARMS:
            path = outdir / "haplotype_cds" / f"Pun-bt-{gene}_{farm}_CDS.fasta"
            if not path.exists():
                continue
            got = read_fasta(path)
            want = read_fasta(DATA / "haplotype_cds" / f"Pun-bt-{gene}_{farm}_haplotype_CDS.fasta")
            report("cds_nt", f"Pun-bt-{gene}_{farm}", {k: v.upper() for k, v in got.items()} == {k: v.upper() for k, v in want.items()},
                   f"{len(got)} vs {len(want)} published")
            same, diff = 0, []
            for name, seq in got.items():
                seq = seq.upper()
                prot = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
                if curated.get((farm == "Resistant",) + cluster_key(name), "").rstrip("*") == prot.rstrip("*"):
                    same += 1
                else:
                    diff.append(name.split()[0])
            # One haplotype (Pun-bt-B Susceptible cluster-24) carries a 1-nt deletion in its
            # pbaa consensus; its frame was restored by hand in the curated alignment.
            known = [d for d in diff if gene == "B" and farm == "Susceptible" and "cluster-24_" in d]
            report("cds", f"Pun-bt-{gene}_{farm}", len(diff) == len(known),
                   f"{same} identical, {len(known)} known curated exception, {len(diff) - len(known)} unexpected")


def vcf_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sorted(tuple(l.split("\t")[i] for i in (0, 1, 3, 4)) for l in fh if not l.startswith("#"))


def check_vcfs(outdir):
    for path in sorted((outdir / "reads_vcf").glob("*.reads.filtered.vcf.gz")):
        got, want = vcf_records(path), vcf_records(DATA / "variants" / "vcf" / path.name)
        report("vcf", path.name.split(".")[0], got == want, f"{len(got)} vs {len(want)} records")


def check_chimeras(outdir):
    for path in sorted((outdir / "chimera_qc").glob("*.uchime.txt")):
        flagged = sum(1 for l in open(path) if l.rstrip("\n").split("\t")[-1] == "Y")
        print(f"INFO\tchimera_qc\t{path.name.split('.')[0]}\t{flagged} reads flagged chimeric")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    args = ap.parse_args()
    check_haplotypes(args.outdir)
    check_matrices(args.outdir)
    check_cds(args.outdir)
    check_vcfs(args.outdir)
    check_chimeras(args.outdir)
    print(f"\n{len(failures)} failure(s)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/awk -f
# Usage: orient_haplotypes.awk records.sam haplotypes.fasta > oriented.fasta
# For each SAM record (coordinate order), print the haplotype; reverse-complement
# it when it maps to the minus strand (flag 16). Mirrors the published job, which
# looped over `samtools view` output of the sorted BAM.
function revcomp(s,   i, c, out) {
    out = ""
    for (i = length(s); i > 0; i--) {
        c = substr(s, i, 1)
        out = out (c in comp ? comp[c] : c)
    }
    return out
}
BEGIN {
    split("A T G C a t g c", f, " "); split("T A C G t a c g", t, " ")
    for (i in f) comp[f[i]] = t[i]
}
FNR == NR { name[++n] = $1; flag[n] = $2; next }
/^>/ { split(substr($0, 2), h, " "); id = h[1]; header[id] = $0; next }
{ seq[id] = seq[id] $0 }
END {
    for (i = 1; i <= n; i++) {
        id = name[i]
        if (!(id in header)) { print "missing haplotype " id > "/dev/stderr"; exit 1 }
        print header[id]
        print (int(flag[i] / 16) % 2 ? revcomp(seq[id]) : seq[id])
    }
}

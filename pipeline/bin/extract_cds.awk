#!/usr/bin/awk -f
# Usage: extract_cds.awk aligned.fasta > cds.fasta
# First record = reference CDS. For every other record, keep the bases at columns
# where the reference has no gap and drop the record's own gaps (as in
# run_build_cluster_cds.sh, the script behind the published CDS).
/^>/ { order[++n] = substr($0, 2); next }
{ seq[n] = seq[n] $0 }
END {
    ref = seq[1]
    for (s = 2; s <= n; s++) {
        out = ""
        for (i = 1; i <= length(ref); i++) {
            if (substr(ref, i, 1) == "-") continue
            c = substr(seq[s], i, 1)
            if (c != "-" && c != "") out = out c
        }
        print ">" order[s]
        print out
    }
}

#!/usr/bin/env python3
"""Group the downloaded GenBank contigs into one FASTA per published GVMAG.

The 74 accessions listed by Booker et al. are contigs. Each carries an `isolate` in its
GenBank title (e.g. "MAG: Nucleocytoviricota sp. isolate vOTU6_CE2_022B_bin14r genomic
sequence") naming the GVMAG it belongs to. This writes one multi-FASTA per isolate, so the
published genomes have the same shape as ours and fastANI compares genome to genome.
"""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--esummary", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--manifest", required=True)
    a = p.parse_args()

    # accession (without version) -> isolate
    d = json.load(open(a.esummary))["result"]
    iso_of, len_of = {}, {}
    for u in d.get("uids", []):
        r = d[u]
        acc = r.get("accessionversion", "")
        base = acc.split(".")[0]
        m = re.search(r"isolate ([^,]+?)(?: genomic sequence| sequence|,|$)", r.get("title", ""))
        iso_of[base] = m.group(1).strip() if m else "UNPARSED"
        len_of[base] = int(r.get("slen", 0))

    # read the FASTA, assigning each record to its isolate
    seqs = defaultdict(list)          # isolate -> [(header, seq)]
    name, buf = None, []
    def flush():
        if name is None:
            return
        acc = name.split()[0].split(".")[0]
        iso = iso_of.get(acc, "UNKNOWN")
        seqs[iso].append((name, "".join(buf)))
    with open(a.fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    flush()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("*.fna"):
        f.unlink()

    with open(a.manifest, "w") as man:
        man.write("booker_gvmag\tbooker_votu\tn_contigs\ttotal_bp\taccessions\tfasta\n")
        for iso in sorted(seqs):
            recs = seqs[iso]
            fn = out / f"{safe(iso)}.fna"
            with open(fn, "w") as fo:
                for h, s in recs:
                    fo.write(f">{h}\n")
                    for i in range(0, len(s), 80):
                        fo.write(s[i:i + 80] + "\n")
            total = sum(len(s) for _, s in recs)
            accs = ",".join(sorted(h.split()[0] for h, _ in recs))
            m = re.match(r"(vOTU\d+)", iso)
            man.write(f"{iso}\t{m.group(1) if m else 'NA'}\t{len(recs)}\t{total}\t{accs}\t{fn}\n")
            print(f"  {iso:34s} {len(recs):3d} contigs  {total:9d} bp -> {fn.name}")

    print(f"\nwrote {len(seqs)} published GVMAGs to {out}")
    unk = [i for i in seqs if i in ("UNKNOWN", "UNPARSED")]
    if unk:
        print(f"WARNING: {sum(len(seqs[i]) for i in unk)} contigs could not be assigned "
              f"to a GVMAG ({unk}); they are written as their own file and should be "
              f"checked before the ANI results are trusted.")


if __name__ == "__main__":
    main()

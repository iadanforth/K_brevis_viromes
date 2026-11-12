#!/usr/bin/env python3
"""Merge metaviralSPAdes and metaSPAdes contigs from one sample.

Two contigs (one from each assembler) are treated as the same sequence when
they share >= MIN_ANI % average nucleotide identity over >= MIN_COV % of the
SHORTER contig (MIUViG-style 95% ANI / 85% AF by default). ANI and coverage
come from `anicalc` (CheckV) run on blastn output.

Clustering is greedy and centroid-based (like CheckV's `aniclust`): contigs are
ranked, and each unassigned contig in rank order becomes a representative and
absorbs every still-unassigned contig from the OTHER assembler that matches it.
Ranking:
    1. metaviralSPAdes circular contigs (likely complete genomes)
    2. longer contigs first
    3. on equal length, metaviralSPAdes before metaSPAdes
Contigs are only compared across assemblers; nothing is deduplicated within
one assembler's output.

Outputs (in OUT_DIR):
    SAMPLE_merged_contigs.fasta   one sequence per cluster (union of both)
    SAMPLE_shared_contigs.fasta   only clusters found by BOTH assemblers
    SAMPLE_clusters.tsv           cluster membership and match statistics
"""
import argparse
import sys


def read_fasta(path):
    seqs = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def topology(name):
    # metaviralSPAdes headers end in _type_circular / _type_linear / _type_linearrepeat
    if "_type_" in name:
        return name.rsplit("_type_", 1)[1]
    return "NA"


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--mv", required=True, help="metaviralSPAdes contigs, headers prefixed MV_")
    p.add_argument("--ms", required=True, help="metaSPAdes contigs, headers prefixed MS_")
    p.add_argument("--ani", required=True, help="anicalc output (MV queries vs MS subjects)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-cov", type=float, default=85.0)
    a = p.parse_args()

    mv = read_fasta(a.mv)
    ms = read_fasta(a.ms)
    seqs = {**mv, **ms}
    if len(seqs) != len(mv) + len(ms):
        sys.exit("Duplicate contig names between assemblies; expected MV_/MS_ prefixes")

    # Cross-assembler matches: edges[contig] -> {other: (ani, cov_of_shorter)}
    edges = {n: {} for n in seqs}
    with open(a.ani) as fh:
        next(fh)  # header: qname tname num_alns pid qcov tcov
        for line in fh:
            q, t, _, pid, qcov, tcov = line.rstrip("\n").split("\t")
            if q not in mv or t not in ms:
                continue
            pid, qcov, tcov = float(pid), float(qcov), float(tcov)
            cov = qcov if len(seqs[q]) <= len(seqs[t]) else tcov
            if pid >= a.min_ani and cov >= a.min_cov:
                prev = edges[q].get(t)
                if prev is None or (pid, cov) > prev:
                    edges[q][t] = edges[t][q] = (pid, cov)

    def rank(n):
        is_mv = n in mv
        return (not (is_mv and topology(n) == "circular"), -len(seqs[n]), not is_mv, n)

    assigned = set()
    clusters = []
    for rep in sorted(seqs, key=rank):
        if rep in assigned:
            continue
        assigned.add(rep)
        members = []
        for other, (pid, cov) in sorted(edges[rep].items(), key=lambda kv: rank(kv[0])):
            if other not in assigned:
                assigned.add(other)
                members.append((other, pid, cov))
        clusters.append((rep, members))

    s = a.sample
    out = f"{a.out_dir}/{s}"
    counts = {"both": 0, "metaviralspades": 0, "metaspades": 0}
    with open(f"{out}_merged_contigs.fasta", "w") as merged, \
         open(f"{out}_shared_contigs.fasta", "w") as shared, \
         open(f"{out}_clusters.tsv", "w") as tsv:
        tsv.write("\t".join([
            "contig_id", "representative", "rep_assembler", "rep_length",
            "mv_topology", "found_by", "n_matched", "matched_contigs",
            "matched_ani", "matched_cov_shorter"]) + "\n")
        for rep, members in clusters:
            rep_is_mv = rep in mv
            if members:
                found_by = "both"
            else:
                found_by = "metaviralspades" if rep_is_mv else "metaspades"
            counts[found_by] += 1

            mv_names = [rep] if rep_is_mv else []
            mv_names += [m for m, _, _ in members if m in mv]
            topos = sorted({topology(n) for n in mv_names})
            mv_topo = ",".join(topos) if topos else "NA"

            cid = f"{s}_{rep}"
            header = (f">{cid} sample={s} found_by={found_by} "
                      f"rep_assembler={'metaviralspades' if rep_is_mv else 'metaspades'} "
                      f"mv_topology={mv_topo} n_matched={len(members)}")
            rec = f"{header}\n{wrap(seqs[rep])}\n"
            merged.write(rec)
            if members:
                shared.write(rec)
            tsv.write("\t".join([
                cid, rep,
                "metaviralspades" if rep_is_mv else "metaspades",
                str(len(seqs[rep])), mv_topo, found_by, str(len(members)),
                ",".join(m for m, _, _ in members) or "NA",
                ",".join(f"{pid:.2f}" for _, pid, _ in members) or "NA",
                ",".join(f"{cov:.2f}" for _, _, cov in members) or "NA",
            ]) + "\n")

    mv_shared = sum(1 for n in mv if edges[n])
    print(f"[{s}] metaviralSPAdes contigs: {len(mv)} "
          f"({mv_shared} matched a metaSPAdes contig)")
    print(f"[{s}] metaSPAdes contigs (after length filter): {len(ms)}")
    print(f"[{s}] merged clusters: {len(clusters)} "
          f"(both={counts['both']}, metaviralspades_only={counts['metaviralspades']}, "
          f"metaspades_only={counts['metaspades']})")


if __name__ == "__main__":
    main()

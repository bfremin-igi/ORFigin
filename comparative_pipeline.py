#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cmd(cmd, log=True):
    """Run a shell command with optional logging."""
    if log:
        print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def combine_and_dedup(genomes_list, combined_fasta, nodup_fasta):
    """Combine FASTA files and remove duplicates using seqkit."""
    with open(combined_fasta, "w") as out:
        for fasta_file in Path(genomes_list).read_text().splitlines():
            with open(fasta_file) as f:
                out.write(f.read())
    # Deduplicate
    run_cmd(f"seqkit rmdup -s {combined_fasta} > {nodup_fasta}")

def run_mmseqs(nodup_fasta, threads, outdir):
    inputDB = Path(outdir, "inputDB")
    clusterDB = Path(outdir, "clusterDB")
    tmp = Path(outdir, "tmp")
    tmp.mkdir(exist_ok=True)

    run_cmd(f"mmseqs createdb {nodup_fasta} {inputDB}")
    run_cmd(f"mmseqs cluster {inputDB} {clusterDB} {tmp} "
            f"--min-seq-id 0.5 -c 0.5 --cov-mode 1 --cluster-mode 2 -e 0.001 "
            f"--threads {threads}")
    tsv = Path(outdir, "cluster_report.tsv")
    run_cmd(f"mmseqs createtsv {inputDB} {inputDB} {clusterDB} {tsv}")
    return tsv

def parse_clusters(tsv_file, min_size=8):
    """Parse MMseqs2 tsv output and filter clusters by minimum size."""
    clusters = {}
    with open(tsv_file) as f:
        for line in f:
            q, t = line.strip().split("\t")[:2]
            clusters.setdefault(q, set()).add(t)
    # Filter by cluster size
    return {cid: members for cid, members in clusters.items() if len(members) >= min_size}

def process_cluster(cluster_id, cluster_members, nodup_fasta, cluster_dir, results_file, threads=1):
    """Process a single cluster: extract sequences, align, run RNAcode, and parse results."""
    cluster_fasta = cluster_dir / f"{cluster_id}.fasta"
    aln_file = cluster_dir / f"{cluster_id}.aln"
    txt_file = cluster_dir / f"{cluster_id}.rnacode.txt"

    tmp_members = cluster_dir / f"{cluster_id}_members.txt"
    with open(tmp_members, "w") as tmp:
        tmp.write("\n".join(cluster_members))

    # seqkit grep: auto-detect sequence type, multi-threaded
    run_cmd(f"seqkit grep -f {tmp_members} {nodup_fasta} -t auto -j {threads} > {cluster_fasta}")

    # mafft alignment: multi-threaded
    run_cmd(f"mafft --clustalout --auto --thread {threads} {cluster_fasta} > {aln_file}")

    # RNAcode: single-threaded
    run_cmd(f"RNAcode {aln_file} > {txt_file}")

    # parse RNAcode results
    with open(txt_file) as f, open(results_file, "a") as out:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 3 and cols[1] == "+1":
                pval = cols[-1]
                try:
                    if float(pval) < 0.05:
                        out.write(f"{cluster_id}\t{line}\n")
                except ValueError:
                    if pval == "<1e-16":
                        out.write(f"{cluster_id}\t{line}\n")

def main(genomes_list, threads=8, outdir="results"):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    cluster_dir = outdir / "clusters"
    cluster_dir.mkdir(exist_ok=True)

    combined_fasta = outdir / "combined.fasta"
    nodup_fasta = outdir / "nodup_combined.fasta"
    results_file = outdir / "rnacode_significant.txt"

    # Clear previous results file
    if results_file.exists():
        results_file.unlink()

    print("🔹 Combining and deduplicating FASTA files...")
    combine_and_dedup(genomes_list, combined_fasta, nodup_fasta)

    print("🔹 Running MMseqs2 clustering...")
    tsv_file = run_mmseqs(nodup_fasta, threads, outdir)

    print("🔹 Parsing clusters...")
    clusters = parse_clusters(tsv_file)
    n_clusters = len(clusters)
    print(f"Found {n_clusters} clusters with ≥8 members.")

    # Determine threads per cluster safely
    threads_per_task = max(1, threads // n_clusters) if n_clusters > 0 else 1

    print("🔹 Processing clusters in parallel...")
    with ThreadPoolExecutor(max_workers=min(threads, n_clusters)) as executor:
        futures = [
            executor.submit(
                process_cluster, cid, members, nodup_fasta, cluster_dir, results_file, threads_per_task
            )
            for cid, members in clusters.items()
        ]
        for f in as_completed(futures):
            f.result()  # Raises exceptions if any

    print(f"\n✅ Pipeline complete. Results saved to {results_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Cluster → Align → RNAcode → Filter pipeline")
    parser.add_argument("genomes_list", help="Text file with paths to FASTA files")
    parser.add_argument("--threads", type=int, default=8, help="Total threads to use")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()
    main(args.genomes_list, args.threads, args.outdir)

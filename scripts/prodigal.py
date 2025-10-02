#!/usr/bin/env python3
import subprocess
from pathlib import Path

def run_prodigal(fasta_file, output_prefix, mode='meta'):
    fasta_file = Path(fasta_file)
    output_prefix = Path(output_prefix)

    # outputs
    ffn_out = output_prefix.with_suffix(".ffn")        # nucleotide sequences
    s_out = output_prefix.with_suffix(".scores")       # start score file

    cmd = [
        "prodigal",
        "-i", str(fasta_file),
        "-d", str(ffn_out),
        "-s", str(s_out),
        "-p", "meta" if mode=='meta' else "single",
        "-q"
    ]

    print(f"Running Prodigal: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"Prodigal completed. Outputs: {ffn_out}, {s_out}")
    return ffn_out, s_out

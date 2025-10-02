#!/usr/bin/env python3
import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import csv
from pathlib import Path
import sys
import subprocess

# import utilities
from utils import extract_strong_starts

# -----------------------------
# Prodigal wrapper
# -----------------------------
def run_prodigal(fasta_file, output_prefix, mode="meta"):
    output_prefix = Path(output_prefix)
    ffn_out = output_prefix.with_suffix(".ffn")
    faa_out = output_prefix.with_suffix(".faa")
    scores_out = output_prefix.with_suffix(".scores")

    cmd = [
        "prodigal",
        "-i", str(fasta_file),
        "-d", str(ffn_out),
        "-a", str(faa_out),
        "-s", str(scores_out),
        "-p", mode,
        "-q"
    ]
    print(f"Running Prodigal: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"Prodigal complete. Outputs: {ffn_out}, {faa_out}, {scores_out}")
    return ffn_out, faa_out, scores_out

# -----------------------------
# DNABERT inference
# -----------------------------
def run_dnabert_inference(ffn_file, output_file, device="cuda"):
    model_dir = Path("model/")  # expects model/ folder with HF model files

    print(f"Loading model from {model_dir} ...")
    tokenizer = AutoTokenizer.from_pretrained(model_dir)
    model = AutoModelForSequenceClassification.from_pretrained(
        model_dir,
        trust_remote_code=True
    )

    # Force PyTorch attention (disable FlashAttention/Triton)
    if hasattr(model.config, "use_flash_attn"):
        model.config.use_flash_attn = False
    for module in model.modules():
        if hasattr(module, "flash_attn_forward"):
            module.flash_attn_forward = None

    model.to(device)
    model.eval()

    # read sequences from FASTA
    from Bio import SeqIO
    fragments = [(record.id, str(record.seq)) for record in SeqIO.parse(ffn_file, "fasta")]

    results = []
    print(f"Running inference on {len(fragments)} sequences ...")
    with torch.no_grad():
        for header, seq in fragments:
            inputs = tokenizer(seq, return_tensors="pt", padding=True, truncation=True).to(device)
            logits = model(**inputs).logits
            prob = torch.softmax(logits, dim=1)[0, 1].item()
            results.append((header, prob))

    # save results
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ORF_ID", "Real_probability"])
        writer.writerows(results)

    print(f"Inference complete. Results saved to {output_file}")

# -----------------------------
# Main CLI
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python main.py prodigal <input.fasta> <output_prefix> [meta|single]")
        print("  python main.py extract <scores_file> <genome.fasta> <output.fasta> [threshold] [length]")
        print("  python main.py inference <input.fasta> <output.csv> [device]")
        sys.exit(1)

    mode = sys.argv[1]

    if mode == "prodigal":
        if len(sys.argv) < 4:
            print("Usage: python main.py prodigal <input.fasta> <output_prefix> [meta|single]")
            sys.exit(1)
        fasta_file = sys.argv[2]
        output_prefix = sys.argv[3]
        prod_mode = sys.argv[4] if len(sys.argv) > 4 else "meta"
        run_prodigal(fasta_file, output_prefix, prod_mode)

    elif mode == "extract":
        if len(sys.argv) < 5:
            print("Usage: python main.py extract <scores_file> <genome.fasta> <output.fasta> [threshold] [length]")
            sys.exit(1)
        scores_file = sys.argv[2]
        genome_fasta = sys.argv[3]
        output_file = sys.argv[4]
        threshold = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
        length = int(sys.argv[6]) if len(sys.argv) > 6 else 150
        # automatically pick up corresponding .faa file from scores_file prefix
        faa_file = Path(scores_file).with_suffix(".faa")
        extract_strong_starts(scores_file, genome_fasta, threshold, length, output_file, faa_file=faa_file)

    elif mode == "inference":
        if len(sys.argv) < 4:
            print("Usage: python main.py inference <input.fasta> <output.csv> [device]")
            sys.exit(1)
        ffn_file = sys.argv[2]
        output_file = sys.argv[3]
        device = sys.argv[4] if len(sys.argv) > 4 else "cuda"
        run_dnabert_inference(ffn_file, output_file, device)

    else:
        print(f"Unknown mode: {mode}")
        sys.exit(1)

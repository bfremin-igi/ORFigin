#!/usr/bin/env python3
import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from utils import extract_top_orf_fragments
import csv
from pathlib import Path

def run_dnabert_inference(ffn_file, output_file, device='cuda'):
    model_dir = Path("model/")

    # Load model & tokenizer
    tokenizer = AutoTokenizer.from_pretrained(model_dir)
    model = AutoModelForSequenceClassification.from_pretrained(model_dir)
    model.to(device)
    model.eval()

    # Extract ORFs
    fragments = extract_top_orf_fragments(ffn_file)

    results = []
    with torch.no_grad():
        for header, seq in fragments:
            inputs = tokenizer(seq, return_tensors="pt", padding=True, truncation=True).to(device)
            logits = model(**inputs).logits
            prob = torch.softmax(logits, dim=1)[0,1].item()
            results.append((header, prob))

    # Save results
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ORF_ID", "Real_probability"])
        writer.writerows(results)

    print(f"Inference complete. Results saved to {output_file}")

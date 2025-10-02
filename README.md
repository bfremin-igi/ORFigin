# 🧬 Comparative ORF Prediction Pipeline

This repository provides a pipeline for **Prodigal-based ORF prediction** combined with a **BERT classifier** for sequence scoring. It includes scripts for inference, utilities, and an end-to-end workflow.

---

## ⚙️ Installation

```bash
# Clone the repository
git clone https://github.com/bfremin-igi/ORFigin.git
cd ORFigin

# Create and activate environment
conda create -n BERT python=3.10 -y
conda activate BERT

# Install dependencies
pip install -r requirements.txt
```


## 🚀 Usage

### Example Pipeline Script

`run_pipeline.sh`

```bash
#!/bin/bash
conda activate BERT

# Step 1: Run Prodigal to generate ORF scores (only requires CPU)
python main.py prodigal $1.fa $1 meta

# Step 2: Extract candidate ORFs and output truncated FASTA (only requires CPU)
python main.py extract $1.scores $1.fa $1_150.fasta 0 150

#Step 3: Runs ORFigin to assign probability of "Real" vs "Not Real" (requires GPU)
python main.py inference $1_150.fasta $1_150out

#Step 4: (Optional) - if you give it a list of 150 bp start of genes from Step 2 (ie. $1_150.fasta), it will cluster, align, and build comparative genomics support (only requires CPU)
python comparative_pipeline.py startslist.txt --threads 16 --outdir my_output



#!/usr/bin/env python3
from Bio import SeqIO
import re
import csv
import pandas as pd
from pathlib import Path


def extract_strong_starts(scores_file, genome_fasta, threshold=0.0, length=150,
                          output_fasta=None, faa_file=None):
    """
    Extract strong start sequences from a Prodigal .scores file, ignoring 'Edge' genes.
    Also skip any starts that Prodigal already predicted (from .faa file).

    Filtering logic:
      - Parse all ORFs from the .faa file and collect (start,end,strand) per contig.
      - When scanning the .scores candidates, skip any candidate whose 'Beg' lies
        within any predicted ORF interval on the same contig (inclusive).
    """

    # Load all contigs (seq records keyed by contig id)
    contigs = {record.id.split()[0]: record.seq for record in SeqIO.parse(genome_fasta, "fasta")}

    # Parse prodigal .faa to collect predicted ORFs per contig
    predicted_orfs = {}  # contig -> list of (start, end, strand_int)
    if faa_file and Path(faa_file).exists():
        for record in SeqIO.parse(faa_file, "fasta"):
            desc = record.description

            # Try a robust regex to capture "# <start> # <end> # <strand> #"
            m = re.search(r'#\s*(\d+)\s*#\s*(\d+)\s*#\s*([+-]?\d+)\s*#', desc)
            if m:
                start = int(m.group(1))
                end = int(m.group(2))
                # strand in .faa is often "-1" or "1"
                try:
                    strand_int = int(m.group(3))
                    if strand_int not in (-1, 1):
                        strand_int = 1 if strand_int > 0 else -1
                except ValueError:
                    strand_int = 1
            else:
                # fallback to splitting on "#"
                parts = desc.split("#")
                if len(parts) >= 4:
                    try:
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand_field = parts[3].strip()
                        try:
                            strand_int = int(strand_field)
                            if strand_int not in (-1, 1):
                                strand_int = 1 if strand_int > 0 else -1
                        except ValueError:
                            # e.g. '+' or '-'
                            strand_int = 1 if strand_field.startswith("+") else -1
                    except Exception:
                        continue
                else:
                    continue

            # record.id is the protein id like "k119_73647_1" -> derive contig base "k119_73647"
            prot_id = record.id.split()[0]
            m2 = re.match(r'^(.+?)_(\d+)$', prot_id)
            contig_base = m2.group(1) if m2 else prot_id

            predicted_orfs.setdefault(contig_base, []).append((start, end, strand_int))

    # Counters / results
    strong_starts = []
    filtered_by_prodigal = 0
    current_contig = None

    # Read scores file and extract candidates
    with open(scores_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# Sequence Data:"):
                m = re.search('seqhdr="([^"]+)"', line)
                if m:
                    # seqhdr may contain additional info after a space; take first token
                    current_contig = m.group(1).split()[0]
                    if current_contig not in contigs:
                        raise ValueError(f"Contig {current_contig} not found in genome FASTA")
                continue
            if line.startswith("#"):
                continue

            cols = line.split()
            if len(cols) < 7 or current_contig is None:
                continue

            try:
                beg = int(cols[0])
                end = int(cols[1])
                strand_field = cols[2]
                strt_score = float(cols[5])
                codon = cols[6]
            except ValueError:
                continue

            # Skip edge cases
            if codon == "Edge":
                continue

            # convert strand to integer oriented like Prodigal: '+' -> 1, '-' -> -1
            if strand_field == "+":
                candidate_strand = 1
            elif strand_field == "-":
                candidate_strand = -1
            else:
                # sometimes it's numeric (e.g. "1" or "-1")
                try:
                    candidate_strand = int(strand_field)
                except ValueError:
                    candidate_strand = 1

            # If this contig has predicted ORFs, skip if the candidate 'beg' lies within any ORF
            skip = False
            for (ps, pe, pstrand) in predicted_orfs.get(current_contig, []):
                # match same strand OR ignore strand? currently require same strand
                if candidate_strand == pstrand and (ps <= beg <= pe):
                    skip = True
                    break
            if skip:
                filtered_by_prodigal += 1
                continue

            if strt_score >= threshold:
                seq_obj = contigs[current_contig]
                if candidate_strand == 1:
                    start_idx = max(0, beg - 1)
                    end_idx = start_idx + length
                    seq = seq_obj[start_idx:end_idx]
                else:
                    # negative strand: take window ending at 'end'
                    end_idx = min(len(seq_obj), end)
                    start_idx = max(0, end_idx - length)
                    seq = seq_obj[start_idx:end_idx].reverse_complement()
                strong_starts.append((f"{current_contig}_start_{beg}_{end}_{'+' if candidate_strand==1 else '-'}", str(seq)))

    if output_fasta:
        with open(output_fasta, "w") as out:
            for header, seq in strong_starts:
                out.write(f">{header}\n{seq}\n")
        print(f"Extracted {len(strong_starts)} strong start sequences to {output_fasta} "
              f"(skipped {filtered_by_prodigal} candidates because Prodigal already predicted an ORF here)")

    return strong_starts


def deduplicate_by_stop(inference_csv, scores_file, output_csv, 
                        weight_start=0.5, weight_inference=0.5):
    """
    Group ORFs by stop codon and keep the best one based on combined scoring.
    
    This handles the case where multiple alternative start codons lead to the same
    stop codon. We want to pick the best one based on both the Prodigal start score
    and the DNABERT inference score.
    
    IMPORTANT: Prodigal coordinates are always in genomic order (lower first):
      - Forward strand (+): beg=start_codon, end=stop_codon
      - Reverse strand (-): beg=stop_codon, end=start_codon
    
    So the actual stop position is:
      - Forward strand: end position
      - Reverse strand: beg position (the lower coordinate)
    
    Args:
        inference_csv: CSV with ORF_ID and Real_probability from DNABERT inference
        scores_file: Original Prodigal scores file (to get start scores)
        output_csv: Path for deduplicated output CSV
        weight_start: Weight for start score in combined metric (0-1)
        weight_inference: Weight for inference score in combined metric (0-1)
    
    Returns:
        DataFrame with deduplicated results
    """
    
    # Load inference results
    df = pd.read_csv(inference_csv)
    
    # Parse ORF_ID to extract contig, start, end, strand
    # Expected format: "contig_start_100_500_+" or "contig_start_100_500_-"
    def parse_orf_id(orf_id):
        match = re.match(r'(.+)_start_(\d+)_(\d+)_([+-])$', orf_id)
        if match:
            return {
                'contig': match.group(1),
                'beg': int(match.group(2)),      # Lower coordinate
                'end': int(match.group(3)),      # Higher coordinate
                'strand': match.group(4)
            }
        return None
    
    df['parsed'] = df['ORF_ID'].apply(parse_orf_id)
    df = df[df['parsed'].notna()].copy()
    
    if len(df) == 0:
        print("Warning: No valid ORF IDs found in inference CSV")
        # Create empty output file
        pd.DataFrame(columns=['ORF_ID', 'Real_probability', 'start_score', 
                              'combined_score', 'num_alternatives']).to_csv(output_csv, index=False)
        return df
    
    df['contig'] = df['parsed'].apply(lambda x: x['contig'])
    df['beg'] = df['parsed'].apply(lambda x: x['beg'])
    df['end'] = df['parsed'].apply(lambda x: x['end'])
    df['strand'] = df['parsed'].apply(lambda x: x['strand'])
    
    # Calculate actual stop position based on strand
    # Forward strand (+): stop is at 'end' (higher coordinate)
    # Reverse strand (-): stop is at 'beg' (lower coordinate)
    df['stop_pos'] = df.apply(
        lambda row: row['end'] if row['strand'] == '+' else row['beg'], 
        axis=1
    )
    
    # Calculate actual start position based on strand
    # Forward strand (+): start is at 'beg' (lower coordinate)
    # Reverse strand (-): start is at 'end' (higher coordinate)
    df['start_pos'] = df.apply(
        lambda row: row['beg'] if row['strand'] == '+' else row['end'], 
        axis=1
    )
    
    # Load start scores from Prodigal scores file
    # Key by (contig, beg) since 'beg' is always the first coordinate in the scores file
    start_scores = {}
    current_contig = None
    
    with open(scores_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Sequence Data:"):
                m = re.search('seqhdr="([^"]+)"', line)
                if m:
                    current_contig = m.group(1).split()[0]
                continue
            if line.startswith("#") or not line:
                continue
            
            cols = line.split()
            if len(cols) >= 6 and current_contig:
                try:
                    beg = int(cols[0])
                    strt_score = float(cols[5])
                    start_scores[(current_contig, beg)] = strt_score
                except (ValueError, IndexError):
                    continue
    
    # Add start scores to dataframe
    df['start_score'] = df.apply(
        lambda row: start_scores.get((row['contig'], row['beg']), 0.0), 
        axis=1
    )
    
    # Normalize scores to 0-1 range
    if df['start_score'].max() > df['start_score'].min():
        df['start_score_norm'] = ((df['start_score'] - df['start_score'].min()) / 
                                  (df['start_score'].max() - df['start_score'].min()))
    else:
        df['start_score_norm'] = 1.0
    
    df['inference_score_norm'] = df['Real_probability']  # Already 0-1
    
    # Combined score
    df['combined_score'] = (weight_start * df['start_score_norm'] + 
                            weight_inference * df['inference_score_norm'])
    
    # Group by (contig, stop_pos, strand) - this correctly identifies shared stop codons
    df['stop_key'] = df.apply(
        lambda row: (row['contig'], row['stop_pos'], row['strand']), 
        axis=1
    )
    
    # Count alternatives per stop
    stop_counts = df.groupby('stop_key').size()
    df['num_alternatives'] = df['stop_key'].map(stop_counts)
    
    # Keep the row with highest combined score per stop position
    best_orfs = df.loc[df.groupby('stop_key')['combined_score'].idxmax()].copy()
    
    # Save results with additional columns for transparency
    output_df = best_orfs[['ORF_ID', 'Real_probability', 'start_score', 
                           'combined_score', 'num_alternatives', 
                           'start_pos', 'stop_pos']].copy()
    output_df.to_csv(output_csv, index=False)
    
    # Summary statistics
    num_duplicates = len(df) - len(best_orfs)
    num_with_alternatives = (best_orfs['num_alternatives'] > 1).sum()
    
    print(f"\nDeduplication complete:")
    print(f"  Input candidates: {len(df)}")
    print(f"  Unique stop codons: {len(best_orfs)}")
    print(f"  Removed duplicates: {num_duplicates}")
    print(f"  Stops with multiple starts: {num_with_alternatives}")
    print(f"  Results saved to {output_csv}")
    
    # Show example of deduplicated group if any exist
    if num_with_alternatives > 0:
        example = df[df['stop_key'] == best_orfs[best_orfs['num_alternatives'] > 1].iloc[0]['stop_key']]
        print(f"\nExample of deduplication:")
        print(f"  Stop position: {example.iloc[0]['stop_pos']} on strand {example.iloc[0]['strand']}")
        print(f"  Number of alternative starts: {len(example)}")
        for _, row in example.iterrows():
            print(f"    Start {row['start_pos']}: score={row['start_score']:.3f}, "
                  f"inference={row['Real_probability']:.3f}, "
                  f"combined={row['combined_score']:.3f} "
                  f"{'← SELECTED' if row['ORF_ID'] == best_orfs[best_orfs['stop_key'] == row['stop_key']].iloc[0]['ORF_ID'] else ''}")
    
    return best_orfs

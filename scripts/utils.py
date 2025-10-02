#!/usr/bin/env python3
from Bio import SeqIO
import re
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

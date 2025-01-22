#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import json
from pathlib import Path

DOCS = """
Transfer annotations between sequence sets.
Two input multi fasta are aligned using a global gapped alignment algorithm.
Queries and targers input are expectec to have 1-to-1 relationship, i.e. the 
same number of sequences in each file.

Queries are processed 1-by-1 to match a query to a target:
  - Each query is aligned vs all targets
  - Both query seq and reverse complement is considered for the alignment
  - The target with the best alignment score is selected as the match
  - This target is no longer included for the next queries

Coordinates mapping from query to target are defined based on the alignment and
written in json format.

Outputs:
  - json files with all mappings <out_dir>/coords.json
    mappings are a list of match with ids a list of matched sequence.
    "ids": [query_id, target_id].
    "coords": coordinates frmo the alignment. 
  - reoriented queries in <out_dir>/<queries.basename>.reoriented.fasta
  - aln files for each match: 
    <out_dir>/<queries.basename>_<query_id>_<targets.basename>_<target.id>.aln

"""

def parse_arguments():
  parser = argparse.ArgumentParser(description=DOCS, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('--queries', required=True, help='Input multifasta file 1')
  parser.add_argument('--targets', required=True, help='Input multifasta file 2')
  parser.add_argument('--min_identity', type=float, default=0.85, help='Minimum identity threshold')
  parser.add_argument('--out_dir', default=".", help='Output directory')
  return parser.parse_args()

def setup_aligner():
  """Configure and return a PairwiseAligner object for local gapped alignment"""
  aligner = PairwiseAligner()
  
  # Set alignment mode to local
  aligner.mode = 'global'
  
  # Scoring parameters
  aligner.match_score = 2.0
  aligner.mismatch_score = -1.0
  
  # Gap parameters
  aligner.open_gap_score = -10.0
  aligner.extend_gap_score = -0.5
  
  # Target only the best alignment
  aligner.target_end_gap_score = 0.0
  aligner.query_end_gap_score = 0.0
  
  return aligner

def get_alignment_stats(alignment):
  """Get alignment statistics using built-in methods"""
  
  # Calculate identity
  matches = alignment.counts().identities
  aligned_length = len(str(alignment.query).replace('-', ''))
  identity = (matches / aligned_length) * 100
  
  return {
      'identity': identity,
      'score': alignment.score,
      'matches': matches,
      'aligned_length': aligned_length
  }

def write_coords(coords_dict, yaml_file):
    """
    Writes coordinate mappings to a YAML file using explicit Python tuples as keys.
    
    Args:
        coords_dict: Dictionary mapping sequence ID pairs to coordinate arrays
        yaml_file: Path to output YAML file
    """
    # Convert numpy arrays to lists for JSON serialization
    output_dict = [
        {
          "ids": key,
          "coords": (arr[0].tolist(), arr[1].tolist())
        }
        for key, arr in coords_dict.items()
    ]
    
    with open(yaml_file, 'w') as f:
      json.dump(output_dict, f)

def write_alignment(alignment, query_id, target_id, stats, fh):
  """Format alignment output similar to BLAST format 0"""
  fh.write(f"Query: {query_id}\n")
  fh.write(f"Subject: {target_id}\n")
  fh.write(f"Score: {stats['score']:.1f}\n")
  fh.write(f"Identity: {stats['identity']:.2f}% ({stats['matches']}/{stats['aligned_length']})\n\n")
  fh.write("## Alignment:\n\n")
  fh.write(str(alignment))
  
def main():
  args = parse_arguments()
  
  # Create output directory
  os.makedirs(args.out_dir, exist_ok=True)
  
  # Read sequences
  queries = list(SeqIO.parse(args.queries, "fasta"))
  targets = list(SeqIO.parse(args.targets, "fasta"))
  
  # Check if both files have same number of records
  if len(queries) != len(targets):
    sys.exit(f"Error: Different number of sequences in input files: {len(queries)} vs {len(targets)}")
  
  # Setup aligner
  aligner = setup_aligner()
  
  # Initialize coordinate mapping dataframe
  coord_maps = {}
  
  # Get basenames for output files
  queries_basename = Path(args.queries).stem
  targets_basename = Path(args.targets).stem
  
  # Process each sequence from seqs1
  available_targets = targets.copy()
  os.makedirs(args.out_dir, exist_ok=True)

  reoriented_queries_file = open(f"{args.out_dir}/{queries_basename}.reoriented.fasta", "w")
  
  for query in queries:
    best_score = float('-inf')
    best_alignment = None
    best_target = None
    best_target_idx = None
    best_stats = None
    
    # Find best matching sequence in seqs2
    for query_seq, is_rev in [ (query.seq, False), (query.seq.reverse_complement(), True) ]:
      for idx, target in enumerate(available_targets):
        if target is None:
          continue
          
        # Perform alignment
        alignments = aligner.align(query_seq, target.seq)
        if alignments:
          # Get the best scoring alignment
          alignment = alignments[0]
          stats = get_alignment_stats(alignment)
          
          if alignment.score > best_score:
            best_score = alignment.score
            best_alignment = alignment
            best_target = target
            best_target_idx = idx
            best_stats = stats
            best_is_rev = is_rev
    
    if best_alignment is None:
      sys.exit(f"Error: No alignment found for sequence {query.id}")
    
    if best_is_rev:
      query.seq = query.seq.reverse_complement()
      query.id = query_seq.id + "_rev"
      SeqIO.write(query, reoriented_queries_file, "fasta")
    
    # Check thresholds
    if (best_stats['identity'] < args.min_identity):
      sys.exit(f"Error: Identity below threshold for {query.id}\n" +
        f"Identity: {best_stats['identity']:.2f}%")
    
    # Write alignment to file
    aln_filename = f"{targets_basename}_{query.id}_{targets_basename}_{best_target.id}.aln"
    aln_path = os.path.join(args.out_dir, aln_filename)
    
    with open(aln_path, 'w') as f:
      write_alignment(
        best_alignment,
        query.id,
        best_target.id,
        best_stats,
        f
      )
    
    # Create coordinate mapping
    coord_maps[(query.id, best_target.id)] = best_alignment.coordinates
    
    # Remove used target from available sequences
    available_targets[best_target_idx] = None
  reoriented_queries_file.close()
    
  # Create and save coordinate mapping DataFrame
  write_coords(coord_maps, f"{args.out_dir}/coords.json")

if __name__ == "__main__":
  main()

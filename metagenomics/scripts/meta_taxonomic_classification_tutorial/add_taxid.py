#!/usr/bin/env python3
"""
This script modifies the headers of each FASTA record by appending the appropriate Kraken-style 
taxonomic ID tag (|kraken:taxid|<taxid>|) for each Sequin ID based on the names.dmp file.
"""
import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Add Kraken2 taxid to FASTA headers"
)
parser.add_argument("--names", type=argparse.FileType("r"), required=True,
                    help="Path to custom names.dmp file")
parser.add_argument("--fasta", type=argparse.FileType("r"), required=True,
                    help="Path to Sequin FASTA file")
args = parser.parse_args()


sequin_to_taxid = {}
for line in args.names:
    if "SQN" not in line:
        continue
    items = [x.strip() for x in line.split("|")]
    if "abundance" in items[1]:
        continue
    sequin_to_taxid[items[1]] = items[0]


for record in SeqIO.parse(args.fasta, "fasta"):
    id = record.id
    
    # Ensure a taxid exists for this sequence
    if id not in sequin_to_taxid:
        sys.stderr.write(f"[WARNING] No taxid found for sequence: {id}\n")
        continue

    tax_id = sequin_to_taxid[record.id]
    record.id = f"{id}|kraken:taxid|{tax_id}|"
    record.name = ""
    record.description = ""
    SeqIO.write(record, sys.stdout, "fasta")

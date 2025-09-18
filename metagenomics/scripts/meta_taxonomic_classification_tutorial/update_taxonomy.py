#!/usr/bin/env python3
"""Create custom `nodes.dmp` and `names.dmp` files for Kraken2 with Sequin entries.

To add Sequins to a Kraken2 database, we need to assign taxonomic IDs from the
NCBI taxonomy database. As Sequins do not have entries in the NCBI database,
this script generates custom taxonomic entries and writes them to `nodes.dmp` and `names.dmp`.
"""
import csv
import os
import sys
import argparse
from typing import Any, TextIO


field_term = "\t|\t"
line_term = "\t|"


class Taxid:
    def __init__(self, start_id: int = 999000001):
        self.current_id = start_id

    def next(self):
        id = self.current_id
        self.current_id += 1
        return str(id)


def print_node_line(
    taxid: str, parent_id: str, rank: str, division_id: str, handle=sys.stdout
) -> None:
    print(
        field_term.join(
            [
                taxid,
                parent_id,
                rank,
                "",
                division_id,
                "0",
                "1",
                "0",
                "0",
                "0",
                "0",
                "0",
                "",
            ]
        )
        + line_term,
        file=handle,
    )


def print_nodes(entries: list[dict[str, Any]], handle=sys.stdout, sequin_parent_node_id="999000000") -> None:
    # This creates a "sequins" node in the taxonomy whos parent is "artificial sequences" (Taxonomy ID: 81077).
    print_node_line(sequin_parent_node_id, "81077", "no rank", "8", handle)
    # Create the sequin nodes. Each sequin is represented by a species level
    # node. For each abundance level in the ladder a genus level node is
    # created. Each sequin has the abundance/genus node as its parent.
    for entry in entries:
        taxid = entry["abundance"][1]
        parent_id = taxid
        print_node_line(taxid, sequin_parent_node_id, "genus", "7", handle)
        for sequin, taxid in entry["sequins"]:
            print_node_line(taxid, parent_id, "species", "7", handle)


def print_name_line(
    taxid: str, name_txt: str, name_class: str, handle=sys.stdout
) -> None:
    print(field_term.join([taxid, name_txt, "", name_class]) + line_term, file=handle)


def print_names(entries: list[dict[str, Any]], handle=sys.stdout, sequin_parent_node_id="999000000") -> None:
    print_name_line(sequin_parent_node_id, "sequin", "scientific name", handle)
    for entry in entries:
        abundance, taxid = entry["abundance"]
        print_name_line(taxid, f"SQN abundance {abundance}", "scientific name", handle)
        for sequin, taxid in entry["sequins"]:
            print_name_line(taxid, sequin, "scientific name", handle)


def read_manifest(abundance_handle: TextIO, taxid_generator: Taxid) -> list[dict[str, Any]]:
    xs: dict[float, list[str]] = {}
    abx = set()
    reader = csv.DictReader(abundance_handle)
    for row in reader:
        id = row["sequin_id"]
        abundance = float(row["abundance"])
        abx.add(abundance)
        xs.setdefault(abundance, []).append(id)
    ys = []
    abundances = sorted(list(abx))
    for idx, abundance in enumerate(sorted(xs.keys())):
        abundance_taxid = taxid_generator.next()
        xx = []
        for id in xs[abundance]:
            taxid = taxid_generator.next()
            xx.append((id, taxid))
        ys.append({
            "abundance": (abundances[idx], abundance_taxid), 
            "sequins": xx
        })
    return ys


def main() -> None:
    p = argparse.ArgumentParser(
        description="Generate `nodes.dmp` and `names.dmp` files for Kraken2 database building, extending the NCBI taxonomy to include synthetic Sequin entries."
    )
    p.add_argument("--nodes", type=argparse.FileType("r"), required=True, 
                   help="Path to the original NCBI `nodes.dmp` file.")
    p.add_argument("--names", type=argparse.FileType("r"), required=True, 
                   help="Path to the original NCBI `names.dmp` file.")
    p.add_argument("--abundances", type=argparse.FileType("r"), required=True, 
                   help="Path to the Sequin abundances manifest file (e.g. `metasequins_abundances.csv`) provided in the resource bundle (CSV format with columns: `sequin_id`, `abundance`).")
    p.add_argument("--sequin-parent-id", type=int, default=999000000,
               help="Taxonomic ID to assign as the parent node for all sequins (must be > max tax ID in nodes.dmp). Default: 999000000")
    args = p.parse_args()

    max_tax_id = 0
    for line in args.nodes:
        tax_id = int(line.split()[0])
        max_tax_id = max(max_tax_id, tax_id)
    print(f"Highest taxonomic ID: {max_tax_id}")

    # Ensure the sequin parent ID is greater than the highest tax ID in the input file
    if max_tax_id >= args.sequin_parent_id:
        print(
            f"The highest taxid in the input file ({max_tax_id}) is greater than or equal to the specified Sequin parent ID ({args.sequin_parent_id}). Please choose a higher Sequin parent ID.",
            file=sys.stderr,
        )
        sys.exit(1)

    if os.path.isfile("nodes.dmp"):
        print(
            "The file `nodes.dmp` already exists in the current directory, please delete or rename it",
            file=sys.stderr,
        )
        sys.exit(1)

    if os.path.isfile("names.dmp"):
        print(
            "The file `names.dmp` already exists in the current directory, please delete or rename it",
            file=sys.stderr,
        )
        sys.exit(1)

    args.nodes.seek(0)
    taxid_generator = Taxid(start_id=args.sequin_parent_id + 1)
    entries = read_manifest(args.abundances, taxid_generator)

    with open("nodes.dmp", "wt") as handle:
        for line in args.nodes:
            print(line, end="", file=handle)
        print_nodes(entries, handle, sequin_parent_node_id=str(args.sequin_parent_id))
    with open("names.dmp", "wt") as handle:
        for line in args.names:
            print(line, end="", file=handle)
        print_names(entries, handle,sequin_parent_node_id=str(args.sequin_parent_id))
    print("Successfully created `nodes.dmp` and `names.dmp` with Sequin entries.")

if __name__ == "__main__":
    main()

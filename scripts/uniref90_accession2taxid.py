#!/usr/bin/env python3
import sys
import gzip


def open_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_uniref90_fasta(fasta_handle, out_handle):
    for line in fasta_handle:
        if not line.startswith(">"):
            continue

        header = line[1:].strip()

        accession = header.split()[0]

        taxid = None
        for field in header.split():
            if field.startswith("TaxID="):
                taxid = field.split("=", 1)[1]
                break

        if taxid is not None:
            out_handle.write(
                f"{accession}\t{accession}\t{taxid}\t{accession}\n"
            )


def uniref90_accession2taxid(infile, outfile):
    with open_file(infile) as fasta, open(outfile, "w") as out:
        out.write("accession\taccession.version\ttaxid\tgi\n")
        parse_uniref90_fasta(fasta, out)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <uniref90.fasta(.gz)> <output.tsv>")

    infile, outfile = sys.argv[1], sys.argv[2]
    uniref90_accession2taxid(infile, outfile)

if __name__ == "__main__":
    main()

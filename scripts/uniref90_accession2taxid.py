#!/usr/bin/env python3
import sys
import gzip


def open_file(path):
    """根据后缀自动打开普通文件或 gzip 文件"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_uniref90_fasta(fasta_handle, out_handle):
    """
    解析 UniRef90 FASTA header，写出 accession → taxid 映射
    """
    for line in fasta_handle:
        if not line.startswith(">"):
            continue

        header = line[1:].strip()

        # UniRef cluster ID 作为 accession
        accession = header.split()[0]

        # 提取 TaxID
        taxid = None
        for field in header.split():
            if field.startswith("TaxID="):
                taxid = field.split("=", 1)[1]
                break

        if taxid is not None:
            # accession.version 和 gi 使用 accession 填充
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
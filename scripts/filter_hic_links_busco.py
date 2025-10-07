# /usr/bin/python3


"""
script to parse a BUSCO full_table.tsv into an expression for samtools bam filtering
"""

import sys
import argparse
import itertools


def parse_table(inpath: str) -> dict[list]:
    """
    parse a BUSCO full_table.tsv into a dict of lists keyed on busco name
    """
    out = {}
    with open(inpath, "r", encoding="utf-8") as inf:
        for line in inf.readlines():
            line = line.strip().split("\t")
            try:
                out[line[0]].append(line[1:])
            except KeyError:
                out[line[0]] = [line[1:]]
    return out


def get_tets(table: dict, ploidy=4) -> list[tuple]:
    """
    get sets of sequences with complete BUSCOs occurring on ploidy chromosomes
    """
    out = []
    for buscos in table.values():
        if len(buscos) < ploidy:
            continue
        seqs = set(l[1] for l in buscos)
        if len(seqs) == ploidy:  # candidates occur on p separate sequences - this rescues tandems
            out.append(tuple(seqs))
    # print(len(out))
    out = list(set(out))
    # print(len(out))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to parse a BUSCO full_table.tsv into an \
                                     expression for samtools bam filtering",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("table", help="BUSCO full_table.tsv")
    parser.add_argument("input", help="input SAM or stdin", nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument("output", help="output SAM or stdout", nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument("-p", "--ploidy", help="target species ploidy", type=int, default=4)
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    table_dict = parse_table(args.table)
    tets = get_tets(table_dict, args.ploidy)
    tet_pairs = set()
    for tet in tets:
        for pair in itertools.combinations(tet, 2):
            tet_pairs.add(pair)
    tet_pairs = list(tet_pairs)
    # print(args.input.name)
    for line in args.input:
        if line.startswith("@"):
            args.output.write(line)
            continue
        line = line.strip().split("\t")
        seq, pair = line[2], line[6]
        if (seq, pair) in tet_pairs or (pair, seq) in tet_pairs:
            continue
        args.output.write("\t".join(line) + "\n")

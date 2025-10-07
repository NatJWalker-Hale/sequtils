# /usr/bin/python3


"""
script to parse a BUSCO full_table.tsv into an expression for samtools bam filtering
"""

import sys
import argparse
import itertools


def parse_table(inpath: str) -> dict:
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
    out = []
    for busco, buscos in table.items():
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
    parser.add_argument("-p", "--ploidy", help="target species ploidy", type=int, default=4)
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    with open("bamfilter.sh", "w", encoding="utf-8") as outf:
        outf.write("""samtools view -h HiC.filtered.subsample.bam | perl -ne '
if (/^@/) { print; next }
@F = split(/\\t/);
$map = $F[2];
$pair = $F[6];
"""
                   )


        table_dict = parse_table(args.table)
        tets = get_tets(table_dict, args.ploidy)
        for tet in tets:
            for p in itertools.combinations(tet, 2):
                outf.write(f"if ((($map eq \"{p[0]}\") && ($pair eq \"{p[1]}\")) || "
                           f"(($map eq \"{p[1]}\") && ($pair eq \"{p[0]}\"))) {{\nnext;\n" + "}\n")
                
        outf.write("print;' \\\n| samtools view -b -o HiC.filtered.subsample.filtered.bam\n")

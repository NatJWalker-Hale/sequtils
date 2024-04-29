#! /usr/bin/python3


import re
import sys
import argparse


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("GFF", help="input GFF3 file to rename")
    parser.add_argument("prefix", help="string to append in front on new id")
    args = parser.parse_args()

    # relies on accurate GFF ordering - will fail otherwise
    corres = {}
    with open(args.GFF, "r") as f:
        line = f.readline().strip()
        count = 0
        while line:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            idstr = fields[-1]
            if fields[2] == "gene":
                count += 1
                ids = idstr.split(";")
                for i in ids:
                    if i.startswith("ID="):
                        id = i.split("=")[1]
            scaff_int = fields[0].split("_")[-1]
            newid = f"{args.prefix}{scaff_int}g{count}"
            idstr = re.sub(id, newid, idstr)
            # outf.write(f"{id}\t{newid}\n")
            corres[id] = newid # will overwrite but always same value so fine
            print("\t".join(fields[:-1] + [idstr]))
            line = f.readline().strip()

    with open("correspondence.tsv", "w") as outf:
        for k, v in corres.items():
            outf.write(f"{k}\t{v}\n")

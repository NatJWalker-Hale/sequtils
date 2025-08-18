#! /usr/bin/python3


import os
import sys
import argparse


class BED:
    def __init__(self):
        self.records = {"": []}
    
    def add(self, BED_rec):
        self.records[BED_rec.name] = BED_rec


class BED_rec:
    def __init__(self, line=[]):
        self.raw = line
        self.chrom = line[0]
        self.start = line[1]
        self.end = line[2]
        self.name = line[3]
        self.score = line[4]
        self.strand = line[5]

    def __str__(self):
        return "\t".join(self.raw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to convert python MCScanX simple anchors \
                                     to links files for circos")
    parser.add_argument("blocks", help="python MCScanX .anchors.simple file")
    parser.add_argument("bed1", help="reference BED")
    parser.add_argument("bed2", help="query BED")
    args = parser.parse_args()

    bed1 = BED()
    with open(args.bed1) as inf:
        for line in inf:
            line = line.strip().split("\t")
            new_rec = BED_rec(line)
            bed1.add(new_rec)

    bed2 = BED()
    with open(args.bed2) as inf:
        for line in inf:
            line = line.strip().split("\t")
            new_rec = BED_rec(line)
            bed2.add(new_rec)

    with open(args.blocks) as inf:
        for line in inf:
            line = line.strip().split("\t")
            # print(line)
            id1 = bed1.records[line[0]].chrom
            id2 = bed2.records[line[2]].chrom
            start1 = bed1.records[line[0]].start
            end1 = bed1.records[line[1]].end
            if line[-1] == "+":
                start2 = bed2.records[line[2]].start
                end2 = bed2.records[line[3]].end
            else:  # reverse when phase is negative
                start2 = bed2.records[line[3]].start
                end2 = bed2.records[line[2]].end
            print(" ".join([id1, start1, end1, id2, start2, end2]))


    
            
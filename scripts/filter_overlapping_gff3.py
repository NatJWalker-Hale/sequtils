#! /usr/bin/python3


import os
import sys
import argparse
import gffutils


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to detect overlapping genes and report \
                                     their attributes")
    parser.add_argument("gff3", help="input gff3 file.")

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    dbname = args.gff3.replace(".gff3", ".db")
    if os.path.isfile(dbname) and os.path.getsize(dbname) > 0:
        db = gffutils.FeatureDB(dbname)
    else:
        db = gffutils.create_db(args.gff3, dbname, id_spec={'gene': 'ID',
                                                            'mRNA': 'ID'},
                                                            keep_order=True
        )

    for gene in db.features_of_type(featuretype="gene", order_by="featuretype"):
        for feature in db.region(gene, featuretype="gene"):
            if gene.id != feature.id:
                print(f"{gene.id} with strand {gene.strand} overlaps with {feature.id} with strand {feature.strand}")

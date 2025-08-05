#! /usr/bin/python3


import os
import re
import sys
import argparse
import gffutils


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to rename GFF3s to format SpStr01G0000100")
    parser.add_argument("gff3", help="input gff3 file. By default expects sequence IDs to look \
                        like CHR01, any other sequence IDs will lead to SpStrUnG0000100")
    parser.add_argument("sp", help="string to use as SpStr")
    args = parser.parse_args()


    dbname = args.gff3.rstrip(".gff3") + ".db"
    if os.path.isfile(dbname) and os.path.getsize(dbname) > 0:
        db = gffutils.FeatureDB(dbname)
    else:
        db = gffutils.create_db(args.gff3, dbname, id_spec={'gene': 'ID',
                                                            'mRNA': 'ID'},
                                                            keep_order=True
                               )

    print("##gff-version 3")
    print(f"# {" ".join(sys.argv)}")
    num = 100
    chrom = ""
    for gene in db.features_of_type(featuretype="gene", order_by="featuretype"):
        if gene.seqid != chrom:  # new chr
            chrom = gene.seqid  # set chr
            if "CHR" in gene.seqid:
                num = 100  # reset number if we go to a new chrom
            # otherwise for unplaced, we just keep counting
        if "CHR" in chrom:  # expects seqids to be CHR01 etc.
            chrom_num = chrom.lstrip("CHR")  # leave just number
        else:
            chrom_num = "Un"  # for unplaced scaffolds/contigs
        new_id = f"{args.sp}{chrom_num}G{str(num).zfill(7)}"  # pad to 7 digits
        num += 100  # iterate num

        # now for attributes we need specific formatting
        gene.attributes['ID'] = new_id
        gene.attributes['Name'] = new_id
        gene.attributes['Alias'] = gene.id
        print(gene)

        trans_n = 0
        transcripts = []
        for mrna in db.children(gene, featuretype="mRNA"):
            cds_length = sum(cds.end - cds.start for cds in db.children(mrna, featuretype="CDS"))
            transcripts.append((mrna, cds_length))
            # print(transcripts)

        for mrna, _ in sorted(transcripts, key=lambda x:x[1], reverse=True):
            # so .1 is longest transcript
            trans_n += 1
            trans_id = new_id + f".{trans_n}"  # .1 .2 etc. formatting
            mrna.attributes['Parent'] = new_id  # keep accurate linking
            mrna.attributes['ID'] = trans_id
            mrna.attributes['Name'] = trans_id  # update name
            mrna.attributes['Alias'] = mrna.id
            print(mrna)

            exon_n = 0
            for exon in db.children(mrna, featuretype="exon", order_by="start"):
                exon_n += 1
                exon.attributes['Parent'] = trans_id
                exon.attributes['ID'] = f"{trans_id}.exon{exon_n}"
                print(exon)

            cds_n = 0
            for cds in db.children(mrna, featuretype="CDS", order_by="start"):
                cds_n += 1
                cds.attributes['Parent'] = trans_id
                cds.attributes['ID'] = f"{trans_id}.cds{cds_n}"
                print(cds)

            fiveputr_n = 0
            for fiveputr in db.children(mrna, featuretype="five_prime_UTR", order_by="start"):
                fiveputr_n += 1
                fiveputr.attributes['Parent'] = trans_id
                fiveputr.attributes['ID'] = f"{trans_id}.utr5p{fiveputr_n}"
                print(fiveputr)

            threeputr_n = 0
            for threeputr in db.children(mrna, featuretype="three_prime_UTR", order_by="start"):
                threeputr_n += 1
                threeputr.attributes['Parent'] = trans_id
                threeputr.attributes['ID'] = f"{trans_id}.utr3p{threeputr_n}"
                print(threeputr)

            




        


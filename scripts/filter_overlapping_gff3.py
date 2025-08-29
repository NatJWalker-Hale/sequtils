#! /usr/bin/python3


import os
import sys
import argparse
import gffutils


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to detect overlapping genes and report \
                                     their attributes")
    subparser = parser.add_subparsers(title="mode", dest="mode", required=True)

    parser_d = subparser.add_parser("detect", help="detect overlapping genes")
    parser_d.add_argument("query_gff3", help="input gff3 file")
    
    parser_f = subparser.add_parser("filter", help="filter overlapping genes")
    parser_f.add_argument("query_gff3", help="input gff3 file")

    parser_c = subparser.add_parser("compare", help="detect overlapping genes and compare to \
                                    reference")
    parser_c.add_argument("query_gff3", help="input query gff3 file")
    parser_c.add_argument("reference_gff3", help="input reference gff3 file")
    parser_c.add_argument("blastp_results", help="blastp results of query vs reference, tabular \
                                    format, qseqid sseqid")
    
    parser_m = subparser.add_parser("merge", help="merge adjacent genes on same strand with \
                                    overlapping CDS features")
    parser_m.add_argument("query_gff3", help="input query gff3 file")

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    # load query into gffutils db
    qdbname = args.query_gff3.replace(".gff3", ".db")
    if os.path.isfile(qdbname) and os.path.getsize(qdbname) > 0:
        qdb = gffutils.FeatureDB(qdbname)
    else:
        qdb = gffutils.create_db(args.query_gff3, qdbname, id_spec={'gene': 'ID',
                                                                    'mRNA': 'ID'},
                                                                    keep_order=True
        )

    if args.mode == "detect":
        visited = []
        for gene in qdb.features_of_type(featuretype="gene", order_by="featuretype"):
            for feature in qdb.region(f"{gene.seqid}:{gene.start}-{gene.end}", featuretype="gene"):
                if feature.id in visited:
                    continue
                if gene.id != feature.id:
                    print(f"{gene.id} with strand {gene.strand} overlaps with {feature.id} "
                        f"with strand {feature.strand}")
            visited.append(gene.id)

    if args.mode == "filter":
        # do something here to filter out overlapping genes based on a criterion
        sys.exit()

    if args.mode == "merge":
        # do something here to merge adjacent genes on same strand with overlapping CDS features
        sys.exit()

    # parse blastp results
    if args.mode == "compare":
        with open(args.blastp_results, "r", encoding="utf-8") as infile:
            blastp_hits = {line.split("\t")[0]: line.split("\t")[1] for line in infile if
                           line.strip()}
        # print(blastp_hits)

        # load reference into gffutils db
        rdbname = args.reference_gff3.replace(".gff3", ".db")
        if os.path.isfile(rdbname) and os.path.getsize(rdbname) > 0:
            rdb = gffutils.FeatureDB(rdbname)
        else:
            rdb = gffutils.create_db(args.reference_gff3, rdbname, id_spec={'gene': 'ID',
                                                                            'mRNA': 'ID'},
                                                                            keep_order=True
            )


        visited = []
        for gene in qdb.features_of_type(featuretype="gene", order_by="featuretype"):
            for feature in qdb.region(f"{gene.seqid}:{gene.start}-{gene.end}", featuretype="gene"):
                if feature.id in visited:
                    continue
                if gene.id != feature.id:
                    print(f"{gene.id} with strand {gene.strand} overlaps with {feature.id} "
                        f"with strand {feature.strand}")
                    qtranscripts = list(qdb.children(gene, featuretype="mRNA", order_by="start"))
                    qmax_exons_transcript = max(qtranscripts,
                                                key=lambda t: len(list(
                                                    qdb.children(t, featuretype="CDS"))))
                    # have to use CDS throughout because some exons are utr only
                    qmax_n_exons = len(list(qdb.children(qmax_exons_transcript, featuretype="CDS")))
                    if qmax_exons_transcript.id in blastp_hits:
                        ref_transcript_id = blastp_hits[qmax_exons_transcript.id]
                        ref_gene_id = list(rdb.parents(ref_transcript_id, featuretype="gene"))[0].id
                        rtranscripts = list(rdb.children(ref_gene_id,
                                                            featuretype="mRNA",
                                                            order_by="start"))
                        rmax_exons_transcript = max(rtranscripts,
                                                    key=lambda t: len(list(
                                                    rdb.children(t, featuretype="CDS"))))
                        rmax_n_exons = len(list(rdb.children(rmax_exons_transcript,
                                                                featuretype="CDS")))
                        print(f"Query gene {gene.id} with {qmax_n_exons} exons "
                                f"matches reference gene {ref_gene_id} with {rmax_n_exons} exons")
                        exons_diff = (qmax_n_exons != rmax_n_exons)
                            # print(f"WARNING: different number of exons between query and reference. "
                            #       f"Check {gene.id} for fixing.")
                        # need a secondary check because some genes have same number of exons but are 
                        # still obviously wrong
                        q_gene_length = gene.end - gene.start + 1
                        r_gene_length = rdb[ref_gene_id].end - rdb[ref_gene_id].start + 1
                        length_diff =  (q_gene_length > 3 * r_gene_length or
                                        r_gene_length > 3 * q_gene_length)
                            # print(f"WARNING: different gene lengths between query and reference. "
                            #       f"Check {gene.id} for fixing.")
    
                        # need also to add a strand check
                        print((gene.strand == rdb[ref_gene_id].strand))
                        strand_same = (gene.strand == feature.strand)
                            # print(f"WARNING: different strand between query and reference. "
                            #       f"Check {gene.id} for fixing.")
                        if exons_diff and length_diff and strand_same:
                            print(f"WARNING: {gene.id} has overlapping genes on same strand "
                                f"different exon count and > 3x length, relative to {ref_gene_id}. "
                                f"Definite fix.")
                        elif any([exons_diff, length_diff, strand_same]):
                            print(f"WARNING: {gene.id} has ones of overlapping genes on same strand "
                                f"different exon count and > 3x length, relative to {ref_gene_id}. "
                                f"Check fix.")

            visited.append(gene.id)

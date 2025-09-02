#! /usr/bin/python3


import os
import sys
import argparse
import gffutils
import gffutils.interface


def write_gene(dbobj: gffutils.interface.FeatureDB, geneid: str) -> str:
    out_str = ""
    geneobj = dbobj[geneid]
    out_str += f"{str(geneobj)}\n"
    for transcript in dbobj.children(geneobj, featuretype="mRNA"):
        out_str += f"{str(transcript)}\n"
        for cds in dbobj.children(transcript, featuretype="CDS"):
            out_str += f"{str(cds)}\n"
        for exon in dbobj.children(transcript, featuretype="exon"):
            out_str += f"{str(exon)}\n"
        for fputr in dbobj.children(transcript, featuretype="five_prime_UTR"):
            out_str += f"{str(fputr)}\n"
        for tputr in dbobj.children(transcript, featuretype="three_prime_UTR"):
            out_str += f"{str(tputr)}\n"
    return out_str


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to detect genes with isoforms that don't \
                                     share any exons, indicative of bad PASA merge")
    parser.add_argument("pasa_gff3", help="input PASA-updated gff3 file.")
    parser.add_argument("evm_gff3", help="original EVM gff3")
    parser.add_argument("-m", "--mode", choices=["report", "replace"], help="report: report \
                        genes with isoforms that don't share exons, replace: replace merged \
                        prediction with original EVM prediction", default="report")
    parser.add_argument("--replace_list", help="file with IDs of genes to fix, one per line. \
                        Only used in fix mode. If not provided, all genes with isoforms that don't \
                        share exons will be fixed.")

    args = parser.parse_args()


    dbname = args.pasa_gff3.replace(".gff3", ".db")
    if os.path.isfile(dbname) and os.path.getsize(dbname) > 0:
        db = gffutils.FeatureDB(dbname)
    else:
        db = gffutils.create_db(args.pasa_gff3, dbname, id_spec={'gene': 'ID',
                                                                 'mRNA': 'ID'},
                                                                 keep_order=True
        )

    evmdbname = args.evm_gff3.replace(".gff3", ".db")
    if os.path.isfile(evmdbname) and os.path.getsize(evmdbname) > 0:
        evmdb = gffutils.FeatureDB(evmdbname)
    else:
        evmdb = gffutils.create_db(args.evm_gff3, evmdbname, id_spec={'gene': 'ID',
                                                                      'mRNA': 'ID'},
                                                                       keep_order=True
        )

    if args.mode == "replace":
        if args.replace_list:
            with open(args.replace_list, encoding="utf-8") as f:
                fix_ids = set(line.strip() for line in f)
        print("##gff-version 3", file=sys.stdout)
    
    for gene in db.features_of_type(featuretype="gene", order_by="featuretype"):  # iterate genes
        print(f"working on {gene.id}", file=sys.stderr)
        transcripts = list(db.children(gene, featuretype="mRNA"))
        if len(transcripts) < 2:
            # only interested in genes with multiple transcripts
            if args.mode == "replace":
                if gene.id not in fix_ids:
                    sys.stderr.write(f"Processing gene {gene.id} - no fix needed\n")
                    print(write_gene(db, gene.id))
                else:
                    sys.stderr.write(f"Processing gene {gene.id}\n")
                    evm_genes = gene.id.split("_")
                    for g in evm_genes:
                        print(evmdb[g])
                        count = 1
                        for c in evmdb.children(g, order_by="featuretype"):
                            if c.featuretype == "CDS":
                                c.attributes['ID'][0] += f".cds{count}"
                                count += 1
                                print(c)
                            else:
                                print(c)
            continue
            # also no way to tell if genes with single transcript are bad merges except by reads
            # or homology - we do this elsewhere
        exons = [list(db.children(mrna, featuretype="CDS")) for mrna in transcripts]
        # note that I'm calling these exons but using the CDS span to avoid differing UTR lengths
        starts = [set(exon.start for exon in exon_list) for exon_list in exons]
        ends = [set(exon.end for exon in exon_list) for exon_list in exons]
        shared_starts = set.intersection(*starts)
        # print(f"Gene {gene.id}", shared_starts)
        shared_ends = set.intersection(*ends)
        # if gene.id == "evm.TU.CHR02.1895":
        #     print(starts)
        #     print(ends)
        #     print(shared_starts)
        #     print(shared_ends)
        # print(f"Gene {gene.id}", shared_ends)
        # if at least one exon start and exon end is shared then the isoforms share at least one
        # exon
        if not shared_starts or not shared_ends:
            # no shared exons, so this is a bad merge
            if args.mode == "report":
                print(f"Gene {gene.id} has isoforms that don't share exons")
            elif args.mode == "replace":
                if gene.id not in fix_ids:
                    # we've decided not to fix this one for some reason, so we need to write it
                    print(write_gene(db, gene.id))
                    continue
                sys.stderr.write(f"Processing gene {gene.id}\n")
                evm_genes = gene.id.split("_")
                for g in evm_genes:
                    print(evmdb[g])
                    count = 1
                    for c in evmdb.children(g, order_by="featuretype"):
                        if c.featuretype == "CDS":
                            c.attributes['ID'][0] += f".cds{count}"
                            count += 1
                            print(c)
                        else:
                            print(c)
        else:
            if args.mode == "replace":
                if gene.id not in fix_ids:
                    sys.stderr.write(f"Processing gene {gene.id} - no fix needed\n")
                    print(write_gene(db, gene.id))
                else:
                    sys.stderr.write(f"Processing gene {gene.id}\n")
                    evm_genes = gene.id.split("_")
                    for g in evm_genes:
                        print(evmdb[g])
                        count = 1
                        for c in evmdb.children(g, order_by="featuretype"):
                            if c.featuretype == "CDS":
                                c.attributes['ID'][0] += f".cds{count}"
                                count += 1
                                print(c)
                            else:
                                print(c)


                # old fixing, too hard
                # # first handle cases with two non-overlapping transcripts
                # if len(transcripts) == 2:
                #     print(f"splitting {gene.id} with two transcripts")
                #     gene_start = min(t.start for t in transcripts)
                #     gene_end = max(t.end for t in transcripts)
                #     for t in transcripts:
                #         new_gene_line = "\t".join(
                #             [
                #                 f"{gene.seqid.split("_")[0]}",
                #                 "EVM",
                #                 "gene",
                #                 f"{gene.start}",
                #                 f"{gene.end}",
                #                 f"{gene.score}",
                #                 f"{gene.strand}",
                #                 f"{gene.frame}",
                #                 f"{";".join(f"{x}={gene.attributes[x][0]}"
                #                             for x in gene.attributes)}"
                #             ]
                #         )
                #         new_gene = gffutils.feature.feature_from_line(
                #             new_gene_line
                #         )
                #         print(new_gene)         
                # # get transcripts that do share exons
                # else:
                #     transcript_clusters = []
                #     for transcript in transcripts:
                #         current_exons = list(db.children(transcript, featuretype="CDS"))
                #         current_starts = set(exon.start for exon in current_exons)
                #         current_ends = set(exon.end for exon in current_exons)
                #         # find other transcripts that share at least one exon start and end
                #         shared_transcripts = []
                #         for other_transcript in transcripts:
                #             if other_transcript == transcript:
                #                 continue
                #             other_exons = list(db.children(other_transcript, featuretype="CDS"))
                #             other_starts = set(exon.start for exon in other_exons)
                #             other_ends = set(exon.end for exon in other_exons)
                #             if (
                #                 current_starts.intersection(other_starts) and
                #                 current_ends.intersection(other_ends)
                #             ):
                #                 shared_transcripts.append(other_transcript)
                #         if shared_transcripts:
                #             transcript_clusters.append([transcript] + shared_transcripts)
                #     # for cluster in transcript_clusters:
                #     #     print(" ".join(t.id for t in cluster))
                #     remove = set.intersection(*[set(t) for t in transcript_clusters])
                #     # this is the transcript or transcripts that we are going to bin to split the
                #     # remaining transcripts into new genes
                #     if remove:
                #         for cluster in transcript_clusters:
                #             if cluster[0] == list(remove)[0]:
                #                 continue
                #             cluster.remove(list(remove)[0])
                #             # for

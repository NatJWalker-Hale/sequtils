#! /usr/bin/python3


import os
import argparse
import iranges
import gffutils


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to detect genes with isoforms that don't \
                                     share any exons, indicative of bad PASA merge")
    parser.add_argument("gff3", help="input gff3 file.")
    parser.add_argument("-m", "--mode", choices=["report", "fix"], help="report: report genes with \
                         isoforms that don't share exons, fix: split such genes into multiple \
                        genes with isoforms that share exons", default="report")
    parser.add_argument("--fix_list", help="file with IDs of genes to fix, one per line. \
                        Only used in fix mode. If not provided, all genes with isoforms that don't \
                        share exons will be fixed.")
    
    args = parser.parse_args()


    dbname = args.gff3.rstrip(".gff3") + ".db"
    if os.path.isfile(dbname) and os.path.getsize(dbname) > 0:
        db = gffutils.FeatureDB(dbname)
    else:
        db = gffutils.create_db(args.gff3, dbname, id_spec={'gene': 'ID',
                                                            'mRNA': 'ID'},
                                                            keep_order=True
        )

    
    for gene in db.features_of_type(featuretype="gene", order_by="featuretype"):  # iterate genes
        transcripts = list(db.children(gene, featuretype="mRNA"))
        if len(transcripts) < 2:
            continue  # only interested in genes with multiple transcripts
            # also no way to tell if genes with single transcript are bad merges except by reads
            # or homology - we do this elsewhere
        exons = [list(db.children(mrna, featuretype="CDS")) for mrna in transcripts]
        # note that I'm calling these exons but using the CDS span to avoid differing UTR lengthss
        starts = [set(exon.start for exon in exon_list) for exon_list in exons]
        ends = [set(exon.end for exon in exon_list) for exon_list in exons]
        shared_starts = set.intersection(*starts)
        # print(f"Gene {gene.id}", shared_starts)
        shared_ends = set.intersection(*ends)
        # print(f"Gene {gene.id}", shared_ends)
        # if at least one exon start and exon end is shared then the isoforms share at least one exon
        if not shared_starts or not shared_ends:
            # no shared exons, so this is a bad merge
            if args.mode == "report":
                print(f"Gene {gene.id} has isoforms that don't share exons.")
            elif args.mode == "fix":
                print(f"Processing gene {gene.id}")
                if args.fix_list:
                    with open(args.fix_list) as f:
                        fix_ids = set(line.strip() for line in f)
                    if gene.id not in fix_ids:
                        continue
                # get transcripts that do share exons
                transcript_clusters = []
                for transcript in transcripts:
                    current_exons = list(db.children(transcript, featuretype="CDS"))
                    current_starts = set(exon.start for exon in current_exons)
                    current_ends = set(exon.end for exon in current_exons)
                    # find other transcripts that share at least one exon start and end
                    shared_transcripts = []
                    for other_transcript in transcripts:
                        if other_transcript == transcript:
                            continue
                        other_exons = list(db.children(other_transcript, featuretype="CDS"))
                        other_starts = set(exon.start for exon in other_exons)
                        other_ends = set(exon.end for exon in other_exons)
                        if current_starts.intersection(other_starts) and current_ends.intersection(other_ends):
                            shared_transcripts.append(other_transcript)
                    if shared_transcripts:
                        transcript_clusters.append([transcript] + shared_transcripts)
                # for cluster in transcript_clusters:
                #     print(" ".join(t.id for t in cluster))
                remove = set.intersection(*[set(t) for t in transcript_clusters])
                # this is the transcript or transcripts that we are going to bin to split the remaining
                # transcripts into new genes
                if remove:
                    for cluster in transcript_clusters:
                        if cluster[0] == list(remove)[0]:
                            continue
                        cluster.remove(list(remove)[0])
                        for 

                
        
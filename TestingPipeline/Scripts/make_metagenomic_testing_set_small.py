#!/usr/bin/python

#usage: python make_metagenomic_testing_set_small.py [TESTING SET INDEX]

import Bio
from Bio import SeqIO
import random
import sys

testing_set_index = sys.argv[1]

#function for getting a subset of sequences
def get_fasta_subset(infile, totalnumseqs, keepnum, seqtype):
    keeplist=random.sample(range(1,totalnumseqs), k=keepnum)
    keeplist.sort()
    #print(keeplist)
    i=0
    for seqrecord in SeqIO.parse(infile, "fasta"):
        print("check seq")
        i+=1
        if (i in keeplist):
            index=keeplist.index(i)
            #print("index: %s\n" % (index))
            #get a random subset of the contig that is ≤ 2 Mb
            limit=len(seqrecord.seq)
            #print(limit)
            if (limit>2000000):
                start=random.randint(0, limit - 2000000)
                end=start+2000000
                trimmedseq = seqrecord.seq[start:end]
                outfile.write(">%s--%s\n%s\n" % (seqtype, seqrecord.id, trimmedseq))
                #print("trimmed seq")
            elif (limit<3000):
                #print("too short")
                #print("keeplist index: %s" % (keeplist[index]))
                keeplistindex=keeplist[index]
                while (keeplistindex in keeplist):
                #   print("in while loop")
                    keeplistindex=keeplistindex+1
                #    print("keeplist index: %s" % (keeplist[index]))
                #    print("new keeplist: %s" % (keeplistindex))
                #    print("keeplist: %s" % (keeplist))
                keeplist.append(keeplistindex)
                #print("new keeplist: %s" % (keeplist))
                keeplist.sort()
            else:
                #print("no issues")
                outfile.write(">%s--%s\n%s\n" % (seqtype, seqrecord.id, seqrecord.seq))
            keeplist.pop(index)
        if (keeplist==[]):
            break

outname="metagenomic_testing_set_%s.fna" % (testing_set_index)

with open(outname, 'at') as outfile:
    #viruses
    print("viruses")
    get_fasta_subset(infile="/nfs/turbo/lsa-duhaimem/VSTE/TestingSet/Viruses/refseq_and_vs2_viruses.fna", totalnumseqs=49500, keepnum=1000, seqtype="virus")
    #bacteria
    print("bacteria")
    get_fasta_subset(infile="/nfs/turbo/lsa-dudelabs/databases/genomes/refseq/bacteria/bacteria_refseq_NOV2019.fna", totalnumseqs=5300000, keepnum=6500, seqtype="bacteria")
    #archaea
    print("archaea")
    get_fasta_subset(infile="/nfs/turbo/lsa-dudelabs/databases/genomes/refseq/archaea/archaea_refseq_NOV2019.fna", totalnumseqs=55100, keepnum=1000, seqtype="archaea")
    #plasmids
    print("plasmids")
    get_fasta_subset(infile="/nfs/turbo/lsa-duhaimem/VSTE/All_plasmids.fna", totalnumseqs=6642, keepnum=500, seqtype="plasmid")
    #protists
    print("protists")
    get_fasta_subset(infile="/nfs/turbo/lsa-dudelabs/databases/genomes/refseq/protists/combined_protists_refseq.fna", totalnumseqs=69600, keepnum=500, seqtype="protist")
    #fungi
    print("fungi")
    get_fasta_subset(infile="/nfs/turbo/lsa-dudelabs/databases/genomes/refseq/fungi/combined_fungi_refseq.fna", totalnumseqs=216500, keepnum=500, seqtype="fungi")
    print("done")

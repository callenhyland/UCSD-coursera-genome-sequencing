# -*- coding: utf-8 -*-
"""
Functions for Week 1 of UCSD Genome Sequencing Course
Callen Hyland, August 2018
"""

from collections import defaultdict


def kmer_sort(Text, k):
    kmers = []
    for i in range(len(Text)-k+1):
        kmer = Text[i:i+k]
        kmers.append(kmer)
    return(sorted(kmers))


def composition(text, k):
    kmers = [text[i:i+k] for i in range(len(text)-k+1)]
    return(sorted(kmers))
    

def assemble(kmers):
    genome = kmers[0]
    for i in range(1,len(kmers)):
        genome = genome + kmers[i][-1]
    return(genome)


def overlap_graph(kmers):
    graph = []
    l = len(kmers[0])-1
    for kmer in kmers:
        suffix = kmer[-l:]
        prefix_matches = []
        for next_kmer in kmers:
            if suffix == next_kmer[0:l]:
                prefix_matches.append(next_kmer)
        if prefix_matches != []:
            graph.append(kmer + " -> " + ",".join(prefix_matches))
    return(graph)


def debruijn_from_string(text, k):
    m = k-1
    adj = defaultdict(list)
    for i in range(len(text)-m):
        adj[text[i:i+m]].append(text[i+1:i+m+1])
    nodes = sorted(adj.keys())
    return([(node + " -> " + ",".join(adj[node])) for node in nodes])


def debruijn_from_kmers(kmers):
    db = defaultdict(list)
    for kmer in kmers:
        db[kmer[0:-1]].append(kmer[1:])
    nodes = sorted(db.keys())
    return([(node + " -> " + ",".join(db[node])) for node in nodes])


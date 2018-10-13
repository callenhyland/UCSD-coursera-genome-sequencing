# Code-graded problem 2: Create overlap graph for a list of kmers
# Input: list of kmers read from a text file
# Return: Overlap graph of (k-1)-mers separated by "->"
# Created July 2018, Callen Hyland


import sys


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


lines = sys.stdin.read().splitlines()
print(*overlap_graph(lines), sep='\n')


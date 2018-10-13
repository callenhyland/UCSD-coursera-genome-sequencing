# Solve debruin graph from a string problem, code-graded question #3
# Input: DNA sequence (string)
# Return: Debruijn Graph with (k-1)-mers separated by "->"
# Created July 2018, Callen Hyland

from collections import defaultdict
import sys


def debruijn_from_string(text, k):
    m = k-1
    adj = defaultdict(list)
    for i in range(len(dna)-m):
        adj[dna[i:i+m]].append(dna[i+1:i+m+1])
    nodes = sorted(adj.keys())
    return([(node + " -> " + ",".join(adj[node])) for node in nodes])


k = int(input())
dna = input().strip()

print(*debruijn_from_string(dna), sep='\n')

# Code graded problem #4, Create debruin graph from list of kmers
# Input: List of kmers
# Return: Debruijn Graph with (k-1)-mers separated by "->"
# Created July 2018, Callen Hyland


from collections import defaultdict
import sys


def debruijn_from_kmers(kmers):
    
    """
    Parameters: list of kmers
    Returns: a Debruijn graph
    """
    
    db = defaultdict(list)
    for kmer in kmers:
        db[kmer[0:-1]].append(kmer[1:])
    nodes = sorted(db.keys())
    return([(node + " -> " + ",".join(db[node])) for node in nodes])


lines = sys.stdin.read().splitlines()
print(*debruijn_from_kmers(lines), sep='\n')

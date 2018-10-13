# Solve the string composition problem (code-graded problem 1)
# INPUT: DNA sequence (string)
# RETURN: sorted list of kmers in the DNA string


def composition(text, k):
    return(sorted([text[i:i+k] for i in range(len(text)-k+1)]))


k = int(input())
text = input().strip()

print(*composition(text, k), sep='\n')

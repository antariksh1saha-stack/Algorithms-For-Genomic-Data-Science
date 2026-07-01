from bm_preproc import BoyerMoore
from kmer_index import Index

def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    alingments=0
    character_comparisons=0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            character_comparisons +=1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
        alingments +=1
    return occurrences#,alingments,character_comparisons

def naive_with_counts(p,t):
    occurrences = []
    alingments=0
    character_comparisons=0
    for i in range(len(t)-len(p)+1):
        match  = True
        for j in range(len(p)):
            character_comparisons +=1
            if not t[i+j]==p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
        alingments +=1
    return occurrences,alingments,character_comparisons

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def naive_2mm(p,t):
    occurrences = []
    for i in range(len(t)-len(p)+1):
        miss_match  = 0
        for j in range(len(p)):
            if not t[i+j]==p[j]:
                miss_match += 1
                if miss_match > 2:
                    break
        if miss_match <= 2:
            occurrences.append(i) 
    return occurrences

import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
def queryIndex(p,t,index):
        k=index.k
        offsets=[]
        for i in index.query(p):
            if p[k:]==t[i+k:i+len(p)]:
                offsets.append(i)
        return offsets
def approximate_matching(p,t,n):
    length_partition=int(round(len(p)/(n+1)))
    all_matches=set()
    index_hits=0
    for i in range(n+1):
        start = i*length_partition
        end = (i+1)*length_partition
        p_bm=BoyerMoore(p[start:end],alphabet='ACGT')
        matches = boyer_moore_with_counts(p[start:end],p_bm,t)
        index_hits +=len(matches)
        for offset in matches:
            if offset < start or (offset - start) +len(p)>len(t):
                continue
            mismatches=0
            for j in range (0,start):
                if not p[j]==t[offset-start + j]:
                    mismatches +=1
                    if mismatches > n:
                        break
            for j in range (end,len(p)):
                if not p[j]==t[offset-start + j]:
                    mismatches +=1
                    if mismatches > n:
                        break
            if mismatches <=n:
                all_matches.add(offset-start)
    return list(all_matches),index_hits


lowercase_alphabet = 'ATGC'
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = readGenome('chr1.GRCh38.excerpt.fasta')
p_bm = BoyerMoore(p, lowercase_alphabet)

print("Naive: ",naive_with_counts(p,t))
print("Boyer Moore: ",boyer_moore_with_counts(p,p_bm,t))
index = Index(t,8)
print("Index hits: ",len(index.query(p)))
print("K_mer: ", queryIndex(p,t,index))
print("Naive 2mm: ",len(naive_2mm(p,t)))
apmatches,indexHits= approximate_matching(p,t,n=2)
print("k_mer approx matching ",len(apmatches))
print("Index hits with kmer: ",indexHits)

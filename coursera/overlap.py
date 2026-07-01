from itertools import permutations
def overlap(a,b,min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0 
        if b.startswith(a[start:]):
            return len(a)-start
        start +=1

def naive_overlap_map(reads,k):
    olaps={}
    for a,b in permutations(reads,2):
        olen=overlap(a,b,min_length=k)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def pick_maximal_overlap(reads,k):
    reada,readb= None , None
    best_olen=0
    for a,b in permutations(reads,2):
        olen=overlap(a,b,min_length=k)
        if olen > best_olen:
            reada,readb = a,b
            best_olen=olen
    return reada,readb,best_olen

def greedy_scs(reads,k):
    read_a,read_b,olen=pick_maximal_overlap(reads,k)
    while olen>0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a,read_b,olen=pick_maximal_overlap(reads,k)
    return ''.join(reads)


#print("Greedy: ",len(greedy_scs(reads,1)))

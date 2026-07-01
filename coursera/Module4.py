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
def scs_modified(ss):
    shortest_sup= ['0']
    for ssperm in permutations(ss):
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i],ssperm[i+1],min_length=1)
            sup +=ssperm[i+1][olen:]
        if shortest_sup[0]=='0':
            shortest_sup[0]=sup
        if len(sup)<len(shortest_sup[0]):
            shortest_sup.append(sup)
    return shortest_sup

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences
def pick_maximal_overlap(reads,k):
    reada,readb= None , None
    best_olen=0
    for a,b in permutations(reads,2):
        olen=overlap(a,b,min_length=k)
        if olen > best_olen:
            reada,readb = a,b
            best_olen=olen
    return reada,readb,best_olen

def overlap_graph(reads, k):
    # kmer -> set(reads containing that kmer)
    all_kmers = {}
    # build dictionary
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in all_kmers:
                all_kmers[kmer] = set()
            all_kmers[kmer].add(read)
    reada,readb=None,None
    best_olen=0
    # find overlaps
    for a in reads:
        suffix = a[-k:]
        # only reads containing suffix kmer
        candidates = all_kmers[suffix]
        for b in candidates:
            # do not overlap read with itself
            if a == b:
                continue
            olen = overlap(a, b, min_length=k)
            if olen > best_olen:
                reada,readb = a,b
                best_olen=olen
    return reada,readb,best_olen

def greedy_scs(reads,k):
    read_a,read_b,olen=overlap_graph(reads,k)
    while olen>0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a,read_b,olen=overlap_graph(reads,k)
    return ''.join(reads)

reads=readFastq('ads1_week4_reads.fq')
assembled_genome=greedy_scs(reads,30)
length=len(assembled_genome)
a_count=0
t_count=0
print(length)
for i in range(length):
    if assembled_genome[i]=='A':
        a_count+=1
    if assembled_genome[i]=='T':
        t_count+=1
print("A counts: ",a_count)
print("T counts: ",t_count)



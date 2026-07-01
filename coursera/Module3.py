def readFastq(filename):
    sequences = []

    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()

            if len(seq) == 0:
                break

            sequences.append(seq)

    return sequences


def overlap(a, b, min_length=30):

    start = 0

    while True:

        start = a.find(b[:min_length], start)

        if start == -1:
            return 0

        if b.startswith(a[start:]):
            return len(a) - start

        start += 1


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

    overlaps = {}
    edges = 0

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

            if olen > 0:
                overlaps[(a, b)] = olen
                edges += 1

    return overlaps, edges


reads = readFastq('ERR266411_1.for_asm.fastq')

overlaps, edges = overlap_graph(reads, 30)

nodes_with_outgoing = set()

for a, b in overlaps:
    nodes_with_outgoing.add(a)

print(len(nodes_with_outgoing))

print(edges)
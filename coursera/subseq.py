import bisect

class SubseqIndex(object):

    def __init__(self, t, k, ival):
        self.k = k
        self.ival = ival
        self.index = []
        self.span = 1 + ival * (k - 1)

        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i+self.span:ival], i))

        self.index.sort()

    def query(self, p):
        subseq = p[:self.span:self.ival]

        i = bisect.bisect_left(self.index, (subseq, -1))

        hits = []

        while i < len(self.index):
            if self.index[i][0] != subseq:
                break

            hits.append(self.index[i][1])
            i += 1

        return hits


def readGenome(filename):
    genome = ''

    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()

    return genome


def approximate_match_subseq(p, t, index, n):
    occurrences = set()
    hits = 0

    for start in range(index.ival):

        # query subsequence partition
        sub_p = p[start:]

        matches = index.query(sub_p)

        hits += len(matches)

        for m in matches:

            # align full pattern position
            pos = m - start

            if pos < 0:
                continue

            if pos + len(p) > len(t):
                continue

            mismatches = 0

            for i in range(len(p)):
                if p[i] != t[pos + i]:
                    mismatches += 1

                    if mismatches > n:
                        break

            if mismatches <= n:
                occurrences.add(pos)

    return list(occurrences), hits


p = 'GGCGCGGTGGCTCACGCCTGTAAT'

t = readGenome('chr1.GRCh38.excerpt.fasta')

index = SubseqIndex(t, 8, 3)

occurrences, hits = approximate_match_subseq(p, t, index, 2)

print("Occurrences:", occurrences)
print("Total index hits:", hits)
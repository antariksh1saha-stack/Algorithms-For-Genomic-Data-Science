import matplotlib.pyplot as plt
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
    return sequences, qualities

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive(p,t):
    occurrences = []
    for i in range(len(t)-len(p)+1):
        match  = True
        for j in range(len(p)):
            if not t[i+j]==p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    p_prime = reverseComplement(p)

    for i in range(len(t)-len(p_prime)+1):
        match  = True
        for j in range(len(p_prime)):
            if not t[i+j]==p_prime[j]:
                match = False
                break
        if match and i not in occurrences:
            occurrences.append(i)
    return occurrences

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

def Phredd33toQ(qual):
    return ord(qual)-33



seqs,quals=readFastq('ERR037900_1.first1000.fastq')
hist = [0]*100
for qual in quals:
    i=0
    for q in qual:
        q=Phredd33toQ(q)
        hist[i] +=q
        i+=1



lowest = hist[0]
i=0
for val in hist:
    if val < lowest:
        lowest = val

for index in range(len(hist)):
    if hist[index]==lowest:
        break
print(index)
plt.bar(range(len(hist)),hist)
plt.show()
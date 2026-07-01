import matplotlib.pyplot as plt
import collections
def readFastQ(filename):
    sequences = []
    qualities = []
    with open(filename,'r') as f:
        while True:
            f.readline()
            seq= f.readline().rstrip()
            if len(seq) == 0:
                break
            f.readline()
            qual=f.readline().rstrip()
            sequences.append(seq)
            qualities.append(qual)
    return sequences,qualities

seq,qual=readFastQ('SRR835775_1.first1000.fastq')

def Phred33toQ(qual):
    return ord(qual)-33

def createHist(qualities):
    hist=[0]*50
    for qual in qualities:
        for phred in qual:
            q=Phred33toQ(phred)
            hist[q] +=1
    return hist




def findGC(reads):
    gc=[0]*100
    totals = [0]*100
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] +=1
            totals[i]+=1
    for i in range(len(gc)):
        if totals[i]>0:
            gc[i]=gc[i]/float(totals[i])
    return gc

count = collections.Counter()
for s in seq:
    count.update(s)
print(count)
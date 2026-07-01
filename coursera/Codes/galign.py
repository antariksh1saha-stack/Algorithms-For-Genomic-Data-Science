alphabet = ['A','C','G','T']
score = [[0,4,2,4,8],\
         [4,0,4,2,8],\
         [2,4,0,4,8],\
         [4,2,4,0,8],\
         [8,8,8,8,8]]

def editdist(x,y):
    D=[]
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(1,len(x)+1):
        D[i][0]=D[i-1][0] + score[alphabet.index(x[i-1])][-1]
    for i in range(1,len(y)+1):
        D[0][i]=D[0][i-1] + score[-1][alphabet.index(y[i-1])]
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            distHor=D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            distVer=D[i-1][j] + score[alphabet.index(x[j-1])][-1]
            if x[i-1]==y[j-1]:
                distDiag=D[i-1][j-1]
            else:
                distDiag=D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            D[i][j]=min(distDiag,distHor,distVer)
    return D[-1][-1]

x='TACGGATGC'
y='GCTAGCTAG'
D=editdist(x,y)
for i in range(len(D)):
    print(D[i])




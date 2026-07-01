def editdist(x,y):
    D=[]
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(len(x)+1):
        D[i][0]=i
    for i in range(len(y)+1):
        D[0][i]=0
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            distHor=D[i][j-1] + 1
            distVer=D[i-1][j] + 1
            if x[i-1]==y[j-1]:
                distDiag=D[i-1][j-1]
            else:
                distDiag=D[i-1][j-1] + 1
            D[i][j]=min(distDiag,distHor,distVer)
    return D

x='shake spea'
y='Shakespear'
D=editdist(x,y)
for i in range(len(x)+1):
    print(D[i])

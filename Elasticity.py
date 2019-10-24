import numpy as np


ni = 2
nj = 2

dx = 1.0 / ni
dy = 1.0 / nj

E = 1.0
v = 0.2

area = dx*dy*0.5

qtdNodes = (ni+1)*(nj+1)

neq = 3*qtdNodes

coords = np.zeros((3, qtdNodes))

border = np.zeros(qtdNodes, dtype="bool")

inode = 0
for j in range(nj+1):
    for i in range(ni+1):

        coords[0, inode] = i * dx
        coords[1, inode] = j * dy
        coords[2, inode] = 0.0

        if i == 0 or j == 0 or i == ni or j == nj:
            border[inode] = True

        inode += 1

triangles = []



iel = 0
for j in range(nj):
    for i in range(ni):

        inode = iel + j

        triangles.append((inode, inode+ni+1, inode+1))
        triangles.append((inode+1, inode+ni+1, inode+ni+2))


        iel+=1

M = np.zeros((3, 3))
deriv = np.zeros((3, 3))
B = np.zeros((6, 9))
D = np.zeros((6, 6))

Kg = np.zeros((neq, neq))


for i in range(3):
    for j in range(3):
        D[i,j] = 1 - v if i == j else v

for i in range(3, 6):
    D[i, i] = (1-2*v)/2

D = (E / ((1+v)*(1-2*v))) * D

np.savetxt("D.txt", D, delimiter="\t")

for t in triangles:
    M[0, 0] = 1.0
    M[1, 0] = 1.0
    M[2, 0] = 1.0

    M[0, 1] = coords[0, t[0]]
    M[1, 1] = coords[0, t[1]]
    M[2, 1] = coords[0, t[2]]

    M[0, 2] = coords[1, t[0]]
    M[1, 2] = coords[1, t[1]]
    M[2, 2] = coords[1, t[2]]

    print(M)

    Minv = np.linalg.inv(M)

    deriv[0, 0] = Minv[1, 0]
    deriv[1, 0] = Minv[2, 0]

    deriv[0, 1] = Minv[1, 1]
    deriv[1, 1] = Minv[2, 1]

    deriv[0, 2] = Minv[1, 2]
    deriv[1, 2] = Minv[2, 2]


    deriv[2, 0] = 0.0
    deriv[2, 1] = 0.0
    deriv[2, 2] = 0.0

    for i in range(3):

        B[0,0 + i*3] = deriv[0, i]
        B[1,1 + i*3] = deriv[1, i]
        B[2,2 + i*3] = deriv[2, i]

        B[3,0 + i*3] = deriv[1, i]
        B[3,1 + i*3] = deriv[0, i]

        B[4,0 + i*3] = deriv[2, i]
        B[4,2 + i*3] = deriv[0, i]

        B[5,1 + i*3] = deriv[2, i]
        B[5,2 + i*3] = deriv[1, i]

    np.savetxt("B.txt", B, delimiter="\t")

    Bt = np.transpose(B)
    K = np.dot(Bt, D)
    K = np.dot(K, B) * area

    np.savetxt("K.txt", K, delimiter="\t")

    g = np.zeros(9, dtype="int")

    for i in range(3):
        g[3*i+0] = 3*t[i]+0
        g[3*i+1] = 3*t[i]+1
        g[3*i+2] = 3*t[i]+2


    for j in range(9):
        for i in range(9):
            Kg[g[i], g[j]] += K[i, j]


for i in range(len(border)):

    if border[i]:

        Kg[3*i+0, :] = 0
        Kg[3*i+1, :] = 0
        Kg[3*i+2, :] = 0

        Kg[3*i+0, 3*i+0] = 1
        Kg[3*i+1, 3*i+1] = 1
        Kg[3*i+2, 3*i+2] = 1

np.savetxt("Kg.txt", Kg, delimiter="\t")

rhs = np.zeros(neq)

inode = 0
for j in range(nj+1):
    for i in range(ni+1):
        if border[inode]:
            rhs[3*inode] = coords[1, inode]

        inode += 1


sol = np.linalg.solve(Kg, rhs)

np.savetxt("solx.txt", rhs[0::3])

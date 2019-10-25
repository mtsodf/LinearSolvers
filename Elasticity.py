import matplotlib.pyplot as plt
import numpy as np
import argparse




# Definicoes dos parametros para submissao
parser = argparse.ArgumentParser(description='Problema de Elasticidade com triangulos')
parser.add_argument('--plot', action='store_true', help='Plota os arquivos de saida')
parser.add_argument('--show', action='store_true', help='Mostra plot')
parser.add_argument('--txt', action='store_true', help='Guarda matriz em txt')

parser.add_argument('--ni', default=2, type=int, help="Tamanho do grid em elementos ni")
parser.add_argument('--nj', default=2, type=int, help="Tamanho do grid em elementos nj")



solanaliticfunction = lambda x, y, z: [y, 0.0, 0.0]
f = lambda x, y, z: [0.0, 0.0, 0.0]

args = parser.parse_args()

ni = args.ni
nj = args.nj

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


if args.plot:
    fig, ax = plt.subplots(1, 1)

iel = 0
for j in range(nj):
    for i in range(ni):

        inode = iel + j

        triangles.append((inode, inode+ni+1, inode+1))
        triangles.append((inode+1, inode+ni+1, inode+ni+2))

        if args.plot:
            closed = list(triangles[-2])
            closed.append(triangles[-2][0])
            ax.plot(coords[0, closed], coords[1, closed])

            closed = list(triangles[-1])
            closed.append(triangles[-1][0])
            ax.plot(coords[0, closed], coords[1, closed])


        iel+=1

if args.plot:
    plt.show()
    plt.close()

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

if args.txt:
    np.savetxt("D.txt", D, delimiter="\t")

for it, t in enumerate(triangles):
    M[0, 0] = 1.0
    M[1, 0] = 1.0
    M[2, 0] = 1.0

    M[0, 1] = coords[0, t[0]]
    M[1, 1] = coords[0, t[1]]
    M[2, 1] = coords[0, t[2]]

    M[0, 2] = coords[1, t[0]]
    M[1, 2] = coords[1, t[1]]
    M[2, 2] = coords[1, t[2]]

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

    if args.txt:
        np.savetxt("B_%d.txt" % it, B, delimiter="\t")

    Bt = np.transpose(B)
    K = np.dot(Bt, D)
    K = np.dot(K, B) * area

    if args.txt:
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

if args.txt:
    np.savetxt("Kg.txt", Kg, delimiter="\t")

rhs = np.zeros(neq)

inode = 0
for j in range(nj+1):
    for i in range(ni+1):

        if border[inode]:
            rhs[3*inode] = coords[1, inode]

        else:
            x, y, z = coords[0, inode], coords[1, inode], coords[2, inode]
            rhs[3*inode], rhs[3*inode+1], rhs[3*inode+2] = f(x, y, z)

        inode += 1


if args.txt:
    np.savetxt("rhsx.txt", rhs[0::3])

sol = np.linalg.solve(Kg, rhs)

if args.txt:
    np.savetxt("solx.txt", sol[0::3])



solanalitic = np.zeros(neq)

for inode in range(qtdNodes):
    x, y, z = coords[0, inode], coords[1, inode], coords[2, inode]

    solanalitic[3*inode], solanalitic[3*inode+1], solanalitic[3*inode+2] = solanaliticfunction(x, y, z)



print("Diferenca entre encontrada e analitica = %e" % (np.linalg.norm(sol - solanalitic)))

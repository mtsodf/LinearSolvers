import matplotlib.pyplot as plt
import numpy as np
import argparse
from math import *
from matplotlib.colors import LinearSegmentedColormap
from scipy.sparse.linalg import gmres
from scipy.sparse import coo_matrix

colors = [(0.000, 0.000, 0.559), (0.000, 0.000, 0.766), (0.000, 0.000, 0.973), (0.000, 0.184, 0.996), (0.000, 0.391, 0.996), (0.000, 0.598, 0.996), (0.000, 0.805, 0.996), (0.012, 0.996, 0.984), (0.219, 0.996, 0.777), (0.426, 0.996, 0.570), (0.633, 0.996, 0.363), (0.840, 0.996, 0.156), (0.996, 0.945, 0.000), (0.996, 0.742, 0.000), (0.996, 0.535, 0.000), (0.996, 0.328, 0.000), (0.996, 0.121, 0.000), (0.910, 0.000, 0.000), (0.703, 0.000, 0.000), (0.500, 0.000, 0.000)]
cm_castelletto = LinearSegmentedColormap.from_list('mateus', colors, N=20)


def ElasticMatrix(E, v):
    D = np.zeros((6, 6))

    for i in range(3):
        for j in range(3):
            D[i,j] = 1 - v if i == j else v

    for i in range(3, 6):
        D[i, i] = (1-2*v)/2

    D = (E / ((1+v)*(1-2*v))) * D

    return D

def RunCase(ni, nj, E, v, verbose=True):

    if args.case == 0:
        solfunc = lambda x, y, z: [y, 0.0, 0.0]
        ffunc = lambda x, y, z: [0.0, 0.0, 0.0]

    elif args.case == 1:

        def ffunc(x, y, z):
            f = [None, None, None]
            f[0] = pi*pi*E*(-4*v + 3)*sin(pi*x)*sin(pi*y)/(2*(v + 1)*(2*v - 1))
            f[1] = -pi*pi*E*cos(pi*x)*cos(pi*y)/(2*(v + 1)*(2*v - 1))
            f[2] = 0

            return f

        def solfunc(x, y, z):
            return [sin(pi*x)*sin(pi*y), 0, 0]

    elif args.case == 2:

        def ffunc(x, y, z):
            f = [None, None, None]
            f[0] = -pi**2*E*cos(pi*x)*cos(pi*y)/(2*(v + 1)*(2*v - 1))
            f[1] =  pi**2*E*(3 - 4*v)*sin(pi*x)*sin(pi*y)/(2*(v + 1)*(2*v - 1))
            f[2] =  0

            return f

        def solfunc(x, y, z):
            return [0, sin(pi*x)*sin(pi*y), 0]

    elif args.case == 3:

        def ffunc(x, y, z):
            f = [None, None, None]
            f[0] = -pi**2*E*(v - 1)*sin(pi*x)/((v + 1)*(2*v - 1))
            f[1] = 0
            f[2] = 0

            return f

        def solfunc(x, y, z):
            return [sin(pi*x), 0, 0]

    elif args.case == 4:

        def ffunc(x, y, z):
            f = [None, None, None]
            f[0] = 0
            f[1] = 0
            f[2] = 0

            return f

        def solfunc(x, y, z):
            if y <= 1e-15:
                return [1 - x/args.li, 0, 0]
            if x <= 1e-15:
                return [1 - y/args.lj, 0, 0]
            return [0, 0, 0]

        E = None
        v = None

        youngElem = np.zeros(ni*nj) + 1
        poissonElem = np.zeros(ni*nj) + 0.2

    elif args.case == 5:

        def ffunc(x, y, z):
            f = [None, None, None]
            f[0] = 0
            f[1] = 0
            f[2] = 0

            return f

        def solfunc(x, y, z):
            if y <= 1e-15:
                return [1 - x/args.li, 0, 0]
            if x <= 1e-15:
                return [1 - y/args.lj, 0, 0]
            return [0, 0, 0]

        E = None
        v = None

        youngElem = np.loadtxt("young.txt")
        poissonElem = np.zeros(ni*nj) + 0.2

        ni = 8
        nj = 8


    dx = args.li / ni
    dy = args.lj / nj

    area = dx*dy*0.5

    nelem = ni*nj
    qtdNodes = (ni+1)*(nj+1)

    neq = 3*qtdNodes

    if verbose:
        print("------------------------------------------")
        print("            Dados Utilizados              ")
        print("------------------------------------------")
        print("      ni = %6d  nj = %6d" % (ni, nj))
        print("      dx = %6.3f  dy = %6.3f" % (dx, dy))
        print("------------------------------------------")
        print("")


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


    if args.plot and False:
        fig, ax = plt.subplots(1, 1)

    iel = 0

    young = []
    poisson = []

    for j in range(nj):
        for i in range(ni):

            inode = iel + j

            triangles.append((inode, inode+ni+1, inode+1))
            triangles.append((inode+1, inode+ni+1, inode+ni+2))

            young.append(youngElem[iel] if E is None else E)
            young.append(youngElem[iel] if E is None else E)

            poisson.append(poissonElem[iel] if v is None else v)
            poisson.append(poissonElem[iel] if v is None else v)

            if args.plot and False:
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


    rhs = np.zeros(neq)

    # Vetores da matriz de rigidez
    row = []
    col = []
    data = []

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

        xc = coords[0, t[0]] + coords[0, t[1]] + coords[0, t[2]]
        yc = coords[1, t[0]] + coords[1, t[1]] + coords[1, t[2]]
        zc = coords[2, t[0]] + coords[2, t[1]] + coords[2, t[2]]

        xc /= 3
        yc /= 3
        zc /= 3

        fvalue = ffunc(xc, yc, zc)

        for inode in t:
            rhs[3*inode]   -= area * fvalue[0] / 3.0
            rhs[3*inode+1] -= area * fvalue[1] / 3.0
            rhs[3*inode+2] -= area * fvalue[2] / 3.0


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

        D = ElasticMatrix(young[it], poisson[it])

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
                inode = g[i] // 3
                if not border[inode]:
                    #Kg[g[i], g[j]] += K[i, j]
                    row.append(g[i])
                    col.append(g[j])
                    data.append(K[i, j])


    for i in range(len(border)):

        if border[i]:

            for offset in range(3):
                row.append(3*i+offset)
                col.append(3*i+offset)
                data.append(1)


    Kg = coo_matrix((data, (row, col)), shape=(neq, neq))

    if args.spy:
        fig, ax = plt.subplots(1, 1)
        ax.spy(Kg)
        plt.savefig("Kg_spy_%d_%d.png" % (ni, nj))
        plt.close()


    if args.txt:
        np.savetxt("Kg.txt", Kg.toarray(), delimiter="\t")


    # Montagem do lado direito para condicoes de Dirichlet
    inode = 0
    for j in range(nj+1):
        for i in range(ni+1):

            x, y, z = coords[0, inode], coords[1, inode], coords[2, inode]

            if border[inode]:
                rhs[3*inode], rhs[3*inode+1], rhs[3*inode+2] = solfunc(x, y, z)

            inode += 1


    if args.txt:
        np.savetxt("rhsx.txt", rhs[0::3])
        np.savetxt("rhs.txt", rhs)


    sol, exitCode = gmres(Kg, rhs)

    if exitCode != 0:
        print("Gmres não convergiu.")

    if args.txt:
        np.savetxt("solx.txt", sol[0::3])
        np.savetxt("sol.txt", sol)



    solanalitic = np.zeros(neq)

    for inode in range(qtdNodes):
        x, y, z = coords[0, inode], coords[1, inode], coords[2, inode]

        solanalitic[3*inode], solanalitic[3*inode+1], solanalitic[3*inode+2] = solfunc(x, y, z)

    errorAbs = np.linalg.norm(sol - solanalitic)
    errorRel = np.linalg.norm(sol - solanalitic) / np.linalg.norm(solanalitic)

    if verbose:
        print("Diferenca entre encontrada e analitica = %e" % (errorAbs))
        print("Diferenca entre encontrada e analitica = %e" % (errorRel))


    if args.plot:
        fig, axs = plt.subplots(3, 2)

        for i in range(3):
            im = axs[i,0].imshow(sol[i::3].reshape((nj+1, ni+1)), interpolation="bilinear", origin='lower', cmap=cm_castelletto)
            fig.colorbar(im, ax=axs[i,0]).ax.tick_params(labelsize=12)

            im = axs[i,1].imshow(solanalitic[i::3].reshape((nj+1, ni+1)), interpolation="bilinear", origin='lower', cmap=cm_castelletto)
            fig.colorbar(im, ax=axs[i,1]).ax.tick_params(labelsize=12)

        plt.show()
        plt.close()

    return errorAbs, errorRel



# Definicoes dos parametros para submissao
parser = argparse.ArgumentParser(description='Problema de Elasticidade com triangulos')
parser.add_argument('--plot', action='store_true', help='Plota os arquivos de saida')
parser.add_argument('--show', action='store_true', help='Mostra plot')
parser.add_argument('--case', '-c', default=0, type=int, help='Numero do caso')
parser.add_argument('--txt', action='store_true', help='Guarda matriz em txt')
parser.add_argument('--spy', action='store_true', help='Salvar Arquivo Spy')
parser.add_argument('--error', action='store_true', help='Roda o caso para diferentes tamanhos de malha')

parser.add_argument('--ni', default=2, type=int, help="Tamanho do grid em elementos ni")
parser.add_argument('--nj', default=2, type=int, help="Tamanho do grid em elementos nj")

parser.add_argument('--li', default=1.0, type=int, help="Tamanho em metros do grid na direcao i")
parser.add_argument('--lj', default=1.0, type=int, help="Tamanho em metros do grid na direcao j")
args = parser.parse_args()


E = 1.0
v = 0.2


ni = args.ni
nj = args.nj

if args.error:
    print("")
    print("-"*42)
    print("| %10s | %10s  | %10s  |" % ("N", "erro abs.", "erro rel."))
    print("-"*42)

    ns = [5, 10, 100, 120]
    errors = []

    for n in ns:
        errorAbs, errorRel = RunCase(n, n , E, v, verbose=False)
        errors.append(errorRel)
        print("| %10d | %10.5e | %10.5e |" % (n, errorAbs, errorRel))
    
    print("-"*42)

    fig, ax = plt.subplots(1, 1)

    ax.loglog(ns, errors)
    ax.loglog([5, 50], [0.3, 0.3/100], ls="--", c="black")
    plt.show()

else:
    RunCase(ni, nj, E, v)
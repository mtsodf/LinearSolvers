import numpy as np


def jabobi(A, B, x0):

    size = A.shape[0]

    D = np.zeros(shape=[size, size])
    Dinv = np.zeros(shape=[size, size])
    for i in range(size):
        D[i][i] = A[i][i]
        Dinv[i][i] = 1/A[i][i]

    LU = A - D
    erro = 10.0
    i = 0
    LD = np.dot(Dinv, B)
    update = np.dot(Dinv, LU)
    while erro > 1e-6 and i < 10000:
        x = np.dot(update, x0)
        x0 = LD - x
        Bcalc = np.dot(A, x0)
        erro_v = Bcalc - B
        erro = np.dot(np.transpose(erro_v), erro_v)
        i += 1

    return x0


if __name__ == '__main__':
    A = np.array([1, 1, 1, 1, 2, 1, 3, 0, 10]).reshape(3, 3)

    A = A / 1.0
    print A
    B = np.array([10.0, 15.0, 36.0]).reshape(3, 1)
    print B
    x0 = np.zeros(shape=[3, 1])

    jabobi(A, B, x0)

    inv = np.linalg.inv(A)
    print A
    print inv
    res = np.dot(inv, B)
    print "Resposta "
    print res
    print np.dot(A, res)
    print np.dot(A, inv)

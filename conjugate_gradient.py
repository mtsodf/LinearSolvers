
import numpy as np
import matplotlib.pyplot as plt


def conjugate_method(A, B, x0):
    r = B - np.dot(A, x0)

    error = np.dot(np.transpose(r), r)
    i = 0
    xs = [x0[0]]
    ys = [x0[1]]

    while error > 1e-15 and i < 100000:

        lambda_ = error

        aux = np.dot(A, r)
        aux = np.dot(np.transpose(r), aux)

        lambda_ /= aux

        x0 = x0 + lambda_ * r
        xs.append(x0[0])
        ys.append(x0[1])

        r = B - np.dot(A, x0)

        error = np.dot(np.transpose(r), r)
        i += 1

    print "Numero de iteracoes: %d" % i

    xs = np.array(xs)
    ys = np.array(ys)
    fig, ax = plt.subplots()
    plt.plot(xs, ys, 'r.-',)
    ax.grid()
    ax.axis('equal')
    plt.show()

    return x0


def conjugate_gradient(A, B, x0):
    r = B - np.dot(A, x0)

    error = np.dot(np.transpose(r), r)
    d = r
    i = 0
    xs = [x0[0]]
    ys = [x0[1]]

    while error > 1e-15 and i < 100000:

        lambda_ = error

        aux = np.dot(A, d)
        aux = np.dot(np.transpose(d), aux)

        lambda_ /= aux

        x0 = x0 + lambda_ * d
        xs.append(x0[0])
        ys.append(x0[1])

        rnext = B - np.dot(A, x0)

        beta = np.dot(np.transpose(rnext), rnext) / np.dot(np.transpose(r), r)

        r = rnext
        d = r + beta * d
        error = np.dot(np.transpose(r), r)

        i += 1

    print "Numero de iteracoes: %d" % i

    xs = np.array(xs)
    ys = np.array(ys)
    fig, ax = plt.subplots()
    plt.plot(xs, ys, 'r.-',)
    ax.grid()
    ax.axis('equal')
    plt.show()

    return x0


if __name__ == '__main__':
    A = np.array([4, 1, 1, 3]).reshape(2, 2)

    A = A / 1.0
    print A
    B = np.array([5.0, 4.0]).reshape(2, 1)
    print B
    x0 = np.array([10.0, 20.0]).reshape(2, 1)

    res = conjugate_gradient(A, B, x0)

    print res

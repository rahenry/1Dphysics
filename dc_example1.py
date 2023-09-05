#!/usr/bin/env python3
import numpy, scipy

# helper function to find local minima in a 2D array
def find_minima(z):
    X, Y = z.shape
    res = []
    for i in range(X-2):
        x = i+1
        for j in range(Y-2):
            y = j+1
            u = z[x][y+1]
            d = z[x][y-1]
            l = z[x-1][y]
            r = z[x+1][y]
            t = z[x][y]
            things = [-1, 0, 1]
            accepted = 1
            for ox in things:
                for oy in things:
                    if oy == ox == 0:
                        continue
                    if t > z[x+ox][y+oy]:
                        accepted = 0

            if accepted:
                res.append((x, y))
    return res

# helper function to make the product of two arrays
def cartesian_product(*arrays):
    la = len(arrays)
    dtype = numpy.result_type(*arrays)
    arr = numpy.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(numpy.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

# this is the function that must be minimised to find the values of k which correspond to EPs
def k_expression(k, L, D=1):
    return (L+D) * numpy.cos((L+D)*k) * numpy.sin(L*k) - L * numpy.sin((L+D)*k)*numpy.cos(L*k)

def abs_k_expression(k, L, D=1):
    return abs(k_expression(k, L, D))


# D=1 for Ising/FP, D=2 for XY
# Get the EPs for some values of L and N
# Searches between the specified k values
def get_eps(L, N, k_real_min, k_real_max, k_imag_min, k_imag_max, npoints=501, D=1):
    m1 = (k_real_max - k_real_min) / (npoints-1)
    m2 = (k_imag_max - k_imag_min) / (npoints-1)
    def kfunc(i, j):
        # divide the imaginary part by L because the EPs approach the real axis as L increases
        return k_real_min + m1 * i + 1.j * (k_imag_min + j * m2) / L

    k_grid = numpy.fromfunction(kfunc, (npoints, npoints))
    def z_func(k):
        return abs_k_expression(k, L, D)

    z_grid = z_func(k_grid)
    minima = find_minima(z_grid)

    # the minimize function needs a real function of two variables, rather than a complex function of one variable
    def min_func(x):
        return z_func(x[0] + 1.j * x[1])
    dx = 2 * m1
    dy = 2 * m2
    res = []
    for i, j in minima:
        k = kfunc(i, j)
        bounds = [(k.real - dx, k.real + dx), (k.imag - dy, k.imag + dy)]
        meth = "Nelder-Mead"
        m = scipy.optimize.minimize(min_func, (k.real, k.imag), method=meth, bounds=bounds, tol=1E-26)
        kimproved = m.x[0] + 1.j * m.x[1]
        res.append(kimproved)
        print(f"Found an EP:\nk={kimproved}\n|k_exp|={z_func(kimproved)}\n")
    return res


if __name__ == "__main__":
    eps = get_eps(5, 3, -numpy.pi, numpy.pi, -numpy.pi, numpy.pi, npoints=500)

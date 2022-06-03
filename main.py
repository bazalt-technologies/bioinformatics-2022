import numpy as np
import time
import matplotlib.pyplot as plt


# Note: mathematical article is not written yet.
# "Formula (x.x)" expression means that this formula is numbered in the article with an x.x number
# See README.md for mathematics


# Returns random matrix W of {-1, 0, 1} | Formula (x.x)
# Probability of element to be 1 or -1: k/(2*N), to be 0: 1 - k/N
# RETURNS NOT A NUMPY.MATRIX
def W_calc(k, N, M):
    W = []
    for i in range(N):
        # Making a new row
        W.append([])
        # Fill a new row
        for j in range(M):
            # Take a random number [0; 100] for emulating probability
            p = np.random.randint(0, 10001) / 100
            if p < 100 * k / N:
                if p < (k / (2 * N)) * 100:
                    W[i].append(1)
                else:
                    W[i].append(-1)
            else:
                W[i].append(0)
    return W


# Sigma function for phenotype calculation | Formula (x.x)
def sigma_calc(z):
    return 1 / (1 + np.exp(z))


# Calculate phenotype: vector of [0; 1] floats | Formula (x.x)
def f_calc(M, N, W, S, h):
    f_wektor = []
    # Filling vector
    for i in range(M):
        temp = 0
        for j in range(N):
            # W.item is a W_ji element of matrix, because W is a numpy.matrix type here
            temp += W.item(j, i) * S[j]
        temp -= h
        f_wektor.append(sigma_calc(temp))
    return f_wektor


# Calculate W_tr (Fitness Potential) | Formula (x.x)
def w_tr_calc(C, f, B, M):
    # Calculate lineal part of W_tr
    lineal = 0
    for i in range(M):
        lineal += C[i] * f[i]
    # Calculate double sum part
    double_sum = 0
    for i in range(M):
        for j in range(M):
            double_sum += B[i][j] * f[i] * f[j]
    return lineal + double_sum / 2


# mutate takes a genotype and mutates it with the following mechanism
# Every coordinate of a vector changes to opposite with probability of p_mut
def mutate(S, p_mut):
    for i in range(len(S)):
        # Emulating probability
        p = np.random.randint(0, 101)
        if p < p_mut:
            if S[i]:
                S[i] = 0
            else:
                S[i] = 1


def d_alg(M, N, h, p_mut, S, C, B, W, f, iterations):
    S_new = []
    for i in range(N):
        S_new.append(S[i])

    for t in range(iterations):
        # Mutate S_new
        mutate(S_new, p_mut)
        # Calculates a phenotype of mutated generation
        f_new = f_calc(M, N, W, S_new, h)
        # If Fitness Potential of mutated generation is more than of previous ...
        if w_tr_calc(C, f, B, M) < w_tr_calc(C, f_new, B, M):
            # ... replace old genotype and phenotype with a new one ...
            for i in range(N):
                S[i] = S_new[i]
            f = f_new
        else:
            # ... else keep the old generation
            S_new = S
    return w_tr_calc(C, f, B, M)


def from_p_mut():
    # M - Size of phenotype vector
    # N - Size of genotype vector

    iterations = 500
    # -- Looking for dependence of Fitness Potential from p_mut --
    # Consts in this case
    M_main, N_main, k_main, h_main = 20, 40, 4, 2
    # Calculate S vector | Formula (x.x)
    # It is genotype of a first generation
    S = []
    for i in range(N_main):
        S.append(np.random.randint(0, 2))
    # Calculate C vector | Formula (x.x)
    C = []
    for i in range(M_main):
        C.append(np.random.randint(0, 2))
    # Calculate **symmetric** B matrix | Formula (x.x)
    B_asymmetric = np.random.randint(0, 2, size=(M_main, M_main), dtype='int')
    B = (B_asymmetric + B_asymmetric.T) % 2
    # Calculate W matrix | Formula (x.x)
    W = np.matrix(W_calc(k_main, N_main, M_main))
    # Calculate f vector | Formula (x.x)
    # It is phenotype of a first generation
    f = f_calc(M_main, N_main, W, S, h_main)
    # S_new is a genotype of generation
    # At first equal to a current generation, then mutates and
    # Either becomes a genotype of a new generation
    # Or dies and turns back into previous
    # -- p_mut, alpha, t --
    # w_arr is the array of tuples (w_tr, p_mut)
    # where w_tr is a Fitness Potential achieved with probability p_mut
    w0 = w_tr_calc(C, f, B, M_main)
    fitness_c = 0.1/w0
    w_arr = []
    p_mut_main = 1
    print(f'{M_main = }, {N_main = }, {k_main = }, {h_main = }')
    while p_mut_main <= 40:
        start = time.time()
        output = []
        for kk in range(5):
            output.append(d_alg(M_main, N_main, h_main, p_mut_main, S, C, B, W, f, iterations))
        end = time.time()
        ou = 0
        for i in range(len(output)):
            ou += output[i]
        ou /= len(output)
        print(f'w_tr = {ou} {p_mut_main = } time: {(end - start)}s')
        w_arr.append((ou, p_mut_main))
        p_mut_main += 3
    for lol in range(1, len(w_arr) - 1):
        try:
            if w_arr[lol - 1][0] - w_arr[lol][0] > 0.7 and w_arr[lol + 1][0] - w_arr[lol][0] > 0.7 \
                    or w_arr[lol - 1][0] - w_arr[lol][0] < -0.7 and w_arr[lol + 1][0] - w_arr[lol][0] < -0.7:
                w_arr.pop(lol)
        except IndexError:
            break
    figure, ax = plt.subplots(nrows=1, ncols=1)
    x = [a[1] for a in w_arr]
    y = [np.exp(fitness_c * a[0]) for a in w_arr]
    ax.plot(x, y)
    ax.set_ylabel("Fitness")
    ax.set_xlabel("Mutation Probability")
    titl = "F(p_mut) calculation\nF(0) = " + str(0.1)
    ax.set(title=titl)
    plt.show()

    print("Fitness Potentials:")
    w_arr.sort(reverse=True)
    print(w_arr)


def from_iterations():
    M_main, N_main, k_main, h_main = 10, 20, 4, 2
    p_mut_main = 35
    # Calculate S vector | Formula (x.x)
    # It is genotype of a first generation
    S = []
    for i in range(N_main):
        S.append(np.random.randint(0, 2))
    # Calculate C vector | Formula (x.x)
    C = []
    for i in range(M_main):
        C.append(np.random.randint(0, 2))
    # Calculate **symmetric** B matrix | Formula (x.x)
    B_asymmetric = np.random.randint(0, 2, size=(M_main, M_main), dtype='int')
    B = (B_asymmetric + B_asymmetric.T) % 2
    # Calculate W matrix | Formula (x.x)
    W = np.matrix(W_calc(k_main, N_main, M_main))
    # Calculate f vector | Formula (x.x)
    # It is phenotype of a first generation
    f = f_calc(M_main, N_main, W, S, h_main)
    # S_new is a genotype of generation
    # At first equal to a current generation, then mutates and
    # Either becomes a genotype of a new generation
    # Or dies and turns back into previous
    # -- p_mut, alpha, t --
    # w_arr is the array of tuples (w_tr, p_mut)
    # where w_tr is a Fitness Potential achieved with probability p_mut
    w0 = w_tr_calc(C, f, B, M_main)
    fitness_c = 0.1/w0
    w_arr = []
    iterations = 100
    print(f'{M_main = }, {N_main = }, {k_main = }, {h_main = }')
    global_start = time.time()
    while iterations <= 3500:
        start = time.time()
        output = []
        for kk in range(5):
            output.append(d_alg(M_main, N_main, h_main, p_mut_main, S, C, B, W, f, iterations))
        end = time.time()
        ou = 0
        for i in output:
            ou += i
        ou /= len(output)
        print(f'w_tr = {ou} {iterations = } time: {(end - start)}s')
        w_arr.append((ou, iterations))
        iterations += 300
    global_end = time.time()

    figure, ax = plt.subplots(nrows=1, ncols=1)
    x = [a[1] for a in w_arr]
    y = [np.exp(fitness_c*a[0]) for a in w_arr]
    ax.plot(x, y)
    ax.set_ylabel("Fitness")
    ax.set_xlabel("Iterations")
    titl = "F(iterations) calculation\nF(0) = " + str(0.1)
    ax.set(title=titl)
    plt.show()

    print("Fitness Potentials:")
    w_arr.sort(reverse=True)
    print(w_arr)
    print(f"Calculated in {global_end - global_start}s")


if __name__ == '__main__':
    from_p_mut()
    from_iterations()

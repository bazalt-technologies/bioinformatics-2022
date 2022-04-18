import numpy as np
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
            p = np.random.randint(0, 101)
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


def d_alg(M, N, k, h, p_mut):
    # Calculate S vector | Formula (x.x)
    # It is genotype of a first generation
    S = []
    for i in range(N):
        S.append(np.random.randint(0, 2))
    # Calculate C vector | Formula (x.x)
    C = []
    for i in range(M):
        C.append(np.random.randint(0, 2))
    # Calculate **symmetric** B matrix | Formula (x.x)
    B_asymmetric = np.random.randint(0, 2, size=(M, M), dtype='int')
    B = (B_asymmetric + B_asymmetric.T) % 2
    # Calculate W matrix | Formula (x.x)
    W = np.matrix(W_calc(k, N, M))
    # Calculate f vector | Formula (x.x)
    # It is phenotype of a first generation
    f = f_calc(M, N, W, S, h)
    # S_new is a genotype of generation
    # At first equal to a current generation, then mutates and
    # Either becomes a genotype of a new generation
    # Or dies and turns back into previous
    S_new = []
    for i in range(N):
        S_new.append(S[i])

    for t in range(1000):
        # Mutate S_new
        mutate(S_new, p_mut)
        # Calculates a phenotype of mutated generation
        f_new = f_calc(M, N, W, S_new, h)
        # If Fitness Potential of mutated generation is more than of previous ...
        if w_tr_calc(C, f, B, M) > w_tr_calc(C, f_new, B, M):
            # ... replace old genotype and phenotype with a new one ...
            for i in range(N):
                S[i] = S_new[i]
            f = f_new
        else:
            # ... else keep the old generation
            for i in range(N):
                S_new[i] = S[i]
    return w_tr_calc(C, f, B, M)


if __name__ == '__main__':
    # M - Size of phenotype vector
    # N - Size of genotype vector

    # -- Looking for dependence of Fitness Potential from p_mut --
    # Expecting p_mut to be high on average.

    # Consts in this case
    M_main, N_main, k_main, h_main = 5, 10, 5, 2

    # w_arr is the array of tuples (w_tr, p_mut)
    # where w_tr is a Fitness Potential achieved with probability p_mut
    w_arr = []

    p_mut_main = 5
    print(f'{M_main = }, {N_main = }, {k_main = }, {h_main = }')
    while p_mut_main <= 100:
        output = d_alg(M_main, N_main, k_main, h_main, p_mut_main)
        w_arr.append((output, p_mut_main))
        p_mut_main += 5

    # Showing 5 most successful runs
    w_arr.sort(reverse=True)
    w_arr = w_arr[0:5]
    w_arr.sort(key=lambda x: x[1], reverse=True)

    print("Biggest Fitness Potentials:")
    print(w_arr)

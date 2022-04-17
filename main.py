import numpy as np


def W_calc(k, N, M):
    W = []
    for i in range(N):
        W.append([])
        for j in range(M):
            p = np.random.randint(0, 101)
            if p < 100 * k / N:
                if p < (k / (2 * N)) * 100:
                    W[i].append(1)
                else:
                    W[i].append(-1)
            else:
                W[i].append(0)
    return W


def sigma_calc(z):
    return 1 / (1 + np.exp(z))


def f_calc(M, N, W, S, h):
    f_wektor = []
    for i in range(M):
        temp = 0
        for j in range(N):
            temp += W.item(j, i) * S[j]
        temp -= h
        f_wektor.append(sigma_calc(temp))
    return f_wektor


def w_tr_calc(C, f, B, M):
    lineal = 0
    for i in range(M):
        lineal += C[i] * f[i]
    double_sum = 0
    for i in range(M):
        for j in range(M):
            double_sum += B[i][j] * f[i] * f[j]
    return lineal + double_sum / 2


def mutate(S, p_mut):
    for i in range(len(S)):
        p = np.random.randint(0, 101)
        # Чем больше p_mut тем больше шансов мутировать
        if p < p_mut:
            if S[i]:
                S[i] = 0
            else:
                S[i] = 1


def d_alg(M, N, k, h, p_mut):
    # Высчитываем рандомный вектор S
    S = []
    for i in range(N):
        S.append(np.random.randint(0, 2))
    # Высчитываем рандомный вектор C
    C = []
    for i in range(M):
        C.append(np.random.randint(0, 2))
    # Высчитываем рандомную симметричную матрицу B
    B_asymmetric = np.random.randint(0, 2, size=(M, M), dtype='int')
    B = (B_asymmetric + B_asymmetric.T) % 2
    # Высчитываем рандомную матрицу W
    W = np.matrix(W_calc(k, N, M))
    # Высчитываем изначальный вектор фенотипа f
    f = f_calc(M, N, W, S, h)
    # S_new будет мутированный генотип с которым будем сравнивать старый
    S_new = []
    for i in range(N):
        S_new.append(S[i])

    for t in range(250):
        # Мутируем гены
        mutate(S_new, p_mut)
        # Высчитываем новый фенотип
        f_new = f_calc(M, N, W, S_new, h)
        # Сравниваем нужно ли мутировать
        if w_tr_calc(C, f, B, M) > w_tr_calc(C, f_new, B, M):
            for i in range(N):
                S[i] = S_new[i]
        else:
            for i in range(N):
                S_new[i] = S[i]
    return w_tr_calc(C, f, B, M)


if __name__ == '__main__':
    # Looking for dependence of F from p_mut

    # Consts in this case
    M_main, N_main, k_main, h_main = 5, 10, 5, 2
    p_mut_main = 10
    print(f'{M_main=}, {N_main=}, {k_main=}, {h_main=}')
    while p_mut_main <= 100:
        print(f'{p_mut_main=} Got:',
              d_alg(M_main, N_main, k_main, h_main, p_mut_main))
        p_mut_main += 10

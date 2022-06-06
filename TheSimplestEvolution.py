# Визуализация простейшего алгоритма эволюции

from random import random
import matplotlib.pyplot as plt

# Генерации случайной матрицы 2 на 2
def choice22Matr(K, N, M, q):
    p = K / N
    W = [[-1 if random() < q else 1 if random() < p else 0 for i in range(M)] for j in range(N)]
    barW = [] # Вспоиогательная величина
    for j in range(N):
        ss = 0
        for i in range(M):
            ss += W[i][j]
        barW.append(ss / M)

    return W, barW

# Вычисление начального вектора
def initGenePool(N):
    s0 = [1 if random() > 0.5 else -1 for _ in range(N)]
    return s0

# Мутация
def randFlip(gamma, N, s):
    snew = s
    rc = gamma / N
    snew = [-s[k] if random() < rc else s[k] for k in range(N)]
    return snew

# Вычисление фитнеса
def fitPotent(barW, N, s):
    Wpot = 0
    for j in range(N):
        Wpot += barW[j] * s[j]
    return Wpot

# Вычисление максимального фитнеса
def maxPotent(barW, N):
    maxW = 0
    for j in range(N):
        maxW += abs(barW[j])
    return maxW

# Простейшая эволюция
def verySimpleEvol(gamma, Tevol, M, N, barW, s0):
    maxW = maxPotent(barW, N)
    s = s0
    sold = s
    Wpotcur = []
    Time = []
    for t in range(Tevol):
        Wpot = fitPotent(barW, N, sold)
        Wpotcur.append(Wpot)
        Time.append(t)
        Wpotold = Wpot
        snew = randFlip(gamma, N, sold)
        Wpot = fitPotent(barW, N, snew)
        Wpotnew = Wpot
        if Wpotnew > Wpotold - 0.0001:
            sold = snew
    sfinal = sold

    return maxW, Wpotcur, Time, sfinal



def main():
    gamma = 1
    Tevol = 2000 # Граница исследования (Число шагов)
    M = 100 # Число признаков
    N = M # Число генов
    q = 0.1 # Вероятность мутации
    K = 10 # Генетическая регуляция
    W, barW = choice22Matr(K, N, M, q)
    s0 = initGenePool(N)
    maxW, Wpotcur, Time, sfinal, = verySimpleEvol(gamma, Tevol, M, N, barW, s0)
    print("MaxW:", maxW)
    figure, ax = plt.subplots(nrows=2, ncols=1)
    ax[0].plot(Wpotcur)
    textLabel = "Изменение фитнеса при простейшей эволюции при gamma = " + str(gamma) + "\nМаксимальный фитнес = " + str(maxW)
    ax[0].set(title = textLabel)
    
    ax[1].pie([s0.count(1), s0.count(-1)], radius=1, labels=["Изначальное кол-во векторов 1", "Изначальное кол-во векторов -1"])
    ax[1].pie([sfinal.count(1), sfinal.count(-1)], radius=0.5, labels=["Конечное кол-во векторов 1", "Конечное кол-во векторов -1"])
    textLabel = "Соотношение начального и конечного вектора при q = " + str(q)
    ax[1].set(title = textLabel)

    figure.set_size_inches(10, 10)
    plt.show()

main()
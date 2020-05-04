import cmath
import math
import random
import matplotlib.pyplot as plt
import time
import numpy as np
from cmath import exp, pi

N = 1024
n = 6
W = 2100

minimal_range = 0
maximum_range = 1


# Функція генератора стаціонарного випадкового сигналу
def stat_random_signal(n, W):
    A = [random.random() for _ in range(n)]
    phi = [random.random() for _ in range(n)]

    def f(t):
        x = 0
        for i in range(n):
            x += A[i]*math.sin(W / n * t * i + phi[i])
            #x *= (math.cos((2*math.pi/N)*i*t))
        return x
    return f


# Функція повертає значення мат. очікування
def get_m(x):
    return sum(x)/len(x)


# Функція повертає значення Дисперсії
def get_D(x, m=None):
    if m is None:
        m = get_m(x)
    return sum([(i - m) ** 2 for i in x]) / (len(x) - 1)


s_gen = stat_random_signal(n, W)
s = [s_gen(i) for i in range(N)]
dft = np.fft.fft(s)

# Функція для обчислення Дискретного перетворення Фур'є
def dft(x):
    t = []
    N_ = len(x)
    for k in range(N_):
        a = 0
        for n_ in range(N):
            a += x[n_]*cmath.exp(-2j*cmath.pi*k*n_*(1/N_))
        t.append(a)
    return t

# Функція для обчислення Дискретного перетворення Фур'є
def fft(x):
    N_ = len(x)
    if N_ <= 1: return x
    even = fft(x[0::2])
    odd =  fft(x[1::2])
    T= [exp(-2j * pi * k / N_) * odd[k] for k in range(N_ // 2)]
    return [even[k] + T[k] for k in range(N_//2)] + \
           [even[k] - T[k] for k in range(N_//2)]



#ddd = dft(s)
ddd = fft(s)
print(s)

m_start_time = time.time()
m = get_m(s)
m_end_time = time.time()
m_time = m_end_time - m_start_time


D_start = time.time()
D = get_D(s, m)
D_end = time.time()
D_time = D_end - D_start


print("Час розрахунку мат. очікування: ", m_time)
print("Час розрахунку дисперсії: ", D_time)
print("m =", m)
print('D =', D)

plt.plot(range(N), s)
plt.plot(range(N), np.array(ddd).real)


plt.show()


from math import pi
import numpy as np


def Sphere(x):
    x = np.asarray_chkfinite(x)
    return sum(x * np.transpose(x))


def Ackley(x, d=20, e=0.2, f=2 * pi):
    x = np.asarray_chkfinite(x)
    n = len(x)
    s1 = sum(x ** 2)
    s2 = sum(np.cos(f * x))
    return -d * np.exp(-e * np.sqrt(s1 / n)) - np.exp(s2 / n) + d + np.exp(1)


def SumSquare(x):
    x = np.asarray_chkfinite(x)
    for i in range(0, 3):
        SUM = sum(x * np.transpose(x) + i ** 2)
    return SUM

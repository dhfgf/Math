from methods import *


def f1(x):
    return x * x * x * x - 2 * x


def f2(x):
    return math.sin(x) * x * x


class Input:
    def __init__(self, f, s, a):
        self.f = f
        self.s = s
        self.a = a


in1 = Input(f2, Section(-10, 10), 0.001)

bisection(in1.f, in1.s, in1.a).print()
golden_section(in1.f, in1.s, in1.a).print()
fibonacci(in1.f, in1.s, in1.a).print()
parabola(in1.f, in1.s, in1.a).print()
brent(in1.f, in1.s, in1.a).print()

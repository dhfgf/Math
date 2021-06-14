import math
from math import sqrt, copysign


class Section:
    def __init__(self, a, b):
        self.begin = a
        self.end = b

    def length(self):
        return self.end - self.begin

    def print(self):
        print("[" + str(self.begin) + "; " + str(self.end) + "]")
        return

    def include(self, u):
        return (u <= self.end) and (u >= self.begin)


class Result:
    def __init__(self, iter, func, sect: Section, meth):
        self.iter_count = iter
        self.func_count = func
        self.sect_res = sect
        self.meth_name = meth

    def print(self):
        print(str(self.meth_name).upper())
        print("Section is [" + str(self.sect_res.begin) + "; " + str(self.sect_res.end) + "]")
        print("Number of iterations: " + str(self.iter_count))
        print("Number of function evaluations: " + str(self.func_count))
        print("Extremum: " + str((self.sect_res.begin + self.sect_res.end) / 2))
        print()


G_CONST = (3 - math.sqrt(5)) / 2


def bisection(function, section: Section, accuracy) -> Result:
    def bisection_iter(_function, _section: Section, _delta):
        mid = (_section.begin + _section.end) / 2
        x1 = mid - _delta
        x2 = mid + _delta
        f1 = _function(x1)
        f2 = _function(x2)
        if f1 < f2:
            return Section(_section.begin, x2)
        else:
            return Section(x1, _section.end)

    iter_c = 0
    func_c = 0
    delta = accuracy / 2 - accuracy * 0.1
    sect = section
    while sect.length() > accuracy:
        sect = bisection_iter(function, sect, delta)
        iter_c += 1
        func_c += 2
    return Result(iter_c, func_c, sect, "bisection")


def golden_section(function, section: Section, accuracy) -> Result:
    class ResConv:
        def __init__(self, s: Section, f, b):
            self.s = s
            self.f = f
            self.b = b

        def print(self):
            print(self.s.begin, self.s.end)
            print("\n")

    def golden_section_iter(_function, sect, f, f_was_left) -> ResConv:
        x1 = sect.begin + G_CONST * (sect.end - sect.begin)
        x2 = sect.end - G_CONST * (sect.end - sect.begin)
        if f_was_left:
            f2 = f
            f1 = _function(x1)
        else:
            f1 = f
            f2 = _function(x2)

        if f1 > f2:
            z = ResConv(Section(x1, sect.end), f2, False)
            return z
        else:
            z = ResConv(Section(sect.begin, x2), f1, True)
            return z

    iter_c = 0
    func_c = 0
    l = ResConv(section, 0, True)
    if section.length() > accuracy:
        x1 = section.begin + G_CONST * (section.end - section.begin)
        x2 = section.end - G_CONST * (section.end - section.begin)
        f1 = function(x1)
        f2 = function(x2)

        iter_c += 1
        func_c += 2

        if f1 > f2:
            l = ResConv(Section(x1, section.end), f2, False)
        else:
            l = ResConv(Section(section.begin, x2), f1, True)

        while l.s.length() > accuracy:
            l = golden_section_iter(function, l.s, l.f, l.b)
            iter_c += 1
            func_c += 1

    return Result(iter_c, func_c, l.s, "golden section")


def parabola(function, section: Section, accuracy) -> Result:
    class ret:
        def __init__(self, sect: Section, f1, f3):
            self.sect = sect
            self.f1 = f1
            self.f3 = f3

    def parabola_iter(function, _section: Section, f1, f3) -> ret:
        x1 = _section.begin
        len = (_section.end + _section.begin) / 2
        x2 = len
        x3 = _section.end
        if x1 + x3 == 0:
            x2 += (x3 - x2) / 2
        f2 = function(x2)

        u = x2 - ((x2 - x1) * (x2 - x1) * (f2 - f3) - (x2 - x3) * (x2 - x3) * (f3 - f1)) / 2 / \
            ((x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1))
        fu = function(u)

        if f2 > fu:
            if x2 < u:
                return ret(Section(x2, x3), f2, f3)
            else:
                return ret(Section(x1, x2), f1, f2)
        else:
            if x2 < u:
                return ret(Section(x1, u), f1, fu)
            else:
                return ret(Section(u, x3), fu, f3)

    iter_c = 0
    func_c = 2
    f1 = function(section.begin)
    f3 = function(section.end)
    while section.length() > accuracy:
        r = parabola_iter(function, section, f1, f3)
        section = r.sect
        f1 = r.f1
        f3 = r.f3
        iter_c += 1
        func_c += 1
    return Result(iter_c, func_c, section, "parabola")


def brent(function, section: Section, accuracy) -> Result:
    a = section.begin
    c = section.end
    x = a + G_CONST * (c - a)
    w = x
    v = x
    f_x = function(x)
    f_w = f_x
    f_v = f_x
    d = c - a
    e = d
    iter_c = 1
    func_c = 1
    while True:
        g = e
        e = d
        u_check = False
        if (x != w) and (x != v) and (w != v):
            u_check = True
            lst = [(x, f_x), (w, f_w), (v, f_v)]
            lst.sort(key=lambda z: z[0])
            x1 = lst[0][0]
            x2 = lst[1][0]
            x3 = lst[2][0]
            f1 = lst[0][1]
            f2 = lst[1][1]
            f3 = lst[2][1]
            u_mb = x2 - ((x2 - x1) * (x2 - x1) * (f2 - f3) - (x2 - x3) * (x2 - x3) * (f3 - f1)) / 2 / (
                        (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1))
        if u_check and Section(a + accuracy, c - accuracy).include(u_mb) and (abs(u_mb - x) < g / 2):
            u = u_mb                            # парабола
            d = abs(u - x)
        else:
            if x < (c - a) / 2:                 # золотое сечение х,с
                u = x + G_CONST * (c - x)
                d = c - x
            else:                               # золотое сечение а,х
                u = x - G_CONST * (c - x)
                d = x - a
        if abs(u - x) < accuracy:
            u = x + copysign(accuracy, u - x)

        f_u = function(u)
        func_c += 1
        if f_u <= f_x:
            if u >= x:
                a = x
            else:
                c = x
            v = w
            w = x
            x = u
            f_v = f_w
            f_w = f_x
            f_x = f_u
        else:
            if u >= x:
                c = u
            else:
                a = u
            if f_u <= f_w or w == x:
                v = w
                w = u
                f_v = f_w
                f_w = f_u
            elif f_u <= f_v or v == x or v == w:
                v = u
                f_v = f_u
        iter_c += 1
        if max(x - a - accuracy * 0.01, c - x - accuracy * 0.01) <= accuracy:
            return Result(iter_c, func_c, Section(x - accuracy/2, x + accuracy/2), "brent")


def fibonacci(function, section: Section, accuracy) -> Result:
    def F(n):
        z = 1 / sqrt(5) * (((1 + sqrt(5)) / 2) ** n - ((1 - sqrt(5)) / 2) ** n)
        return round(z)

    c_for_n = section.length() / accuracy
    n = 1
    while c_for_n >= F(n):
        n += 1

    a = section.begin
    b = section.end
    x1 = a + (F(n - 2) / F(n)) * (b - a) - ((-1) ** n) / F(n) * accuracy
    x2 = a + (F(n - 1) / F(n)) * (b - a) - ((-1) ** n) / F(n) * accuracy
    f1 = function(x1)
    f2 = function(x2)
    func_c = 2
    iter_c = 0
    j = 1
    res = 0
    while j <= n - 1:
        if (n-j+1) % 2 == 0:
            sign = 1
        else:
            sign = -1
        if f1 <= f2:
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + (F(n-j-1)/F(n-j+1))*(b-a) - sign * (accuracy / F(n-j+1))
            f1 = function(x1)
            res = x2
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + (F(n-j)/F(n-j+1))*(b-a) - sign * (accuracy / F(n-j+1))
            f2 = function(x2)
            res = x1
        func_c += 1
        iter_c += 1
        j += 1
    return Result(iter_c, func_c, Section(res-accuracy/2, res+accuracy/2), "fibonacci")

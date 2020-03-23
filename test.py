def t0(x):
    return 1

def t1(x):
    return x

def t2(x):
    # t2 = 2 * x * t1 - t0
    # t2 = 2 * x * x - 1
    return 2*x*x - 1

def t3(x):
    # t3 = 2 * x * t2 - t1
    # t3 = 2 * x * (2 * x * x - 1) - x
    # t3 = 4 * x**3 - 2x -x
    # t3 = 4 * x**3 - 3x
    return 4 * x ** 3 - 3 * x

def t4(x):
    return 8 * x ** 4 - 8 * x**2 + 1

def tn(x, n, num, denom):
    l = []
    for i in range(0, n + 1):
        if i == 0:
            l.append(1)
        elif i == 1:
            l.append(x)
        else:
            l.append(2 * x * l[-1] - l[-2])
    return (l[-1] / denom) * num


# tn(cos(x), n) = cos(n * x)

# cos(x) = a0 * t0(x) + a1 * t1(x) + a2 * t2(x) + a3 * t3(x) ...
# cos(x) = a0 * 1 + a1 * x + a2 * (2*x*x - 1) + a3 * (4 * x * x * x - 3 * x)
#
# Solve for a0, a1, a2, and a3 given known values of x:
# Notes for wolfram alpha:
def make_into_string(answers):
    print "linel:4096;"
    print "solve(["
    for x in answers:
        x2 = x[1] * x[1]
        x3 = x[1] * x[1] * x[1]
        x4 = x2 * x2
        x5 = x[1] * x4
        x6 = x3 * x3
        x7 = x4 * x3
        x8 = x4 * x4
        x10 = x8 * x2
        x12 = x6 * x6
        x14 = x7 * x7
        x16 = x8 * x8

        b1 = 1<<62
        b2 = b1*b1
        b3 = b1 * b2
        b4 = b2 * b2
        b5 = b3 * b2
        b6 = b3 * b3
        b7 = b3 * b4
        b8 = b4 * b4
        b10 = b8 * b2
        b12 = b6 * b6
        b14 = b7 * b7
        b16 = b8 * b8

        print "%s = a" % (x[0]),
        #print " + b * %i" % (x[1]),
        print " + c * (2 * ({0}/{1}) - 1)".format(x2, b2),
        #print " + d * (4 * {0} - 3 * {1})".format(x3, x[1]),
        print " + e * (8 * ({0}/{2}) - 8 * ({1}/{3}) + 1)".format(x4, x2, b4, b2),
        #print " + f * (16 * {0} - 20 * {1} + 5 * {2})".format(x5, x3, x[1]),
        print " + g * (32 * ({0}/{3}) - 48 * ({1}/{4}) + 18 * ({2}/{5}) - 1)".format(x6, x4, x2, b6, b4, b2),
        #print " + h * (64 * {0} - 112 * {1} + 56 * {2} - 7 * {3})".format(x7, x5, x3, x[1]),
        print " + j * (128 * ({0}/{4}) - 256 * ({1}/{5}) + 160 * ({2}/{6}) - 32 * ({3}/{7}) + 1)".format(x8, x6, x4, x2, b8, b6, b4, b2),
        print " + b * (512 * ({0}/{5}) - 1280 * ({1}/{6}) + 1120 * ({2}/{7}) - 400 * ({3}/{8}) + 50 * ({4}/{9}) - 1)".format(x10, x8, x6, x4, x2, b10, b8, b6, b4, b2),
        print " + d * (2048 * ({0}/{6}) - 6144 * ({1}/{7}) + 6192 * ({2}/{8}) - 3584 * ({3}/{9}) + 840 * ({4}/{10}) - 72 * ({5}/{11}) + 1)".format(x12, x10, x8, x6, x4, x2, b12, b10, b8, b6, b4, b2),
        print " + f * (8192 * ({0}/{7}) - 28672 * ({1}/{8}) + 39424 * ({2}/{9}) - 26880 * ({3}/{10}) + 9408 * ({4}/{11}) - 1568 * ({5}/{12}) + 98 * ({6}/{13}) - 1)".format(x14, x12, x10, x8, x6, x4, x2, b14, b12, b10, b8, b6, b4, b2),
        print " + h * (32768 * ({0}/{8}) - 131072 * ({1}/{9}) + 212992 * ({2}/{10}) - 180224 * ({3}/{11}) + 84480 * ({4}/{12}) - 21054 * ({5}/{13}) + 2688 * ({6}/{14}) - 128 * ({7}/{15}) + 1)".format(
                x16, x14, x12, x10, x8, x6, x4, x2, b16, b14, b12, b10, b8, b6, b4, b2)
        if x != answers[-1]:
            print ","
        else:
            print ""
    print "], [a, b, c, d, e, f, g, h, j])$"
    for i in range(len(answers)):
        print "string(%%th(%i)[1][%i]);" % (i + 1, i + 1)

pi = 1<<63
make_into_string([
    (1.0, 0),
    ("sqrt(2 + sqrt(2 + sqrt(2))) / 2", pi / 16), # sqrt(2 + sqrt(2 + sqrt(2))) / 2.0
    ("sqrt(2 + sqrt(2)) / 2", 2 * pi / 16), # sqrt(2 + sqrt(2)) / 2.0
    ("sqrt(2 + sqrt(2 - sqrt(2))) / 2", 3 * pi / 16), # sqrt(2 + sqrt(2 - sqrt(2))) / 2.0
    ("sqrt(2) / 2", 4 * pi / 16), # sqrt(2) / 2
    ("sqrt(2 - sqrt(2 - sqrt(2))) / 2", 5 * pi / 16), # sqrt(2 - sqrt(2 - sqrt(2))) / 2
    ("sqrt(2 - sqrt(2)) / 2", 6 * pi / 16), # sqrt(2 - sqrt(2)) / 2
    ("sqrt(2 - sqrt(2 + sqrt(2))) / 2", 7 * pi / 16), # sqrt(2 - sqrt(2 + sqrt(2))) / 2
    ("0", pi / 2),
])

def cos(x):
    coeffs = [0.7651967488801196,
           -1.672470895113307e-5,
           -0.2297973984541904,
           -1.325934365379417e-5,
           0.004966809308769983,
           -7.298967585728185e-6,
           -3.904335691434468e-5]

    ts = [1.0, float(x)]
    s = 0.0
    for i in range(len(coeffs)):
        if i == 0:
            s = ts[i] * coeffs[i]
        elif i == 1:
            s += ts[i] * coeffs[i]
        else:
            ts.append(2 * x * ts[-1] - ts[-2])
            s += ts[i] * coeffs[i]
    return s


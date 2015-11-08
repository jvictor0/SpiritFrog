

def c_add (x,y):
    return (x[0] + y[0], x[1] + y[1])

def c_sub (x,y):
    return (x[0] - y[0], x[1] - y[1])

def c_mult (x, y):
    a,b = x
    c,d = y
    return (a*c - b*d, b*c + a*d)

def c_abs2(x):
    return x[0] * x[0] + x[1] * x[1]

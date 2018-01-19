def f(x,num):
    return x**2 - num

def df(x):
    return 2*x

def sqrt(f, df, num):
    a = 1000
    for i in range(1000):
        a = a - f(a,num)/df(a)

    print(a)

while True:
    num = float(input("What number do you want to take the sqrt of \n: "))
    sqrt(f,df,num)

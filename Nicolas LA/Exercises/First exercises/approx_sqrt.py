from time import perf_counter
import matplotlib.pyplot as plt
def f(x,num):
    return x**2 - num

def df(x):
    return 2*x

def sqrt_newton(f, df, num, initial_value):
    root = initial_value
    for i in range(1000):
        root = root - f(root,num)/df(root)
    return root

def sqrt_search(num,error):
    root = num/2
    limit_upper = num
    limit_lower = 0
    while abs(root**2 - num ) > error:
        if root**2 > num:
            limit_upper = root
            root = (limit_upper-limit_lower)/2 + limit_lower
        if root**2 < num:
            limit_lower = root
            root = (limit_upper-limit_lower)/2 + limit_lower
    return(root)



initial_values = [1+1j, (-5) , -1-1j, 2]


for num in initial_values:
    print("Initial value", num,  )
    print("Result ",sqrt_newton(f,df ,-1 ,num),"\n")

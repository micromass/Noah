
import matplotlib.pyplot as plt
from numpy import linspace, zeros, loadtxt
def f(x,num):
    return x**2 - num

def df(x):
    return 2*x

def sqrt_newton(f, df, num, initial_value):
    root = initial_value
    for i in range(1000):
        root = root - f(root,num)/df(root)
    return root




 # Initial set of variables
cells = 20
rows = cells * 2 + 1
columns = cells * 2 + 1
num = -1

#Creates table with complex numbers and then evaluates the root of such table

initial_values = zeros((rows , columns), complex)
for x  in range(-cells,cells + 1):
    for y in range(-cells, cells + 1):
        initial_values[x + cells][y + cells] = complex(x,y)
initial_values = sqrt_newton(f, df, num, initial_values)


#Creates output talbe
plot_values = zeros((rows,columns))
for x in range(-cells,cells + 1):
    for y in range(-cells, cells +1):
        if initial_values[x + cells][y+ cells] ==1j:
            plot_values[x + cells][y + cells] = 2
        elif initial_values[x + cells][y + cells] ==-1j:
            plot_values[x + cells][y + cells] = 1
        else:
            plot_values[x + cells][y + cells] = 0

text_file = open("complexArray.txt","w")
for x in range(-cells,cells+1):
    for y in range(-cells,cells+1):
        text_file.write(f"{initial_values[x][y]}  ")
    text_file.write("\n")
text_file.close()



# x equals the real value y equals the imaginary value
plt.imshow(plot_values, extent =[-cells, cells, cells, -cells])
plt.gray()
plt.xlabel("Imaginary Number")
plt.ylabel("Real Number")
plt.show()

"""
for num in initial_values:
    print("Initial value", num,  )
    print("Result ",sqrt_newton(f,df ,-1 ,num),"\n")
"""

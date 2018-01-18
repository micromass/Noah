from operations import calcbinomial_form, calcbinomial_tail, calcbinomial_multiply, calcbinomial_reduce
from pylab import plot, show
import time


x = range(1,15)
y_form = []
y_tail = []
y_multiply = []
y_reduce =[]

for n in x:
    a_1 = time.perf_counter()
    calcbinomial_form(2*n,n)
    a_2 = time.perf_counter()

    y_form.append(a_2 - a_1)


for n in x:
    a_1 = time.perf_counter()
    calcbinomial_tail(2*n,n)
    a_2 = time.perf_counter()

    y_tail.append(a_2 - a_1)


for n in x:
    a_1 = time.perf_counter()
    calcbinomial_multiply(2*n,n)
    a_2 = time.perf_counter()

    y_multiply.append(a_2 - a_1)

for n in x:
    a_1 = time.perf_counter()
    calcbinomial_reduce(2*n,n)
    a_2 = time.perf_counter()

    y_reduce.append(a_2 - a_1)


plot(x, y_form,"r")
plot(x, y_multiply,"g")
plot(x, y_tail,"b")
plot(x, y_reduce,"y")

show()

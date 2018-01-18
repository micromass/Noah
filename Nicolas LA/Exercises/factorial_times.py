from operations import calcfactorial
from pylab import plot, show
import time

y =[]
x = range(1,500)

for i in x:
    a_1 = time.perf_counter()
    calcfactorial(i)
    a_2 = time.perf_counter()
    y.append(a_2 - a_1)

print(y)
plot(x,y)
show()

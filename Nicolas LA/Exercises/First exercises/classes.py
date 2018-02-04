
def f(x,num):
    return x**2 - num

def df(x):
    return 2*x

def sqrt_newton(f, df, num, error):
    root = num/2
    while abs(root** - num ) > error:
        root = root - f(a,num)/df(a)

class Complex():
    def __init__(self, real, imaginary):
        if not isinstance(real,float):
            raise TypeError('Input should be real numbers')
        if not isinstance(imaginary, float):
            raise TypeError("Input should be real numbers")
        self.real = real
        self.imaginary = imaginary

    def __repr__(self):
        if self.imaginary < 0:
            return f"{self.real} - {abs(self.imaginary)}i"
        else:
            return f"{self.real} + {self.imaginary}i"

    def __add__(self,other):
        if isinstance( other, Complex):
            return Complex(self.real + other.real, self.imaginary + other.imaginary)
        elif isinstance(other, float):
            return Complex(self.real + other, self.imaginary)
        else:
            raise TypeError("You can only add two complex numbers")
    def __radd__(self,other):
        other + self

    def sqrt_complex(self):

c = Complex(3,4)
a = Complex(3,5)

print(c + a)

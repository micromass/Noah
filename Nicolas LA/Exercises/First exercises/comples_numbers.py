"""This program computes complex numbers without using OOP"""
from numpy import empty


def array_to_complex(c_r):
    return f"{c_r[0]} + {c_r[1]}i"

def complex_sum(c_1,c_2):
    """Computes the sum of two complex numbers"""
    return c_1 + c_2

def complex_difference(c_1,c_2):
    """Computes the difference of
     two complex numbers"""
    return c_1 - c_2

def complex_multiplication(c1,c2,cr):
    """Computes the multiplication of two complex numbers"""
    cr[0] = c1[0]*c2[0] - c1[1]*c2[1]
    cr[1] = c1[0]*c2[1] + c1[1]*c2[0]
    return cr

def complex_reader(operation,input):
    """Reads a string and extracts complex numbers. Return, values of list each list per index and operation """
    input = input[1:len(input)-1]

    #extracts the first complex number and arranges the real part and imaginary part in list
    complex_1 = input[: input.find(')')].split('+') #Creates a list with the imaginary and real part of the first complex number
    complex_1[1] = complex_1[1][:complex_1[1].find('i')] # removes the i from the second part


    #extracts the second complex number and arranges the real part and imaginary part  in list
    complex_2 = input[input.find('(') + 1 : ].split('+')# Creates a list with the imaginary and real part of the complex numbers
    complex_2[1] = complex_2[1][:complex_2[1].find('i')] # removes the i from the second part



    #extracts operation
    operation  = input[ input.find(')') + 1 : input.find('(') ] # finds the operation and assigns it to the variable

    return complex_1[0], complex_1[1] , complex_2[0], complex_2[1], operation

def complex_inverse(c1,cr):
    """Gives the inverse of a complex number"""
    

#Defines variables to work with
operation = None
c1 = empty(2,float)
c2 = empty(2,float)
cr = empty(2,float)
c1[0], c1[1], c2[0], c2[1], operation = complex_reader(operation,"(3+4i)*(2+5i)")

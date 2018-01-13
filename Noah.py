import math

def is_prime(n):
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d = d+1
    return n > 1

def Modulo(n):
    class aux(object):
        def __init__(self,input):
            if isinstance(input,str):
                self.value = int(input) % n
            elif isinstance(input,int):
                self.value = input % n
            else:
                raise TypeError('Type not recognized')
        def __add__(self,other):
            if isinstance(other,aux):
                return aux(self.value + other.value)
            else:
                raise TypeError('Arguments are of wrong types')
        def __radd__(self,other):
            if isinstance(other,aux):
                return aux(self.value + other.value)
            else:
                raise TypeError('Arguments are of wrong types')
        def __sub__(self,other):
            if isinstance(other,aux):
                return aux(self.value - other.value)
            else:
                raise TypeError('Arguments are of wrong types')
        def __mul__(self,other):
            if isinstance(other,aux):
                return aux(self.value * other.value)
            elif isinstance(other,int):
                return aux(self.value * other)
            else:
                raise TypeError('Arguments are of wrong types')
        def __rmul__(self,other):
            if isinstance(other,aux):
                return aux(self.value * other.value)
            elif isinstance(other,int):
                return aux(self.value * other)
            else:
                raise TypeError('Arguments are of wrong types')
        def __pow__(self,other):
            if isinstance(other, int):
                if other < 0:
                    return NotImplemented
                else:
                    a,b,A = self, aux(1), other
                    while A>0:
                        if A % 2 == 1:
                            b = b*a
                        a = a*a
                        A = math.floor(A/2)
                    return b                
            else:
                raise TypeError('Arguments are of wrong types')
        def __floordiv__(self,other):
            if isinstance(other, aux):
                raise TypeError('Arguments are of wrong types')
            elif is_prime(n):
                return self / other
            else:
                raise TypeError('Ring is not a Euclidean Domain')
                
        def __mod__(self,other):
            if isinstance(other, aux):
                raise TypeError('Arguments are of wrong types')
            elif is_prime(n):
                return self.ReturnZero()
            else:
                raise TypeError('Ring is not a Euclidean Domain')
        def __neg__(self):
            return aux(-self.value)
        def __eq__(self,other):
            if isinstance(other,aux):
                return self.value == other.value
            else: 
                return False
        def divides(self,b):
            if not isinstance(b,aux):
                raise TypeError('Arguments are of wrong types')
            return Divides(ComputeGCD(self.value, n), b.value)
        def __truediv__(self,a):
            if not isinstance(a,aux):
                raise TypeError('Arguments are of wrong types')
            else:
                bez = Bezout(a.value, n)
                if not Divides(bez['gcd'], n):
                    raise ValueError('Arguments do not divide')
                else:
                    sol = int(bez['coefficientfirst']*self.value/bez['gcd']) % int(n/bez['gcd'])
                    if bez['gcd'] == 1:
                        return aux(sol)
                    else:
                        L = list(range(0, bez['gcd']-1))
                        return [aux(sol+int(n/bez['gcd'])*x) for x in L]
        def Invert(self):
            return self.ReturnIdentity() / self
        def __repr__(self):
            return str(self.value)
        def ReturnZero(self):
            return aux(0)
        def ReturnIdentity(self):
            return aux(1)
        def Commutative(self):
            return True
        def Identity(self):
            return True
        def IsRng(self):
            return True
        def HasZeroDivisors(self):
            return is_prime(n)
        def HasInverses(self):
            return is_prime(n)
        def HasDivAlgo(self):
            return is_prime(n)
    return aux

def Divides(a,b):
    if not type(a)==type(b):
        raise TypeError('Arguments dont have the same type')
    if IsEuclideanDomain(a):
       return IsZero(b % a)
    else:
        return a.divides(b)
    
def IsInvertible(a):
    if not IsRing(a):
        raise TypeError('Argument are not in ring')
    else:
        return Divides(a, Identity(a))
    
def Invert(a):
    if not IsRing:
        raise TypeError('Argument are not in ring')
    else:
        return a.Invert()
        
    

def Zero(typ):
    if isinstance(typ,int):
        return 0
    elif isinstance(typ,float):
        return 0
    else:
        return typ.ReturnZero()

def Identity(typ):
    if isinstance(typ,int):
        return 1
    elif isinstance(typ,float):
        return 1
    else:
        return typ.ReturnIdentity()

def IsRng(typ):
    if isinstance(typ,int):
        return True
    elif isinstance(typ,float):
        return True
    else: 
        return typ.IsRng()

def IsCommutativeRng(typ):
    if isinstance(typ,int):
        return True
    elif isinstance(typ,float):
        return True
    elif IsRng(typ):
        return typ.Commutative()
    else:
        return False


def IsRing(typ):
    if isinstance(typ,int):
        return True
    elif isinstance(typ,float):
        return True
    elif IsRng(typ):
        return typ.Identity()
    else:
        return False

def IsCommutativeRing(typ):
    return IsRing(typ) and IsCommutativeRng(typ)

def IsIntegralDomain(typ):
    if isinstance(typ,int):
        return True
    elif isinstance(typ,float):
        return True
    elif IsCommutativeRing(typ):
        return typ.HasZeroDivisors()
    else:
        return False
    
def IsDivisionRing(typ):
    if isinstance(typ,int):
        return False
    elif isinstance(typ,float):
        return True
    elif IsRing(typ):
        return typ.HasInverses()
    else:
        return False

def IsField(typ):
    return IsDivisionRing(typ) and IsCommutativeRng(typ)

def IsEuclideanDomain(typ):
    if isinstance(typ,int):
        return True
    elif isinstance(typ,float):
        return False
    elif IsField(typ):
        return True
    elif IsIntegralDomain(typ):
        return typ.HasDivAlgo()
    else:
        return False
    
def IsZero(n):
    if not IsRng(n):
        raise ValueError('Inputs do not belong to a Rng')
    else:
        return n == Zero(n)
    
def IsIdentity(n):
    if not IsRing(n):
        raise ValueError('Inputs do not belong to a Ring')
    else:
        return n == Identity(n)
    
def ComputeGCD(n,m):
    if not type(n) == type(m):
        raise ValueError('Inputs are not of the same type')
    elif not IsEuclideanDomain(n):
        raise ValueError('Inputs are not elements of a Euclidean Domain')
    elif IsZero(m):
        return n
    else:
        return ComputeGCD(m, n % m)

def Bezout(a,b):
    if not type(a) == type(b):
        raise ValueError('Inputs are not of the same type')
    elif not IsEuclideanDomain(a):
        raise ValueError('Inputs are not elements of a Euclidean Domain')
    else:
        if IsZero(b):
            return {'gcd' : a, 'coefficientfirst' : Identity(a), 'coefficientsecond' : Zero(a)}
        else:
            u,g,x,y = Identity(a),a,Zero(a),b
            while not IsZero(y):
                u,g,x,y = x,y,u-x*(g // y), g % y
            v = (g-a*u)/b
            if isinstance(a,int):
                aux = math.ceil(-u*g/b)
                u = int(u + aux*b/g)
                v = int(v - aux*a/g)
            return {'gcd' : g, 'coefficientfirst' : u, 'coefficientsecond' : v}
            

    
    


    
        
    
        

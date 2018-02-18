import math
import re
import copy
import time

def is_prime(n):
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d = d+1
    return n > 1

def Modulo():
    AlreadyDone = []
    def CreateObject(n):
        class aux(object):
            def __init__(self,input):
                if input == "default":
                    self.value = 0
                elif isinstance(input,int):
                    self.value = input % n
                elif isinstance(input,str):
                    input = re.sub(' ', '',input)
                    input = StringReplace(input)
                    k = FindMainField(input)
                    L = len(input)
                    if k == "Primitive":
                        self.value = int(input) % n
                    elif input[k] == "(":
                        hulp = aux(input[1:L-1])
                        self.value = hulp.value
                    elif input[k] == "-":
                        if k == 0:
                            hulp = - aux(input[1:])
                            self.value = hulp.value
                        else:
                            hulp = aux(input[:k]) - aux(input[(k+1):])
                            self.value = hulp.value
                    elif input[k] == "+":
                        hulp = aux(input[:k]) + aux(input[(k+1):])
                        self.value = hulp.value
                    elif input[k] == "*":
                        hulp = aux(input[:k]) * aux(input[(k+1):])
                        self.value = hulp.value
                    elif input[k] == "/":
                        hulp = aux(input[:k]) + aux(input[(k+1):])
                        self.value = hulp.value
                    elif input[k] == "^":
                        hulp = aux(input[:k]) ** Z(input[(k+1):])
                        self.value = hulp.value
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
                        return Inverse(self ** (-other))
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
                    if IsZero(other):
                        return self
                    else: 
                        return aux(0)
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
                if self.divides(aux(1)):
                    return aux(1) / self
                else:
                    return("Element does not have an inverse")
            def __repr__(self):
                return str(self.value)
            def ElementOf(self):
                return internal
        def internal(input):
            if input == "name":
                return "Z_n"
            elif input == "Zero":
                return aux(0)
            elif input == "Identity":
                return aux(1)
            elif input == "HasIdentity":
                return True
            elif input == "IsRng":
                return True
            elif input == "HasZeroDivisors":
                return not is_prime(n)
            elif input == "HasInverses":
                return not is_prime(n)
            elif input == "HasDivAlgo":
                return not is_prime(n)
            elif input == "IsCommutative":
                return True
            elif input == "Ordered":
                return False
            elif input == "GroundRing":
                return Z
            elif input == "Structure":
                return "Quotient"
            else:
                return aux(input)
        if not isinstance(n,int):
            raise TypeError("Input not an integer")
        if n<1:
            raise TypeError("Input not a positive integer")
        if AlreadyDone == []:
            AlreadyDone.append([n,internal])
            return internal
        else:
            for iter in range(0,len(AlreadyDone)):
                if AlreadyDone[iter][0] == n:
                    return AlreadyDone[iter][1]
            AlreadyDone.append([n,internal])
            return internal
    return CreateObject

Mod = Modulo()
            
    
def Z(input):
        if input == "name":
            return "Z"
        elif input == "Zero":
            return int(0)
        elif input == "Identity":
            return int(1)
        elif input == "IsRng":
            return True
        elif input == "IsCommutative":
            return True
        elif input == "HasZeroDivisors":
            return False
        elif input == "HasInverses":
            return False
        elif input == "HasDivAlgo":
            return True
        elif input == "HasIdentity":
            return True
        elif input == "Ordered":
            return True
        elif input == "Structure":
            return "Basic"
        elif isinstance(input,str):
            input = re.sub(' ', '',input)
            input = StringReplace(input)
            k = FindMainField(input)
            L = len(input)
            if k == "Primitive":
                return int(input)
            elif input[k] == "(":
                return Z(input[1:L-1])
            elif input[k] == "-":
                if k == 0:
                    return - Z(input[1:])
                else:
                    return Z(input[:k]) - Z(input[(k+1):])
            elif input[k] == "+":
                return Z(input[:k]) + Z(input[(k+1):])
            elif input[k] == "*":
                hulp = Z(input[:k]) * Z(input[(k+1):])
            elif input[k] == "/":
                hulp1 = Z(input[:k]) 
                hulp2 = Z(input[(k+1):])
                if Divides(hulp2, hulp1):
                    return hulp1 // hulp2
                else:
                    raise TypeError("Division not an integer")
            elif input[k] == "^":
                hulp1 = Z(input[:k])
                hulp2 = Z(input[(k+1):])
                if hulp2 >= 0:
                    return hulp1 ** hulp2
                else:
                    if hulp1 == 1 or hulp1 == -1:
                        return hulp1 ** (- hulp2)
                    else:
                        raise ValueError("Exponentiation not possible")
        else:
            return int(input)
    
def R(input):
        if input == "name":
            return "R"
        elif input == "Zero":
            return int(0)
        elif input == "Identity":
            return int(1)
        elif input == "IsRng":
            return True
        elif input == "IsCommutative":
            return True
        elif input == "HasZeroDivisors":
            return False
        elif input == "HasInverses":
            return True
        elif input == "HasDivAlgo":
            return True
        elif input == "Ordered":
            return True
        elif input == "Structure":
            return "Basic"
        elif isinstance(input,str):
            input = re.sub(' ', '',input)
            input = StringReplace(input)
            k = FindMainField(input)
            L = len(input)
            if k == "Primitive":
                return float(input)
            elif input[k] == "(":
                return R(input[1:L-1])
            elif input[k] == "-":
                if k == 0:
                    return - R(input[1:])
                else:
                    return R(input[:k]) - R(input[(k+1):])
            elif input[k] == "+":
                return R(input[:k]) + R(input[(k+1):])
            elif input[k] == "*":
                hulp = R(input[:k]) * R(input[(k+1):])
            elif input[k] == "/":
                hulp1 = R(input[:k]) 
                hulp2 = R(input[(k+1):])
                if Divides(hulp2, hulp1):
                    return hulp1 // hulp2
                else:
                    raise ValueError("Division not possible")
            elif input[k] == "^":
                return R(input[:k]) ** R(input[(k+1):])
        else:
            return float(input)

def ElementOf(a):
    if isinstance(a,int):
        return Z
    elif isinstance(a,float):
        return R
    else: 
        return a.ElementOf()
    
def IsElementOf(a,Rng):
        return ElementOf(a) == Rng

def Divides(a,b):
    if not type(a)==type(b):
        raise TypeError('Arguments dont have the same type')
    if IsEuclideanDomain(ElementOf(a)):
       return IsZero(b % a)
    else:
        return a.divides(b)
    
def IsInvertible(a):
    if not IsRing(ElementOf(a)):
        raise TypeError('Argument are not in ring')
    else:
        return Divides(a, Identity(ElementOf(a)))
    
def Invert(a):
    if not IsRing(ElementOf(a)):
        raise TypeError('Argument are not in ring')
    else:
        return a.Invert()
        
    

def Zero(typ):
    return typ("Zero")

def Identity(typ):
    return typ("Identity")

def IsRng(typ):
    return typ("IsRng")

def name(typ):
    return typ("name")

def IsCommutativeRng(typ):
    return IsRng(typ) and typ("IsCommutative")


def IsRing(typ):
    return IsRng(typ) and typ("HasIdentity")
    
def IsOrderedRing(typ):
    return IsRing(typ) and typ("Ordered")

def IsCommutativeRing(typ):
    return IsRing(typ) and IsCommutativeRng(typ)

def IsIntegralDomain(typ):
    return IsRing(typ) and not typ("HasZeroDivisors")
    
def IsDivisionRing(typ):
    return IsRing(typ) and typ("HasInverses")

def IsField(typ):
    return IsDivisionRing(typ) and IsCommutativeRng(typ)

def IsEuclideanDomain(typ):
   return IsRing(typ) and typ("HasDivAlgo")
    
def IsZero(n):
    if not IsRng(ElementOf(n)):
        raise ValueError('Inputs do not belong to a Rng')
    else:
        return n == Zero(ElementOf(n))
    
def IsIdentity(n):
    if not IsRing(ElementOf(n)):
        raise ValueError('Inputs do not belong to a Ring')
    else:
        return n == Identity(ElementOf(n))
    
def HasSameCategory(n,m):
    return ElementOf(n) == ElementOf(m)
    
def ComputeGCD(n,m):
    if not HasSameCategory(n,m):
        raise ValueError('Inputs are not of the same type')
    elif not IsEuclideanDomain(ElementOf(n)):
        raise ValueError('Inputs are not elements of a Euclidean Domain')
    elif IsZero(m):
        return n
    else:
        return ComputeGCD(m, n % m)

def Bezout(a,b):
    if not HasSameCategory(a,b):
        raise ValueError('Inputs are not of the same type')
    elif not IsEuclideanDomain(ElementOf(a)):
        raise ValueError('Inputs are not elements of a Euclidean Domain')
    else:
        category = ElementOf(a)
        if IsZero(b):
            return {'gcd' : a, 'coefficientfirst' : Identity(category), 'coefficientsecond' : Zero(category)}
        else:
            u,g,x,y = Identity(category),a,Zero(category),b
            while not IsZero(y):
                u,g,x,y = x,y,u-x*(g // y), g % y
            v = (g-a*u)/b
            if category == Z:
                aux = math.ceil(-u*g // b)
                u = u + aux*b//g
                v = v - aux*a//g
            return {'gcd' : g, 'coefficientfirst' : u, 'coefficientsecond' : v}
        
def StringReplace(string):
    L = len(string)
    if L == 0:
        return string
    Res = string[0]
    prev = string[0]
    for iter in range(1,L):
        if string[iter] == '-':
            if prev in ["(", "+", "*","-","/","^"]:
                Res = Res + string[iter]
            else:
                Res = Res + "+" + string[iter]
        else:
            Res = Res + string[iter]
        prev = string[iter]
    return Res
            
    

def FindMainField(string):
    L = len(string)
    Depth = 0
    var = False
    bracketoccuring = False
    for iter in range(0,L-1):
        if iter>0 and Depth == 0 and string[iter] in ["+","*","-","/","^"]:
            var=True
        if string[iter] == "(":
            Depth = Depth + 1
            bracketoccuring = True
        elif string[iter] == ")":
            Depth = Depth - 1
        elif string[iter] == "+":
            if Depth == 0:
                return iter
        elif string[iter] == "-":
            if iter == 0:
                k = FindMainField(string[1:])
                if k == "Primitive":
                    return 0
                elif string[k+1] == "(":
                    return 0
            elif Depth == 0:
                return iter
    if not var and bracketoccuring:
        if string[L-1] == ")":
                    return 0
        else:
            raise ValueError("bracket mismatch")
    Depth = 0
    for iter in range(0,L-1):
        if string[iter] == "(":
            Depth = Depth + 1
        elif string[iter] == ")":
            Depth = Depth - 1
        elif string[iter] == "*":
            if Depth == 0:
                return iter
        elif string[iter] == "/":
            if Depth == 0:
                return iter
    for iter in range(0,L-1):
        if string[iter] == "(":
            Depth = Depth + 1
        elif string[iter] == ")":
            Depth = Depth - 1
        elif string[iter] == "^":
            if Depth == 0:
                return iter
    return("Primitive")
            

def RingOfFractionsCreator():
    AlreadyDone = []
    def RingOfFractions(IntDom):
        if not IsIntegralDomain(IntDom):
             raise TypeError('Argument must be an integral domain')
        else:
            A = IsEuclideanDomain(IntDom)
            class aux(object):
                def __init__(self,input):
                    if input == "default":
                        self.numerator = Zero(IntDom)
                        self.denominator = Identity(IntDom)
                    elif isinstance(input,int):
                        hulp = Identity(IntDom) * input
                        self.numerator = hulp.numerator
                        self.denominator = hulp.denominator
                    elif isinstance(input,str):
                        input = re.sub(' ', '',input)
                        input = StringReplace(input)
                        k = FindMainField(input)
                        L = len(input)
                        if k == "Primitive":
                            self.numerator = IntDom(input)
                            self.denominator = Identity(IntDom)
                        elif input[k] == "(":
                            hulp = aux(input[1:L-1])
                            self.numerator = hulp.numerator
                            self.denominator = hulp.denominator
                        elif input[k] == "-":
                            if k == 0:
                                hulp = - aux(input[1:])
                                self.numerator = hulp.numerator
                                self.denominator = hulp.denominator
                            else:
                                hulp = aux(input[:k]) - aux(input[(k+1):])
                                self.numerator = hulp.numerator
                                self.denominator = hulp.denominator
                        elif input[k] == "+":
                            hulp = aux(input[:k]) + aux(input[(k+1):])
                            self.numerator = hulp.numerator
                            self.denominator = hulp.denominator
                        elif input[k] == "*":
                            hulp = aux(input[:k]) * aux(input[(k+1):])
                            self.numerator = hulp.numerator
                            self.denominator = hulp.denominator
                        elif input[k] == "/":
                            hulp = aux(input[:k]) / aux(input[(k+1):])
                            self.numerator = hulp.numerator
                            self.denominator = hulp.denominator
                        elif input[k] == "^":
                            hulp = aux(input[:k]) ** int(input[(k+1):])
                            self.numerator = hulp.numerator
                            self.denominator = hulp.denominator
                    elif ElementOf(input) == IntDom:
                        self.numerator = input
                        self.denominator = Identity(IntDom)
                    else:
                        raise TypeError('Type not recognized')
                
                def __add__(self,other):
                    if isinstance(other,aux):
                        hulp = aux("default")
                        num = (self.numerator * other.denominator) + (self.denominator * other.numerator)
                        denom = self.denominator*other.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom
                        return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __sub__(self,other):
                    if isinstance(other,aux):
                        hulp = aux("default")
                        num = (self.numerator * other.denominator) - (self.denominator * other.numerator)
                        denom = self.numerator*other.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom
                        return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __mul__(self,other):
                    if isinstance(other,aux):
                        hulp = aux("default")
                        num = self.numerator * other.numerator
                        denom = self.denominator*other.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom
                        return(hulp)
                    elif isinstance(other,int):
                        hulp = aux("default")
                        num = self.numerator * other
                        denom = self.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom 
                        return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __rmul__(self,other):
                    if isinstance(other,aux):
                        hulp = aux("default")
                        num = self.numerator * other.numerator
                        denom = self.denominator*other.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom
                        return(hulp)
                    elif isinstance(other,int): 
                        hulp = aux("default")
                        num = self.numerator * other
                        denom = self.denominator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        hulp.numerator = num
                        hulp.denominator = denom 
                        return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __pow__(self,other):
                    if isinstance(other,int):
                        if(other<0):
                            return(Inverse(self ** (-other)))
                        else:
                            hulp = aux("default")
                            num = self.numerator ** other
                            denom = self.denominator ** other
                            if A:
                                G = ComputeGCD(num,denom)
                                num = num // G
                                denom = denom // G
                            hulp.numerator = num
                            hulp.denominator = denom 
                            return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __truediv__(self,other):
                    if isinstance(other,aux):
                        hulp = aux("default")
                        num = self.numerator * other.denominator
                        denom = self.denominator*other.numerator
                        if A:
                            G = ComputeGCD(num,denom)
                            num = num // G
                            denom = denom // G
                        if IsZero(denom):
                            raise ValueError('Division by Zero')
                        if IsOrderedRing(IntDom):
                            if denom < Zero(IntDom):
                                num = - num
                                denom = -denom
                        hulp.numerator = num
                        hulp.denominator = denom
                        return(hulp)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __mod__(self,other):
                    if isinstance(other,aux):
                        if IsZero(other):
                            return self
                        else:
                            return aux(Zero(IntDom))
                def __floordiv__(self,other):
                    return self / other
                def divides(self,b):
                    return IsZero(b)
                def Invert(self):
                    return self.ReturnIdentity() / self
                def __neg__(self):
                    hulp = aux("default")
                    hulp.numerator = -self.numerator
                    hulp.denominator = self.denominator
                    return hulp
                def __le__(self,other):
                    if not IsOrderedRing(IntDom):
                        return TypeError("Comparison not implemented")
                    if isinstance(other,aux):
                        return (self.numerator*otherdenominator) <= (self.denominator * other.numerator)
                    else:
                        return TypeError("Comparison not implemented")
                def __lt__(self,other):
                    return self <= other and not self == other
                def __gt__(self,other):
                    return other < self
                def __ge__(self,other):
                    return other <= self
                def __repr__(self):
                    if IsInvertible(self.denominator):
                        self.numerator = self.numerator // self.denominator
                        self.denominator = 1
                        return str(self.numerator)
                    else:
                        return "(" + str(self.numerator) + "/" + str(self.denominator) + ")"
                def ElementOf(self):
                    return internal
                def __eq__(self,other):
                    if isinstance(other,aux):
                        return (self.numerator * other.denominator) == (self.denominator * other.numerator)
                    elif isinstance(other,str):
                        return False
                    elif ElementOf(other) == IntDom:
                        return self.numerator == (self.denominator*other)
                    else:
                        return False
        def internal(input):
            if input == "name":
                if IntDom == Z:
                    return "Q"
                else:
                    return "Q(" + Rng("name") + ")"
            elif input == "Zero":
                return aux(Zero(IntDom))
            elif input == "Identity":
                return aux(Identity(IntDom))
            elif input == "HasIdentity":
                return True
            elif input == "IsRng":
                return True
            elif input == "HasZeroDivisors":
                return False
            elif input == "HasInverses":
                return True
            elif input == "HasDivAlgo":
                return True
            elif input == "IsCommutative":
                return True
            elif input == "GroundRing":
                return IntDom
            elif input == "Ordered":
                return IsOrderedRing(IntDom)
            elif input == "Structure":
                return "Localization"
            else:
                return aux(input)
        if AlreadyDone == []:
            AlreadyDone.append([IntDom,internal])
            return internal
        else:
            for iter in range(0,len(AlreadyDone)):
                if AlreadyDone[iter][0] == IntDom:
                    return AlreadyDone[iter][1]
            AlreadyDone.append([IntDom,internal])
            return internal
    return RingOfFractions

RingOfFractions = RingOfFractionsCreator()
Q = RingOfFractions(Z)


def MatrixRingCreator():
    AlreadyDone = []
    def CreateCat(nrow, ncol,Rng):
        if not IsRng(Rng):
            raise TypeError('Must be presented with rng')
        if not isinstance(nrow, int):
            raise TypeError('Number of rows must be an integer')
        if not isinstance(ncol, int):
            raise TypeError('Number of columns must be an integer')
        if nrow < 1:
            raise TypeError('Numbers of rows must be positive')
        if ncol <1:
            raise TypeError('Number of columns must be positive')
        class aux(object):
            def __init__(self,input):
                if input == "default":
                    hulp = [None for i in range(0,nrow)]
                    for rw in range(0, nrow):
                        hulp[rw] = [None for i in range(0,ncol)]
                        for cl in range(0, ncol):
                            hulp[rw][cl] = Zero(Rng)
                    self.value = hulp
                else:
                    if not len(input) == nrow:
                        raise TypeError('Insufficient number of rows')
                    else:
                        hulp = [None for i in range(0,nrow)]
                        for rw in range(0, nrow):
                            if not len(input[rw]) == ncol:
                                raise TypeError("Insufficient number of columns in row" + str(rw+1))
                            hulp[rw] = [None for i in range(0,ncol)]
                            for cl in range(0, ncol):
                                entry = input[rw][cl]
                                if isinstance(entry, T):
                                    hulp[rw][cl] = entry
                                elif isinstance(entry,str) or isinstance(entry, int):
                                    hulp[rw][cl] = Rng(entry)
                                else:
                                    raise TypeError("Unknown type in entry" + str(rw) + str(cl))
                        self.value = hulp
            def getitem(self,rw,cl):
                if not isinstance(rw, int):
                     raise TypeError('Row number must be an integer')
                if not isinstance(cl, int):
                     raise TypeError('Column number must be an integer')
                if rw < 1 or rw > nrow:
                     raise TypeError('Row number not accessible')
                if cl< 1 or cl > ncol:
                     raise TypeError('Column number not accessible')
                return self.value[rw-1][cl-1]
            def __repr__(self):
                hulp = ""
                for rw in range(0,nrow):
                    for cl in range(0,ncol):
                        if cl == 0:
                            hulp = hulp + str(self.value[rw][cl])
                        else:
                            hulp = hulp + "\t & \t" + str(self.value[rw][cl])
                    hulp = hulp + "\n"
                return(hulp)
            def ElementOf(self):
                    return internal
        callnumber = 0
        def internal(input):
            if input == "callnumber":
                callnumber = callnumber + 1
                return callnumber
            if input == "name":
                return "M_{" + str(nrow) + "," + str(ncol) + "}(" + Rng("name") + ")"
            elif input == "GroundRing":
                return Rng
            elif input == "Dimensions":
                return [nrow,ncol]
            elif input == "Structure":
                return "MatrixRing"
            else:
                return aux(input)
        if AlreadyDone == []:
            AlreadyDone.append([nrow,ncol,Rng,internal])
            return internal
        else:
            for iter in range(0,len(AlreadyDone)):
                if AlreadyDone[iter][0:3] == [nrow, ncol,Rng]:
                    return AlreadyDone[iter][2]
            AlreadyDone.append([nrow,ncol,Rng,internal])
            return internal
    return CreateCat

def GetEntry(mat, rw, cl):
    return mat.getitem(rw,cl)
    
MatrixRing = MatrixRingCreator()
 



    
def LEX(alpha, beta):
    if not len(alpha) == len(beta):
        raise ValueError("Arguments not of same length")
    for iter in range(0,len(alpha)):
        if alpha[iter] < beta[iter]:
            return True
        if alpha[iter] > beta[iter]:
            return False
    return False

def GRLEX(alpha, beta):
    if not len(alpha) == len(beta):
        raise ValueError("Arguments not of same length")
    if sum(alpha) < sum(beta):
        return True
    if sum(alpha) > sum(beta):
        return False
    for iter in range(0,len(alpha)):
        if alpha[iter] < beta[iter]:
            return True
        if alpha[iter] > beta[iter]:
            return False
    return False

def GREVLEX(alpha, beta):
    if not len(alpha) == len(beta):
        raise ValueError("Arguments not of same length")
    if sum(alpha) < sum(beta):
        return True
    if sum(alpha) > sum(beta):
        return False
    n = len(alpha)
    for iter in range(0,n):
        if alpha[n-1-iter] < beta[n-1-iter]:
            return False
        if alpha[n-1-iter] > beta[n-1-iter]:
            return True
    return False

def INVLEX(alpha, beta):
    if not len(alpha) == len(beta):
        raise ValueError("Arguments not of same length")
    n = len(alpha)
    for iter in range(0,n):
        if alpha[n-1-iter] < beta[n-1-iter]:
            return True
        if alpha[n-1-iter] > beta[n-1-iter]:
            return False
    return False

def CREATE_MIXED_ORDER(args):
    if len(args) == 0:
        raise ValueError("No order is given")
    TOT = 0
    for iter in range(0, len(args)):
        hlp = args[iter]
        TOT = TOT + hlp[1]
    def AUX(alpha, beta):
        if not len(alpha) == TOT:
            return ValueError("Argument not of proper length")
        if not len(beta) == TOT:
            return ValueError("Argument not of proper length")
        TOTa = 0
        for iter in range(0,len(args)):
            hlp = args[iter]
            fun = hlp[0]
            if iter < len(args)-1:
                if alpha[TOTa:(TOTa + hlp[1])] == beta[TOTa:(TOTa + hlp[1])]:
                    TOTa = TOTa + hlp[1]
                    continue
            if not fun(alpha[TOTa:(TOTa + hlp[1])], beta[TOTa:(TOTa + hlp[1])]):
                return False
            if fun(alpha[TOTa:(TOTa + hlp[1])], beta[TOTa:(TOTa + hlp[1])]):
                return True
            TOTa = TOTa + hlp[1]
        return True
    return AUX


def WEIGHT_ORDER(alpha, beta, u, alt):
    N = len(u)
    if not len(alpha) == len(beta):
        raise ValueError("Arguments not of same length")
    if len(alpha) < len(u):
        raise ValueError("Arguments not of proper length")
    res1 = 0
    res2 = 0
    if not len(u) == 0:
        for iter in range(0,N):
            res1 = res1 + u[iter] * alpha[iter]
            res2 = res2 + u[iter] * beta[iter]
    if res1 < res2:
        return True
    if res1 > res2:
        return False
    if res1 == res2:
        return alt(alpha,beta)

def CREATE_ELIMINATION_ORDER(l):
    def AUX(alpha,beta):
        return WEIGHT_ORDER(alpha,beta, [1 for iter in range(0,l)], GREVLEX)
    return AUX


def PolynomialRingCreator():  
    AlreadyDone = []
    def check(lst):
        Res = []
        for l in lst:
            if l in Res:
                raise ValueError("Variable list contains duplicates")
            Res.append(l)
    def CreateRing(Variables,Rng):
        if not IsCommutativeRing(Rng):
            raise TypeError('Argument must be a commutative ring')
        if len(Variables) < 1:
            return Rng
        else:
            check(Variables)
            B = IsEuclideanDomain(Rng)
            C = IsField(Rng)
            class aux(object):
                VarlistInternal = copy.deepcopy(Variables)
                VarlistInternalTime = None
                OrderInternal = None
                OrderInternalTime=None
                def __init__(self,input):
                    if input == "default":
                        hulp = [0 for i in range(len(Variables)+1)]
                        hulp[len(Variables)] = Zero(Rng)
                        self.value = [hulp]
                    elif isinstance(input,int):
                        hulp = [0 for i in range(len(Variables)+1)]
                        hulp[len(Variables)] = Identity(Rng) * input
                        self.value = [hulp]
                    elif isinstance(input,str):
                        input = re.sub(' ', '',input)
                        input = StringReplace(input)
                        k = FindMainField(input)
                        L = len(input)
                        if k == "Primitive":
                            if input in Variables:
                                idx = Variables.index(input)
                                hulp = [0 for i in range(len(Variables)+1)]
                                hulp[idx] = 1
                                hulp[len(Variables)] = Identity(Rng)
                                self.value = [hulp]
                            else:
                                val = Rng(input)
                                hulp = [0 for i in range(len(Variables)+1)]
                                hulp[len(Variables)] = val
                                self.value = [hulp]
                        elif input[k] == "(":
                            hulp = aux(input[1:L-1])
                            self.value = hulp.value
                        elif input[k] == "-":
                            if k == 0:
                                hulp = - aux(input[1:])
                                self.value = hulp.value
                            else:
                                hulp = aux(input[:k]) - aux(input[(k+1):])
                                self.value = hulp.value
                        elif input[k] == "+":
                            hulp = aux(input[:k]) + aux(input[(k+1):])
                            self.value = hulp.value
                        elif input[k] == "*":
                            hulp = aux(input[:k]) * aux(input[(k+1):])
                            self.value = hulp.value
                        elif input[k] == "/":
                            hulp = aux(input[:k]) / aux(input[(k+1):])
                            self.value = hulp.value
                        elif input[k] == "^":
                            hulp = aux(input[:k]) ** Z(input[(k+1):])
                            self.value = hulp.value
                    elif ElementOf(input) == Rng:
                        hulp = [0 for i in range(len(Variables)+1)]
                        hulp[len(Variables)] = input
                        self.value = [hulp]
                    else:
                        raise TypeError('Type not recognized')
                
                def clean(self):
                    L = len(self.value)
                    hulp = []
                    for iter in range(0,L):
                            if not IsZero(self.value[iter][len(Variables)]):
                                hulp.append(self.value[iter])
                    if len(hulp) == 0:
                        hulp.append([0 for i in range(len(Variables)+1)])
                        hulp[0][len(Variables)] = Zero(Rng)
                    Varlist = self.CurrentVarList()
                    fun = self.CurrentOrder()

                    def CorrectOrder(alpha):
                        res = []
                        for iter in range(0,len(Varlist)):
                            RightVar = Varlist[iter]
                            pos = 0
                            for iter2 in range(0,len(Variables)):
                                if(Variables[iter2] == RightVar):
                                    pos = iter2
                                    break
                            res.append(alpha[pos])
                        return res
                    
                    def sortpoly(arg):
                        k = len(arg)
                        if (k == 0) or (k == 1):
                            return copy.deepcopy(arg)
                        k1 = math.floor(k/2)
                        k2 = k - k1
                        P1 = []
                        P2 = []
                        for s in range(0,k):
                            if s == k1:
                                continue
                            pivot = CorrectOrder(arg[k1])
                            if fun(CorrectOrder(arg[s]), pivot):
                                P2.append(arg[s])
                            else:
                                P1.append(arg[s])
                        P1 = sortpoly(P1)
                        P2 = sortpoly(P2)
                        res = []
                        res.extend(P1)
                        res.append(arg[k1])
                        res.extend(P2)
                        return res
                    
                    self.value = sortpoly(hulp)    
                
                def __add__(self,other):
                    if isinstance(other,aux):
                        hulp = self.value
                        hulp2 = other.value
                        Final = copy.deepcopy(hulp)
                        L = len(hulp2)
                        LL = len(hulp)
                        for iter in range(0,L):
                            Now = hulp2[iter]
                            Now2 = None
                            for iter2 in range(0,LL):
                                if hulp[iter2][:len(Variables)] == Now[:len(Variables)]:
                                    Now2 = iter2
                                    break
                            if Now2 == None:
                                Final.append(Now)
                            else:
                                Final[iter2][len(Variables)] = hulp[iter2][len(Variables)] + Now[len(Variables)]
                        Res = aux("default")
                        Res.value = copy.deepcopy(Final)
                        Res.clean()
                        return Res
                    elif ElementOf(other) == Rng:
                        return self + aux(other)
                    else:
                        raise TypeError('Arguments are of wrong types')
                    def __radd__(self,other):
                        return self + other
                def __neg__(self):
                    hulp = copy.deepcopy(self.value)
                    L = len(self.value)
                    for iter in range(0,L):
                        hulp[iter][len(Variables)] = - hulp[iter][len(Variables)]
                    Res = aux("default")
                    Res.value = copy.deepcopy(hulp)
                    Res.clean()
                    return Res
                def __sub__(self,other):
                    return self + (-other)
                def __mul__(self,other):
                    if isinstance(other,aux):
                        hulp = copy.deepcopy(self.value)
                        hulp2 = copy.deepcopy(other.value)
                        Final = []
                        for iter in range(0,len(hulp)):
                            Now = hulp[iter]
                            for iter2 in range(0,len(hulp2)):
                                Now2 = hulp2[iter2]
                                Now3 = copy.deepcopy(Now2)
                                for iter3 in range(0,len(Variables)):
                                    Now3[iter3] = Now[iter3] + Now2[iter3]
                                Now3[len(Variables)] = Now[len(Variables)]*Now2[len(Variables)]
                                Test = None
                                for iter4 in range(0,len(Final)):
                                    if Final[iter4][:len(Variables)] == Now3[:len(Variables)]:
                                        Test = iter4
                                        break
                                if Test == None:
                                    Final.append(Now3)
                                else:
                                    Final[iter4][len(Variables)] = Final[iter4][len(Variables)] + Now3[len(Variables)]
                        Res = aux("default")
                        Res.value = copy.deepcopy(Final)
                        Res.clean()
                        return Res
                    elif ElementOf(other) == Rng:
                        return self * aux(other)
                    elif isinstance(other,int):
                        return self * aux(other)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __rmul__(self,other):
                    return other * self
                def __pow__(self,other):
                    if isinstance(other,int):
                        if other < 0:
                            return(Inverse(self ** (-other)))
                        elif other == 0:
                            hulp = [0 for i in range(len(Variables)+1)]
                            hulp[len(Variables)] = Identity(Rng)
                            Res = aux("default")
                            Res.value = copy.deepcopy([hulp])
                            return Res
                        else:
                            if len(self.value) == 1:
                                hulp = copy.deepcopy(self.value[0])
                                for iter in range(0,len(Variables)):
                                    hulp[iter] = other*hulp[iter]
                                hulp[len(Variables)] = hulp[len(Variables)] ** other
                                Res = aux("default")
                                Res.value = copy.deepcopy([hulp])
                                return Res
                            else:
                                a,b,A = self, aux("default") ** 0, other
                                while A>0:
                                    if A % 2 == 1:
                                        b = b*a
                                    a = a*a
                                    A = math.floor(A/2)
                                return b    
                    else:
                        raise TypeError('Arguments are of wrong types')
                def __truediv__(self,other):
                    if isinstance(other,aux):
                        if not len(other.value) == 1:
                            raise ValueError('Division is not well defined')
                        else:
                            hulp = copy.deepcopy(other.value[0])
                            for iter in range(0,len(Variables)):
                                if not hulp[iter] == 0:
                                    raise ValueError('Division is not well defined')
                            divider = hulp[len(Variables)]
                            hulp2 = copy.deepcopy(self.value)
                            for iter in range(0,len(hulp2)):
                                if not Divides(divider, hulp2[iter][len(Variables)]):
                                    raise ValueError('Division is not well defined')
                                hulp2[iter][len(Variables)] = hulp2[iter][len(Variables)] // divider
                            Res = aux("default")
                            Res.value = copy.deepcopy(hulp2)
                            Res.clean()
                            return Res
                    if ElementOf(other) == Rng:
                        return self / aux(other)
                    else:
                        raise TypeError('Arguments are of wrong types')
                def divides(self,b):
                    if isinstance(self,aux):
                        if not len(self.value) == 1:
                            return False
                        hulp = self.value[0]
                        for iter in range(0,len(Variables)):
                            if not hulp[iter] == 0:
                                return False
                        divider = hulp[len(Variables)]
                        hulp2 = b.value
                        for iter in range(0,len(hulp2) - 1):
                            if not Divides(divider, hulp2[iter][len(Variables)]):
                                return False
                        return True
                def Invert(self):
                    if not(self.value) == 1:
                        raise ValueError('Element not invertible')
                    hulp = copy.deepcopy(self.value[0])
                    for iter in range(0,len(Variables)):
                        if not hulp[iter] == 0:
                            raise ValueError('Element not invertible')
                    hulp[len(Variables)] = Invert(hulp[len(Variables)])
                    Res = aux("default")
                    Res.value = copy.deepcopy([hulp])
                    return Res
                def ChangeOrderVar(self, newvar):
                    check(newvar)
                    if not ListEqualAsSet(newvar,Variables):
                        raise ValueError("List of new variables does not coincide with old variables")
                    self.VarlistInternal = copy.deepcopy(newvar)
                    self.VarlistInternalTime = time.monotonic()
                    self.clean()
                def ChangeOrder(self,neworder):
                    self.OrderInternal = neworder
                    self.OrderInternalTime = time.monotonic()
                    self.clean()
                def CurrentVarList(self):
                    if self.VarlistInternalTime == None:
                        Varlist = copy.deepcopy(VarlistExternal)
                    elif VarlistExternalTime <= self.VarlistInternalTime:
                        Varlist = copy.deepcopy(self.VarlistInternal)
                    else:
                        Varlist = copy.deepcopy(self.VarlistExternal)
                    return copy.deepcopy(Varlist)
                def CurrentOrder(self):
                    if self.OrderInternalTime == None:
                        fun = OrderExternal
                    elif OrderExternalTime <= self.OrderInternalTime:
                        fun = self.OrderInternal
                    else:
                        fun = self.OrderExternal
                    return fun
                
                def multideg(self):
                    self.clean()
                    return self.value[0][0:len(Variables)]
                
                def LC(self):
                    self.clean()
                    return self.value[0][len(Variables)]
                
                def LT(self):
                    self.clean()
                    pol = aux("default")
                    pol.value = [copy.deepcopy(self.value[0])]
                    return pol
                
                def LM(self):
                    pol = self.LT()
                    pol.value[0][len(Variables)] = Identity(Rng)
                    return pol
                    
                def ObtainCoefficient(self,mydict):
                    Res = []
                    for V in Variables:
                        if V in mydict:
                            Res.append(mydict[V])
                        else:
                            Res.append(0)
                    for iter in range(0,len(self.value)):
                        if self.value[iter][0:len(Variables)] == Res:
                            return self.value[iter][len(Variables)]
                    return Zero(Rng)
                        
                def __repr__(self):
                    Res = ""
                    self.clean()
                    hulp = self.value
                    if self.VarlistInternalTime == None:
                        Varlist = copy.deepcopy(VarlistExternal)
                    elif VarlistExternalTime <= self.VarlistInternalTime:
                        Varlist = copy.deepcopy(self.VarlistInternal)
                    else:
                        Varlist = copy.deepcopy(self.VarlistExternal)
                    if hulp[0][:len(Variables)] == [0 for i in range(len(Variables))]:
                        Res = Res + str(hulp[0][len(Variables)])
                    else:                   
                        if not IsIdentity(hulp[0][len(Variables)]):
                            if IsOrderedRing(Rng) and hulp[0][len(Variables)] == -Identity(Rng):
                                Res = Res + "-"
                                test1 = True
                            else: 
                                Res = Res + str(hulp[0][len(Variables)])
                                test1 = False
                        else:
                            test1 = True
                        for iter2 in range(0,len(Variables)):
                            iter3 = Variables.index(Varlist[iter2])
                            now = hulp[0][iter3]
                            if not now == 0:
                                if now == 1:
                                    if test1:
                                        Res = Res + Variables[iter3]
                                        test1 = False
                                    else: Res = Res + "*" + Variables[iter3]
                                else:
                                    if test1:
                                        Res = Res + Variables[iter3] + "^" + str(now)
                                        test1 = False
                                    else:
                                        Res = Res + "*" + Variables[iter3] + "^" + str(now)
                    if len(hulp) == 1:
                        if hulp[0][:len(Variables)] == [0 for i in range(len(Variables))]:
                            return Res
                        else:
                            return "(" + Res + ")"
                    else:
                        for iter in range(1,len(hulp)):
                            if hulp[iter][:len(Variables)] == [0 for i in range(len(Variables))]:
                                if IsOrderedRing(Rng):
                                    if hulp[iter][len(Variables)] < 0:
                                        Res = Res + str(hulp[iter][len(Variables)])
                                    else:
                                        Res = Res + "+" + str(hulp[iter][len(Variables)])
                                else:
                                    Res = Res + "+" + str(hulp[iter][len(Variables)])
                            else:
                                if not IsIdentity(hulp[iter][len(Variables)]):
                                    if (Rng == Z or Rng == R) and hulp[iter][len(Variables)] == -1:
                                        Res = Res + "-"
                                        test1 = True
                                    elif IsOrderedRing(Rng):
                                        if hulp[iter][len(Variables)] < 0:
                                            Res = Res + str(hulp[iter][len(Variables)])
                                        else:
                                            Res = Res + "+" + str(hulp[iter][len(Variables)])
                                            
                                        test1 = False
                                    else: 
                                        Res = Res + "+" + str(hulp[iter][len(Variables)])
                                        test1 = False
                                else:
                                    test1 = True
                                    Res = Res + "+"
                                for iter2 in range(0,len(Variables)):
                                    iter3 = Variables.index(Varlist[iter2])
                                    now = hulp[iter][iter3]
                                    if not now == 0:
                                        if now == 1:
                                            if test1:
                                                Res = Res + Variables[iter3]
                                                test1 = False
                                            else: Res = Res + "*" + Variables[iter3]
                                        else:
                                            if test1:
                                                Res = Res + Variables[iter3] + "^" + str(now)
                                                test1 = False
                                            else:
                                                Res = Res + "*" + Variables[iter3] + "^" + str(now)
                        return "(" + Res + ")"
                            
                        
                def ElementOf(self):
                    return internal
                def __eq__(self, other):
                    if isinstance(other, aux):
                        hulp = copy.deepcopy(self.value)
                        hulp2 = copy.deepcopy(other.value)
                        if not len(hulp) == len(hulp2):
                            return False
                        else:
                            for iter in range(0,len(hulp)):
                                test = None
                                for iter2 in range(0,len(hulp2)):
                                    if hulp[iter] == hulp2[iter2]:
                                        test = iter2
                                        break
                                if test == None:
                                    return False
                                
                            return True
                    if isinstance(other,str):
                        return False
                    if ElementOf(other) == Rng:
                        return self == aux(other)
                    else:
                        return False
                    
            VarlistExternal = copy.deepcopy(Variables)
            VarlistExternalTime = time.monotonic()  
            OrderExternal = GRLEX
            OrderExternalTime = time.monotonic()
            def internal(input):
                if input == "name":
                    res = Rng("name") + "["
                    for iter in range(0,len(Variables)):
                        res = res + Variables[iter]
                        if not iter == len(Variables)-1:
                            res = res + ","
                        else:
                            res = res + "]"
                    return(res)
                elif input == "Zero":
                    return aux(Zero(Rng))
                elif input == "Identity":
                    return aux(Identity(Rng))
                elif input == "HasIdentity":
                    return True
                elif input == "IsRng":
                    return True
                elif input == "HasZeroDivisors":
                    return Rng("HasZeroDivisors")
                elif input == "HasInverses":
                    return False
                elif input == "HasDivAlgo":
                    return IsField(Rng) and len(Variables) == 1
                elif input == "IsCommutative":
                    return True
                elif input == "GroundRing":
                    return Rng
                elif input == "Variables":
                    return VarlistExternal
                elif input == "Order":
                    return OrderExternal
                elif input == "Structure":
                    return "PolynomialRing"
                elif input == "ChangeOrderVar":
                    def reordering(NewVar):
                        check(NewVar)
                        VarlistExternal = copy.deepcopy(Variables)
                        varlistExternalTime = time.monotonic()
                    return reordering
                elif input == "ChangeOrder":
                    def reorder(neworder):
                        OrderExternal = neworder
                        OrderExternalTime = time.monotonic()
                    return reorder
                else:
                    return aux(input)
            if AlreadyDone == []:
                AlreadyDone.append([Variables,Rng,internal])
                return internal
            else:
                for iter in range(0,len(AlreadyDone)):
                    if ListEqualAsSet(AlreadyDone[iter][0],Variables) and AlreadyDone[iter][1] == Rng:
                        Res = AlreadyDone[iter][2]
                        Res("Reorder")(Variables)
                AlreadyDone.append([Variables,Rng,internal])
                return internal
    return CreateRing


def RecastPolynomial(PolyRing2, Poly):
    VarList1 = []
    VarList2 = []
    PolyRing1 = ElementOf(Poly)
    Loop1 = True
    RecRing1 = PolyRing1
    while Loop1:
        if RecRing1("Structure") == "PolynomialRing":
            VarList1.extend(RecRing1("Variables"))
            RecRing1 = RecRing1("GroundRing")
        else:
            Loop1 = False
    Loop2 = True
    RecRing2 = PolyRing2
    while Loop2:
        if RecRing2("Structure") == "PolynomialRing":
            VarList2.extend(RecRing2("Variables"))
            RecRing2 = RecRing2("GroundRing")
        else:
            Loop2 = False
    if not RecRing1 == RecRing2:
        raise ValueError("Underlying Ground Rings not the same")
    for var in VarList1:
        Bool = False
        for var2 in VarList2:
            if var == var2:
                Bool = True
                break
        if not Bool:
            raise ValueError("Not all variables in polynomial ring coincide")
    return PolyRing2(str(Poly))


def ListEqualAsSet(list1, list2):
    L = len(list1)
    if not L == len(list1):
        return False
    if L == 0:
        return True
    for iter in list1:
        if not iter in list2:
            return False
    return True

def Substitute(Poly, mydict, NewVar):
    VarList = []
    PolyRing = ElementOf(Poly)
    Loop = True
    RecRing = PolyRing
    while Loop:
        if RecRing("Structure") == "PolynomialRing":
            VarList.extend(RecRing("Variables"))
            RecRing = RecRing("GroundRing")
        else:
            Loop = False
    VarListNew = []
    for V in VarList:
        if not V in mydict:
            VarListNew.append(V)
    VarListNew.extend(NewVar)
    NewRing = PolynomialRing(VarListNew, RecRing)
    string = str(Poly)
    newstring = ""
    for iter in range(0,len(string)):
        if string[iter] in mydict:
            newstring = newstring + "(" + str(mydict[string[iter]]) + ")"
        else:
            newstring = newstring + string[iter]
    return NewRing(newstring)

def MonomialDivides(mon1, mon2):
    if not HasSameCategory(mon1, mon2):
        return TypeError("Monomials are not elements of the same polynomial ring")
    if len(mon1.value) > 1:
        return ValueError("First argument is not a monomial")
    if len(mon2.value) > 1:
        return ValueError("Second argument is not a monomial")
    mon1info = mon1.value[0]
    mon2info = mon2.value[0]
    for iter in range(0,len(mon1info)):
        if iter == len(mon1info) - 1:
            return Divides(mon1info[iter], mon2info[iter])
        if not mon1info[iter] <= mon2info[iter]:
            return False
        
def MonomialDivision(mon2, mon1):
    if not HasSameCategory(mon1, mon2):
        return TypeError("Monomials are not elements of the same polynomial ring")
    if len(mon1.value) > 1:
        return ValueError("First argument is not a monomial")
    if len(mon2.value) > 1:
        return ValueError("Second argument is not a monomial") 
    PolRng = ElementOf(mon1)
    RES = PolRng("default")
    mon1info = mon1.value[0]
    mon2info = mon2.value[0]
    for iter in range(0,len(mon1info)):
        if iter == len(mon1info) - 1:
            RES.value[0][iter] = mon1info[iter] // mon2info[iter]
        else:
            RES.value[0][iter] = mon1info[iter]-mon2info[iter]
    return RES

def DivisionAlgorithm(poly, polylist):
    PolRng = ElementOf(poly)
    if len(polylist) == 0:
        return poly
    for pol in polylist:
        if not IsElementOf(pol, PolRng):
            ValueError("Polynomials not of correct ring")
    Quolist = [Zero(PolRng) for i in range(0,len(polylist))]
    Rem = Zero(PolRng)
    p = poly
    while not p == Zero(PolRng):
        i = 0
        divisionoccured = False
        while i<=len(polylist)-1 and divisionoccured == False:
            if MonomialDivides(polylist[i].LT(), p.LT()):
                div = MonomialDivision(polylist[i].LT(), p.LT())
                Quolist[i] = Quolist[i] + div
                p = p - div*polylist[i]
                divisionoccured = True
            else:
                i = i+1
        if divisionoccured == False:
            Rem = Rem + p.LT()
            p = p - p.LT()
    return {"Quotient":Quolist, "Remainder":Rem}


PolynomialRing = PolynomialRingCreator()
A = PolynomialRing(["X","Y","Z"],Q)
A("ChangeOrder")(LEX)
MonomialDivides(A("X*Y*Z"), A("X^2*Y^2"))
DivisionAlgorithm(A("X^7*Y^2 + X^3*Y^2 - Y + 1"), [A("X*Y^2 - X"), A("X-Y^3")])

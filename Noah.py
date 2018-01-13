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
            if input == "default":
                self.value = 0
            elif isinstance(input,int):
                self.value = input % n
            elif isinstance(input,str):
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
                    hulp = aux(input[:k]) ** int(input[(k+1):])
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
            return not is_prime(n)
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
        return not typ.HasZeroDivisors()
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
            

def RingOfFractions(IntDom):  
    T = type(IntDom)
    if not IsIntegralDomain(IntDom):
         raise TypeError('Argument must be an element of an integral domain')
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
                    k = FindMainField(input)
                    L = len(input)
                    if k == "Primitive":
                        self.numerator = T(input)
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
                    hulp.denomonator = denom
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
                    hulp.denomonator = denom
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
                    hulp.numerator = num
                    hulp.denominator = denom
                    return(hulp)
                else:
                    raise TypeError('Arguments are of wrong types')
            def divides(self,b):
                return IsZero(b)
            def Invert(self):
                return self.ReturnIdentity() / self
            def __neg__(self):
                hulp = aux("default")
                hulp.numerator = -self.numerator
                hulp.denominator = self.denominator
                return hulp
            def __repr__(self):
                return str(self.value)
            def ReturnZero(self):
                hulp = aux("default")
                hulp.numerator = Zero(IntDom)
                hulp.denominator = Identity(IntDom)
                return(hulp)
            def ReturnIdentity(self):
                hulp = aux("default")
                hulp.numerator = Identity(IntDom)
                hulp.denominator = Identity(IntDom)
                return(hulp)
            def Commutative(self):
                return True
            def Identity(self):
                return True
            def IsRng(self):
                return True
            def HasZeroDivisors(self):
                False
            def HasInverses(self):
                True
            def HasDivAlgo(self):
                True
            def __repr__(self):
                if IsInvertible(self.denominator):
                    self.numerator = self.numerator // self.denominator
                    self.denominator = 1
                    return(str(self.numerator))
                else:
                    return "(" + str(self.numerator) + "/" + str(self.denominator) + ")"
    return aux


def MatrixRing(nrow, ncol, Rng):
    T = type(Rng)
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
                                hulp[rw][cl] = T(entry)
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
    return aux



Q = RingOfFractions(0)



def GetEntry(mat, rw, cl):
    return mat.getitem(rw,cl)
    
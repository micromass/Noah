
def calcfactorial(n):
    """Calculates the factorial"""
    if n == 1:
        return 1
    else:
        return n * calcfactorial(n-1)


def calcbinomial_form(n,r):
      """Calculates the binomial with formula"""
      return calcfactorial(n)/(calcfactorial(n-r)*calcfactorial(r))


def calcbinomial_multiply(n,r):
    """ Using method of tail recursion with multiplication """
    if r == 1:
        return n
    else:
        return (n-r+1)*calcbinomial_multiply(n,r-1)/r

def calcbinomial_tail_efficienct(n,r):
    list_old = [1] + [0] * (n)
    list_new = [0] * (n + 1)


    for i in range(1,n + 1):
        for k in range(i+1):
            if k == 1:
                list_new[k] = i
            elif k == 0  or k == i:
                list_new[k] = 1
            else:
                list_new[k] = list_old [k] + list_old  [k-1]
        for step in range(i+1):
            list_old [step] = list_new[step]

    return(list_new[r])




def calcbinomial_tail(n,r):
    """Calculates binomial with tail recurssion"""
    if n == r:
        return 1
    elif r == 1:
        return n
    else:
       return calcbinomial_tail(n - 1,r) + calcbinomial_tail(n - 1, r - 1)


def calcbinomial_reduce(n,r):
    """n*(n-1) * ... * (n-k+1)/ k! """
    denominator = calcfactorial(r)
    numerator = 1
    i = n - r + 1
    while n >= i:
        numerator *= n
        n = n - 1
    return(numerator/denominator)


def prime_list(n):
    """ Gives a list with the primes of the number"""
    primes = []
    k = 2
    while k < n:
        if n % k == 0:
            primes.append(k)
            n /= k
        else:
            k += 1
    return primes


def is_prime(n):
    """ Checks a number for primality"""
    if len(prime_list(n)) == 0:
        return True
    else:
        return False


def average_list(list):
    """Computes the average of a list"""
    sum = 0
    for number in list:
        sum += number

    average = sum/len(list)

    return(average)

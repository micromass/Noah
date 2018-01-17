
def calcfactorial(n):
    """Calculates the factorial"""
    if n == 1:
        return 1
    else:
        return n * calcfactorial(n-1)

def calcbinomial(n,r):
    """Calculates the binomial"""
    return calcfactorial(n)/(calcfactorial(n-r)*calcfactorial(r))

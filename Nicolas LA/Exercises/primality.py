from operations import is_prime
while True:
    n = int(input("Which number do you want to check for primality: "))

    if is_prime(n) == True:
        print(f"{n} is a prime number")
    elif is_prime(n) == False:
        print(f"{n} is not a prime number")

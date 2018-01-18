from operations import calcbinomial_tail, calcbinomial_multiply, calcbinomial_form, calcbinomial_reduce

while True:
    bi_1 = input("Input binomial coefficient in the form (n,r): ")
    bi_2 = bi_1.replace("(","")
    bi_def= bi_2.replace(")","")
    data = bi_def.split(",")

    n = int(data[0])
    r = int(data[1])

    print("Tail ", calcbinomial_tail(n,r))
    print("Multiply ", calcbinomial_multiply(n,r))
    print("Form ", calcbinomial_form(n,r))
    print("Reduce ",calcbinomial_reduce(n,r))

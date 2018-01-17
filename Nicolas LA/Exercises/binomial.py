from operations import calcbinomial

while True:
    bi_1 = input("Input binomial coefficient in the form (n,r): ")
    bi_2 = bi_1.replace("(","")
    bi_def= bi_2.replace(")","")
    data = bi_def.split(",")

    n = int(data[0])
    r = int(data[1])

    print(calcbinomial(n,r))

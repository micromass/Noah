from operations import average_list

n = 4
list = []
while True:
    print("""Please type all the numbers you want to take the average of.
    Follow each number with an enter. After you are done type 'ok'
    """)

    while True:
        n = input("> ")
        if n == "ok":
            break
        else:
            list.append(int(n))


    print(average_list(list))

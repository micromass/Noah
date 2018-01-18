
from numpy import *



def g(g_raw,num):
    x = num
    return eval(g_raw)

def op_list(g, g_string,old_list, new_list):
    length = len(old_list)
    for i in range(length):
        new_list.append(g(g_string,old_list[i]))
    print(old_list, "-------->",new_list)



old_list = []
new_list = []
print("Please type all the numbers you want to perform the function on.")
print("Follow each number with an enter. After you are done type 'ok'")
while True:
    n = input("> ")
    if n == "ok":
        break
    else:
        old_list.append(float(n))

g_string = input("Please input your function using standard python notation and x as a variable: ")
op_list(g,g_string,old_list, new_list)

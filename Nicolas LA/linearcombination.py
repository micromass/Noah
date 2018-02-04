#This program verifies whether one vector in R^3 is a linear combination of combination of two other vectors

def is_linear_combination(v1, v2, v3):
    try:
        scalar_3 = (v1[1]*v2[0] - v2[1]*v1[0])/(v3[1]*v2[0]-v3[0]*v2[1])
        scalar_2 =(v1[0]/v2[0]) - (v3[0]/v2[0])*scalar_3
        print(f"The first vector is a combination of the other two and the scalars are {scalar_2} {scalar_3}")
    except ZeroDivisionError:
        print("The first vector is not a combination of the other two")



v1 = [-2 , 0 , 3]
v2 = [1 , 3 , 0]
v3 = [2 , 4 , -1]

is_linear_combination(v1,v2,v3)

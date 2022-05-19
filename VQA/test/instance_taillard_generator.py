
import random


def generator(min,max):
    x = random.randint(0,2**31-1)
    a,b,c = 16807,127773,2836
    m = 2**31-1
    k = x/b
    x = a*(x % b) - k*c
    if x < 0 : x =x + m
    gen = x/m
    return int(min + gen * (max+1-min))

def instance_generator(N,M,min,max):
    pm = [[0 for i in range(N)] for j in range(M)]
    for i in range(M):
        for j in range(N):
            pm[i][j] = generator(min,max)
    return pm 


print(instance_generator(2,3,1,5))  
    





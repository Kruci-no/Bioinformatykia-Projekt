import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62

str1 = "AAAABBBBCCCCMMMMAAABBBWWW"
str2 = "AAAABZAABCPANMMMAAABBBWWW"
import time

start_time = datetime.now()
start_time = time.time()
print("Iteration test for linear meamory algorytm")
x = "%s" % (time.time() - start_time)
next_iter = 1
while(next_iter):
    start_time = time.time()
    go =  0
    ge = -1
    alignment([str1, str2],go,ge,matrix,ret_max = 1,linear_memory=True)
    x = "%s" % (time.time() - start_time)
    x =float(x)
    print("time taken",x)
    print(len(str1),len(str2))
    next_iter = input("Do you want to try next iteration")
    if(x > 10):
        break
    str1 = str2 + str1
    str2 = str1 + str2

   
"""
output:
Iteration test for linear meamory algorytm
time taken 0.07812833786010742
25 25

Do you want to try next iteration1
time taken 0.0781254768371582
50 75

Do you want to try next iteration1
time taken 0.26563596725463867
125 200

Do you want to try next iteration1
time taken 0.9219188690185547
325 525

Do you want to try next iteration1
time taken 5.103905439376831
850 1375

Do you want to try next iteration1
time taken 26.65145969390869
2225 3600

The change of memory used by program was not so big after each iteration
""""
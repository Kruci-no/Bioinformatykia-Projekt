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
print("Iteration test for afinic gap algorytm")
x = "%s" % (time.time() - start_time)
next_iter = 1
while(next_iter):
    start_time = time.time()
    go = -1
    ge = -1
    alignment([str1, str2],go,ge,matrix,ret_max = 1,linear_memory=False)
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
output
Iteration test for afinic gap algorytm
time taken 0.0
25 25

Do you want to try next iteration1
time taken 0.0312497615814209
50 75

Do you want to try next iteration1
time taken 0.18750476837158203
125 200

Do you want to try next iteration1
time taken 1.250060796737671
325 525

Do you want to try next iteration1
time taken 8.582213163375854
850 1375

Do you want to try next iteration1
time taken 58.51440906524658
2225 3600



"""
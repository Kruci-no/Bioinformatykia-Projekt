import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
str1 = "BBBWWW"
str2 = "AABWWW"
str3 = "ABBWCC"
import time

start_time = datetime.now()
start_time = time.time()
print("Iteration test for 3 seq ")
x = "%s" % (time.time() - start_time)
next_iter = 1
while(next_iter):
    start_time = time.time()
    go = 0
    ge = -1
    alignment([str1, str2, str3],go,ge,matrix,ret_max = 1,linear_memory=False)
    x = "%s" % (time.time() - start_time)
    x =float(x)
    print("time taken",x)
    print(len(str1),len(str2),len(str3))
    next_iter = input("Do you want to try next iteration")
    if(x > 10 or next_iter==0):
        break
    str1 =str1 + str2
    str2 =str2 + str3
    str3 =str3 + str3 
    
    """
    Iteration test for 3 seq 
    time taken 0.07812714576721191
    6 6 6
    
    Do you want to try next iteration1
    time taken 0.4687633514404297
    12 12 12
    
    Do you want to try next iteration1
    time taken 3.354745864868164
    24 24 24
    
    Do you want to try next iteration 1
    time taken 25.800235986709595
    48 48 48  
    """
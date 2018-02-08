import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
go= -1
ge= -1
print("afinic_gap" ,"simple_cost_function") 
start_time = datetime.now() 
pairwise2.align.globalxs(MHC.Homo_sapiens, MHC.Mus_musculus,-2,-1)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Homo_sapiens, MHC.Mus_musculus],go,ge,[1,0],ret_max = 2,local = False,linear_memory= False)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

for a in pairwise2.align.globalxs(MHC.Homo_sapiens, MHC.Mus_musculus,-2,-1):
    print(format_alignment(*a))
    break
print(X)
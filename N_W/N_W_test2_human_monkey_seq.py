import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
go= -1
ge= -1
print("afinic_gap") 
start_time = datetime.now() 
pairwise2.align.globalds(MHC.Homo_sapiens, MHC.Pan_troglodytes, matrix,-2,-1)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Homo_sapiens, MHC.Pan_troglodytes],go,ge,matrix,ret_max = 20)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

for a in pairwise2.align.globalds(MHC.Homo_sapiens, MHC.Pan_troglodytes, matrix,-2,-1):
    print(format_alignment(*a))
print(X)
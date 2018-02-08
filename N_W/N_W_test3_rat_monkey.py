import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
go= -1
ge= -1
print("afinic_local_gap") 
start_time = datetime.now() 
pairwise2.align.localds(MHC.Rattus_norvegicus[:100], MHC.Pan_troglodytes[:100], matrix,-2,-1)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Rattus_norvegicus[:100], MHC.Pan_troglodytes[:100]],go,ge,matrix,ret_max = 20,local = True)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

for a in pairwise2.align.localds(MHC.Rattus_norvegicus[:100], MHC.Pan_troglodytes[:100], matrix,-2,-1):
    print(format_alignment(*a))
print(X)

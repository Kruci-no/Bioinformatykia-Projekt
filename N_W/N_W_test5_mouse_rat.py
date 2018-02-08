import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
go= 0
ge= 0
print("A match score is the score of identical chars, otherwise mismatchscore") 
start_time = datetime.now() 
pairwise2.align.globalmx(MHC.Rattus_norvegicus, MHC.Mus_musculus, 0.1,-10)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Rattus_norvegicus, MHC.Mus_musculus],go,ge,[0.1,-10],ret_max = 2,local = False,linear_memory= False)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

for a in pairwise2.align.globalmx(MHC.Rattus_norvegicus, MHC.Mus_musculus, 0.1,-10):
    print(format_alignment(*a))
print(X)
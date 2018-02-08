import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from datetime import datetime 
matrix = mi.blosum62
go= 0
ge= -1
print("multiseq" ,"simple_cost_function") 
start_time = datetime.now() 
X=alignment([MHC.Homo_sapiens[0:40], MHC.Mus_musculus[0:40] ,MHC.Rattus_norvegicus[0:40] ],go,ge,[1,0],ret_max = 2,local = False,linear_memory= False)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))


print(X)
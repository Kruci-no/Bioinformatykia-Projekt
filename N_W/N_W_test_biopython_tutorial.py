import MHC_test as MHC
from N_W import alignment
from Bio.SubsMat import MatrixInfo as mi
from Bio.pairwise2 import format_alignment
import time
from datetime import datetime 
s=mi.blosum62


#X=alignment([MHC.Homo_sapiens,MHC.Pan_troglodytes],go,ge,s,ret_max = 1)
#print(X)
from Bio import pairwise2
#Bio python documentation example comparing
alignments = pairwise2.align.globalxx("ACCGT", "ACG")
print(format_alignment(*alignments[0]))
go=0
ge=0
X=alignment(["ACCGT","ACG"],go,ge,[1,0],ret_max = 1)
print(X)
for a in pairwise2.align.globalmx("ACCGT", "ACG", 2, -1):
    print(format_alignment(*a))
go=0
ge=0
X=alignment(["ACCGT","ACG"],go,ge,[2,-1],ret_max = 20)
print(X)

for a in pairwise2.align.globalms("ACCGT", "ACG", 2, -1, -.5, -.1):
    print(format_alignment(*a))

go= -.4
ge= -.1    
X=alignment(["ACCGT","ACG"],go,ge,[2,-1],ret_max = 20)
print(X)


matrix = mi.blosum62
for a in pairwise2.align.globaldx("KEVLA", "EVL", matrix):
    print(format_alignment(*a))

go=0
ge=0    
X=alignment(["KEVLA", "EVL"],go,ge,matrix,ret_max = 20)
print(X)

matrix = mi.blosum62

go=0
ge=0 

start_time = datetime.now() 
pairwise2.align.globaldx(MHC.Homo_sapiens, MHC.Pan_troglodytes, matrix)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Homo_sapiens, MHC.Pan_troglodytes],go,ge,matrix,ret_max = 20)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

go= -2
ge= -1
print("afinic_gap") 
pairwise2.align.globalds(MHC.Homo_sapiens, MHC.Pan_troglodytes, matrix,go,ge)
time_elapsed = datetime.now() - start_time 
print('biopython','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
start_time = datetime.now() 
X=alignment([MHC.Homo_sapiens, MHC.Pan_troglodytes],go,ge,matrix,ret_max = 20)
time_elapsed = datetime.now() - start_time 
print('our implemation','Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))



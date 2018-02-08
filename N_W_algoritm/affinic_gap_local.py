import numpy as np
from matching_function import create_matching_function

def init_matrix_H(seqs_lenght,go ,ge):
    """
    seqs lenght - array with contain lenght of all sequence
    go - gap open penalyty
    ge - gap extextion penalyty
    """
    H = np.zeros(tuple(seqs_lenght))   
    return H
def init_matrix_E(seqs_lenght, go ,ge):
    """
    seqs lenght - array with contain lenght of all sequence
    go - gap open penalyty
    ge - gap extextion penalyty
    """
    E = np.zeros(tuple(seqs_lenght))
    for i in range(1, seqs_lenght[0]):
        E[i,0] = ge * i  + go
    for i in range( seqs_lenght[1]):
        E[0,i] = - float("inf")
    return E

def init_matrix_F(seqs_lenght, go ,ge ):
    """
    seqs lenght - array with contain lenght of all sequence
    go - gap open penalyty
    ge - gap extextion penalyty
    """
    F = np.zeros(tuple(seqs_lenght))
    #for i in range(1, seqs_lenght[0]):
    #    F[0,i] = ge * i  + go
    for i in range( seqs_lenght[0]):
        F[i,0] = - float("inf")
    for i in range(1, seqs_lenght[1]):
        F[0,i] = ge * i  + go

    return F

def fill_H(H ,E ,F ,seq1 , seq2, x, y, d):
    H[x,y] = max([H[x-1][y-1] + d[(seq1[x-1],seq2[y-1])],E[x,y], F[x,y] , 0])
    
def fill_E(H ,E , x, y, go, ge):
    E[x,y] = max([H[x-1][y] + go + ge ,E[x-1][y] + ge]) 
    
def fill_F(H ,F , x, y, go, ge):
    F[x,y] = max([H[x][y - 1] + go + ge ,F[x][y - 1] + ge])
    
def fill_all(H ,E ,F ,seqs, go, ge, d):
    
    for x in range(1,H.shape[0]):
        for y in range(1,H.shape[1]):
            fill_E(H ,E , x, y, go, ge)
            fill_F(H ,F , x, y, go, ge)
            fill_H(H ,E ,F ,seqs[0] , seqs[1], x, y, d)
    
def matching(H,E,F,seqs,d,max_number_of_matching):
    seq1 = seqs[0]
    seq2 = seqs[1]
    H_arg_maxs = np.argwhere(H == np.max(H))
    array=[(["" for seq in seqs] ,tuple(arg_max) ) for arg_max in H_arg_maxs]
    ended_matching_array=[]
    while(array):
        next_array=[]
        for record in array:
            index = record[1]
            matching = record[0]
            if(H[index] == 0):
                ended_matching_array.append( ([matching[i]   for i in range(len(seqs))],index ) )
            else:
                x = index[0]
                y = index[1]
                if(H[x,y] == H[x-1,y-1] + d[(seq1[x-1],seq2[y-1])] ):
                    symbols = (seq1[x-1] , seq2[y-1])
                    next_matching = [symbols[i] + matching[i]   for i in range(len(seqs))]
                    next_array.append((next_matching,(x-1,y-1)))
                if(H[x,y] == E[x,y]):
                    codition = True
                    symbols = ("","")
                    while(codition):
                        symbols = (seq1[x-1] + symbols[0], "_" + symbols[1])
                        if(H[x-1][y] + go + ge == E[x,y]):
                            next_matching = [symbols[i] + matching[i]   for i in range(len(seqs))]
                            next_array.append((next_matching,(x - 1,y)))
                        codition = (E[x,y] == E[x-1][y] + ge)
                        x = x - 1 
                    pass
                if(H[x,y] == F[x,y]):
                    codition = True
                    symbols = ("","")
                    while(codition):
                        symbols = ("_" + symbols[0], seq2[y - 1] + symbols[1])
                        if(H[x][y - 1] + go + ge == F[x ,y]):
                            next_matching = [symbols[i] + matching[i]   for i in range(len(seqs))]
                            next_array.append((next_matching,(x ,y - 1)))
                        codition = (F[x,y] == F[x][y - 1] + ge)
                        y = y - 1
        array = next_array[0:max_number_of_matching - len(ended_matching_array)]
    return ended_matching_array



def linear_gap_local_algorytm(seqs, go, ge , s ,max_number_of_matching):
    H = init_matrix_H([len(s) + 1 for s in seqs],go , ge)
    E = init_matrix_E([len(s) + 1 for s in seqs],go , ge)
    F = init_matrix_F([len(s) + 1 for s in seqs],go , ge)
    d = create_matching_function(s , ge)
    fill_all(H, E ,F ,seqs , go, ge ,d)
    X = matching(H,E,F,seqs,d,max_number_of_matching)
    return X

if __name__ == "__main__":
    print()
    s1= "AAAA"
    s2= "AA"
    seqs = [s1, s2 ]
    ge = - 1
    go = - 3
    X = linear_gap_local_algorytm(seqs, go, ge , [1,-1],100)
    print(X)
    #H = init_matrix_F([len(s) + 1 for s in seqs],go , ge)
    #print(H)
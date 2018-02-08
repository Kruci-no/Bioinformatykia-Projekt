import numpy as np
from matching_function import create_matching_function
from multidimensional_N_W import multidimesional_N_W_algoritm
def init_matrix_H(seqs_lenght ,ge):
    """
    seqs lenght - array with contain lenght of all sequence
    go - gap open penalyty
    ge - gap extextion penalyty
    """
    H ={}
    for k, seq_leght in enumerate(seqs_lenght):
        for i in range(1, seq_leght):
            A = np.zeros(len(seqs_lenght))
            A[k] = i
            H[tuple(A.astype(int))] = ge * i 
    H[(0,0)] = 0        
    return H


def fill_H(H ,seq1 , seq2,ge ,x, y, d):
    H[(x,y)] = max([H[(x-1,y-1)] + d[(seq1[x-1],seq2[y-1])], H[(x - 1,y)] + ge, H[(x,y - 1)] + ge])
    del H[(x-1,y-1)]
    
    
def fill_all(H ,seqs,  ge, d):
    for x in range(1, len(seqs[0])//2 + 1):
        for y in range(1,len(seqs[1]) + 1):
            fill_H(H ,seqs[0] , seqs[1], ge , x, y, d)

    return np.array([H[(len(seqs[0])//2 ,i)] for i in range(len(seqs[1]) + 1)])
            
def linear_seqs_division(seqs,  ge , s ):
    H = init_matrix_H([len(s) + 1 for s in seqs], ge)
    d = create_matching_function(s , ge)
    P = fill_all(H ,seqs, ge, d)
    print(P)
    a = len(seqs[0])//2

    B = np.argwhere(P == np.max(P))

    array = []
    for b in B:
        b = b.item()
        seqs_0_begin = seqs[0][0:a]
        seqs_1_begin = seqs[1][0:b]
        seqs_0_end = seqs[0][a:]
        seqs_1_end = seqs[1][b:]
        array.append(([seqs_0_begin,seqs_1_begin],[seqs_0_end[::-1],seqs_1_end[::-1]]))
    return array



def merge(matches_begin, matches_end, max_number_of_matching):
    matches = []
    for match_begin in matches_begin:
        for match_end in matches_end:
            matches.append( [match_begin[0] + match_end[0][::-1],match_begin[1] + match_end[1][::-1] ] )
            if(len(matches)== max_number_of_matching):
                return matches
    return matches
def linear_algorytm(seqs,  ge , s, max_number_of_matching, L):
    if((len(seqs[0]) + 1) * (len(seqs [1]) + 1) <= L):
        return multidimesional_N_W_algoritm(seqs,  ge , s, max_number_of_matching)
    else:
        
        array = linear_seqs_division(seqs,  ge , s)
        maximum = -float("inf")
        print(array)
        for seqs1, seqs_rev1 in array:
            score_begin, matches_begin = linear_algorytm(seqs1,  ge , s, max_number_of_matching, L)
            score_end, matches_end = linear_algorytm(seqs_rev1,  ge , s, max_number_of_matching, L)
            if (maximum < score_begin +score_end):
                maximum = score_begin +score_end
                max_matches_begin = matches_begin
                max_matches_end = matches_end
        return maximum, merge(max_matches_begin, max_matches_end, max_number_of_matching)
        
if __name__ == "__main__":
    s1= "AAABB"
    s2= "AAA"
    seqs = [s1, s2 ]
    ge = - 1
    d = create_matching_function([1,-1] , ge)
    H = init_matrix_H([len(seq) + 1 for seq in seqs] ,ge)
    X = linear_algorytm(seqs,  ge , [1,-1], 100, 7*2)
    print(X)

    
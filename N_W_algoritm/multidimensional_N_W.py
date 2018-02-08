
from Bio.SubsMat import MatrixInfo as mi
import numpy as np
def init_matrix(seqs_lenght, g):
    """
    seqs lenght - array with contain lenght of all sequence
    g - gap penalyty
    """
    H = np.zeros(tuple(seqs_lenght))
    for k, seq_leght in enumerate(seqs_lenght):
        for i in range(1, seq_leght):
            A = np.zeros(len(seqs_lenght))
            A[k] = i
            H[tuple(A.astype(int))] = g * i * (len (seqs_lenght ) - 1)
    return H
        
def save_number_as_tuple(k,n):
    """
    k lenght of tuple
    n number to save in tuple
    """
    return tuple([(n % 2**(i+1))//2**(i) for i in range(0,k)])
    
def save_sybols_array(seqs, index, tuple_number):
    """
    seqs is array of seq that we want to match
    index is index of element H array that we want to compute value
    tuple_number indicate what case we are considering
    """
    g = lambda i : seqs[i][index[i] - 1] if tuple_number[i] % 2 == 1 else "_"
    return [g(i) for i in range(len(index))]

def fill_H(H, seqs, index, d):
    """
    d is dictionary with compute matching function for two symbols
    
    """
    maximum = (-1) * float("inf")
    for i in range(1,  2 ** len(seqs) ) :
         tuple_number = save_number_as_tuple(len(seqs),i)
         H_previous_index = tuple([index[i] - tuple_number[i] for i in range(len(seqs))])
         if(np.min(H_previous_index) > - 1):
             symbols = save_sybols_array(seqs,index,tuple_number)             
             x = q(symbols, d) + H[H_previous_index]
             if maximum < x :
                maximum = x
    H[index] = maximum 
        
def fill_H_all(H, seqs , d):
    seqs_lenght = np.array([len(seq) for seq in seqs])
    index = np.zeros(len(seqs_lenght))
    orgin = np.zeros(len(seqs_lenght))
    for i in range(np.prod(seqs_lenght + 1)  - 1):
        next_index(seqs_lenght,index)
        if np.sum(index == orgin) < len(seqs_lenght) -1:
            fill_H(H,seqs,tuple(index.astype(int)),d)
        
            
        
    
    
def next_index(seqs_lenght,index):
    for i in range(len(seqs_lenght)):
        if index[i]  != seqs_lenght[i]:
            index[i] = index[i] + 1
            return
        else :
            index[i] = 0
       
    


def q(symbols,d):
    """
    symbols is contener with symbols that we want to compute 
    matching function
    d is dictionary with compute matching function for two symbols
    return value of matching function for symbols array
    """
    S=[]
    for i in range(len(symbols)):
        for j in range(0,i):
            S.append((symbols[i], symbols[j]))      
    return sum([d[s] for s in S])
    
def create_matching_function(s, g):
    """
    s is matrix-dictioniany like blosum or 
    it is array with contain two elements score for
    matching or unmatching
    g gap penalyty
    """
    d={}
    d[('_','_')] = 0
    if isinstance(s,dict):
        
        for key in s.keys():
            d[key]=s[key]
            d[(key[1],key[0])]=d[key]
            d[('_',key[0])] = g
            d[(key[0],'_')] = g
            d[('_',key[1])] = g
            d[(key[1],'_')] = g
            return d
    else:
        
        for key in mi.blosum62.keys():
            if key[1] == key[0]:
                d[key] = s[0]
                d[(key[1],key[0])] = s[0]
            else :
                d[key] = s[1]
                d[(key[1],key[0])] = s[1]
            d[('_',key[0])] = g
            d[(key[0],'_')] = g
            d[('_',key[1])] = g
            d[(key[1],'_')] = g
        return d

def matching(H,seqs,d,max_number_of_matching):
    starting_index = tuple([len(seq)   for seq in seqs])
    starting_matching = ["" for seq in seqs]
    orgin = np.zeros(len(seqs))
    array=[(starting_matching ,starting_index)]
    ended_matching_array=[]
    while(array):
        next_array=[]
        for record in array:
            index = record[1]
            matching = record[0]
            if(np.sum(index) == 0):
                ended_matching_array.append([matching[i]   for i in range(len(seqs))])
            elif(np.sum(index == orgin)  >=  len(seqs) -1):
                n = np.max(index)
                for i in range(len(seqs)):
                    if(index[i] == 0):
                        matching[i] = n * "_" + matching[i]
                    else:
                        matching[i] = seqs[i][0:n] + matching[i]
                ended_matching_array.append([matching[i]   for i in range(len(seqs))])
                
            else:
                for i in range(1,  2 ** len(seqs) ) :
                    tuple_number = save_number_as_tuple(len(seqs),i)       
                    symbols = save_sybols_array(seqs,index,tuple_number)
                    H_previous_index = tuple([index[i] - tuple_number[i] for i in range(len(seqs))])
                    if(np.min(H_previous_index)>-1):
                        x = q(symbols, d) + H[H_previous_index]
                        if(x == H[tuple(index)]):
                            next_matching = [symbols[i] + matching[i]   for i in range(len(seqs))]
                            next_array.append((next_matching,H_previous_index))
        array = next_array[0:max_number_of_matching - len(ended_matching_array)]
    return ended_matching_array


def multidimesional_N_W_algoritm(seqs, g , s , max_mathing):
    H = init_matrix([len(s) + 1 for s in seqs],g)
    d = create_matching_function(s,g)
    fill_H_all(H, seqs , d)
    return H[tuple([len(seq) for seq in seqs])] ,matching(H,seqs,d,max_mathing)


if __name__ == "__main__":
    s1= "AAAAAA"
    s2= "AAA"
    s3= "A"
    tab = [s1, s2 ,s3]
    ge = -0.5
    go = 0
    s =[1,-3]
    max_mathing = 3
    X = multidimesional_N_W_algoritm(tab, g , s ,max_mathing)
    print(X)
    """H = init_matrix([len(s) + 1 for s in tab],g)

    d = create_matching_function([1,-3],g)

    fill_H_all(H, tab , d)
    print(H)
    print(matching(H,tab,d,100))"""

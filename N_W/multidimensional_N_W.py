from matching_function import create_matching_function
import numpy as np
def init_matrix(seqs_lenght, g):
    """
    seqs lenght - array with contain lenght of all sequences
    g - gap penalty
    return initial H matrix
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
    return n save in binar form in tuple of lenght k
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
    H is matrix that we are filling
    seqs is array of seq that we want to match
    index is index of element H matrix that we want to compute value
    d is dictionary that contain cost of matching between symbols
    function compute H[index]
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
    """
    H is matrix that we are filling
    seqs is array of seq that we want to match
    d is dictionary that contain cost of matching between symbols
    function compute whole H
    """
    seqs_lenght = np.array([len(seq) for seq in seqs])
    index = np.zeros(len(seqs_lenght))
    orgin = np.zeros(len(seqs_lenght))
    for i in range(np.prod(seqs_lenght + 1)  - 1):
        next_index(seqs_lenght,index)
        if np.sum(index == orgin) < len(seqs_lenght) -1:
            fill_H(H,seqs,tuple(index.astype(int)),d)
        
            
        
    
    
def next_index(seqs_lenght,index):
    """
    seqs_lenght is array containing lengths of sequences
    index is array that can be tranform into index of H matrix
    function add to index one on the first position if it posible 
    if not then first position is set to be zero and second position 
    is increase and so on
    """

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
    
def matching(H,seqs,d,max_number_of_matching):
    """
    H is filled matrix
    seqs is array of seq that we want to match
    d is dictionary with compute matching function for two symbols
    max_number_of_matching maximal number of matching
    function return array with contain optimal matching.
    this array is not bigger that max_number_of_matching
    """
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
    """
    seqs is array of seq that we want to match
    g is gap extention penalty
    s is substitution matrix or array [x,y] x is cost of matching 
    when symbols that we are matching are the same and y if they are not
    max_mathing is maximal number of matching that are returned
    function return pair containg cost of matching and matches
    """
    H = init_matrix([len(seq) + 1 for seq in seqs],g)
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
    s =[10,-3]
    max_mathing = 3
    X = multidimesional_N_W_algoritm(tab, ge , s ,max_mathing)
    print(X)

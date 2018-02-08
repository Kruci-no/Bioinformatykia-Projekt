import numpy as np
def create_matching_function(s, g):
    """
    s is matrix-dictioniany like blosum or 
    it is array with contain two elements score for
    matching or unmatching
    g gap penalyty
    return d dictionary containing containd cost of matching for pair of symbols
    """
    blosum6 = np.load('blosum62.npy').item()
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
        
        for key in blosum6.keys():
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
if __name__ == "__main__":
    s=[1,-1]
    g=-2

    d = create_matching_function(s,g)
    print(d[('_','M')])
    if(not d[('_','A')] == -2):
        raise AssertionError 
    if(not d[('_','_')] == 0):
        raise AssertionError 
    if(not d[('A','A')] == 1):
        raise AssertionError
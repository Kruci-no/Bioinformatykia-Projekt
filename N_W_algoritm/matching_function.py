# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 10:48:33 2018

@author: kuba
"""
from Bio.SubsMat import MatrixInfo as mi
def create_matching_function(s, g):
    """
    s is matrix-dictioniany like blosum or 
    it is array with contain two elements score for
    matching or unmatching
    g gap penalyty
    """
    d={}
    d[('_','_')] = g
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
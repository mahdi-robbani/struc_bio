#!/usr/bin/env python
# coding: utf-8

# In[1]:


##Adapted from Stefan E Seemann and Giulia I. Corsi
import numpy as np

def init_matrix(seq):
    '''Initializes empty score matrix'''
    l = len(seq)
    m = np.matrix([[0 for i in range(l)] for j in range(l)])
    return m

def get_score(i, j, seq, weighted):
    '''Given a position i and j, returns either a weighted or unweighted score'''
    score = 0
    pair = seq[i] + seq[j]
    #choose scoring system
    if weighted:
        scores = {'AU': 2, 'UA': 2, 'GU': 1, 'UG': 1, 'GC': 3, 'CG': 3}
    else:
        scores = {'AU': 1, 'UA': 1, 'GU': 1, 'UG': 1, 'GC': 1, 'CG': 1}
    #get score
    if pair in scores:
        score = scores[pair]
    return score

def max_score(m, i, j, loop_size, weighted):
    '''Returns the max value for a cell in a matrix
    Adapted from Giulia I. Corsi'''
    max_no_bf = max(m[i, j-1], m[i+1, j],
                    m[i+1, j-1] + get_score(i, j, seq, weighted))
    bf_scores = [m[i, k] + m[k+1, j] for k in range(i+1+loop_size, j-1-loop_size)]
    if len(bf_scores) > 0:
        return max(max_no_bf, max(bf_scores))
    else:
        return max_no_bf

def score_matrix(seq, con, loop_size, weighted):
    '''Returns the filled nussinov matrix for a given sequence
    Adapted from Stefan E Seemann and Giulia I. Corsi'''
    l = len(seq)
    m = init_matrix(seq) # initialzie matrix
    pairs = []
    
    if con == None: #if no constraint create a string of size l
        con = "." * l 
    
    for diag in range(loop_size+1, l): # the diagonal of the matrix to loop over
        for i in range(0, l-diag): # the entry on the diagonal to fill
            j = i + diag
            if con[i] == "(" and con[j] == ")":
                m[i, j] = m[i+1, j-1] + get_score(i, j, seq, weighted)
            elif con[i] == "x" or con[j] == "x": #check for unpaired sites and force skip
                m[i, j] = max(m[i, j-1], m[i+1, j])
            elif con[i] == "." and con[j] == ".":
                 m[i, j] = max_score(m, i, j, loop_size, weighted)
    return m
    
def backtrack(seq, matrix, loop_size, weighted):
    '''Perform backtracking on a matrix and return the dot bracket structure
    Adapted from Stefan E Seemann'''
    l = len(seq)
    m = matrix
    h = loop_size
    structure = ['.' for i in range(l)]
    stack = []
    stack.append((0,l-1))
    while len(stack) > 0:
        top = stack.pop(),
        i = top[0][0]
        j = top[0][1]
        if i >= j:
            continue
        elif m[i+1, j] == m[i, j]:
            stack.append((i+1,j))
        elif m[i, j-1] == m[i, j]:
            stack.append((i,j-1))
        elif m[i+1, j-1] + get_score(i, j, seq, weighted) == m[i, j]:
            structure[i] = "("
            structure[j] = ")"
            stack.append((i+1,j-1))
        else:
            for k in range(i+1+h, j-1-h):
                if m[i, k]+ m[k+1, j] == m[i, j]:
                    stack.append((k+1,j))
                    stack.append((i,k))
                    break
    return structure

def print_struc(mat, seq, backtrack):
    '''Prints the matrix, sequence, dot bracket structure and index'''
    l = len(seq)
    ind = []
    for i in range(l):
        ind.append(str(i % 10))
    
    print(mat)
    print("Seq:", seq, "Length:", l)
    print("DBS:", backtrack, "Score:", mat[0, l-1])
    print("Ind:", "".join(ind))
    
def nussinov(sequence, constraint = None, minimum_loop_size = 3, weighted = True):
    '''Takes in a sequence and optional constraint. Prints filled nussniov matrix,
    and prints and returns dotbracket structure
    '''
    seq = sequence
    con = constraint
    min_ls = minimum_loop_size
    
    mat = score_matrix(seq, con, min_ls, weighted)
    bac = backtrack(seq, mat, min_ls, weighted)
    bac = "".join(bac)
    
    print_struc(mat, seq, bac)
    return bac
        
def BP_set(string):
    '''creates a set of base pair positions'''
    stack = []
    pairs = []
    for pos, elem in enumerate(string):
        if elem == "(":
            stack.append(pos)
        if elem == ")":
            Open = stack.pop() #capitalize because open is python keyword
            Close = pos
            pairs.append([Open, Close])
    pairs_set = set(tuple(row) for row in pairs)
    return pairs_set

def BP_distance(string1, string2):
    '''calculates Base pair distance'''
    set1 = BP_set(string1)
    set2 = BP_set(string2)
    difference = len(set1-set2) + len(set2-set1)
    return difference


# In[2]:


if __name__=="__main__":
    
    seq = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
    con = ".........................(((((xxxxxxx))))).............................."
    
    print("====== Unconstrained sequeunce ======")
    seq_DBS = nussinov(seq, minimum_loop_size=3, weighted=True)
    
    print("\n" + "====== Constrained sequence ======")
    con_DBS = nussinov(seq, con, minimum_loop_size=3, weighted=True)
    
    BP_dist = BP_distance(seq_DBS, con_DBS)
    print("\n" + "Base Pair Distance:", BP_dist)


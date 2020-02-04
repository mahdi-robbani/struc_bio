#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Bio.PDB as PDB
import random
import numpy as np
import os
import matplotlib.pyplot as plt

def get_sidechain(res):
    '''Get a list of side chain atoms from a residue'''
    assert PDB.is_aa(res) #make sure residue is actually an amino acid
    sidechain = []
    exclude = ["N", "C", "O", "OXT"] #exclude nitrogen, carbonyl carbon, oxygen and special case
    for atom in res.get_atoms():
        if atom.get_id() in exclude or atom.element == "H": #ignore exclusion list and hydrogen atoms
            continue
        else:
            sidechain.append(atom)
    return sidechain #either returns list of side chain atoms or an empty list

def select_res(aa_name, prot_name, prot_AA_dict):
    """
    This function takes an AA name, a protein name and a nested dict:
    nested_dict = {Prot1: 
                        {AA1: [res1, res2], AA2: [res7, res8]}}
    It returns a randomly chosen res from the nested dict. Return empty list if no AA exists
    """
    inner_dict = prot_AA_dict[prot_name]
    res_list = inner_dict[aa_name] #get a list of residues
    while True: #loop forever
        if len(res_list) > 0: #if residue list is not empty, choose random residue
            res = random.choice(res_list)
            sc_list = get_sidechain(res) #get list of sidechain atoms
            if len(sc_list) == 0: 
                res_list.remove(res) #if no sidechain exists, remove residue from list and try again
            else:
                return sc_list #if side chain exists, return side chain
                break
        else: #if all residues removed, return empty list so that this protein is removed (not searched again)
            return [] 
    
def get_pair(aa_name, prot_AA_dict):
    """Takes in an AA name and dictionary and returns a pair of side chain lists"""
    prot_list = list(prot_AA_dict.keys())
    while True: #loop forever
        prot1, prot2 = random.sample(prot_list, k = 2) #sample 2 proteins without replacement
        SC1 = select_res(aa_name, prot1, prot_AA_dict) #get side chain
        SC2 = select_res(aa_name, prot2, prot_AA_dict) #get side chain
        if len(SC1) == 0: 
            prot_list.remove(prot1) #if side chain is empty, remove that protein from future searches
        elif len(SC2) == 0:
            prot_list.remove(prot2) #if side chain is empty, remove that protein from future searches
        elif len(SC1) != len(SC2):
            pass #if lengths are unequal, ignore both and sample again
        else:
            break #side chains are both non zero and the same length so end loop
    #assert (len(SC1) == len(SC2)), f'SC1 = {SC1}, SC2 = {SC2}' #ensure side chains are equal elngth
    SC_pair = [SC1, SC2]
    return(SC_pair)

def get_coord(sidechain):
    '''Returns a matrix of coordinates for a list of atoms'''
    coord_list = []
    for atom in sidechain:
        coord_list.append(atom.get_coord())
    coord_matrix = np.matrix(coord_list).T #Turn the list into a matrix and transpose it
    return coord_matrix

def get_center(coord_matrix):
    '''Returns the center of mass for a matrix of coordinates'''
    center_of_mass = coord_matrix.sum(1)/coord_matrix.shape[1] #calcualte center of mass of side chain
    centered_matrix = coord_matrix - center_of_mass #center the matrix
    return centered_matrix

def RMSD(a, b):
    '''This function takes in two matrices, each containing all the
    coordinates of a protein and returns the transformation matrix
    and the RMSD value'''
    #center a and b
    centered_a = get_center(a)
    centered_b = get_center(b)

    #calculate n
    n = np.size(centered_a)/len(centered_a)

    #calculate R
    R = centered_b * centered_a.T

    #SVD
    V, S, W_t = np.linalg.svd(R)
    
    #Calculate U
    U = W_t.T * V.T

    #check if U is a reflection (-1)
    if round(np.linalg.det(U)) == -1 :
        Z = np.diag([1,1,-1])
        U = W_t.T * Z * V.T #reflect U so that it is no longer a reflection
        S[2] = -S[2]

    #calculate E_0 and n
    E_0 = np.linalg.norm(centered_a)**2 + np.linalg.norm(centered_b)**2

    #calculate RMSD using the formula
    RMSD = np.sqrt(1/n * (E_0 - 2 * sum(S)))
    return(RMSD)

def get_AAlist(aa_name, prot):
    """This functions takes in an amino acids string name and a protein pdb structure and 
    returns a list of the AAs inside that structure"""
    aa_list = [] #create empty list
    assert PDB.is_aa(aa_name) # Make sure aa_name is amino acid (works for string or residue object)
    for res in prot.get_residues(): #loop through all residues in the protein
        if res.get_resname() == aa_name: #check the right AA is selected
            aa_list.append(res)
    return aa_list

def extract_AA(AA_list, directory = "top500H"):
    '''This function takes a list of AA names (3 letter code) and a directory
    in the same folder as this file    than contains all the PDB files'''
    #create parser
    parser = PDB.PDBParser(QUIET = True) #Quiet ignores warnings
    
    prot_dict = {} #dictionary of proteins
    for filename in os.listdir(directory): #loop through all files in directory
        try:
            structure = parser.get_structure(filename, os.path.join(directory, filename))
            AA_dict = {} #dictionary of amino acids
            for AA in AA_list:
                residues =  get_AAlist(AA, structure) #get list of residues
                AA_dict[AA] = residues #insert residues in amino acid list
            prot_dict[filename] = AA_dict #insert dictionary of amino acids into protein dictionary
        except (TypeError, ValueError) as e: #check for type error and value error
            print("Could not parse", filename)
    return prot_dict

def generate_data(AA_list, prot_AA_dict, N = 1000):
    '''This function takes in a list of AAs and a nested dictionary. For everu AA in the list,
    it creates N pairs of AAs from two different proteins'''
    RMSD_dict = {} #dict = {AA1 : {RMSD1, RMSD2}, AA2 : {RMSD1, RMSD2} ..}
    for AA in AA_list:
        RMSD_list = []
        for i in range(N): #Do N times
            pair = get_pair(AA, prot_AA_dict) #create 2 lists of side chain atoms
            C1 = get_coord(pair[0]) #get matrix of coodrinates for first list
            C2 = get_coord(pair[1]) #get matrix of coodrinates for second list
            R = RMSD(C1, C2) #calculate RMSD
            RMSD_list.append(R) #append to RMSD list
        RMSD_dict[AA] = RMSD_list # Assign RMSD list to corresponding protein
    return RMSD_dict

def save_histogram(rlist, aa):
    """
    Make histogram and save as png file.
    Adapted from Thomas
    """
    # The histogram of the data
    n, bins, patches = plt.hist(rlist, 50, density=True, 
        facecolor='b', alpha=0.75)

    # Labels
    plt.xlabel('RMSD (Angstrom)')
    plt.ylabel('Probability')
    plt.title('Histogram of RMSD for random pairs of %s' % aa)
    # x and y limits of plot
    plt.xlim(0, 1.5)
    plt.ylim(0, 35)
    plt.grid(True)
    # Save
    plt.savefig("RMSD_dist"+aa+".png")
    # Clear canvas for next plot
    plt.clf()


# In[2]:


if __name__=="__main__":
    
    #create list of AAs
    AA_list = ["CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", 
               "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    #create nested dictionary
    Prot_dict = extract_AA(AA_list)
    
    #create dictionary of RMSD
    AA_dict = generate_data(AA_list, Prot_dict)
    
    #save graphs
    for aa_name in AA_list:
        save_histogram(AA_dict[aa_name], aa_name)
        m = np.mean(AA_dict[aa_name])
        s = np.std(AA_dict[aa_name])
        print(aa_name+":", "Mean =", m, "SD =", s)


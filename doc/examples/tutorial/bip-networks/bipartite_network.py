import networkx as nx
import numpy as np
import random
from mpmath import *

#--------------------------------------------------------#
# Module for Generating Bipartite Networks
#--------------------------------------------------------#
def sf_sequence(size,exponent,max_deg):
  '''
  Generate powerlaw sequence 
  Input: size, exponent, max degree
  '''
  seq=[random.paretovariate(exponent-1) for i in range(size)]
  #round to desired range
  dseq = [min(max_deg, max(int(round(s)),0)) for s in seq]
  return dseq


def sf_cut_sequence(size,exponent,cutoff,max_deg):
  '''
  Generate a powerlaw sequence with exponential cutoff 
  p_k~ k^(-exponent)exp(-k/cutoff)
  Input: size, eponent, cutoff, max degree
  '''  
  seq = []
  # add to sequence using a rejection sampling method
  while len(seq) < size:
    r = random.uniform(0.0,1.0)
    x = 1-cutoff*ln(1-r)
    u = random.uniform(0.0,1.0)
    if u <= power(x,-exponent):
       seq.append(float(x))
  # round to desired range  
  dseq = [min(max_deg, max(int(round(s)),0)) for s in seq]
  return dseq  



def poi_sequence(size, avg_deg, max_deg):
    '''
    Returns discrete poisson degree sequence in range [0, max_deg]
    Input: size, average degree of distribution, maximum degree
    '''
    seq = np.random.poisson(avg_deg,size)
    # round to desired range 
    dseq = [min(max_deg, max(int(round(s)),0)) for s in seq]

    return dseq            



def bipartite_sequence_generator(n,m,maxtries,flag, **kwargs):  
    '''
    Build sequences for bipartite random graph from two degree sequences
    in at most maxtries 
    INPUT: m - size of large pop
           n - size of small pop
           flag - dist of large pop ('poi' or 'sfc' or 'sf')
           **kwargs - parameters for sequence of large pop
                    - avg_deg, (for poi)
                    - exp, cut (for sfc), exp (for sf)
    '''
    tries=0
    while tries < maxtries:
        # First build C_seq
        if flag == 'poi':
            C_seq = poi_sequence(m,kwargs['avg_deg'],n)
        elif flag == 'sfc':
            C_seq = sf_cut_sequence(m,kwargs['exp'], kwargs['cut'],n)
        elif flag == 'sf':
	        C_seq = sf_sequence(m,kwargs['exp'],n)
				# Now build S_seq
        mean_deg_S = np.divide(np.sum(C_seq), n)
        S_seq = poi_sequence(n,mean_deg_S,m)
        # Check that num stubs match
        sumS = np.sum(S_seq)
        sumC = np.sum(C_seq)
       # pdb.set_trace()
        while not sumS == sumC:    
            if sumS > sumC:
                p = np.random.randint(0,n)
                # reduce num stubs by 1
                S_seq[p] = S_seq[p]-1
            else:
                # reduce num stubs by 1
                p = np.random.randint(0,m)
                #p=C_seq.index(max(C_seq))
                C_seq[p] = C_seq[p]-1        
            sumS = np.sum(S_seq)
            sumC = np.sum(C_seq)
        if flag == 'poi':
          if np.abs(np.divide(sumC,m)-kwargs['avg_deg']) < 1:           
            return C_seq,S_seq
        else:
            return C_seq,S_seq
        tries+=1
    raise ValueError('Max tries exceeded')



def bipartite_network_generator(n,m,name,max_tries,flag,**kwargs):
    """
    Generate bipartite network with two degree sequences
    ------------------------------------------------------
    INPUT: 
        n - size of smaller population
        m - size of larger population
        name - desired name of output file
        max_tries - max num of attempts
        flag - deg dist of larger pop 
       
        **kwargs - parameters for deg dist
    OUTPUT:
        name.gpickle.bz2 - bipartite network


    """
    # generate sequences
    [S_seq, C_seq] = bipartite_sequence_generator(n,m,max_tries, flag, **kwargs)  

    #----------------------- Build Graph --------------------------------------

  
    # Match stubs randomly with configuration model
    G = nx.bipartite.configuration_model(S_seq, C_seq, create_using = nx.Graph())

    # Ensure no self-loops 
    G.remove_edges_from(G.selfloop_edges())

    # Save 
    nx.write_gpickle(G, name + '.gpickle.bz2')
 


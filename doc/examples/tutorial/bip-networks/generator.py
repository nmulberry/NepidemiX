from bipartite_network import bipartite_network_generator
import sys
import numpy as np
'''
Generating bipartite networks with specified degree distributions
'''

#-------------SFC & POI-------------------------------#
#n = int(sys.argv[1]) #size of poi seq
#m = int(sys.argv[2]) #size of other pop

#a = float(sys.argv[3]) #desired exponent
#k = float(sys.argv[4]) #desired cutoff


#for net in range(40):
#    bipartite_network_generator(n,m,str(net),10,'sfc',exp=a, cut=k)


#------------POI & POI------------------------------#

n = int(sys.argv[1]) #size of smaller pop 
m = int(sys.argv[2]) #size of larger pop 
a = float(sys.argv[3]) #desired average degree
for net in range(40):

	bipartite_network_generator(n,m,str(net),10,'poi',avg_deg=a)

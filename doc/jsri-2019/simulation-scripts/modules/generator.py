from bipartite_network import bipartite_network_generator
import sys
import numpy as np
#-------------SFC & POI-------------------------------#
n = int(sys.argv[1]) #size of poi seq
m = int(sys.argv[2]) #size of other pop

a = float(sys.argv[3]) #desired exponent
k = float(sys.argv[4]) #desired cutoff


for net in range(1):
    bipartite_network_generator(n,m,'smalltest',10,'sfc',exp=a, cut=k, min_deg =0)



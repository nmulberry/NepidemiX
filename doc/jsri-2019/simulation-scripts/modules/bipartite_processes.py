from nepidemix.process import BipartiteStateProcess

from nepidemix.utilities.networkxtra import attributeCount, neighbors_data_iter
import datetime
import numpy

import pdb
from bipartite_network import bipartite_network_generator

class BipartiteSIRProcess(BipartiteStateProcess):
    """
    S I R process on Bipartite network

    Attributes
    ----------
    beta_ab - Infection rate from nodes in set A to nodes in set B.
    beta_ba - Infection rate from nodes in set B to nodes in set A.
    gamma - Recovery rate.
    """
    def __init__(self, beta_ab, beta_ba, gamma):

        super(BipartiteSIRProcess, self).__init__(['SA', 'IA', 'RA','SB', 'IB', 'RB'],
                                         [],
                                         runNodeUpdate = True,
                                         runEdgeUpdate = False,
                                         runNetworkUpdate = False,
                                         constantTopology = True)
        self.beta_ab = float(beta_ab)
        self.beta_ba = float(beta_ba)
        self.gamma = float(gamma)	


    def nodeUpdateRule(self, node, srcNetwork, dt):
        # Read original node state.
        srcState = node[1][self.STATE_ATTR_NAME]
        # By default we have not changed states.
        dstState = srcState
        # Start out with a dictionary of zero neighbors in each state.
        nNSt = dict(zip(self.nodeStateIds,[0]*len(self.nodeStateIds)))
        # Calculate the actual numbers and update dictionary.
        nNSt.update(attributeCount(neighbors_data_iter(srcNetwork, node[0]),
                                   self.STATE_ATTR_NAME))
        # Pick a random number.
        eventp = numpy.random.random_sample()
        # Go through each state name, and chose an action.
        if srcState == 'SA':
            # Allow for asymmetric transmission.        
            if eventp < self.beta_ba*nNSt['IB']*dt:
               dstState = 'IA'
        elif srcState == 'SB':
            if eventp < self.beta_ab*nNSt['IA']*dt:
               dstState = 'IB'
        elif srcState == 'IA':
            if eventp < self.gamma*dt:
                dstState = 'RA'
        elif srcState == 'IB':
            if eventp < self.gamma*dt:
                dstState = 'RB'

        node[1][self.STATE_ATTR_NAME] = dstState

        return node


class BipartiteSALTSProcess(BipartiteStateProcess):
    """
    S A L T S process on Bipartite network with Vaccination State

    Attributes
    ----------
    beta    - Infectivity of latent nodes
    eta     - Transition rate from acute to latent 
    gamma_b - Treatment rate for nodes in set B.
    gamma_a - Treatment rate for nodes in set A.
    mu_a    - Migration rate (return to susceptible state).
    mu_b    - Migration rate (return to susceptible state).
    delta   - Natural death rate
    tau     - Death due to AIDS (return to susceptible state w/o entering treated state).
    alpha   - Backround HIV rate
    c       - Parameter which controls how much more infectious nodes in acute stage are.
    nu_a    - Vaccination rate
    nu_b    - Vaccination rate
    """
    def __init__(self, beta, eta, gamma_a, gamma_b, mu_a, mu_b, tau, c, delta, alpha, nu_a, nu_b):

        super(BipartiteSALTSProcess, self).__init__(['VA','VB','SA', 'AA', 'LA', 'TA', 'SB', 'AB', 'LB', 'TB'],
                                         [],
                                         runNodeUpdate = True,
                                         runEdgeUpdate = False,
                                         runNetworkUpdate = False,
                                         constantTopology = True)
        self.beta       = float(beta)
        self.delta      = float(delta)
        self.eta        = float(eta)
        self.gamma_a    = float(gamma_a)
        self.gamma_b    = float(gamma_b)
        self.mu_a       = float(mu_a)
        self.mu_b       = float(mu_b)	
        self.tau        = float(tau)
        self.c          = float(c)
        self.alpha      = float(alpha)
        self.nu_a       = float(nu_a)
        self.nu_b       = float(nu_b)
    def nodeUpdateRule(self, node, srcNetwork, dt):
        # Read original node state.
        srcState = node[1][self.STATE_ATTR_NAME]
        # By default we have not changed states.
        dstState = srcState
        # Start out with a dictionary of zero neighbors in each state.
        nNSt = dict(zip(self.nodeStateIds,[0]*len(self.nodeStateIds)))
        # Calculate the actual numbers and update dictionary.
        nNSt.update(attributeCount(neighbors_data_iter(srcNetwork, node[0]),
                                   self.STATE_ATTR_NAME))
        # Pick a random number.
        eventp = numpy.random.random_sample()
        # Chose an action:
        #-------------------------------------------------#
            

        # Move from Susceptible to Acute
        if srcState == 'SA': 
            eventq = numpy.random.random_sample() 
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.nu_a*dt:
                dstState = 'VA'
            elif eventq < self.beta*(self.c*nNSt['AB']+nNSt['LB'])*dt:
               dstState = 'AA'
        elif srcState == 'SB':
            eventq = numpy.random.random_sample() 
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.nu_b*dt:
                dstState = 'VB'
            elif eventq < self.beta*(self.c*nNSt['AA']+nNSt['LA'])*dt:
               dstState = 'AB'
        # Move from Acute to Latent (or migrate)
        elif srcState == 'AA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.eta*dt:
                dstState = 'LA'
        elif srcState == 'AB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.eta*dt:
                dstState = 'LB'
        # Move from Latent to either Susceptible or Treated
        elif srcState == 'LA':
            eventq = numpy.random.random_sample()
            if eventp < (self.tau+self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.gamma_a*dt:
                dstState = 'TA'
        elif srcState == 'LB':
            eventq = numpy.random.random_sample()
            if eventp < (self.tau+self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.gamma_b*dt:
                dstState = 'TB'
        # Move from Treated to Susceptible
        elif srcState == 'TA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
        elif srcState == 'TB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
        # Vaccinated state
        elif srcState == 'VA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
        elif srcState == 'VB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'

        # update node with chosen state
        node[1][self.STATE_ATTR_NAME] = dstState

        return node



class DynamicBipartiteSALTSProcess(BipartiteStateProcess):
    """
    S A L T S process on Dynamic Bipartite network

    Attributes
    ----------
    beta    - Infectivity of latent nodes
    eta     - Transition rate from acute to latent 
    gamma_b - Treatment rate for nodes in set B.
    gamma_a - Treatment rate for nodes in set A.
    mu_a    - Migration rate (return to susceptible state).
    mu_b    - Migration rate (return to susceptible state).
    delta   - Natural death rate
    tau     - Death due to AIDS (return to susceptible state w/o entering treated state).
    alpha   - Backround HIV rate
    c       - Parameter which controls how much more infectious nodes in acute stage are.
    k       - cut-off for sfc distribution
    a       - exonent for sfc distribution
    """
   def __init__(self, beta, eta, gamma_a, gamma_b, mu_a, mu_b, tau, c, delta, alpha, nu_a, nu_b, a,k):

        super(DyanamicBipartiteSALTSProcess, self).__init__(['VA','VB','SA', 'AA', 'LA', 'TA', 'SB', 'AB', 'LB', 'TB'],
                                         [],
                                         runNodeUpdate = True,
                                         runEdgeUpdate = False,
                                         runNetworkUpdate = True)
        self.beta       = float(beta)
        self.delta      = float(delta)
        self.eta        = float(eta)
        self.gamma_a    = float(gamma_a)
        self.gamma_b    = float(gamma_b)
        self.mu_a       = float(mu_a)

        self.mu_b       = float(mu_b)
        self.tau        = float(tau)
        self.c          = float(c)
        self.alpha      = float(alpha)
        self.nu_a       = float(nu_a)
        self.nu_b       = float(nu_b)

        self.k          = float(k)
        self.a          = float(a)







    def networkUpdateRule(self, network, dt):
        # remove all current edges
        network.remove_edges_from(network.edges())        
        # rewire
        nodes = network.nodes(data=True)
        sizeA = sum(1 for n in nodes if n[-1]['bipartite']==0)
        sizeB = sum(1 for n in nodes if n[-1]['bipartite']==1)
        network=bipartite_network_generator(sizeB, sizeA,'sfc', network, exp=self.a, cut=self.k, min_deg = 0) 
        return network

    def nodeUpdateRule(self, node, srcNetwork, dt):
        # Read original node state.
        srcState = node[1][self.STATE_ATTR_NAME]
        # By default we have not changed states.
        dstState = srcState
        # Start out with a dictionary of zero neighbors in each state.
        nNSt = dict(zip(self.nodeStateIds,[0]*len(self.nodeStateIds)))
        # Calculate the actual numbers and update dictionary.
        nNSt.update(attributeCount(neighbors_data_iter(srcNetwork, node[0]),
                                   self.STATE_ATTR_NAME))
        # Pick a random number.
        eventp = numpy.random.random_sample()
        # Chose an action:
        #-------------------------------------------------#
            

        # Move from Susceptible to Acute
        if srcState == 'SA': 
            eventq = numpy.random.random_sample() 
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.nu_a*dt:
                dstState = 'VA'
            elif eventq < self.beta*(self.c*nNSt['AB']+nNSt['LB'])*dt:
               dstState = 'AA'
        elif srcState == 'SB':
            eventq = numpy.random.random_sample() 
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.nu_b*dt:
                dstState = 'VB'
            elif eventq < self.beta*(self.c*nNSt['AA']+nNSt['LA'])*dt:
               dstState = 'AB'
        # Move from Acute to Latent (or migrate)
        elif srcState == 'AA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.eta*dt:
                dstState = 'LA'
        elif srcState == 'AB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.eta*dt:
                dstState = 'LB'
        # Move from Latent to either Susceptible or Treated
        elif srcState == 'LA':
            eventq = numpy.random.random_sample()
            if eventp < (self.tau+self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
            elif eventq < self.gamma_a*dt:
                dstState = 'TA'
        elif srcState == 'LB':
            eventq = numpy.random.random_sample()
            if eventp < (self.tau+self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
            elif eventq < self.gamma_b*dt:
                dstState = 'TB'
        # Move from Treated to Susceptible
        elif srcState == 'TA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
        elif srcState == 'TB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'
        # Vaccinated state
        elif srcState == 'VA':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_a+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LA'
                else:
                    dstState = 'SA'
        elif srcState == 'VB':
            eventq = numpy.random.random_sample()
            if eventp < (self.mu_b+self.delta)*dt:
                if eventq < self.alpha*dt:
                    dstState = 'LB'
                else:
                    dstState = 'SB'

        # update node with chosen state
        node[1][self.STATE_ATTR_NAME] = dstState

        return node


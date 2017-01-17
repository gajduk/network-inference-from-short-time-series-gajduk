import random
import numpy as np
from adapters import one_pair_at_a_time_wrapper
from correlation_based import partial_corr,cross_correlation,iota,kendall,granger,cross_correlation_correct
from symbolic import symbolSequenceSimilarity, mutualInformationOfSymbols
from srep05030_network_reconstruction import zoranDerivative


@one_pair_at_a_time_wrapper
def random_(x1, x2):
	return random.random()

def holy_grail(i):
	return np.abs(np.array(i.y.reshape((i.n_nodes, i.n_nodes))))


#depracated methods "granger_partial_r": granger_partial_r, ,

methods = {"holy_grail": holy_grail, "random": random_, "cross_correlation": cross_correlation, "iota": iota,
           "kendall": kendall,  "partial_corr": partial_corr, "symbol_similarity": symbolSequenceSimilarity,
           "symbol_mutual": mutualInformationOfSymbols, "derivative_correlation": zoranDerivative, "granger": granger, "cross_correlation_correct":cross_correlation_correct}

good_methods = {f: methods[f] for f in
                ["derivative_correlation", "cross_correlation_correct", "holy_grail", "kendall", "random", "partial_corr",
                 "symbol_similarity", "symbol_mutual","granger", "iota"]}  # ,  ,"granger_partial_r"]}

if __name__ == "__main__":

	#generate some random data
	from src.model.instance import Instance
	n_nodes = 10
	time_series = [np.random.randn(n_nodes,20)]
	A = np.random.randint(0,2,(n_nodes,n_nodes))
	i = Instance(time_series,A,n_nodes)

	#use a method to reconstruct the network
	print methods['cross_correlation'](i)
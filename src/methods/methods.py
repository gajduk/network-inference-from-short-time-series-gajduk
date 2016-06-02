import random
import numpy as np
from adapters import one_pair_at_a_time_wrapper
from correlation_based import partial_corr,cross_correlation,iota,kendall
from symbolic import symbolSequenceSimilarity, mutualInformationOfSymbols
from srep05030_network_reconstruction import zoranDerivative


@one_pair_at_a_time_wrapper
def random_(x1, x2):
	return random.random()

def holy_grail(i):
	return np.abs(np.array(i.y.reshape((i.n_nodes, i.n_nodes))))


#depracated methods "granger_partial_r": granger_partial_r, "granger": granger,

methods = {"holy_grail": holy_grail, "random": random_, "cross_correlation": cross_correlation, "iota": iota,
           "kendall": kendall,  "partial_corr": partial_corr, "symbol_similarity": symbolSequenceSimilarity,
           "symbol_mutual": mutualInformationOfSymbols, "derivative_correlation": zoranDerivative}

good_methods = {f: methods[f] for f in
                ["derivative_correlation", "cross_correlation", "holy_grail", "kendall", "random", "partial_corr",
                 "symbol_similarity", "symbol_mutual"]}  # ,  "granger","granger_partial_r"]}

import numpy as np
import math

from src.utils import load_r_file, deprecated
from adapters import one_pair_at_a_time_wrapper


@deprecated
@one_pair_at_a_time_wrapper
def granger_r(x1, x2):
	d = np.vstack((np.matrix(x1), np.matrix(x2))).T
	granger_a = load_r_file('gc1a.R', "granger_a")
	return granger_a.granger(d, 1)[0]


@deprecated
def granger_partial_r(i):
	'''
	Depracated
	'''
	n = i.n_nodes
	cc = granger_r(i)

	ccc_mean = np.zeros((n, n))
	ccc_max = np.zeros((n, n))
	ccc_min = np.zeros((n, n))
	ccc = np.zeros((n, n, n))
	granger_a = load_r_file('gc1a.R', "granger_a")
	for ii in range(n):
		xii, _ = i.get(ii)
		for jj in range(n):
			xjj, _ = i.get(jj)
			for kk in range(n):
				xkk, _ = i.get(kk)
				d = np.vstack((xii.T, xjj.T, xkk.T)).T
				temp = granger_a.granger_part(d, 1)[0]
				if math.isnan(temp) or math.isinf(temp):
					temp = 0.0
				ccc[ii][jj][kk] = temp
			ccc_mean[ii][jj] = np.mean(ccc[ii][jj][:])
			ccc_max[ii][jj] = np.max(ccc[ii][jj][:])
			ccc_min[ii][jj] = np.min(ccc[ii][jj][:])

	return ccc_max

import random
import math

import numpy as np
from scipy import stats, linalg
import statsmodels.tsa.stattools as sts

from src.utils import load_r_file

def one_pair_at_a_time_wrapper(method):
	def inner(i):
		res = np.zeros((i.n_nodes,i.n_nodes))
		for idx_node1 in range(i.n_nodes):
			x1,_ = i.get(idx_node1)
			for idx_node2 in range(i.n_nodes):
				x2,_ = i.get(idx_node2)
				res[idx_node1][idx_node2] = method(x1,x2)
		return res
	return inner


def time_series_columnwise_wrapper(method):
	def inner(i):
		data = np.zeros((i.n_time_points,i.n_nodes))
		for idx_node in range(i.n_nodes):
			x,_ = i.get(idx_node)
			data[:,idx_node] = x
		return method(data)
	return inner

@one_pair_at_a_time_wrapper
def random_(x1,x2):
    return random.random()

@one_pair_at_a_time_wrapper
def granger(x1,x2):
	x1[0] += .00000000001
	x2[0] += .00000000001
	res = sts.grangercausalitytests(np.vstack((x1,x2)).T,2,verbose=False)[1][0]['params_ftest'][0]
	return res

@one_pair_at_a_time_wrapper
def cross_correlation(x1,x2):
    mean_x1,mean_x2 = np.mean(x1),np.mean(x2)
    std_x1,std_x2 = np.std(x1),np.std(x2)
    return np.dot((x1-mean_x1),(x2-mean_x2))/std_x1/std_x2

@one_pair_at_a_time_wrapper
def iota(x1,x2):
    pi1 = np.argsort(x1)
    g_k_l = x2[pi1]
    n = len(x1)
    res = 0.0
    for i in range(n-2):
        for j in range(i+1,n-1):
            res += 1 if (g_k_l[j+1]-g_k_l[i])*(g_k_l[i]-g_k_l[j]) > 0 else 0
    res /= (n-1)*(n-2)/2
    return 1.0-res

@one_pair_at_a_time_wrapper
def kendall(x,y):
    numer = 1
    d = 1
    for i in range(1,len(x)):
        for j in range(i):
            numer += np.sign(x[i] - x[j]) * np.sign(y[i] - y[j])
            if x[i] != x[j] and y[i] != y[j]:
                d += 1
    return numer*1.0/d

@one_pair_at_a_time_wrapper
def granger_r(x1,x2):
	d = np.vstack((x1,x2)).T
	granger_a = load_r_file('gc1a.R',"granger_a")
	return granger_a.granger(d,1)[0]

def granger_partial_r(i):
	n = i.n_nodes
	cc = granger_r(i)

	ccc_mean = np.zeros((n,n))
	ccc_max = np.zeros((n,n))
	ccc_min = np.zeros((n,n))
	ccc = np.zeros((n,n,n))
	granger_a = load_r_file('gc1a.R',"granger_a")
	for ii in range(n):
		xii,_ = i.get(ii)
		for jj in range(n):
			xjj,_ = i.get(jj)
			for kk in range(n):
				xkk,_ = i.get(kk)
				d = np.vstack((xii,xjj,xkk)).T
				temp = granger_a.granger_part(d,1)[0]
				if math.isnan(temp) or math.isinf(temp):
					temp = 0.0
				ccc[ii][jj][kk] = temp
			ccc_mean[ii][jj] = np.mean(ccc[ii][jj][:])
			ccc_max[ii][jj] = np.max(ccc[ii][jj][:])
			ccc_min[ii][jj] = np.min(ccc[ii][jj][:])

	return ccc_max

def holy_grail(i):
	return i.y.reshape((i.n_nodes,i.n_nodes))

"""
Partial Correlation in Python (clone of Matlab's partialcorr)
This uses the linear regression approach to compute the partial
correlation (might be slow for a huge number of variables). The
algorithm is detailed here:
    http://en.wikipedia.org/wiki/Partial_correlation#Using_linear_regression
Taking X and Y two variables of interest and Z the matrix with all the variable minus {X, Y},
the algorithm can be summarized as
    1) perform a normal linear least-squares regression with X as the target and Z as the predictor
    2) calculate the residuals in Step #1
    3) perform a normal linear least-squares regression with Y as the target and Z as the predictor
    4) calculate the residuals in Step #3
    5) calculate the correlation coefficient between the residuals from Steps #2 and #4;
    The result is the partial correlation between X and Y while controlling for the effect of Z
Date: Nov 2014
Author: Fabian Pedregosa-Izquierdo, f@bianp.net
Testing: Valentina Borghesani, valentinaborghesani@gmail.com
"""
@time_series_columnwise_wrapper
def partial_corr(C):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in C, controlling
    for the remaining variables in C.
    Parameters
    ----------
    C : array-like, shape (n, p)
        Array with the different variables. Each column of C is taken as a variable
    Returns
    -------
    P : array-like, shape (p, p)
        P[i, j] contains the partial correlation of C[:, i] and C[:, j] controlling
        for the remaining variables in C.
    """

    C = np.asarray(C)
    p = C.shape[1]
    P_corr = np.zeros((p, p), dtype=np.float)
    for i in range(p):
        P_corr[i, i] = 1
        for j in range(i+1, p):
            idx = np.ones(p, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j = linalg.lstsq(C[:, idx], C[:, i])[0]

            res_j = C[:, j] - C[:, idx].dot( beta_i)
            res_i = C[:, i] - C[:, idx].dot(beta_j)

            corr = stats.pearsonr(res_i, res_j)[0]
            P_corr[i, j] = corr
            P_corr[j, i] = corr

    return P_corr


methods = {"holy_grail":holy_grail,"granger_partial_r":granger_partial_r,"random":random_,"cross_correlation":cross_correlation,"iota":iota,"kendall":kendall,"granger":granger,"partial_corr":partial_corr,"granger_r":granger_r}
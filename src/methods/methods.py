import math
import random

import numpy as np
import statsmodels.tsa.stattools as sts
from scipy import stats, linalg

from src.utils import load_r_file, deprecated
from utils import one_pair_at_a_time_wrapper, time_series_columnwise_wrapper


@one_pair_at_a_time_wrapper
def random_(x1, x2):
	return random.random()


@one_pair_at_a_time_wrapper
def granger(x1, x2):
	x1[0] += .00000000001
	x2[0] += .00000000001
	res = sts.grangercausalitytests(np.vstack((x1, x2)).T, 2, verbose=False)[1][0]['params_ftest'][0]
	return res


@one_pair_at_a_time_wrapper
def cross_correlation(x1, x2):
	mean_x1, mean_x2 = np.mean(x1), np.mean(x2)
	std_x1, std_x2 = np.std(x1), np.std(x2)
	return np.dot((x1 - mean_x1).T, (x2 - mean_x2)) / std_x1 / std_x2


@one_pair_at_a_time_wrapper
def iota(x1, x2):
	pi1 = np.argsort(x1)
	g_k_l = x2[pi1]
	n = len(x1)
	res = 0.0
	for i in range(n - 2):
		for j in range(i + 1, n - 1):
			res += 1 if (g_k_l[j + 1] - g_k_l[i]) * (g_k_l[i] - g_k_l[j]) > 0 else 0
	res /= (n - 1) * (n - 2) / 2
	return 1.0 - res


@one_pair_at_a_time_wrapper
def kendall(x, y):
	numer = 1
	d = 1
	for i in range(1, len(x)):
		for j in range(i):
			numer += np.sign(x[i] - x[j]) * np.sign(y[i] - y[j]).T
			if x[i] != x[j] and y[i] != y[j]:
				d += 1
	return numer * 1.0 / d


@one_pair_at_a_time_wrapper
def granger_r(x1, x2):
	d = np.vstack((x1.T, x2.T)).T
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


def holy_grail(i):
	return np.absolute(np.asarray(i.y.reshape((i.n_nodes, i.n_nodes))))


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
		for j in range(i + 1, p):
			idx = np.ones(p, dtype=np.bool)
			idx[i] = False
			idx[j] = False
			beta_i = linalg.lstsq(C[:, idx], C[:, j])[0]
			beta_j = linalg.lstsq(C[:, idx], C[:, i])[0]

			res_j = C[:, j] - C[:, idx].dot(beta_i)
			res_i = C[:, i] - C[:, idx].dot(beta_j)

			corr = stats.pearsonr(res_i, res_j)[0]
			P_corr[i, j] = corr
			P_corr[j, i] = corr

	return P_corr


from symbolic import symbolSequenceSimilarity, mutualInformationOfSymbols


def g_generator(n_harmonics=10):
	pi = np.pi
	hw = 10.0
	def coef_generator():
		r = random.random()
		return math.pow(10, r * 2.00432137) - 1.0

	a = np.matrix([coef_generator() for _ in range(n_harmonics)])
	b = np.matrix([coef_generator() for _ in range(n_harmonics)])
	k = np.matrix([i + 1 for i in range(n_harmonics)])
	def g(x):
		x_norm = (x+hw)/(2.0*hw)
		temp = np.sum(np.multiply(a,np.sin(k*pi*x_norm)) + np.multiply(b,np.cos(k*pi*x_norm)))
		temp = temp/np.sum(a+b)*10.0
		return temp

	return g


def example1(instance):
	x = instance.x
	N = instance.n_nodes
	L = instance.n_time_points
	A0 = np.array(instance.y.reshape((N,N)))
	xx,dx,f_F,h_F = np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1))
	tau = 0.2
	R = 1
	import time

	def f(x):
		return - x

	def h(x):
		return math.tanh(x)

	for i in range(N):
		for k in range(L-1):
			xx[i][k] = ( x[i][k] + x[i][k+1] ) / 2.0
			dx[i][k] = ( x[i][k+1] - x[i][k] )
			f_F[i][k] = f(xx[i][k])
			h_F[i][k] = h(xx[i][k])


	min_delta = 100000
	reconstructed_A = -1
	triple_for,inversion,normalization_ = 0.0,0.0,0.0

	for iter_ in range(100000):
		if iter_ % 500 == 0:
			print iter_
			print triple_for,inversion,normalization_
			print min_delta
		t = time.time()
		g = g_generator()
		g_F = np.zeros((N,L-1))

		for i in range(N):
			for k in range(L-1):
				g_F[i][k] = g(xx[i][k])

		B, C, E = np.zeros((N,N)),np.zeros((N,N)),np.zeros((N,N))
		for i in range(N):
			for j in range(N):
				B[i][j] = np.sum(np.multiply(g_F[i][:],dx[j][:])) / ( tau * (L-1) * R)
				C[i][j] = np.sum(np.multiply(g_F[i][:],f_F[j][:]))  / ( 1.0 * (L-1) * R)
				E[i][j] = np.sum(np.multiply(g_F[i][:],h_F[j][:])) / ( 1.0 * (L-1) * R)
		elapsed = time.time() - t
		triple_for += elapsed
		t = time.time()
		A = np.linalg.inv(E)*(B - C)
		elapsed = time.time() - t
		inversion += elapsed
		#A = ((A-np.min(A[:]))/(np.max(A[:])-np.min(A[:]))-0.5)*20

		t = time.time()
		delta_A = 0.0
		normalization = 0
		for i in range(N):
			for j in range(N):
				delta_A += (A[i][j]-A0[i][j])**2
				normalization += A0[i][j]**2
		delta_A = math.sqrt(delta_A/normalization)
		if delta_A < min_delta:
			min_delta = delta_A
			reconstructed_A = A
		elapsed = time.time() - t
		inversion += elapsed
	import matplotlib.pylab as plt
	plt.subplot(1,3,1)
	plt.set_cmap('bwr')
	plt.imshow(A0,interpolation='none', vmin=-10, vmax=10)
	plt.subplot(1,3,2)
	plt.set_cmap('bwr')
	plt.imshow(reconstructed_A,interpolation='none', vmin=-10, vmax=10)
	plt.subplot(1,3,3)
	plt.set_cmap('bwr')
	plt.imshow(np.abs(reconstructed_A-A0),interpolation='none', vmin=-10, vmax=10)
	plt.show()

@deprecated
def g_generator_old(n_harmonics=10):
	def coef_generator():
		r = random.random()
		return math.pow(10, r * 2.00432137) - 1.0

	a = np.matrix([coef_generator() for _ in range(n_harmonics)])
	b = np.matrix([coef_generator() for _ in range(n_harmonics)])
	k = np.matrix([i + 1 for i in range(n_harmonics)])

	def g(x):
		return np.sum(np.multiply(a,np.sin(k * x)) + np.multiply(b,np.cos(k * x)))/(np.sum(a)+np.sum(b))

	return g

@deprecated
def example1_old(instance):
	x = instance.x
	n = instance.n_nodes
	L = instance.n_time_points

	def f(x):
		return - x

	def h(x):
		return math.tanh(x)
	min_delta = 100000
	reconstructed_R = -1
	A = np.array(instance.y.reshape((n,n)))
	for asd in range(10000):
		g = g_generator()
		B, C, E = np.zeros((n,n)),np.zeros((n,n)),np.zeros((n,n))

		for i in range(n):
			for j in range(n):
				B[i][j] = np.sum([g(x[i][k])*x[j][k] for k in range(L)])
				C[i][j] = np.sum([g(x[i][k])*f(x[j][k]) for k in range(L)])
				E[i][j] = np.sum([g(x[i][k])*h(x[j][k]) for k in range(L)])
		R = np.linalg.inv(E)*(B - C)
		R = ((R-np.min(R[:]))/(np.max(R[:])-np.min(R[:]))-0.5)*20

		delta_A = 0
		normalization = 0
		for i in range(n):
			for j in range(n):
				delta_A += (R[i][j]-A[i][j])*(R[i][j]-A[i][j])
				normalization += A[i][j]*A[i][j]
		delta_A = math.sqrt(delta_A)/36
		if delta_A < min_delta:
			min_delta = delta_A
			reconstructed_R = R
	print min_delta
	import matplotlib.pylab as plt
	plt.subplot(1,3,1)
	plt.set_cmap('bwr')
	plt.imshow(A,interpolation='none')
	plt.subplot(1,3,2)
	plt.set_cmap('bwr')
	plt.imshow(reconstructed_R,interpolation='none')
	plt.subplot(1,3,3)
	plt.set_cmap('bwr')
	plt.imshow(reconstructed_R-A,interpolation='none')
	plt.show()



methods = {"holy_grail": holy_grail, "random": random_, "cross_correlation": cross_correlation, "iota": iota,
           "kendall": kendall, "granger": granger, "partial_corr": partial_corr, "granger_partial_r": granger_partial_r,
           "granger_r": granger_r, "symbol_similarity": symbolSequenceSimilarity,
           "symbol_mutual": mutualInformationOfSymbols, "example1": example1}

good_methods = {f: methods[f] for f in
                ["example1", "cross_correlation", "holy_grail", "kendall", "random", "partial_corr",
                 "symbol_similarity", "symbol_mutual"]}  # ,  "granger","granger_partial_r"]}

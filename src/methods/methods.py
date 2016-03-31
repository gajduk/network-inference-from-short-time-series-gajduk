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
	#a = [[ 11.16713159 , 24.50579577 , 19.88689278 ,  2.46239161 , 66.00032627,
    #2.91314177 , 15.71485056 , 24.20857391 ,  1.44829805 , 12.51065127]]
	b = np.matrix([coef_generator() for _ in range(n_harmonics)])
	#b = [[  2.04817081,   1.75099185,   5.97149345,   8.87731769,   5.66473496,
    #4.34743746,  50.58817116,  65.15472992,   0.78322186 , 80.86623373]]
	k = np.matrix([i + 1 for i in range(n_harmonics)])
	def g(x):
		x_norm = k*pi*(x+hw)/(2.0*hw)
		temp = np.sum(np.multiply(a,np.sin(x_norm)) + np.multiply(b,np.cos(x_norm)))
		#temp = temp/np.sum(a+b)*10.0
		return temp

	return g


def example1(instance):
	'''
	example1
	x = [[-0.12389172666, 0.587579759577, -0.311814570479, -0.414858833614, -0.0738960295328, -0.660661139367],
		[-1.41899704578, -1.40545533491, -0.483119910044, -1.11585183184, 1.34934504315, -0.432930242154],
		[-2.5422552576, -4.07820651294, -1.73518916233, -0.800369827126, 2.55187805734, -0.169093037506],
		[-2.72917539045, -5.85888569981, -3.28653016681, 0.0353143965753, 2.71427307399, 0.0529312258093],
		[-1.52180685357, -6.18742188307, -5.68056949745, 1.08736075478, 1.33330005824, 0.230367898216],
		[0.693843044648, -5.03975857841, -7.17032428001, 2.15796803712, -0.592372556191, 0.253360537644],
		[3.46035824222, -2.41488977518, -6.33502991938, 2.80778639595, -2.20915830073, 0.0287826649899],
		[5.80486659771, -0.254529489474, -5.23240510083, 2.95911569661, -3.41619131684, -0.169334361816],
		[7.54743437683, 1.43855144176, -2.71397554699, 2.7727483063, -4.29897526371, -0.331608570583],
		[8.76727149074, 2.94383846414, -0.187925068392, 2.19438353812, -4.93093954588, -0.464468145224],
		[8.85144239896, 4.92394875873, 1.34916061019, 0.54549904057, -5.14283530688, -0.573244284032],
		[6.82311418653, 5.44959059523, 3.94912395845, -1.30853866102, -3.33370698306, -0.662302480009],
		[4.08009264104, 5.10103678313, 6.94836395588, -2.94520180654, -0.764424014411, -0.735195680896],
		[0.859764343595, 4.68625619598, 9.36734272598, -4.35512051998, 1.42347495019, -0.784962449961],
		[-2.29617659039, 2.51411740383, 8.94345480114, -5.42671086302, 3.18936688594, -0.553526228177],
		[-4.86250934292, 0.132308235983, 7.34727098289, -5.90861290221, 4.50129360858, -0.260874090987],
		[-6.77153554682, -1.71396004278, 4.31280270062, -5.86794150524, 5.42792789475, -0.0206202919944],
		[-8.26737528831, -3.14918445592, 1.36120313212, -5.43188407684, 6.05668342342, 0.176087720454],
		[-9.05762671138, -4.86937298424, -0.579915104924, -4.12902634344, 6.46752923251, 0.337138594162],
		[-9.17571334396, -6.95405138423, -1.58461422162, -2.06485181171, 6.7123603964, 0.468995771491]]

example2
	x = [[0.910958300485961,0.270013480278752,1.98294524959933,1.52693713168222,1.44662584747648,1.52400202048153],
[0.887424100918529,0.473056213012262,1.99704196175860,1.30898088285714,1.20720132407867,1.29521920091175],
[0.878949943699696,0.639125477990157,2.02211032526094,1.12202728118220,1.01148130425606,1.10825594665204],
[0.877520963758689,0.775163603251518,2.04953543664256,0.965027240449318,0.851342789681346,0.958056348235911],
[0.879592498693157,0.886732590558362,2.07605124126255,0.834544155138376,0.720216689781918,0.839205078146950],
[0.883400991266060,0.978302844368203,2.10038766554664,0.726816403227460,0.612796119644950,0.746765601909131],
[0.887962317677017,1.05350981656100,2.12206986958830,0.638355847823220,0.524754748485061,0.676391759694401],
[0.892732317308835,1.11530631686476,2.14103046445655,0.566080080019455,0.452574749999146,0.624314495688058],
[0.897404248554033,1.16610770447526,2.15739176381989,0.507321150558368,0.393378111248520,0.587247848417564],
[0.901817174449033,1.20788409491391,2.17137272659791,0.459802258713669,0.344818569519768,0.562347085014882],
[0.905892868401172,1.24225182418278,2.18323061770350,0.421584902983621,0.304973022151196,0.547136430621417],
[0.909603610451264,1.27053236050669,2.19323161761776,0.391032526226688,0.272271399097124,0.539489234780230],
[0.912949147483651,1.29381175843818,2.20163242011013,0.366762329022382,0.245425407720068,0.537584745099830],
[0.915944820245268,1.31297919227412,2.20867019274071,0.347615769480502,0.223382341357178,0.539897947639169],
[0.918613368573787,1.32876621563439,2.21455777975901,0.332621898427703,0.205277850165796,0.545165879025547],
[0.920980869924315,1.34177229452283,2.21948142568765,0.320974982829887,0.190405255614050,0.552373119856687],
[0.923074121753960,1.35249097372435,2.22360160101523,0.312006901441193,0.178184060753811,0.560716868476115],
[0.924919389609139,1.36132692832106,2.22705412846626,0.305169809389877,0.168139464354082,0.569585851971540],
[0.926541680131453,1.36861350020310,2.22995310230459,0.300015063607905,0.159881297418866,0.578525286856564],
[0.927964400971387,1.37462413810646,2.23239315043847,0.296179654536419,0.153090311154041,0.587214165361722],
[0.929209214763485,1.37958406442249,2.23445255646597,0.293370490150977,0.147504212413062,0.595435637461200]]
	L,N = np.matrix(x).shape
	x = np.matrix(x).T.tolist()
	A0 = [[0,7.2,9.3,0,0,-1],[-1,0,7,0,0,0],[-3,4.2,-3.5,-4.4,0,0],[8.5,6.7,-7.7,0,-9,0],[-5,0,0,0,0,0],[0,2.5,0,10,-2.5,0]]
	'''
	x = instance.x
	x = [[0.910958300485961, 0.270013480278752, 1.98294524959933, 1.52693713168222, 1.44662584747648, 1.52400202048153],
	     [0.887424100918529, 0.473056213012262, 1.99704196175860, 1.30898088285714, 1.20720132407867, 1.29521920091175],
	     [0.878949943699696, 0.639125477990157, 2.02211032526094, 1.12202728118220, 1.01148130425606, 1.10825594665204],
	     [0.877520963758689, 0.775163603251518, 2.04953543664256, 0.965027240449318, 0.851342789681346,
	      0.958056348235911],
	     [0.879592498693157, 0.886732590558362, 2.07605124126255, 0.834544155138376, 0.720216689781918,
	      0.839205078146950],
	     [0.883400991266060, 0.978302844368203, 2.10038766554664, 0.726816403227460, 0.612796119644950,
	      0.746765601909131],
	     [0.887962317677017, 1.05350981656100, 2.12206986958830, 0.638355847823220, 0.524754748485061,
	      0.676391759694401],
	     [0.892732317308835, 1.11530631686476, 2.14103046445655, 0.566080080019455, 0.452574749999146,
	      0.624314495688058],
	     [0.897404248554033, 1.16610770447526, 2.15739176381989, 0.507321150558368, 0.393378111248520,
	      0.587247848417564],
	     [0.901817174449033, 1.20788409491391, 2.17137272659791, 0.459802258713669, 0.344818569519768,
	      0.562347085014882],
	     [0.905892868401172, 1.24225182418278, 2.18323061770350, 0.421584902983621, 0.304973022151196,
	      0.547136430621417],
	     [0.909603610451264, 1.27053236050669, 2.19323161761776, 0.391032526226688, 0.272271399097124,
	      0.539489234780230],
	     [0.912949147483651, 1.29381175843818, 2.20163242011013, 0.366762329022382, 0.245425407720068,
	      0.537584745099830],
	     [0.915944820245268, 1.31297919227412, 2.20867019274071, 0.347615769480502, 0.223382341357178,
	      0.539897947639169],
	     [0.918613368573787, 1.32876621563439, 2.21455777975901, 0.332621898427703, 0.205277850165796,
	      0.545165879025547],
	     [0.920980869924315, 1.34177229452283, 2.21948142568765, 0.320974982829887, 0.190405255614050,
	      0.552373119856687],
	     [0.923074121753960, 1.35249097372435, 2.22360160101523, 0.312006901441193, 0.178184060753811,
	      0.560716868476115],
	     [0.924919389609139, 1.36132692832106, 2.22705412846626, 0.305169809389877, 0.168139464354082,
	      0.569585851971540],
	     [0.926541680131453, 1.36861350020310, 2.22995310230459, 0.300015063607905, 0.159881297418866,
	      0.578525286856564],
	     [0.927964400971387, 1.37462413810646, 2.23239315043847, 0.296179654536419, 0.153090311154041,
	      0.587214165361722],
	     [0.929209214763485, 1.37958406442249, 2.23445255646597, 0.293370490150977, 0.147504212413062,
	      0.595435637461200]]

	L = instance.n_time_points
	N = instance.n_nodes

	L, N = np.matrix(x).shape
	x = np.matrix(x).T.tolist()

	A0 = np.array(np.squeeze(instance.y)[0]).reshape((N,N))
	xx,dx,f_F,h_F = np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1))
	tau = 0.25
	R = 1

	def f(x):
		return - x

	def h(x):
		return math.tanh(x)

	for i in range(N):
		for k in range(L-1):
			xx[i][k] = ( x[i][k] + x[i][k+1] ) / 2.0
			dx[i][k] = ( x[i][k+1] - x[i][k] )
			f_F[i][k] = f(xx[i][k])
			if i < 3:
				h_F[i][k] = 1.0*xx[i][k]**5/(1.0+xx[i][k]**5)
			else:
				h_F[i][k] = 1.0 / (1.0 + xx[i][k] ** 5)

	min_delta = 100000
	reconstructed_A = -1

	for iter_ in range(100):
		g = g_generator()
		g_F = np.zeros((N,L-1))

		for i in range(N):
			for k in range(L-1):
				g_F[i][k] = g(xx[i][k])

		B, C, E = np.zeros((N,N)),np.zeros((N,N)),np.zeros((N,N))
		for i in range(N):
			for j in range(N):
				B[i][j] = np.sum(np.multiply(g_F[:][i],dx[:][j])) / ( tau * (L-1) * R)
				C[i][j] = np.sum(np.multiply(g_F[:][i],f_F[:][j]))  / ( 1.0 * (L-1) * R)
				E[i][j] = np.sum(np.multiply(g_F[:][i],h_F[:][j])) / ( 1.0 * (L-1) * R)

		A = np.dot(np.linalg.inv(E),(B - C))

		delta_A = 0.0
		normalization = 0
		for i in range(N):
			for j in range(N):
				delta_A += (A[i][j]-A0[i][j])**2
				normalization += A0[i][j]**2
		delta_A = math.sqrt(delta_A/normalization)
		print delta_A
		if delta_A < min_delta:
			min_delta = delta_A
			reconstructed_A = A
	import matplotlib.pylab as plt
	plt.subplot(1,3,1)
	plt.set_cmap('bwr')
	plt.imshow(A0,interpolation='none', vmin=-1, vmax=1)
	plt.subplot(1,3,2)
	plt.set_cmap('bwr')
	plt.imshow(reconstructed_A,interpolation='none', vmin=-1, vmax=1)
	plt.subplot(1,3,3)
	plt.set_cmap('bwr')
	plt.imshow(np.abs(reconstructed_A-A0),interpolation='none', vmin=-1, vmax=1)
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

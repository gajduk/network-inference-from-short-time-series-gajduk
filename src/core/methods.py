import random
import numpy as np
import statsmodels.tsa.stattools as sts

def random_(x1,x2):
    return random.random()

def granger(x1,x2):
	x1[0] += .00000000001
	x2[0] += .00000000001
	res = sts.grangercausalitytests(np.vstack((x1,x2)).T,2,verbose=False)[1][0]['params_ftest'][0]
	return res

def cross_correlation(x1,x2):
    mean_x1,mean_x2 = np.mean(x1),np.mean(x2)
    std_x1,std_x2 = np.std(x1),np.std(x2)
    return np.dot((x1-mean_x1),(x2-mean_x2))/std_x1/std_x2

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

def kendall(x,y):
    numer = 1
    d = 1
    for i in range(1,len(x)):
        for j in range(i):
            numer += np.sign(x[i] - x[j]) * np.sign(y[i] - y[j])
            if x[i] != x[j] and y[i] != y[j]:
                d += 1
    return numer*1.0/d

methods = {"random":random_,"cross_correlation":cross_correlation,"iota":iota,"kendall":kendall,"granger":granger}
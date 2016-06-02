from sklearn.preprocessing import binarize
import numpy as np
from sklearn.linear_model import LogisticRegression
import matplotlib.pylab as plt


#number of nodes
n = 6
#hill coeficient
h = 5

#time resolution, for ode solving
dt = .001
#when the simulation stops
end_time = 5.0
#number of time points
n_time_points = 10
#fraction of links
frac_plus = 0.3
frac_minus = 0.3


W0_plus = np.multiply(np.random.rand(n,n),np.random.binomial(1, frac_plus, size=(n,n)))
W0_minus = np.multiply(np.random.rand(n,n),np.random.binomial(1, frac_minus, size=(n,n)))



for i in range(n):
	W0_plus[i,i] = 0
	W0_minus[i,i] = 0
	for k in range(n):
		if W0_plus[i][k] > 0:
			W0_minus[i][k] = 0

def ode_model(xx):
	xxh = np.power(xx,h)
	f_plus_xx = np.divide(xxh,1.0+xxh)
	f_minus_xx = np.divide(1.0,1.0+xxh)
	dx = np.dot(W0_plus,f_plus_xx)+np.dot(W0_minus,f_minus_xx)-xx
	return dx

x0 = np.random.rand(n,1)

tspan = np.arange(0,5,dt)

X = np.zeros((len(tspan),n))

for i,t in enumerate(tspan):
	X[i:i+1,:] = x0.T+ode_model(x0).T*dt
	x0 = X[i,:]


time_idxs = [i*end_time/dt/n_time_points for i in range(n_time_points)]
time_idxs[-1] -= 1

observedX = X[time_idxs,:]
ddt = end_time*1.0/n_time_points

X_for_lasso = observedX[:-1,:]

Y_for_lasso = binarize(observedX[1:,:]-observedX[:-1,:])





clf = LogisticRegression(penalty='l1',C=2)

clf.fit(X_for_lasso,Y_for_lasso)


plt.subplot(1,2,1)
plt.plot(time_idxs,observedX)
plt.legend([str(i+1) for i in range(n)])


predicted = [i for e in clf.coef_ for i in e]
W0_plus_flat = [i for e in W0_plus for i in e]

predicted_null = [abs(predicted[i]) for i,v in enumerate(W0_plus_flat) if v==0]
predicted_not_nul = [abs(predicted[i]) for i,v in enumerate(W0_plus_flat) if not v==0]


plt.subplot(1,2,2)
plt.set_cmap('bwr')
plt.plot(range(len(predicted_null)),predicted_null,'o',color='#336688')
plt.plot(range(len(predicted_not_nul)),predicted_not_nul,'o',color='#FF9900')
plt.title('Difference')



plt.show()
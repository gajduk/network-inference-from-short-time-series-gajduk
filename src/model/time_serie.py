from bisect import bisect_right as binary_search

import numpy as np


class TimeSerieSampler:
	'''
	Samples the time series only at certain time points.
	In reality measurements can only be done once several minutes, so you can not use continous time series for prediction.
	The time time points can be specified directly via the times parameter or the n_time_points.
	:param n_time_points: if set then the time series will be sampled at uniformly sized times steps:
	(end_time/n_time_points-1) in the interval [0,end_time], defaults to 10
	:param times - an array of doubles e.g. [0,3,5,8,10,30,100] - if you want to sample the initial dynamics more precisely,
	                 this makes sense since most of the interesting stuff happens in the beginning
	:return: a tuple of t,new_x where new_x is just x sampled at several discrete time points
	'''

	def __init__(self, n_time_points=-1, times=None):
		if times is None:
			if n_time_points < 1 or n_time_points > 100:
				n_time_points = 10
		self.times = times
		self.n_time_points = n_time_points

	def sample(self, time_serie):
		if self.times is None:
			times = np.linspace(0, np.max(time_serie.t), self.n_time_points)
		times_idxs = [min(binary_search(time_serie.t, t), len(time_serie.t) - 1) for t in times]
		new_x = np.zeros((self.n_time_points, time_serie.x.shape[1]))
		for i, time_idx in enumerate(times_idxs):
			if time_idx == 0:
				new_x[i, :] = time_serie.x[time_idx, :]
			else:
				t2 = time_serie.t[time_idx]
				t1 = time_serie.t[time_idx - 1]
				t = times[i]
				x2 = time_serie.x[time_idx, :]
				x1 = time_serie.x[time_idx - 1, :]
				k = np.max((t - t1) * (t2 - t1))
				temp = x1 + (x2 - x1) * k
				#temp[temp < 0] = 0
				new_x[i, :] = temp
		return times, new_x.T

	def get_n_time_points(self):
		if not self.times is None:
			return len(self.times)
		return self.n_time_points


class TimeSerie:
	'''
	Stores a single time serie for an unspecified graph.
	t: time array
	x: protein activities at different time points (matrix  [len(t) X num_nodes])
	'''

	def __init__(self, t, x):
		self.t = np.array(t)
		self.x = x

import matplotlib.pyplot as plt
import numpy as np

from src.graph_plotting import plotDiGraph_, plotDiGraphViaGraphViz
from src.utils import deprecated


class InstanceSingleTimeSeriesView:
	'''
	Defines a view of a PIN with a just a single time series
	t an array of len n_time_points
	x a matrix with shape n_nodes, n_time_points
	y a matrix with shape n_nodes, n_nodes
	'''

	def __init__(self, t, x, graph_adjacency_matrix, n_nodes):
		self.t = t
		self.x = x
		self.y = graph_adjacency_matrix
		self.n_nodes = n_nodes
		self.n_time_points = len(self.t)
		if n_nodes == 6:
			self.labels = ['EGFR', 'Raf', 'MEK', 'ERK', 'Akt', 'Src']
		else:
			self.labels = [str(i + 1) for i in range(n_nodes)]
		self.pos = None

	def get(self, idx_node):
		return self.x[idx_node,:][:].tolist(),self.y[idx_node*self.n_time_points:(idx_node+1)*self.n_time_points]

	def getebunch(self):
		return [(k, i, self.y[self.n_nodes * i + k]) for i in range(self.n_nodes) for k in range(self.n_nodes) if
		        self.y[self.n_nodes * i + k] != 0]

	@deprecated
	def plotGraph_(self, cmap=plt.cm.Accent):
		plotDiGraph_(self.n_nodes, self.getebunch(), cmap, node_labels=self.labels)

	def plotDiGraphViaGraphViz_(self, cmap=plt.cm.Accent):
		self.pos = plotDiGraphViaGraphViz(self.n_nodes, self.getebunch(), cmap, node_labels=self.labels)

	def plotTimeSeries_(self, cmap=plt.cm.Accent):
		plt.xlabel('Time [min]')
		plt.ylabel('Activity [A.U.]')
		colors = [cmap(i) for i in np.linspace(0, 1, self.n_nodes)]
		for i, color in enumerate(colors):
			x, y = self.get(i)
			plt.plot(self.t, x, color=color)
		legend = plt.legend(self.labels, bbox_to_anchor=(1.4, 1.1))
		plt.grid()
		return legend

	def plot(self):
		plt.figure(figsize=(30, 16))
		cmap = plt.cm.Accent
		plt.subplot(1, 2, 1)
		self.plotDiGraphViaGraphViz_(cmap)
		plt.subplot(1, 2, 2)
		plt.subplots_adjust(left=0.001, right=0.8, top=0.9, bottom=0.1)
		legend = self.plotTimeSeries_(cmap)
		legend.get_frame().set_linewidth(3)

	def exportAsCsv(self, filename):
		np.savetxt(filename + "_time_series.csv", np.reshape(self.x, (self.n_nodes, self.n_time_points)), delimiter=",")
		np.savetxt(filename + "_links.csv", np.reshape(self.y, (self.n_nodes, self.n_nodes)), delimiter=",")


class Instance:
	'''
	Defines a PIN with a list of time series, which should correspond to different simulation setups: different input signals, inhibitions of certain proteins etc
	'''

	def __init__(self, time_series, graph_adjacency_matrix, n_nodes):
		self.time_series = time_series
		self.graph_adjecency_matrix = graph_adjacency_matrix
		self.n_nodes = n_nodes
		self.n_time_series = len(self.time_series)

	def getViewForTimeSeries(self, time_series_idx, sampler):
		'''
		Returns a view of the instance that holds the graph adjacency matrix and a sampled view of one of the time series
		:param time_series_idx: which time serie to sample
		:param sampler: defines how to sample the time series see TimeSerieSampler
		:return: a :class:`~src.model.InstanceSingleTimeSeriesView` object
		'''
		if time_series_idx < 0 or time_series_idx >= len(self.time_series):
			raise ValueError(
				"time_series_idx ({0}) is out of range [0,{1}]".format(time_series_idx, len(self.time_series)))
		t, x = sampler.sample(self.time_series[time_series_idx])
		sampler.get_n_time_points()
		return InstanceSingleTimeSeriesView(t, x, self.graph_adjecency_matrix, self.n_nodes)

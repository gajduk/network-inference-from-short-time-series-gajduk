import numpy as np
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from src.utils import PROJECT_DIR,s_timestamp,plotDiGraph_
from random import shuffle

class Instance:

	def __init__(self,xs,y,n_nodes,n_time_points,n_time_series=-1):
		self.x = None
		self.t = None
		self.xs = xs
		self.y = y
		self.n_nodes = n_nodes
		self.n_time_points = n_time_points
		if n_time_points > -1:
			self.n_time_series = n_time_series
		else:
			self.n_time_series = len(xs)
		if n_nodes == 6:
			self.labels = ['EGFR','Raf','MEK','ERK','Akt','Src']
		else:
			self.labels = [str(i+1) for i in range(n_nodes)]

	def get(self,idx_node):
		return self.x[self.n_time_points*idx_node:self.n_time_points*(idx_node+1)],self.y[self.n_nodes*idx_node:self.n_nodes*(idx_node+1)]

	def getebunch(self):
		return [(k,i,self.y[self.n_nodes*i+k]) for i in range(self.n_nodes) for k in range(self.n_nodes) if self.y[self.n_nodes*i+k]!=0]

	def plotGraph_(self,cmap=plt.cm.Accent):
		plotDiGraph_(self.n_nodes,self.getebunch(),cmap,node_labels=self.labels)

	def plotTimeSeries_(self,cmap=plt.cm.Accent):
		plt.xlabel('Time [min]')
		plt.ylabel('Activity [A.U.]')
		colors = [cmap(i) for i in np.linspace(0, 1, self.n_nodes)]
		for i, color in enumerate(colors):
			x,y = self.get(i)
			plt.plot(self.t,x, color=color)
		plt.legend(self.labels,bbox_to_anchor=(-0.2, 1.1))

	def setx(self,s_idx):
		self.t,self.x = self.xs[s_idx].getx(self.n_time_points)

	def plot(self):
		plt.figure(figsize=(9.5,4.5))
		cmap = plt.cm.Accent
		plt.subplot(1,2,1)
		plt.axis('off')
		#self.plotGraph_(cmap)
		plt.subplot(1,2,2)
		self.setx(self.n_time_series-1)
		self.plotTimeSeries_(cmap)

	def exportAsCsv(self,filename):
		np.savetxt(filename+"_time_series.csv", np.reshape(self.x,(self.n_nodes,self.n_time_points)), delimiter=",")
		np.savetxt(filename+"_links.csv", np.reshape(self.y,(self.n_nodes,self.n_nodes)), delimiter=",")

class Dataset:

	def __init__(self,instances,n_time_points,n_time_series,n_nodes):
		self.instances = instances
		self.n_instances = len(self.instances)
		self.n_nodes = n_nodes
		self.n_time_points = n_time_points
		self.n_time_series = n_time_series

	def get(self,idx_instance):
		return self.instances[idx_instance]

	def plotAll(self):
		with PdfPages(join(PROJECT_DIR,'output','visualization_graph_and_time_series',s_timestamp()+'.pdf')) as pdf:
			for i in range(self.n_instances):
				instance = self.get(i)
				instance.plot()
				pdf.savefig()  # saves the current figure into a pdf page
				plt.close()



if __name__ == "__main__":
	pass
	#d.plotAll()
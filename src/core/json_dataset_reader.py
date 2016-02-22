import zipfile
import json
from os.path import join
from src.utils import PROJECT_DIR
from bisect import bisect_right as binary_search
from datasets import Instance,Dataset
import numpy as np
import matplotlib.pylab as plt
from random import shuffle

class TimeSeries:

	def __init__(self,t,x):
		self.t = t
		self.x = x

	def getx(self,n_time_points):
		times = np.linspace(0,self.t[-1],n_time_points)
		times_idxs = [min(binary_search(self.t,t),len(self.t)-1) for t in times]
		return times,np.reshape(self.x[times_idxs,:].T,(-1,1))

class CsvDatasetReader:

	DEFAULT_DIR = join(PROJECT_DIR,'data')

	def getDataset(self,filename,n_nodes=14,n_time_points=13):
		prefix = join(self.DEFAULT_DIR,filename)
		instances = []
		with open(prefix+'_X.csv','r') as pinX:
			with open(prefix+'_Y.csv','r') as pinY:
				while True:
					xs = []
					temp_y = pinY.readline()
					if temp_y is None or len(temp_y) <10:
						break
					y = [float(a) for a in temp_y.split(',')]
					y = np.reshape(np.matrix(y),(n_nodes*n_nodes,1))

					t = np.linspace(0,120,n_time_points)

					x = [float(a) for a in pinX.readline().split(',')]
					x = np.reshape(np.matrix(x),(n_nodes,n_time_points-1)).T
					x = np.vstack((np.zeros((1,n_nodes)),x))

					ts = TimeSeries(t,x)
					xs.append(ts)

					instances.append(Instance(xs,y,n_nodes,n_time_points,1))

		return Dataset(instances,n_time_points,1,n_nodes)

class JsonDatasetReader:

	DEFAULT_DIR = join(PROJECT_DIR,'data')

	def __init__(self,filename='initial_testing1.json.zip'):
		self.filename = filename
		self._json = self._readJson(join(JsonDatasetReader.DEFAULT_DIR,filename))

	def _readJson(self,filename):
		if zipfile.is_zipfile(filename):
			with zipfile.ZipFile(filename, 'r') as pin:
				data = pin.read(pin.namelist()[0])
		else:
			with open(filename, 'r') as pin:
				data = pin.read()
		return json.loads(data)

	def getDataset(self,n_time_points=12,n_time_series=3,n_nodes=6):
		root = self._json['root']
		instances = []
		count = 1
		for instance_i in root:
			instance = instance_i[0]
			y = np.reshape(np.matrix(instance['pin']['f']).T,(n_nodes*n_nodes,1))
			xs = []
			pin_simulation_results = instance['pin_simulation_results']
			from copy import deepcopy
			temp = deepcopy(pin_simulation_results)
			shuffle(temp)
			for pin_simulation_result in temp:
				t = np.array(pin_simulation_result['t'])
				x = np.matrix(pin_simulation_result['y'])
				ts = TimeSeries(t,x)
				xs.append(ts)
			instances.append(Instance(xs,y,n_nodes,n_time_points,n_time_series))
			count += 1
		return Dataset(instances,n_time_points,n_time_series,n_nodes)

def main():
	reader = CsvDatasetReader()
	d = reader.getDataset('out_inhibition')
	d.plotAll()
	i = d.get(5)
	i.setx(0)
	i.exportAsCsv(join(PROJECT_DIR,'output','inhibition'))

if __name__ == "__main__":
	main()
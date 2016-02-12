import zipfile
import json
from os.path import join
from src.utils import PROJECT_DIR
from bisect import bisect_right as binary_search
from datasets import Instance,Dataset
import numpy as np
import matplotlib.pylab as plt

class TimeSeries:

	def __init__(self,t,x):
		self.t = t
		self.x = x

	def getx(self,n_time_points):
		times = np.linspace(0,self.t[-1][0],n_time_points)
		times_idxs = [min(binary_search(self.t,t),len(self.t)-1) for t in times]
		return times,np.reshape(self.x[times_idxs,:].T,(-1,1))



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

	def getDataset(self,n_time_points=12,n_time_series=3,n_nodes=14):
		root = self._json['root']
		instances = []
		count = 1
		for instance_i in root:
			instance = instance_i[0]
			y = np.reshape(np.matrix(instance['pin']['f']).T,(n_nodes*n_nodes,1))
			xs = []
			pin_simulation_results = instance['pin_simulation_results']
			for pin_simulation_result in pin_simulation_results:
				t = np.array(pin_simulation_result['t'])
				x = np.matrix(pin_simulation_result['y'])
				ts = TimeSeries(t,x)
				xs.append(ts)
			instances.append(Instance(xs,y,n_nodes,n_time_points,n_time_series))
			count += 1
			if count > 1:
				break
		return Dataset(instances,n_time_points,n_time_series,n_nodes)

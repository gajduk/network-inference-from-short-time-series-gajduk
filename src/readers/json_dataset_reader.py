import json
import zipfile
from os.path import join
from random import shuffle

import numpy as np

from src.model.datasets import Dataset
from src.model.instance import Instance
from src.model.time_serie import TimeSerie, TimeSerieSampler
from src.utils import PROJECT_DIR

DEFAULT_DIR = join(PROJECT_DIR, 'data')


class JsonDatasetReader:
	def __init__(self, filename='initial_testing1.json.zip'):
		self.filename = filename
		self._json = self._readJson(join(DEFAULT_DIR, filename))

	def _readJson(self, filename):
		if zipfile.is_zipfile(filename):
			with zipfile.ZipFile(filename, 'r') as pin:
				data = pin.read(pin.namelist()[0])
		else:
			with open(filename, 'r') as pin:
				data = pin.read()
		return json.loads(data)

	def getDataset(self, n_instances=100):
		root = self._json['root']
		instances = []
		count = 1
		for instance_i in root:
			instance = instance_i[0]
			n_nodes = instance['pin']['n']
			y = np.reshape(np.matrix(instance['pin']['f']).T, (n_nodes * n_nodes, 1))
			time_series = []
			pin_simulation_results = instance['pin_simulation_results']
			from copy import deepcopy
			temp = deepcopy(pin_simulation_results)
			shuffle(temp)
			for pin_simulation_result in temp:
				t = np.array(pin_simulation_result['t'])
				x = np.matrix(pin_simulation_result['y'])
				inhibit = 'unknown'
				if 'inhibit' in pin_simulation_result:
					inhibit = pin_simulation_result['inhibit']
				ts = TimeSerie(t, x, inhibit)
				time_series.append(ts)
			time_series = sorted(time_series,key=lambda  ts : int(ts.inhibit.split('_')[1]) if len(ts.inhibit.split('_')) > 1 else 0)
			instances.append(Instance(time_series, y, n_nodes))
			count += 1
			if count > n_instances:
				break
		return Dataset(instances)


def main():
	reader = JsonDatasetReader('example2.json.zip')
	d = reader.getDataset(n_instances=100)
	d.plotAll(TimeSerieSampler(n_time_points=21))


if __name__ == "__main__":
	main()

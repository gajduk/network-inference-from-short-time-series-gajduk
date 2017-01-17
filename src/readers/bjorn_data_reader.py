from os.path import join

import numpy as np
import matplotlib

from src.model.instance import InstanceSingleTimeSeriesView
from src.utils import PROJECT_DIR
from src.methods.methods import methods
from src.graph_plotting import plotDiGraphViaGraphVizASD

import src.utils
import matplotlib.pylab as plt

method_symmetrical = {"derivative_correlation":False, "cross_correlation_correct":False, "kendall":True, "partial_corr":True,
                 "symbol_similarity":True,"granger":False}

class BjornDataReader:


	def __init__(self):
		pass

	def load(self,file):
		line_num = 2
		t,x = [],[]
		with open(join(PROJECT_DIR,file),"r") as pin:
			pin.readline()
			for line in pin:
				s_line = line.split(',')
				try:
					t.append(int(s_line[0]))
					x.append([float(e) for e in s_line[1:]])
				except ValueError:
					raise ValueError('Values in row {0} are not a number.'.format(line_num))
				line_num += 1
		y = np.zeros((4, 4))
		edges = [(1, 0, 1), (1,3,1), (3,2,1),(2,1,-1)]
		#edges = [(2,0,1),(2,1,1),(7,2,1),(2,6,1),(6,2,-1),(1,5,1),(5,3,1),(3,1,-1)]
		for from_,to_,v_ in edges:
			y[from_][to_] = v_
		instance = InstanceSingleTimeSeriesView(t, np.array(x).T, y, 4)
		return instance

def reconstruct_networks_for_instance(instance,cell_line):
	cmap = matplotlib.cm.get_cmap('Dark2')
	#pos = [(-1, 1), (0, 1), (0, 0), (0, 3), (-2, 2), (0, 2), (2, 1), (1, 1)]
	pos = [(-1, 1), (0, 0), (0, 2), (0, 1)]
	pos = [(e[0], 3 - e[1]) for e in pos]
	#node_labels = ['Akt', 'bRaf', 'egfr', 'Erk', 'Jnk', 'Mek', 'Shp2', 'Src']
	node_labels = ['Akt',  'egfr', 'Erk', 'Mek']
	n = instance.n_nodes
	instance.labels = node_labels
	for method in method_symmetrical:
		print method
		src.utils.s_timestamp_prefix = cell_line+'_'+method
		A = methods[method](instance)

		temp = sorted([abs(e) for i, a in enumerate(A) for k, e in enumerate(a) if not i == k])
		if method_symmetrical[method]:
			th = temp[-8]
			ebunch = [(i, k, A[i][k]) for i in range(n) for k in range(i+1,n) if abs(A[i][k]) >= th and not i == k]
		else:
			th = temp[-6]
			ebunch = [(i, k, A[i][k]) for i in range(n) for k in range(n) if abs(A[i][k]) >= th and not i == k]
		plotDiGraphViaGraphVizASD(n, ebunch, cmap, pos=pos, node_labels=node_labels,method_symmetrical=method_symmetrical[method])

if __name__ == "__main__":
	instance = BjornDataReader().load('import_for_nr_hkh2_4.csv')
	reconstruct_networks_for_instance(instance,'hkh2_4')

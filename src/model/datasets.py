from os.path import join

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from src.utils import PROJECT_DIR, s_timestamp


class Dataset:
	def __init__(self, instances):
		self.instances = instances
		self.n_instances = len(self.instances)

	def get(self, idx_instance):
		return self.instances[idx_instance]

	def plotAll(self, sampler):
		font = {'size': 38}
		plt.rc('font', **font)
		plt.rcParams['lines.linewidth'] = 3
		plt.rcParams['axes.linewidth'] = 3

		with PdfPages(
				join(PROJECT_DIR, 'output', 'visualization_graph_and_time_series', s_timestamp() + '.pdf')) as pdf:
			for i in range(self.n_instances):
				instance = self.get(i)
				for n_time_idx in range(instance.n_time_series):
					istsv = instance.getViewForTimeSeries(n_time_idx, sampler)
					istsv.plot()
					plt.title(instance.time_series[n_time_idx].inhibit)
					pdf.savefig(dpi=100)
					plt.close()

		font = {'size': 12}
		plt.rc('font', **font)


if __name__ == "__main__":
	pass
# d.plotAll()

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

	def plotAll(self):
		font = {'size': 38}
		plt.rc('font', **font)
		plt.rcParams['lines.linewidth'] = 3
		plt.rcParams['axes.linewidth'] = 3

		with PdfPages(
				join(PROJECT_DIR, 'output', 'visualization_graph_and_time_series', s_timestamp() + '.pdf')) as pdf:
			for i in range(self.n_instances):
				instance = self.get(i)
				for n_time_idx in range(instance.n_time_series):
					instance.setx(n_time_idx)
					instance.plot()
					pdf.savefig(dpi=100)  # saves the current figure into a pdf page
					plt.close()

		font = {'size': 12}
		plt.rc('font', **font)


if __name__ == "__main__":
	pass
	# d.plotAll()

from os.path import join

from src.utils import PROJECT_DIR, deprecated

DEFAULT_DIR = join(PROJECT_DIR, 'data')


class CsvDatasetReader:
	@deprecated
	def getDataset(self, filename):
		'''
		prefix = join(DEFAULT_DIR, filename)
		instances = []
		with open(prefix + '_X.csv', 'r') as pinX:
			with open(prefix + '_Y.csv', 'r') as pinY:
				while True:
					time_series = []
					temp_y = pinY.readline()
					if temp_y is None or len(temp_y) < 10:
						break
					y = [float(a) for a in temp_y.split(',')]

					y = np.reshape(np.matrix(y), (n_nodes * n_nodes, 1))

					t = np.linspace(0, 120, n_time_points)

					x = [float(a) for a in pinX.readline().split(',')]
					x = np.reshape(np.matrix(x), (n_nodes, n_time_points - 1)).T
					x = np.vstack((np.zeros((1, n_nodes)), x))

					ts = TimeSerie(t, x)
					time_series.append(ts)

					instances.append(Instance(time_series, y, n_nodes))

		return Dataset(instances)
		'''
		pass


def main():
	reader = CsvDatasetReader('for_presentation.zip')
	d = reader.getDataset()
	d.plotAll()


if __name__ == "__main__":
	main()

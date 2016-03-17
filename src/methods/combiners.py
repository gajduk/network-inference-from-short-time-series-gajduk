import numpy as np
import matplotlib.pylab as plt

from utils import normalize_rowwise
from src.utils import debug

def mutliple_time_series_combiner(method, i):
	cc = np.zeros((i.n_nodes, i.n_nodes))
	if debug:
		plt.subplot(3, 5, 1)
		plt.imshow(np.absolute(np.reshape(i.y, (14, 14))))
		plt.subplot(3, 5, 6)
		plt.imshow(np.absolute(np.reshape(i.y, (14, 14))))

	for idx_time_series in range(i.n_time_series):
		i.setx(idx_time_series)
		res_method = np.absolute(method(i))
		cc += normalize_rowwise(res_method)
		if debug:
			plt.subplot(3, 5, idx_time_series + 2)
			plt.imshow(res_method)
			plt.subplot(3, 5, 5 + idx_time_series + 2)
			plt.imshow(normalize_rowwise(res_method))
			plt.subplot(3, 5, 10 + idx_time_series + 2)
			plt.imshow(cc)
	cc = normalize_rowwise(cc)
	return cc

import numpy as np

def one_pair_at_a_time_wrapper(method):
	def inner(i):
		res = np.zeros((i.n_nodes, i.n_nodes))
		for idx_node1 in range(i.n_nodes):
			x1, _ = i.get(idx_node1)
			for idx_node2 in range(i.n_nodes):
				x2, _ = i.get(idx_node2)
				res[idx_node1][idx_node2] = method(x1, x2)
		return res

	return inner

def time_series_columnwise_wrapper(method):
	def inner(i):
		data = np.zeros((i.n_time_points, i.n_nodes))
		for idx_node in range(i.n_nodes):
			x, _ = i.get(idx_node)
			data[:, idx_node] = np.array(x).T.reshape((i.n_time_points,))
		return method(data)

	return inner

def normalize_rowwise(x):
	return np.absolute(x) / np.max(np.absolute(x), axis=1, keepdims=True)  # /np.std(x,axis=1,keepdims=True)

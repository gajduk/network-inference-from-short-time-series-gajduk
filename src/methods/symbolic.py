from itertools import permutations

from numpy import argsort, matrix, zeros
from sklearn.metrics import normalized_mutual_info_score


def one_pair_at_a_time_symbol(method):
	def inner(instance):
		x = matrix(instance.x)
		n_nodes, n_time_points = x.shape
		if n_time_points > 12:
			delta = 6
		else:
			delta = int(n_time_points / 2)
		res = zeros((n_nodes, n_nodes))
		s, lenP = convertToSymbolSequence(x, delta)
		for idx_node1 in range(n_nodes):
			for idx_node2 in range(n_nodes):
				res[idx_node1][idx_node2] = method(s[idx_node1, :].tolist()[0], s[idx_node2, :].tolist()[0], lenP)
		return res

	return inner


def chooseK(items, k):
	if k == 0:
		yield []
	else:
		for i in xrange(len(items)):
			for cc in chooseK(items[i + 1:], k - 1):
				yield [items[i]] + cc


def convertToSymbolSequence(x, delta):
	res = []
	l = range(delta)
	P = {tuple(permutation): i for i, permutation in enumerate(permutations(l))}
	for x_row in x:
		symbols = [P[tuple(argsort(temp).tolist())] for temp in chooseK(x_row.tolist()[0], delta)]
		res.append(symbols)
	return matrix(res), len(P)


@one_pair_at_a_time_symbol
def symbolSequenceSimilarity(s1, s2, lenP):
	p1, p2 = 0, 0
	for i in range(len(s1)):
		p1 += 1 if s1[i] == s2[i] else 0
		p2 += 1 if s1[i] == lenP - s2[i] else 0
	return max(p1, p2) * 1.0 / len(s1)


@one_pair_at_a_time_symbol
def mutualInformationOfSymbols(s1, s2, lenP):
	return normalized_mutual_info_score(s1, s2)


if __name__ == "__main__":
	x = matrix([[0.5, 0.7, 0.4, 0.6], [.1, .5, .7, .0]])
	print mutualInformationOfSymbols(x)

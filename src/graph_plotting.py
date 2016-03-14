import networkx as nx
import matplotlib.pylab as plt
import os
import struct
from subprocess import call
import matplotlib.image as mpimg
import time

from utils import PROJECT_DIR, s_timestamp


def plotDiGraph_(n_nodes, ebunch, cmap, node_size=400, font_size=12, node_labels=None):
	G = nx.DiGraph()
	G.add_nodes_from(range(n_nodes))
	G.add_weighted_edges_from(ebunch)
	if node_labels is None:
		node_labels = [str(i + 1) for i in range(G.number_of_nodes())]
	pos = nx.shell_layout(G)
	nx.draw_networkx_nodes(G, pos, node_color=range(G.number_of_nodes()), cmap=cmap, node_size=node_size)
	nx.draw_networkx_edges(G, pos, edgelist=[(i, k) for i, k, w in ebunch if w > 0], edge_color='g', width=1)
	nx.draw_networkx_edges(G, pos, edgelist=[(i, k) for i, k, w in ebunch if w < 0], edge_color='r', width=1)
	nx.draw_networkx_labels(G, pos, {i: node_labels[i] for i in range(G.number_of_nodes())}, font_size=font_size)
	plt.axis('off')


def plotDiGraphViaGraphViz(n_nodes, ebunch, cmap, pos=None, node_labels=None):
	def getForwardLinks():
		return [(from_node, to_node) for from_node, to_node, _ in ebunch if from_node < to_node]

	if node_labels is None:
		node_labels = [str(i + 1) for i in range(n_nodes)]
	file_ = os.path.join(PROJECT_DIR, 'output', 'temp', s_timestamp() + ".dot")
	flag = ''
	#flag = '\tsplines = true;\n'
	if pos is None:
		pos = _posTree(getForwardLinks())
		flag = '\tsplines = true;\n'
	with open(file_, 'w') as pout:
		pout.write('digraph G {{\n'.format(n_nodes))
		pout.write('\tsplines="curved";\n')
		pout.write(flag)
		pout.write('\tInput [pos="0,0!" color=white];\n')
		for node_idx in range(n_nodes):
			r, g, b, a = cmap(node_idx * 1.0 / n_nodes)

			color = struct.pack('BBB', *(r * 255, g * 255, b * 255)).encode('hex')
			pout.write(
				'\t{0} [pos="{2},{3}!" shape=circle style=filled fillcolor="#{1}" width=.5  fixedsize=true];\n'.format(
					node_labels[node_idx], color, pos[node_idx][0], pos[node_idx][1]))
		for from_node, to_node, weight in ebunch:
			color = "black"
			extra = ""
			if from_node > to_node:
				extra = ''
				if weight < 0:
					color = "red"
					extra += ' arrowhead="tee"'
				else:
					color = "green"

			pout.write(
				'\t{0} -> {1} [color="{2}" {3}];\n'.format(node_labels[from_node], node_labels[to_node], color, extra))
		pout.write('\tInput -> {0} [color="black"];\n'.format(node_labels[0]))
		pout.write('}\n')
	out_file = file_ + ".png"
	command = 'neato -Tpng "' + file_ + '" > "' + out_file
	call(command, shell=True, env=os.environ)
	img = -1
	while 1:
		try:
			img = mpimg.imread(out_file)
			break
		except:
			time.sleep(1)
			print 'Sleep'

	plt.imshow(img)
	plt.axis('off')
	return pos


def rankTogether(ebunch):
	def getChildren(node):
		return [to_node for from_node, to_node in ebunch if from_node == node]

	queue = [0]
	queue_len = 1
	same_rank = []
	count = 0
	while len(queue) > 0:
		current = queue.pop(0)
		count += 1
		for child in getChildren(current):
			queue.append(child)
		if count == queue_len:
			if len(queue) > 0:
				same_rank.append([i for i in queue])
			queue_len = len(queue)
			count = 0
	return same_rank


def _posTree(ebunch):
	root = 0
	vspace = .8
	res = {0: (0, -vspace)}

	def getChildren(node):
		return [to_node for from_node, to_node in ebunch if from_node == node]

	def getChildrenPosRec(parent, parent_pos, res):
		children = getChildren(parent)
		n = len(children)
		if n == 0:
			return
		if n == 1:
			res.update({children[0]: (parent_pos[0], parent_pos[1] - vspace)})
		else:
			res.update({child: (parent_pos[0] - 1 + 2.0 / (n - 1) * i, parent_pos[1] - vspace) for i, child in
			            enumerate(children)})
		for child in children:
			getChildrenPosRec(child, res[child], res)

	getChildrenPosRec(root, res[root], res)
	return res

import matplotlib.pylab as plt
import numpy as np

from src.readers.mass_spec_data_reader import InsulinMassSpecAndModelDataReader


def plot_all_alpha():
	dataset,t,_ = InsulinMassSpecAndModelDataReader().loadnmeth_2703()
	x = []
	for protein_name in dataset:
		x.append(dataset[protein_name])
	plt.plot(np.matrix(x).T, 'b', alpha=0.3, color='#6D8AC4')
	plt.xlabel('Time [min]')
	plt.ylabel('Log2 - FoldChange compared to t = 0 [min]')
	plt.xticks([0, 1, 2, 3, 4], t)
	plt.ylim([-1.5, 1.5])
	# plt.legend([e for e in dataset])
	plt.show()


def plot_relevant_proteins():
	dataset,t,relevant_proteins = InsulinMassSpecAndModelDataReader().loadnmeth_2703()
	x = []
	for protein_name in relevant_proteins:
		x.append(dataset[relevant_proteins[protein_name]])
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_color_cycle(["#D24FA2", "#68C738", "#D768DA", "#ABB336", "#A275D5", "#E0A233", "#7187D9", "#DD5425", "#4FB4D2", "#D9464D", "#57C97A", "#E03375", "#51953F", "#AF75AB", "#61752C", "#C55E7F", "#34A08C","#C1624D", "#5C7BA7", "#AC6F28"])
	for xs in x:
		xs[0] = 0
	plt.plot(np.matrix(x).T)
	plt.xlabel('Time [min]')
	plt.ylabel('Log2 - FoldChange compared to t = 0 [min]')
	plt.xticks([0, 1, 2, 3, 4], t)
	plt.ylim([-.7, 1.2])
	ax.legend([e for e in relevant_proteins]).draggable()
	plt.show()

if __name__ == "__main__":
	plot_relevant_proteins()
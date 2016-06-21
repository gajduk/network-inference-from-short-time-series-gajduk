from os import path

import numpy as np
from sklearn.metrics import roc_curve, auc

from src.core.predictor import SingleSeriesPredictor
from src.methods.methods import good_methods
from src.model.time_serie import TimeSerieSampler
from src.readers.json_dataset_reader import JsonDatasetReader
from src.utils import five_pin_setups,OUTPUT_DIR


class Evaluator:

	def __init__(self):
		pass

	def evaluate(self, prediction):
		'''
		:param prediction: an object of type Prediction
		:return: roc_auc: area under the curve
		'''
		fpr, tpr, _ = roc_curve(np.abs(prediction.y_true) > 0.01, prediction.y_pred)
		roc_auc = auc(fpr, tpr)
		return roc_auc,fpr,tpr

def evaluateAllMethods(dataset):
	time_sampler = TimeSerieSampler(n_time_points=15)
	with open(path.join(OUTPUT_DIR,"evaluation","presentation_1A.csv"),"w") as pout:
		pout.write(','+','.join(five_pin_setups)+'\n')
		for method in good_methods:
			print method
			pout.write(method)
			predictor = SingleSeriesPredictor(good_methods[method], time_sampler)
			evaluator = Evaluator()
			for time_series_idx in range(5):
				prediction = predictor.predictAllInstancesCombined(dataset,time_series_idx)
				roc_auc,fpr,tpr = evaluator.evaluate(prediction)
				pout.write(', {0:.3f}'.format(roc_auc))
			pout.write('\n')

def evaluateNumTimepoints(dataset):
	ntp = [e for e in range(4,15)]
	ntp.extend([e for e in range(15,20,2)])
	methods = ['cross_correlation','kendall','symbol_mutual']
	evaluator = Evaluator()
	time_series_idx = 0
	res = [[] for m in methods]
	for n_time_points in ntp:
		print n_time_points
		time_sampler = TimeSerieSampler(n_time_points=n_time_points)
		for i,method in enumerate(methods):
			predictor = SingleSeriesPredictor(good_methods[method], time_sampler)
			prediction = predictor.predictAllInstancesCombined(dataset,time_series_idx)
			roc_auc, fpr, tpr = evaluator.evaluate(prediction)
			res[i].append(roc_auc)
	print res

def plotFirstTacROC(dataset):
	import matplotlib.pylab as plt
	from os.path import join
	from src.utils import PROJECT_DIR
	plt.figure(figsize=(6, 6))
	time_sampler = TimeSerieSampler(n_time_points=12)
	evaluator = Evaluator()
	time_series_idx = 0
	methods = {'cross_correlation':'Cross corr.   ','kendall':'Kendall        ','symbol_mutual':'Symbol MI    ','symbol_similarity':'Symbol sim.'}
	for method in methods:
		print method
		predictor = SingleSeriesPredictor(good_methods[method], time_sampler)
		prediction = predictor.predictAllInstancesCombined(dataset,time_series_idx)
		roc_auc, fpr, tpr = evaluator.evaluate(prediction)
		plt.plot(fpr, tpr, label=methods[method] + ' (auc = %0.3f)' % roc_auc)
	plt.legend(loc="lower right")
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.grid()
	plt.savefig(join(PROJECT_DIR, 'output', 'firstTACROC.pdf'))

def main():
	reader = JsonDatasetReader('presentation_1A.json.zip')
	dataset = reader.getDataset(n_instances=100)
	plotFirstTacROC(dataset)

if __name__ == "__main__":
	main()

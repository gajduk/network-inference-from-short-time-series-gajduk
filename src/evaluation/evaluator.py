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
	time_sampler = TimeSerieSampler(n_time_points=10)
	with open(path.join(OUTPUT_DIR,"evaluation","2_inhibitions_top.json.csv"),"w") as pout:
		pout.write(' ,'+','.join(five_pin_setups)+'\n')
		for method in good_methods:
			print method
			pout.write(method)
			predictor = SingleSeriesPredictor(good_methods[method], time_sampler)
			evaluator = Evaluator()
			for time_series_idx in range(1):
				prediction = predictor.predictAllInstancesCombined(dataset,time_series_idx)
				roc_auc,fpr,tpr = evaluator.evaluate(prediction)
				pout.write(', {0:.2f}'.format(roc_auc))
			pout.write('\n')


def main():
	reader = JsonDatasetReader('2_inhibitions_top.json.zip')
	dataset = reader.getDataset(n_instances=1)
	evaluateAllMethods(dataset)

if __name__ == "__main__":
	main()

import numpy as np

from prediction import Prediction


class SingleSeriesPredictor:
	'''
	Defines how to get from a instance to a prediction using a prespecified correlation metric and a time_sampler - only uses information from a single time series
	'''

	def __init__(self, correlation_metric, time_sampler):
		self.correlation_metric = correlation_metric
		self.time_sampler = time_sampler

	def predictSingleSeries(self, instance, time_series_idx):
		'''
		:param instance: an object of type Instance
		:param time_series_idx: which time_series to use for prediction
		:return: a prediction object
		'''
		istsv = instance.getViewForTimeSeries(time_series_idx, self.time_sampler)
		y_pred = self.correlation_metric(istsv)
		return Prediction(istsv.y, y_pred, self)

	def predictAllSeriesCombined(self, instance):
		y_true, y_pred = [], []
		for time_series_idx in range(instance.n_time_series):
			istsv = instance.getViewForTimeSeries(time_series_idx, self.time_sampler)
			y_true.extend(np.asarray(istsv.y).tolist())
			y_pred.extend(np.asarray(self.correlation_metric(istsv).flatten()).tolist())
		return Prediction(y_true, y_pred, self)

	def predictAllInstancesCombined(self, dataset,time_series_idx):
		y_true, y_pred = [], []
		for instance_idx in range(dataset.n_instances):
			instance = dataset.get(instance_idx)
			istsv = instance.getViewForTimeSeries(time_series_idx, self.time_sampler)
			y_true.extend(np.asarray(istsv.y).tolist())
			temp = self.correlation_metric(istsv)
			n,_ = temp.shape
			for i in range(n):
				temp[i,i] = 0.0
			y_pred.extend(np.asarray(temp.flatten()).tolist())
		if len(y_true) == 1:
			y_true = y_true[0]
		else:
			try:
				y_true = [e[0] for e in y_true]
			except:
				pass
		return Prediction(y_true, y_pred, self)


	def __str__(self):
		return str(self.correlation_metric)


class MultipleSeriesPredictor:
	'''
	Defines how to get from a instance to a prediction using a prespecified correlation metric by combining information from multiple time series
	'''

	def __init__(self, correlation_metric, combiner, time_sampler):
		self.correlation_metric = correlation_metric
		self.combiner = combiner
		self.time_sampler = time_sampler

	def predict(self, instance):
		'''
		:param instance: an object of type Instance
		:return: a prediction object
		'''
		y_predictions = []
		for time_series_idx in range(instance.n_time_series):
			istsv = instance.getViewForTimeSeries(time_series_idx, self.time_sampler)
			y_predictions.append(self.correlation_metric(istsv))
		y_pred = self.combiner(y_predictions)
		return Prediction(instance.graph_adjecency_matrix, y_pred, self)

	def __str__(self):
		return str(self.correlation_metric) + str(self.combiner)

from prediction import Prediction


class SingleSeriesPredictor:
	'''
	Defines how to get from a instance to a prediction using a prespecified correlation metric and a time_sampler - only uses information from a single time series
	'''

	def __init__(self, correlation_metric, time_sampler):
		self.correlation_metric = correlation_metric
		self.time_sampler = time_sampler

	def predict(self, instance, time_series_idx):
		'''
		:param instance: an object of type Instance
		:param time_series_idx: which time_series to use for prediction
		:return: a prediction object
		'''
		istsv = instance.getViewForTimeSeries(time_series_idx, self.time_sampler)
		y_pred = self.correlation_metric(istsv)
		return Prediction(istsv.y, y_pred, self)

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

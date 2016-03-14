class Prediction:
	'''
	Holds a prediction for a single graph: the original graph, the predicted graph, the method used to get the prediction
	'''

	def __init__(self, y_true, y_pred, predictor):
		self.y_true = y_true
		self.y_pred = y_pred
		self.predictor = predictor

from sklearn.metrics import roc_curve,auc

class Evaluator:

	def __init__(self):
		pass

	def evaluate(self,prediction):
		'''
		:param prediction: an object of type Prediction
		:return: roc_auc: area under the curve
		'''
		fpr,tpr,_ = roc_curve(prediction.y_true,prediction.y_pred)
		roc_auc = auc(fpr, tpr)
		return roc_auc
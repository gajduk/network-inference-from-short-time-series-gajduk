import numpy as np
import matplotlib.pylab as plt
from methods import good_methods as methods
from methods import mutliple_time_series_combiner
from sklearn.metrics import roc_curve,auc
from src.utils import PROJECT_DIR,s_timestamp,plotDiGraph_,getFeedbackLinks,getForwardLinks
import src.utils
from os.path import join
from matplotlib.backends.backend_pdf import PdfPages
from json_dataset_reader import JsonDatasetReader

debug = False

def normalize_rowwise(x):
	return np.absolute(x)/np.max(np.absolute(x),axis=1,keepdims=True)#/np.std(x,axis=1,keepdims=True)

def evaluateMethodOnInstance(i,method,normalize=False):
	cc = mutliple_time_series_combiner(method,i)

	for idx_node in range(i.n_nodes):
		if normalize:
			if cc[idx_node][idx_node] != 0:
				cc[idx_node] /= np.absolute(cc[idx_node][idx_node])
		cc[idx_node][idx_node] = 0
	if normalize:
		cc = normalize_rowwise(cc)
	y_true = np.reshape(1.0*(np.absolute(i.y)>0.01),(i.n_nodes*i.n_nodes,1))
	cc = np.absolute(cc)
	cc_flat = cc.flatten()
	if debug:
		plt.figure()
		plt.subplot(1,2,1)
		plt.imshow(np.reshape(y_true,(i.n_nodes,i.n_nodes)))
		plt.subplot(1,2,2)
		plt.imshow(np.reshape(cc_flat,(i.n_nodes,i.n_nodes)))
		plt.title(method.__name__)
		plt.show()
	return np.reshape(y_true,(i.n_nodes*i.n_nodes,)),np.reshape(cc_flat,(i.n_nodes*i.n_nodes,))

def evaluateMethod(dataset,method,normalize=False):
	n = dataset.n_instances
	combined_y_true = np.empty((n,dataset.n_nodes*dataset.n_nodes), dtype=np.float64)
	combined_y_pred = np.empty((n,dataset.n_nodes*dataset.n_nodes), dtype=np.float64)

	result = [evaluateMethodOnInstance(dataset.get(idx_instance),method,normalize) for idx_instance in range(n)]

	for i in range(n):
		y_true,y_pred = result[i]
		combined_y_true[i,:] = y_true
		combined_y_pred[i,:] = y_pred
	return combined_y_true,combined_y_pred

def evaluateCombinedTotalRocCurves(predictions,true,methods):
	res = "{0: <20}          {1: <19} {2: <19} {3: <19}\n".format(" ","auc","auc forward","auc feedbacks")
	plt.figure(figsize=(35,12))
	best_auc = 0;
	for f in methods:
		y_true,y_pred = true[f],predictions[f]
		feedbacks_y_true = np.reshape([getFeedbackLinks(temp) for temp in y_true], (-1, 1))
		feedbacks_y_pred = np.reshape([getFeedbackLinks(temp) for temp in y_pred], (-1, 1))
		forward_y_true = np.reshape([getForwardLinks(temp) for temp in y_true], (-1, 1))
		forward_y_pred = np.reshape([getForwardLinks(temp) for temp in y_pred], (-1, 1))
		combined_y_pred = np.reshape(y_pred,(-1,1))
		combined_y_true = np.reshape(y_true,(-1,1))
		plt.subplot(1,3,1)
		roc_auc = plotROC(combined_y_true,combined_y_pred,f)
		if roc_auc > best_auc and roc_auc < 0.99:
			best_auc = roc_auc
		plt.subplot(1,3,2)
		roc_auc_forward = plotROC(forward_y_true,forward_y_pred,f)
		plt.title('ROC for forward only')
		plt.subplot(1,3,3)
		roc_auc_feedbacks = plotROC(feedbacks_y_true,feedbacks_y_pred,f)
		plt.title('ROC for feedbacks only')
		res += "{0: <20} {1:16.3f} {2:16.3f}  {3:16.3f}\n".format(f,roc_auc,roc_auc_forward,roc_auc_feedbacks)
	plt.savefig(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.pdf'))
	with open(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.txt'),'w') as pout:
		pout.write(res)
	return best_auc

def plotROC(y_true,y_pred,label):
	fpr,tpr,_ = roc_curve(y_true,y_pred)
	roc_auc = auc(fpr, tpr)
	plt.plot(fpr, tpr, label=label+' (auc = %0.2f)' % roc_auc)
	plt.legend(loc="lower right")
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	return roc_auc

def plotPredicted(y_pred,label,predict_n,cmap,n_nodes):
	y_pred[np.argsort(y_pred)[:-predict_n]] = 0
	ebunch = [(k,i,y_pred[n_nodes*i+k]) for i in range(n_nodes) for k in range(n_nodes) if y_pred[n_nodes*i+k]!=0]
	plotDiGraph_(n_nodes,ebunch,cmap,500,14)
	plt.title(label)

def evaluateIndividualRocCurvesAndPredictions(d,predictions,true,predict_n,methods):
	with PdfPages(join(PROJECT_DIR,'output','visualization_graph_predictions',s_timestamp()+'.pdf')) as pdf:
		cmap = plt.cm.Accent
		for idx_instance in range(d.n_instances):
			instance = d.get(idx_instance)
			n_nodes = instance.n_nodes
			plt.figure(figsize=(40,20))
			for subplot_idx,f in enumerate(methods):
				y_true = true[f][idx_instance][:]
				y_pred = predictions[f][idx_instance][:]

				#plot the roc curve for the instance
				plt.subplot(4,4,len(predictions)+1)
				plotROC(y_true,y_pred,f)

				#plot the roc curve for the feedbacks only
				plt.subplot(4,4,len(predictions)+2)
				plotROC( getFeedbackLinks(y_true), getFeedbackLinks(y_pred), f)
				plt.title('ROC for feedbacks only')


				#plot the roc curve for the feedbacks only
				plt.subplot(4,4,len(predictions)+3)
				plotROC( getForwardLinks(y_true), getForwardLinks(y_pred), f)
				plt.title('ROC for forward only')

				#plot the predicted networks
				plt.subplot(4,4,subplot_idx+1)
				plotPredicted(y_pred,f,predict_n,cmap,n_nodes)

			plt.subplot(4,4,len(predictions)+4)
			instance.plotTimeSeries_(cmap)

			pdf.savefig()  # saves the current figure into a pdf page
			plt.close()



def evaluateAll(d,normalize=False,predict_n=20,methods=methods):
	predictions = {}
	true = {}
	for f in methods:
		y_true,y_pred = evaluateMethod(d,methods[f],normalize=normalize)
		predictions[f] = y_pred
		true[f] = y_true

	res = evaluateCombinedTotalRocCurves(predictions,true,methods)

	evaluateIndividualRocCurvesAndPredictions(d,predictions,true,predict_n,methods)
	return res




if __name__ == "__main__":
	reader = JsonDatasetReader('initial_12_02_16.json.zip')
	best_aucs = []
	ns_inhibitions = range(1,10)
	for n_inhibitions in ns_inhibitions:
		print n_inhibitions
		src.utils.s_timestamp_prefix = str(n_inhibitions)+"___"
		d = reader.getDataset(n_time_series=n_inhibitions)
		d.plotAll()
		best_aucs.append(evaluateAll(d,predict_n=20))
		print best_aucs
	'''
	plt.close()
	plt.figure()
	plt.plot(ns_inhibitions,best_aucs)
	plt.xlabel('Number of inhibitions')
	plt.ylabel('auc')
	plt.show()
	'''
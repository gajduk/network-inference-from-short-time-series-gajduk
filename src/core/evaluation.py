import numpy as np
import matplotlib.pylab as plt
from methods import methods
from datasets import load,loadEGFR
from sklearn.metrics import roc_curve,auc
from src.utils import PROJECT_DIR,s_timestamp,plotDiGraph_
from os.path import join
from matplotlib.backends.backend_pdf import PdfPages

debug = False

def normalize_rowwise(x):
	return np.absolute(x)/np.max(np.absolute(x),axis=1,keepdims=True)#/np.std(x,axis=1,keepdims=True)

def evaluateMethodOnInstance(i,method,normalize=False):
	cc = method(i)
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


def evaluateAll(d,normalize=True,predict_n=20):
	res = "{0: <20}          {1: <19}\n".format(" ","auc")
	predictions = {}
	for f in methods:
		y_true,y_pred = evaluateMethod(d,methods[f],normalize=normalize)
		predictions[f] = y_pred
		combined_y_pred = np.reshape(y_pred,(-1,1))
		combined_y_true = np.reshape(y_true,(-1,1))
		fpr,tpr,_ = roc_curve(combined_y_true,combined_y_pred)
		roc_auc = auc(fpr, tpr)
		res += "{0: <20} {1:16.3f}\n".format(f,roc_auc)
		plt.plot(fpr, tpr, label=f+' (auc = %0.2f)' % roc_auc)
	plt.legend(loc="lower right")
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.savefig(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.pdf'))
	with open(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.txt'),'w') as pout:
		pout.write(res)
	with PdfPages(join(PROJECT_DIR,'output','visualization_graph_predictions',s_timestamp()+'.pdf')) as pdf:
		cmap = plt.cm.Accent
		for idx_instance in range(d.n_instances):
			n_nodes = d.get(idx_instance).n_nodes
			plt.figure(figsize=(40,20))
			for subplot_idx,f in enumerate(predictions):
				plt.subplot(3,4,subplot_idx+1)
				y_pred = predictions[f][idx_instance][:]
				y_pred[np.argsort(y_pred)[:-predict_n]] = 0
				ebunch = [(k,i,y_pred[n_nodes*i+k]) for i in range(n_nodes) for k in range(n_nodes) if y_pred[n_nodes*i+k]!=0]
				plotDiGraph_(n_nodes,ebunch,cmap,500,14)
				plt.title(f)
			pdf.savefig()  # saves the current figure into a pdf page
			plt.close()



if __name__ == "__main__":
	d = loadEGFR()
	evaluateAll(d,predict_n=8)
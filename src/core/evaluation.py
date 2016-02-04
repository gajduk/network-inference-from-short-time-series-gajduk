import numpy as np
import matplotlib.pylab as plt
from methods import methods
from datasets import load
from sklearn.metrics import roc_curve,auc
from src.utils import PROJECT_DIR,s_timestamp
from os.path import join

debug = False

def normalize_rowwise(x):
	return np.absolute(x)/np.max(np.absolute(x),axis=1,keepdims=True)#/np.std(x,axis=1,keepdims=True)

def evaluateMethod(dataset,method,normalize=False):
	precisions = []
	recalls = []
	accuracies = []
	combined_y_true = []
	combined_y_pred = []
	for idx_instance in range(dataset.n_instances):
		i = dataset.get(idx_instance)

		cc = np.zeros((i.n_nodes,i.n_nodes))
		for idx_node1 in range(i.n_nodes):
			x1,_ = i.get(idx_node1)
			for idx_node2 in range(i.n_nodes):
				x2,_ = i.get(idx_node2)
				cc[idx_node1][idx_node2] = method(x1,x2)

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
		y_pred = np.zeros(cc_flat.shape)
		combined_y_true.extend(y_true)
		combined_y_pred.extend(cc_flat)
		y_pred[np.argsort(cc_flat)[-20:]] = 1.0
		y_pred = np.reshape(y_pred,(i.n_nodes*i.n_nodes,1))

		t_p = np.sum(np.logical_and(y_true==1,y_pred==1))
		t_n = np.sum(np.logical_and(y_true==0,y_pred==0))
		f_p = np.sum(np.logical_and(y_true==0,y_pred==1))
		f_n = np.sum(np.logical_and(y_true==1,y_pred==0))

		p = t_p*1.0/(t_p+f_p+1)
		r = t_p*1.0/(t_p+f_n+1)
		a = (t_p+t_n)*1.0/len(y_true)
		precisions.append(p)
		recalls.append(r)
		accuracies.append(a)
		if debug:
			print 'Instance {1} Precision   {0:.3f}'.format(p,idx_instance)
			print 'Instance {1} Recall      {0:.3f}'.format(r,idx_instance)
			print 'Instance {1} Accuracy    {0:.3f}'.format(a,idx_instance)
			plt.figure()
			plt.subplot(1,3,1)
			plt.imshow(cc)
			plt.colorbar()
			plt.subplot(1,3,2)
			plt.imshow(y_true.reshape(i.n_nodes,i.n_nodes))
			plt.colorbar()
			plt.show()
	if debug:
		print 'Precision   {0:.3f}'.format(np.mean(precisions))
		print 'Recall      {0:.3f}'.format(np.mean(recalls))
		print 'Accuracy    {0:.3f}'.format(np.mean(accuracies))
	return np.mean(precisions),np.mean(recalls),np.mean(accuracies),combined_y_true,combined_y_pred


def evaluateAll(d,normalize=True):
	res = "{0: <20}     {1: <19} {2: <19} {3: <19} {4: <19}\n".format(" ","Precision@20","Recall@20","f1@20","auc")
	for f in methods:
		p,r,a,combined_y_true,combined_y_pred = evaluateMethod(d,methods[f],normalize=normalize)
		fpr,tpr,_ = roc_curve(combined_y_true,combined_y_pred)
		roc_auc = auc(fpr, tpr)
		res += "{0: <20} {1:16.3f} {2:16.3f} {3:16.3f} {3:16.3f}\n".format(f,p,r,2*p*r/(p+r),auc)
		plt.plot(fpr, tpr, label=f+' (auc = %0.2f)' % roc_auc)
	plt.legend(loc="lower right")
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.savefig(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.pdf'))
	with open(join(PROJECT_DIR,'output','evaluation',s_timestamp()+'.txt'),'w') as pout:
		pout.write(res)

if __name__ == "__main__":
	d = load()
	evaluateAll(d)
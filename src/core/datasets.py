import numpy as np
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from src.utils import PROJECT_DIR,s_timestamp,plotDiGraph_
from random import shuffle

class Instance:

    def __init__(self,x,y,n_nodes,n_time_points):
        self.x = x
        self.y = y
        self.n_nodes = n_nodes
        self.n_time_points = n_time_points

    def get(self,idx_node):
        return self.x[self.n_time_points*idx_node:self.n_time_points*(idx_node+1)],self.y[self.n_nodes*idx_node:self.n_nodes*(idx_node+1)]

    def getebunch(self):
        return [(k,i,self.y[self.n_nodes*i+k]) for i in range(self.n_nodes) for k in range(self.n_nodes) if self.y[self.n_nodes*i+k]!=0]

    def plotGraph_(self,cmap):
        plotDiGraph_(self.n_nodes,self.getebunch(),cmap)

    def plotTimeSeries_(self,cmap):
        plt.xlabel('Time [min]')
        plt.ylabel('Activity [A.U.]')
        colors = [cmap(i) for i in np.linspace(0, 1, self.n_nodes)]
        t = np.linspace(0,120,12)
        for i, color in enumerate(colors):
            x,y = self.get(i)
            plt.plot(t,x, color=color)

    def plot(self):
        plt.figure(figsize=(9,4))
        cmap = plt.cm.Accent
        plt.subplot(1,2,1)
        self.plotGraph_(cmap)
        plt.subplot(1,2,2)
        self.plotTimeSeries_(cmap)


class Dataset:

    def __init__(self,X,Y):
        self.X = X
        self.Y = Y
        self.n_instances,self.n_features = X.shape
        self.n_nodes = int(np.sqrt(Y.shape[1]))
        self.n_time_points = self.n_features/self.n_nodes

    def get(self,idx_instance):
         return Instance(self.X[idx_instance],self.Y[idx_instance],self.n_nodes,self.n_time_points)

    def plotAll(self):
        with PdfPages(join(PROJECT_DIR,'output','visualization_graph_and_time_series',s_timestamp()+'.pdf')) as pdf:
            for i in range(self.n_instances):
                instance = self.get(i)
                instance.plot()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()

def load(n_instances=-1,filename_prefix='out_'):
    path_prefix = join(PROJECT_DIR,'data',filename_prefix)
    X, Y = np.loadtxt(path_prefix+'X.csv', delimiter=","),np.loadtxt(path_prefix+'Y.csv', delimiter=",")
    if n_instances != -1:
        idxs = range(len(X))
        shuffle(idxs)
        choice = idxs[:n_instances]
        X = X[choice]
        Y = Y[choice]
    return Dataset(X,Y)

def loadEGFR():
    return load(filename_prefix='out_egfr_')

def loadOscilatory(n_instances=-1):
    return load(n_instances=n_instances,filename_prefix='out_oscilatory_')

if __name__ == "__main__":
    d = loadEGFR()
    d.plotAll()
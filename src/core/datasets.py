import numpy as np
from os.path import join
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from src.utils import PROJECT_DIR,s_timestamp

class Instance:

    def __init__(self,x,y,n_nodes,n_time_points):
        self.x = x
        self.y = y
        self.n_nodes = n_nodes
        self.n_time_points = n_time_points

    def get(self,idx_node):
        return self.x[self.n_time_points*idx_node:self.n_time_points*(idx_node+1)],self.y[self.n_nodes*idx_node:self.n_nodes*(idx_node+1)]

    def getnxDiGraph(self):
        G = nx.DiGraph()
        ebunch = [(k,i,self.y[self.n_nodes*i+k]) for i in range(self.n_nodes) for k in range(self.n_nodes) if self.y[self.n_nodes*i+k]!=0]
        G.add_weighted_edges_from(ebunch)
        return G,ebunch

    def plotGraph_(self,cmap):
        G,ebunch = self.getnxDiGraph()
        pos=nx.shell_layout(G)
        nx.draw_networkx_nodes(G,pos,node_color=range(G.number_of_nodes()),cmap=cmap,node_size=400)
        nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w > 0],edge_color='g',width=1)
        nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w < 0],edge_color='r',width=1)
        nx.draw_networkx_labels(G,pos,{i:str(i+1) for i in range(G.number_of_nodes())},font_size=12)
        plt.axis('off')

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
            for i in range(d.n_instances):
                instance = d.get(i)
                instance.plot()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()

def load():
    X, Y = np.loadtxt(join(PROJECT_DIR,'data','out_X.csv'), delimiter=","),np.loadtxt(join(PROJECT_DIR,'data','out_Y.csv'), delimiter=",")
    return Dataset(X,Y)

if __name__ == "__main__":
    d = load()
    d.plotAll()
import os
import datetime
import networkx as nx
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np

PROJECT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')


def s_timestamp():
    return datetime.datetime.now().strftime("%H_%M_%S %d_%m_%y")

r_namespaces = {}

def load_r_file(filename,namespace):
    if namespace not in r_namespaces:
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        if PROJECT_DIR not in filename:
            filename = os.path.join(PROJECT_DIR,'r_src','forAndrej',filename)
        with open(filename,'r') as pout:
            source = pout.read()
        res = SignatureTranslatedAnonymousPackage(source, namespace)
        r_namespaces[namespace] = res
    return r_namespaces[namespace]


def plotDiGraph_(n_nodes,ebunch,cmap,node_size=400,font_size=12):
    G = nx.DiGraph()
    G.add_nodes_from(range(n_nodes))
    G.add_weighted_edges_from(ebunch)
    pos=nx.shell_layout(G)
    nx.draw_networkx_nodes(G,pos,node_color=range(G.number_of_nodes()),cmap=cmap,node_size=node_size)
    nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w > 0],edge_color='g',width=1)
    nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w < 0],edge_color='r',width=1)
    nx.draw_networkx_labels(G,pos,{i:str(i+1) for i in range(G.number_of_nodes())},font_size=font_size)
    plt.axis('off')

def getFeedbackAndForwardLinks(M):
    if M.ndim == 2:
        n,m = M.shape
    else:
        n,m = M.shape[0],1
    if n != m:
        if n == 1 and int(sqrt(m)) == sqrt(m):
            M = np.reshape(M,(sqrt(m),sqrt(m)))
            n = int(sqrt(m))
        elif m == 1 and int(sqrt(n)) == sqrt(n):
            M = np.reshape(M,(sqrt(n),sqrt(n)))
            n = int(sqrt(n))
        else:
           raise ValueError('input matrix must be square')
    return [M[i,k] for i in range(n) for k in range(n) if i-k>1],[M[i,k] for i in range(n) for k in range(n) if i-k<=1]

def getFeedbackLinks(M):
    res,_ = getFeedbackAndForwardLinks(M)
    return res


def getForwardLinks(M):
    _,res = getFeedbackAndForwardLinks(M)
    return res
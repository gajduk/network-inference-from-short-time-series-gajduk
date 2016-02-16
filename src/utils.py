import os
import datetime
import networkx as nx
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import warnings

PROJECT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')

s_timestamp_prefix = ""

def s_timestamp():
	return s_timestamp_prefix+datetime.datetime.now().strftime("%H_%M_%S %d_%m_%y")

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


def plotDiGraph_(n_nodes,ebunch,cmap,node_size=400,font_size=12,node_labels=None):
    G = nx.DiGraph()
    G.add_nodes_from(range(n_nodes))
    G.add_weighted_edges_from(ebunch)
    if node_labels is None:
        node_labels = [str(i+1) for i in range(G.number_of_nodes())]
    pos=nx.shell_layout(G)
    nx.draw_networkx_nodes(G,pos,node_color=range(G.number_of_nodes()),cmap=cmap,node_size=node_size)
    nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w > 0],edge_color='g',width=1)
    nx.draw_networkx_edges(G,pos,edgelist=[(i,k) for i,k,w in ebunch if w < 0],edge_color='r',width=1)
    nx.draw_networkx_labels(G,pos,{i:node_labels[i] for i in range(G.number_of_nodes())},font_size=font_size)
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

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""
    def newFunc(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__,
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc

def getSubplots(n):
    temp = {1:(1,1),2:(1,2),3:(1,3),4:(2,2),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3),10:(3,4),11:(3,4),12:(4,4),13:(4,4),14:(4,4),15:(4,4),16:(4,4)}
    if n not in temp:
        smaller = int(sqrt(n))
        if smaller*smaller == n:
            return smaller,smaller
        else:
            return smaller,smaller+1
    return temp[n]
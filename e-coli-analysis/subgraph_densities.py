import argparse, os ,sys
import networkx as nx

"""
Command: 
python Graph_check.py --i sample.txt --o sample_converted.txt [--s file_of_communities_whose_union_induced_subgraph_will_be_reported]
"""

def getSubgraph(G, words):
    subg = nx.Graph()
    for (x,y) in G.edges():
        if set([x,y]) <= set(words):
            subg.add_edges_from([(x,y)])
            #subg.add_nodes_from([x,y])
    return subg


parser = argparse.ArgumentParser()
parser.add_argument('--g', type=str, default=None, help='graph file in edgelist format. This graph cannot contain self-loops/multi-edges.')
parser.add_argument('--s', type=str, default="", help='tab separated file containing communities nodes (space separated list, in the last column) whose induced subgraphs are to be computed')
parser.add_argument('--m', type=int, default=3, help='min mumber (must be > 1) of nodes in a community to be considered for subgraph calculation')
#parser.add_argument('--o', type=str, default=None, help='edgelist output file')
args = parser.parse_args()

G = nx.read_edgelist(args.g)

i=1
with open(args.s) as file:
    for line in file:
        cols = line.rstrip().split("\t")
        words = cols[-1].split()
        subgraph = getSubgraph(G, words)
        #whole = G.subgraph(set(words)).copy()
        #take largest connected component, as we are computing connected pseudo-cliques
        ccs = [cc for cc in nx.connected_components(subgraph)]
        if len(ccs) == 0:
            continue
        gcc = max(ccs, key=len)
        sub = getSubgraph(G, gcc)

        #SG.graph.update(G.graph)
        n = len(sub.nodes())
        if n >= args.m:
            m = len(sub.edges())
            den = 2*m/(n*(n-1))
            print('\t'.join(map(str,cols[:-1])), ' '.join(map(str,sorted(sub.nodes()))), n, m, den, sep="\t")
            #nx.write_edgelist(sub, path= cols[0]+'-'+args.g, delimiter=' ', data=False)
            i += 1
#print("Output file is saved in",args.o)

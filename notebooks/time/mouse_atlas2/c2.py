import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
import networkx as nx
Path="/home/mgander/mouse_atlas/data"

Lit=pd.read_csv(f'{Path}/Utils/Curated_transitions.csv', sep=';')
G= nx.DiGraph()
for n in Lit['Cell_type']:
    G.add_node(n)

for i in range(len(Lit)):
    cell_type=Lit['Cell_type'][i]
    desc=Lit['Known_descendants'][i].split(', ')
    for des in desc:
        G.add_edge(cell_type, des)
        G.add_edge(cell_type, cell_type)


def evaluate_using_curated_transitions(df, cutoff=0.05, G=G):
    child_ct=[a.split(':')[1] for a in list(df.columns)]
    parent_ct=[a.split(':')[1] for a in list(df.index)]
    known_edges=list(G.edges)
    known_transition=0
    unknown_transition=0
    M=df.values
    for i,p in zip(range(len(parent_ct)), parent_ct):
        for j,c in zip(range(len(child_ct)),child_ct):
            edge_weight=M[i,j]
            if edge_weight>cutoff:
                if (p, c) in known_edges or c==p:
                    known_transition+=edge_weight
                else:
                    unknown_transition+=edge_weight
    all_transitions=unknown_transition+known_transition
    accuracy=known_transition/all_transitions
    dfr=pd.DataFrame({'Accuracy':[accuracy], 'Total_weight':all_transitions})
    return(dfr)


def evaluate_using_germ_layers(df, cutoff=0.05, Lit=Lit):
    child_ct=[a.split(':')[1] for a in list(df.columns)]
    parent_ct=[a.split(':')[1] for a in list(df.index)]
     # Construct germ layer dict
    Germ_layer_dict=dict(zip(Lit['Cell_type'], Lit['Germ_layer']))

     # Manually add the cell types that have been sub-clusterd
    Germ_layer_dict['Osteoblast progenitors A']='Mesoderm'
    Germ_layer_dict['Osteoblast progenitors B']='Mesoderm'
    Germ_layer_dict['Paraxial mesoderm A']='Mesoderm'

    Germ_layer_dict['Paraxial mesoderm B']='Mesoderm'
    Germ_layer_dict['Paraxial mesoderm C']='Mesoderm'
    Germ_layer_dict['Amniochorionic mesoderm A']='Mesoderm'
    Germ_layer_dict['Amniochorionic mesoderm B']='Mesoderm'

    known_transition=0
    unknown_transition=0
    M=df.values

    for i,p in zip(range(len(parent_ct)), parent_ct):
        for j,c in zip(range(len(child_ct)),child_ct):
            edge_weight=M[i,j]
            if edge_weight>cutoff:

                p_germ=Germ_layer_dict[p]
                c_germ=Germ_layer_dict[c]
                # Don't consider cell types that are not asigned a germ layer
                if p_germ!='Other' and c_germ!='Other':
                # Special case: Neural crest (neuroectoderm) is known to become osteoblasts (mesoderm) of the head
                    if p=='Neural crest' and c in ['Osteoblast progenitors A', 'Osteoblast progenitors B']:
                        known_transition+=edge_weight

                    else:
                        if p_germ==c_germ:
                            known_transition+=edge_weight
                        else:
                            unknown_transition+=edge_weight
    all_transitions=unknown_transition+known_transition
    if all_transitions==0:
        accuracy=1
    else:
        accuracy=known_transition/all_transitions
    dfr=pd.DataFrame({'Accuracy':[accuracy], 'All_transitions':all_transitions})
    return(dfr)
















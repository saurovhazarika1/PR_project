import sys
import numpy as np
import mdtraj as md
import networkx as nx
import matplotlib.pyplot as plt

a = sys.argv[1]

corr_matrix = np.genfromtxt(a)

threshold = 0.75

G = nx.Graph()

#traj = '2ovm/2ovm_strip.dcd'
#top = '2ovm/2ovm_strip.prmtop'

trajectory = md.load('pcca2-samples_1.pdb')

selection = trajectory.topology.select('name CA')
traj_ca = trajectory.atom_slice(selection)

num_residues = traj_ca.topology.n_residues

print('num_residues=', num_residues)
G.add_nodes_from(range(num_residues))


for i in range(num_residues):
    for j in range(i + 1, num_residues):
        if np.abs(corr_matrix[i, j]) > threshold:
            G.add_edge(i, j, weight=corr_matrix[i, j])

print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")

degree_centrality = nx.degree_centrality(G)

betweenness_centrality = nx.betweenness_centrality(G)

clustering_coefficient = nx.clustering(G)

sorted_degree = sorted(degree_centrality.items(), key=lambda item: item[1], reverse=True)
sorted_degree_2 = sorted(betweenness_centrality.items(), key=lambda item: item[1], reverse=True)
print("Top residues by degree centrality:")
for residue, centrality in sorted_degree[:10]:
    print(f"residue {residue}: {centrality:.3f}")

print("Top residues by betweenness centrality:")
for residue, centrality in sorted_degree_2[:10]:
    print(f"residue {residue}: {centrality:.3f}")

communities = nx.algorithms.community.greedy_modularity_communities(G)

communities = [set(community) for community in communities]

#print(communities)

#print(len(communities))

for i, community in enumerate(communities):
    print(f'Community {i + 1}: {community}')

print(len(communities))


import matplotlib.pyplot 
import seaborn as sns


colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
pos = nx.spring_layout(G)
plt.figure(figsize=(12, 12))
for i, community in enumerate(communities):
    nx.draw_networkx_nodes(G, pos, nodelist=community, node_color=colors[i % len(colors)], node_size=100, alpha=0.8)
nx.draw_networkx_edges(G, pos, alpha = 0.75)
plt.title('Communities in Protein Network')
plt.show()


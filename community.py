import sys
import numpy as np
import mdtraj as md
import networkx as nx
import matplotlib.pyplot as plt

a = sys.argv[1]

corr_matrix = np.genfromtxt(a)

threshold = 0.75

G = nx.Graph()

traj = '2ovm/2ovm_strip.dcd'
top = '2ovm/2ovm_strip.prmtop'

trajectory = md.load(traj, top=top)

selection = trajectory.topology.select('name CA')
traj_ca = trajectory.atom_slice(selection)

num_residues = traj_ca.topology.n_residues
G.add_nodes_from(range(num_residues))


for i in range(num_residues):
    for j in range(i + 1, num_residues):
        if np.abs(corr_matrix[i, j]) > threshold:
            G.add_edge(i, j, weight=corr_matrix[i, j])

print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")

from community import community_louvain

partition = community_louvain.best_partition(G)

community_map = {residue.resid: partition[residue.index] for residue in trajectory.topology.residues}
bfactors = np.zeros(len(trajectory.topology.atoms))
for residue in trajectory.topology.residues:
    for atom in residue.atoms:
        bfactors[atom.index] = community_map[residue.resid]

trajectory.topology.atoms.tempfactors = bfactors

trajectory.topology.atoms.write("protein_with_communities.pdb")


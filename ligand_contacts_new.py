import MDAnalysis as mda
from glob import glob

#u = mda.Universe('topology.pdb', 'trajectory.xtc')

traj = glob('2ovm/2ovm_strip.dcd')
#print(trajectory_files)
topology_file = '2ovm/2ovm_strip.prmtop'

u = mda.Universe(topology_file, traj)

ligand = u.select_atoms('resname ASO and not name H*')  # Replace 'LIG' with your ligand residue name

protein = u.select_atoms('protein and not name H*')

from MDAnalysis.lib.distances import capped_distance
import numpy as np

n_frames = len(u.trajectory)
contact_residues = np.zeros((n_frames, len(protein.residues)), dtype=bool)

for ts in u.trajectory:
    ligand_coords = ligand.positions
    protein_coords = protein.positions

    distances, pairs = capped_distance(ligand_coords, protein_coords, 5.0, box=u.dimensions)

    for _, protein_atom_idx in pairs:
        res_idx = protein[protein_atom_idx].residue.index
        contact_residues[ts.frame, res_idx] = True


contact_frequencies = np.sum(contact_residues, axis=0) / n_frames

frequent_contacts = np.where(contact_frequencies >= 0.75)[0]
contact_residues_list = [protein.residues[i] for i in frequent_contacts]

np.savetxt('residues_ligand_contacts.txt', contact_residues_list)

print("Ligand-contacting residues (75% of time):")
for res in contact_residues_list:
    print(f'{res.resname} {res.resid}')


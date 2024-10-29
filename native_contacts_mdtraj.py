import numpy as np
import mdtraj as md
from itertools import combinations
from glob import glob

trajectory_file = glob('4a2j/4a2j_strip_*.dcd')
print(trajectory_file)
topology_file = '4a2j/4a2j_strip.prmtop'

frac_native_contacts = np.array([])
print(frac_native_contacts.shape)

def native_contacts(traj, native):
    BETA_CONST = 50  
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.60

    heavy = native.topology.select_atom_indices('heavy')
    
    heavy_pairs = np.array([(i,j) for (i,j) in combinations(heavy,2) if abs(native.topology.atom(i).residue.index - native.topology.atom(j).residue.index) > 3])
    heavy_pairs_distance = md.compute_distances(native, heavy_pairs)[0]
    native_contacts = heavy_pairs[heavy_pairs_distance < NATIVE_CUTOFF]

    print('number of native contacts = ', native_contacts)
    
    #heavy_pairs_distance_curr = md.compute_distances(traj, heavy_pairs)[0]
    #current_contacts = heavy_pairs[heavy_pairs_distance_curr < NATIVE_CUTOFF]

    #native_contacts_set = set(native_contacts)
    #current_contacts_set = set(current_contacts)
    #common_contacts = native_contacts_set.intersection(current_contacts_set)
   
    #q = len(common_contacts)/len(native_contacts)

    r = md.compute_distances(traj, native_contacts)

    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)

    return q

native = md.load_pdb('4a2j_ref.pdb')
native_sel_1 = native.topology.select('resid 200 to 240')
native_sel_2 = native.topology.select('resid 28 to 73')
native_sel = np.concatenate((native_sel_1, native_sel_2))
native_traj = native.atom_slice(native_sel)

for i in trajectory_file:
    traj = md.load(i, top=topology_file)
    traj_sel_1 = traj.topology.select('resid 200 to 240')
    traj_sel_2 = traj.topology.select('resid 28 to 73')
    traj_sel = np.concatenate((traj_sel_1, traj_sel_2))
    sub_traj = traj.atom_slice(traj_sel)
    q = np.array(native_contacts(sub_traj, native_traj))
    print(q)
    frac_native_contacts = np.concatenate((frac_native_contacts, q))

print(len(frac_native_contacts))

np.savetxt('Q_frac_4a2j.dat', frac_native_contacts)

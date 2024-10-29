import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from glob import glob
import numpy as np

trajectory_files = glob('2ovm/2ovm_strip_1.dcd')
print(trajectory_files)
topology_file = '2ovm/2ovm_strip.prmtop'

#trajs_4a2j = glob('4a2j/4a2j_strip_1.dcd')
#top_4a2j = '4a2j/4a2j_strip.prmtop'

u_ref = mda.Universe('4a2j_ref.pdb')
protein_ref = u_ref.select_atoms("protein")

#protein = u.select_atoms("protein")

native_contacts = contacts.Contacts(u_ref, select=(protein_ref, protein_ref), refgroup=(protein_ref, protein_ref), radius=4.0)
#native_contacts.run()  

print(native_contacts)

native_contact_pairs = native_contacts.results.contacts[0]

q_over_time = []

for traj in trajectory_files:
        u = mda.Universe(topology_file, traj)
        for ts in u.trajectory:
            protein = u.select_atoms("protein")
            frame_contacts = contacts.Contacts(u, select=(protein, protein), refgroup=(protein_ref, protein_ref), radius=6.0)
            frame_contacts.run(start=ts.frame, stop=ts.frame+1)
    
            current_contacts = frame_contacts.results.contact_pairs[0]
    
            q = len(set(current_contacts).intersection(set(native_contact_pairs))) / len(native_contact_pairs)
            q_over_time.append(q)

import matplotlib.pyplot as plt
plt.plot(u.trajectory.times, q_over_time)
plt.xlabel('Time (ps)')
plt.ylabel('Fraction of Native Contacts (Q)')
plt.title('Fraction of Native Contacts During Folding')
plt.show()


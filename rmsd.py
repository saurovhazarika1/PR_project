import numpy as np
import sys
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib.ticker as tick
from matplotlib.axis import Axis


traj_anta = sys.argv[1]
top_anta = sys.argv[2]

traj_ago = sys.argv[3]
top_ago = sys.argv[4]

traj_anta = md.load(traj_anta, top=top_anta)
traj_ago = md.load(traj_ago, top=top_ago)

protein_sel_anta = traj_anta.topology.select('name CA and resid 223 to 241')
protein_traj_anta = traj_anta.atom_slice(protein_sel_anta)
rmsd_anta = md.rmsd(protein_traj_anta, protein_traj_anta, 0)

rmsd_anta = np.array(rmsd_anta)
rmsd_anta = rmsd_anta * 10

protein_sel_ago = traj_ago.topology.select('name CA and resid 223 to 241')
protein_traj_ago = traj_ago.atom_slice(protein_sel_ago)
rmsd_ago = md.rmsd(protein_traj_ago, protein_traj_ago, 0)

rmsd_ago = np.array(rmsd_ago)
rmsd_ago = rmsd_ago * 10

np.savetxt('RMSD_anta.txt', rmsd_anta)
np.savetxt('RMSD_ago.txt', rmsd_ago)

fig, ax = plt.subplots()

plt.plot(rmsd_anta, label='Antagonist', color='brown')
plt.plot(rmsd_ago, label='Agonist', color='darkgreen')
ax.xaxis.set_minor_locator(tick.AutoMinorLocator())
#plt.title('rmsd of 2ovm-H12')
plt.legend()
plt.show()

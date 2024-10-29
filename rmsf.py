import numpy as np
import sys
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib.ticker as tick 
from matplotlib.axis import Axis

traj_anta = '2ovm/2ovm_strip.dcd'
top_anta = '2ovm/2ovm_strip.prmtop'

traj_ago = '4a2j/4a2j_strip.dcd'
top_ago = '4a2j/4a2j_strip.prmtop'

traj_anta = md.load(traj_anta, top=top_anta)
protein_sel_anta = traj_anta.topology.select('protein and name CA')
protein_traj_anta = traj_anta.atom_slice(protein_sel_anta)

rmsf_anta = md.rmsf(protein_traj_anta, protein_traj_anta, 0)

rmsf_anta = np.array(rmsf_anta)
rmsf_anta = rmsf_anta * 10

traj_ago = md.load(traj_ago, top=top_ago)
protein_sel_ago = traj_ago.topology.select('protein and name CA')
protein_traj_ago = traj_ago.atom_slice(protein_sel_ago)

rmsf_ago = md.rmsf(protein_traj_ago, protein_traj_ago, 0)

rmsf_ago = np.array(rmsf_ago)
rmsf_ago = rmsf_ago * 10

fig, ax = plt.subplots()

plt.plot(rmsf_anta[2:], label='Antagonist', color='brown')
plt.plot(rmsf_ago[2:], label='Agonist', color='darkgreen')
ax.xaxis.set_minor_locator(tick.AutoMinorLocator())
#plt.title('rmsf_2ovm')
plt.legend()
plt.show()
plt.savefig('rmsf.png')



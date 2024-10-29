import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.ticker as tick
from matplotlib.axis import Axis


dcd_2ovm = '2ovm/2ovm_strip.dcd'
top_2ovm = '2ovm/2ovm_strip.prmtop'

dcd_4a2j = '4a2j/4a2j_strip.dcd'
top_4a2j = '4a2j/4a2j_strip.prmtop'


traj_2ovm = md.load(dcd_2ovm, top=top_2ovm)
traj_4a2j = md.load(dcd_4a2j, top=top_4a2j)

protein_sel_2ovm = traj_2ovm.topology.select('resid 226 to 239')

protein_sel_4a2j = traj_4a2j.topology.select('resid 226 to 239')


traj_2ovm = traj_2ovm.atom_slice(protein_sel_2ovm)

traj_4a2j = traj_4a2j.atom_slice(protein_sel_4a2j)


dssp_2ovm = md.compute_dssp(traj_2ovm[50:], simplified=True)

dssp_4a2j = md.compute_dssp(traj_4a2j[50:], simplified=True)


Helix_2ovm = []
Coil_2ovm = []

Helix_4a2j = []
Coil_4a2j = []


for i in range(len(dssp_2ovm)):
    helix = 0.0
    coil = 0.0
    for j in range(len(dssp_2ovm[i])):
        if dssp_2ovm[i][j]=='H':
            helix += 1
        elif dssp_2ovm[i][j]=='C':
            coil += 1
    Helix_2ovm.append(helix)
    Coil_2ovm.append(coil)

for i in range(len(dssp_4a2j)):
    helix = 0.0
    coil = 0.0
    for j in range(len(dssp_4a2j[i])):
        if dssp_4a2j[i][j]=='H':
            helix += 1
        elif dssp_4a2j[i][j]=='C':
            coil += 1
    Helix_4a2j.append(helix)
    Coil_4a2j.append(coil)

Helix_4a2j = np.array(Helix_4a2j)
Helix_4a2j_mean = np.mean(Helix_4a2j)

Helix_4a2j_mean_arr = np.ones(len(Helix_4a2j))
Helix_4a2j_mean_arr = Helix_4a2j_mean_arr*Helix_4a2j_mean


Helix_2ovm = np.array(Helix_2ovm)
Helix_2ovm_mean = np.mean(Helix_2ovm)

Helix_2ovm_mean_arr = np.ones(len(Helix_2ovm))
Helix_2ovm_mean_arr = Helix_2ovm_mean_arr*Helix_2ovm_mean

no_frames = np.array(list(range(len(Helix_2ovm))))

np.savetxt('DSSP_anta.txt', Helix_2ovm)
np.savetxt('DSSP_ago.txt', Helix_4a2j)

#fig, ax = plt.subplots()

plt.plot(Helix_2ovm, color='brown')
#plt.plot(Helix_4a2j, color='darkgreen')
plt.plot(Helix_2ovm_mean_arr, label='mean')
#plt.plot(Helix_4a2j_mean_arr, label='mean')
#plt.title('DBD Helix content')
#ax.xaxis.set_minor_locator(tick.AutoMinorLocator())
plt.legend()
plt.savefig('dssp_anta.png')


#print(helix_1)
#print(helix_400)
# open file for writing
#f = open("dssp_dbd.txt","w")

# write file
#f.write( str(dssp_dbd) )

# close file
#f.close()


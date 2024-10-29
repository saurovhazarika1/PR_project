import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory and topology

traj = '4a2j/4a2j_strip_22.dcd'
top = '4a2j/4a2j_strip.prmtop'
      
traj = md.load(traj, top=top)

sasa_values = md.shrake_rupley(traj, mode='atom')
average_sasa_per_frame = sasa_values.sum(axis=1)

average_sasa = np.mean(average_sasa_per_frame)
print(f'Average SASA: {average_sasa:.2f} nm^2')

plt.plot(average_sasa_per_frame)
plt.xlabel('Frame')
plt.ylabel('SASA (nm^2)')
plt.title('SASA Over Time')
plt.show()

sasa_per_residue = sasa_values.mean(axis=0)

for i, sasa in enumerate(sasa_per_residue):
    print(f'Residue {i}: SASA = {sasa:.2f} nm^2')


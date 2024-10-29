##This code calculates the phi and psi dihedral angles for each residue in a MD trajectory and then use the data to calculate the pearson correlation with TICA components##
## Written by Saurov##

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

traj_files = glob('2ovm/2ovm_strip_*.dcd')
print(traj_files)
top = '2ovm/2ovm_strip.prmtop'

phi_angle_sin = np.array([[] for i in range(249)])
phi_angle_cos = np.array([[] for i in range(249)])

psi_angle_sin = np.array([[] for i in range(249)])
psi_angle_cos = np.array([[] for i in range(249)])

for i in traj_files:
    traj = md.load(i, top=top)

    phi_indices, phi_angles = md.compute_phi(traj)
    psi_indices, psi_angles = md.compute_psi(traj)
    omega_indices, omega_angles = md.compute_omega(traj)

    phi_angles_deg = np.rad2deg(phi_angles)
    psi_angles_deg = np.rad2deg(psi_angles)
    omega_angles_deg = np.rad2deg(omega_angles)

    phi_angle_deg_sin = np.transpose(np.sin(phi_angles_deg))
    phi_angle_deg_cos = np.transpose(np.cos(phi_angles_deg))

    psi_angle_deg_sin = np.transpose(np.sin(psi_angles_deg))
    psi_angle_deg_cos = np.transpose(np.cos(psi_angles_deg))
    print(phi_angle_deg_sin.shape)

    phi_angle_sin = np.hstack((phi_angle_sin, phi_angle_deg_sin))
    phi_angle_cos = np.hstack((phi_angle_cos, phi_angle_deg_cos))
    
    psi_angle_sin = np.hstack((psi_angle_sin, psi_angle_deg_sin))
    psi_angle_cos = np.hstack((psi_angle_cos, psi_angle_deg_cos))

    #phi_angle_sin.append(phi_angle_deg_sin)
    #phi_angle_cos.append(phi_angle_deg_cos)

    #psi_angle_sin.append(psi_angle_deg_sin)
    #psi_angle_sin.append(psi_angle_deg_cos)

print(phi_angle_sin.shape)
#print(np.array(phi_angle_sin)[0].shape)




np.savetxt('phi_angle_deg_sin.dat',np.array(phi_angle_sin))
np.savetxt('phi_angle_deg_cos.dat',np.array(phi_angle_cos))

np.savetxt('psi_angle_deg_sin.dat',np.array(psi_angle_sin))
np.savetxt('psi_angle_deg_cos.dat',np.array(psi_angle_cos))

#plt.plot(phi_angles_deg[:, 0], label="Phi (ϕ) angles")
#plt.plot(psi_angles_deg[:, 0], label="Psi (ψ) angles")
#plt.plot(omega_angles_deg[:, 0], label="Omega (ω) angles")
#plt.xlabel("Frame")
#plt.ylabel("Angle (degrees)")
#plt.legend()
#plt.title("Phi, Psi, and Omega Torsional Angles Over Time")
#plt.show()


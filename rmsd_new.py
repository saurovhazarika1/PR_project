import numpy as np
import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.ticker as tick
from matplotlib.axis import Axis


dcd_1 = '2ovm/2ovm_strip.dcd'
psf_1 = '2ovm/2ovm_strip.prmtop'

dcd_ref_1 = '4a2j/4a2j_strip.dcd'
psf_ref_1 = '4a2j/4a2j_strip.prmtop'

#psf_1 = sys.argv[1]
#dcd_1 = sys.argv[2]

#psf_ref_1 = sys.argv[7]
#dcd_ref_1 = sys.argv[8]

u_1 = mda.Universe(psf_1, dcd_1)
ref_1 = mda.Universe(psf_ref_1, dcd_ref_1)

H12 = 'backbone and resid 226-239'
H11 = 'backbone and resid 181-214'

R_1 = rms.RMSD(u_1, ref_1, select='backbone', weights='mass', groupselections=[H12], ref_frame=0)

R2_1 = rms.RMSD(ref_1, ref_1, select='backbone', weights='mass', groupselections=[H12], ref_frame=0)

R_2 = rms.RMSD(u_1, u_1, select='backbone', weights='mass', groupselections=[H12], ref_frame=0)

R2_2 = rms.RMSD(ref_1, u_1, select='backbone', weights='mass', groupselections=[H12], ref_frame=0)


R_1.run()

R2_1.run()

R_2.run()

R2_2.run()

#print(R_2ovm.rmsd.shape)
#print(R2.rmsd.shape)


np.savetxt('RMSD_anta_refAgo_H12.txt', R_1.rmsd[:,3])

np.savetxt('RMSD_ago_refAgo_H12.txt', R2_1.rmsd[:,3])

np.savetxt('RMSD_anta_refAnta_H12.txt', R_2.rmsd[:,3])

np.savetxt('RMSD_ago_refAnta_H12.txt', R2_2.rmsd[:,3])


#fig, ax = plt.subplots()

#plt.plot(R_1.rmsd[:,3], label='antagonistic conformation', color='brown')
#plt.plot(R_2.rmsd[:,3], color='brown')
#plt.plot(R_3.rmsd[:,3], color='brown')
#plt.plot(R2_1.rmsd[:,3], label='agonistic conformation', color='darkgreen')
#plt.plot(R2_2.rmsd[:,3], color='darkgreen')
#plt.plot(R2_3.rmsd[:,3], color='darkgreen')
#ax.xaxis.set_minor_locator(tick.AutoMinorLocator())
#plt.ylabel('rmsd')
#plt.legend()
#plt.show()

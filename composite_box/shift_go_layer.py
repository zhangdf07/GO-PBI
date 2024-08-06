import mdtraj as md
import numpy as np
import os
import sys


fin = 'layered-5gof.gro'

traj = md.load(fin, top=fin)
topo = traj.topology

bbox = traj.unitcell_lengths[0]

natom = topo.n_atoms
nresd = topo.n_residues
nfram = len(traj)

#print( traj.xyz[:] )
shift_vec = bbox*np.array([0.2,0.2,0])

for resid in range(nresd):
    go = [ at.index for at in topo.atoms if at.residue.index==resid]
    traj.xyz[0,go,:] = traj.xyz[0,go,:] +shift_vec*resid

traj.save('layered-5gof-shifted.gro')


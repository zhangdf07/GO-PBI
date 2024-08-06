import mdtraj as md
import numpy as np
import os
import sys


def know_mdtraj(fin,top):
    traj = md.load(fin, top=top)
    topo = traj.topology

    bbox = traj.unitcell_lengths[0]

    natom = topo.n_atoms
    nresd = topo.n_residues
    nfram = len(traj)

    print(set([ a.residue.name for a in topo.atoms ]))
    #zc = [ traj.xyz[0,at.index,2] for at in topo.atoms if at.element.symbol == "C" ]
    #print(np.amin(zc),np.mean(zc),np.amax(zc))
    #exit()
    #mass = [at.element.mass for at in topo.atoms if at.residue.name not in ['SOL','BDO'] ]
    #density = sum(mass) / np.prod(bbox)
    #print( density )

    #print( mass[:5])
    #print( sum(mass) )
    #print( (5**3 * 1.3) / ( sum(mass)* 1.66054e-27), " \n")

    #a2k = [at.index for at in topo.atoms if at.residue.name not in ['SOL','BDO'] ]
    #traj.restrict_atoms(a2k)
    #density = md.density(traj)
    #print( density )

    atm = [at.index for at in topo.atoms if at.residue.index == 1 and at.element.symbol == 'C']
    for n in range(nfram):
      zatm  = [ traj.xyz[n,a,2] for a in atm ]#if topo.atom(a).element.symbol == 'C' ]
      print(n, np.mean(zatm))
    #alo = [at.index for at in topo.atoms if traj.xyz[0,at.index,2] < 5.0]


    #shf = -np.mean(zlo)
    #for n in range(natom):
    #    traj.xyz[0,n,2] += shf

    #zhi = [ traj.xyz[0,a,2] for a in ahi ]#if topo.atom(a).element.symbol == 'C']
    #shf = 2* np.mean(zhi)
    #for n in ahi:
    #    traj.xyz[0,n,2] *= -1.0
    #for n in ahi:
    #    traj.xyz[0,n,2] += shf
    #print(np.amin(zlo),np.mean(zlo),np.amax(zlo))
    #print(np.amin(zhi),np.mean(zhi),np.amax(zhi))
    exit()
    #a2m = [at.index for at in topo.atoms if at.residue.name not in ['GRA'] ]

    #shf = 0.3
    #shf = -min([ traj.xyz[0,n,2] for n in range(natom) ])+0.5

    #for n in range(natom):
    #for n in a2m:
    #    traj.xyz[0,n,2] += shf
    fout = 'out.gro'
    traj.save(fout)


#fin = input("gro file name:")
#top = input("top file name:")
fin = sys.argv[1]
top = sys.argv[2]
know_mdtraj(fin,top)
exit()

fins = sorted([ l+'/nvt2.gro' for l in os.listdir('./') if 'job-msd01-' in l ])
for fin in fins:
    print(fin)
    know_mdtraj(fin)


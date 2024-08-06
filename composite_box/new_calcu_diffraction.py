import sys
import numpy as np
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import os
from natsort import natsorted


def get_xrd(fpos, fxrd):
  #print("Start XRD calculation: ", fpos) # From a VASP POSCAR file 
  data = Structure.from_file(fpos)
  xrd = XRDCalculator().get_pattern(data,scaled=False,two_theta_range=(0,80))
  #print( xrd.x, xrd.y )  
  output = np.transpose( [xrd.x, xrd.y] )
  np.savetxt( fxrd, output, delimiter=",")


# All gro files.
jdir = np.array( natsorted([ l for l in os.listdir('./') if '.gro' in l ]) )
#jdir = jdir[::2]
jdir = jdir[-100:]
jdir = jdir[::-1]

# All POSCAR that has been generated.
jposcar = [ l for l in os.listdir('./') if 'POSCAR-' in l ]
# All XRD that has been generated.
jxrd = [ l for l in os.listdir('./') if 'xrd-' in l ]

# compute XRD for every X frames
for i in range(0,len(jdir)):
  fin = jdir[i]
  at1 = read(fin)

  # structure as is
  jposcar = [ l for l in os.listdir('./') if 'POSCAR-' in l ]
  fout = 'POSCAR-'+ fin[:-4]+ '.vasp'
  if fout not in jposcar:
    write( fout, at1)

  jxrd = [ l for l in os.listdir('./') if 'xrd-' in l ]
  fxrd = 'xrd-'+fin[:-4]+'.log'
  if fxrd not in jxrd:
    get_xrd(fout,fxrd)

  # structure after removing BDO and SOL
  jposcar = [ l for l in os.listdir('./') if 'POSCAR-' in l ]
  fout = 'POSCAR-'+ fin[:-4]+ '-removeBDOandSOL.vasp'
  if fout not in jposcar:
    bdo = at1.get_array('residuenames')=='BDO'
    sol = at1.get_array('residuenames')=='SOL'
    mask = ~np.any([bdo,sol], axis=0)
    at2 = at1[mask]
    write( fout, at2)

  jxrd = [ l for l in os.listdir('./') if 'xrd-' in l ]
  fxrd = 'xrd-'+fin[:-4]+'-removeBDOandSOL.log'
  if fxrd not in jxrd:
    get_xrd(fout,fxrd)

exit()


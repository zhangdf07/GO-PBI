#!/bin/bash
  echo "BDO_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-bdo-x -tu ns -n index.ndx -type x
  echo "SOL_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-wat-x -tu ns -n index.ndx -type x
  echo "BDO_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-bdo-yz -tu ns -n index.ndx -lateral x
  echo "SOL_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-wat-yz -tu ns -n index.ndx -lateral x

nvt=md-npt-03
  echo "BDO_NAME_SPC" | gmx_d msd -f $nvt.xtc -s $nvt.tpr -o msd-$nvt-bdo-y -tu ns -n index.ndx -type y
  echo "SOL_NAME_SPC" | gmx_d msd -f $nvt.xtc -s $nvt.tpr -o msd-$nvt-wat-y -tu ns -n index.ndx -type y
  echo "BDO_NAME_SPC" | gmx_d msd -f $nvt.xtc -s $nvt.tpr -o msd-$nvt-bdo-z -tu ns -n index.ndx -type z
  echo "SOL_NAME_SPC" | gmx_d msd -f $nvt.xtc -s $nvt.tpr -o msd-$nvt-wat-z -tu ns -n index.ndx -type z



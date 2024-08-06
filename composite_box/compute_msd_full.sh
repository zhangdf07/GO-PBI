#!/bin/bash
  echo "BDO_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-bdo -tu ns -n index.ndx
  echo "SOL_NAME_SPC" | gmx_d msd -f md-npt-03.xtc -s md-npt-03.tpr -o msd-md-npt-03-wat -tu ns -n index.ndx

  echo "BDO_NAME_SPC" | gmx_d msd -f md-npt-02.xtc -s md-npt-02.tpr -o msd-md-npt-02-bdo -tu ns -n index.ndx
  echo "SOL_NAME_SPC" | gmx_d msd -f md-npt-02.xtc -s md-npt-02.tpr -o msd-md-npt-02-wat -tu ns -n index.ndx

  echo "BDO_NAME_SPC" | gmx_d msd -f md-npt-01.xtc -s md-npt-01.tpr -o msd-md-npt-01-bdo -tu ns -n index.ndx
  echo "SOL_NAME_SPC" | gmx_d msd -f md-npt-01.xtc -s md-npt-01.tpr -o msd-md-npt-01-wat -tu ns -n index.ndx


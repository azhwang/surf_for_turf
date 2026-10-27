mol new solute_raw.psf
mol addfile min.coor waitfor all
[atomselect top all] writepdb solute.pdb
quit

package require solvate
package require autoionize

solvate complex.psf complex.pdb -t 10 -o solv
autoionize -psf solv.psf -pdb solv.pdb -neutralize -salt 0.15

quit

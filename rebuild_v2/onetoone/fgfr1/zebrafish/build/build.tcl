
package require psfgen
package require solvate
package require autoionize
resetpsf
readpsf /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/fgfr1/zebrafish/build/solute.psf
coordpdb /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/fgfr1/zebrafish/build/solute_oriented.pdb
writepsf pair.psf
writepdb pair.pdb
set h [expr 125.0/2.0]
solvate pair.psf pair.pdb -minmax [list [list -$h -$h -$h] [list $h $h $h]] -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation POT -anion CLA -seg IONS -o ionized
quit

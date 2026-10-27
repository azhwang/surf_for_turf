
package require psfgen
package require solvate
package require autoionize
resetpsf
topology /Users/alexanderwang/Desktop/Surf_for_Turf/all_metrics_analysis/1fq9_goldfish/charmm-gui-8035857277/toppar/top_all36_prot.rtf
topology /Users/alexanderwang/Desktop/Surf_for_Turf/all_metrics_analysis/1fq9_goldfish/charmm-gui-8035857277/toppar/top_all36_carb.rtf
readpsf /Users/alexanderwang/Desktop/Surf_for_Turf/all_metrics_analysis/1fq9_goldfish/charmm-gui-8035857277/step1_pdbreader.psf
coordpdb /Users/alexanderwang/Desktop/Surf_for_Turf/all_metrics_analysis/1fq9_goldfish/charmm-gui-8035857277/step1_pdbreader.pdb
delatom PROC
delatom PROD
pdbalias residue HIS HSD
foreach seg {PROC PROD} {
    segment $seg { first NTER; last CTER; pdb fgf_$seg.pdb }
    coordpdb fgf_$seg.pdb $seg
}
guesscoord
writepsf solute.psf
writepdb solute.pdb
set m [mol new solute.psf]
mol addfile solute.pdb waitfor all
set all [atomselect $m all]
$all moveby [vecscale -1.0 [measure center $all]]
set mm [measure minmax $all]
puts "SOLUTE_MINMAX $mm"
$all writepdb solute_centered.pdb
mol delete $m
set h [expr 124.0/2.0]
solvate solute.psf solute_centered.pdb -minmax [list [list -$h -$h -$h] [list $h $h $h]] -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation POT -anion CLA -seg IONS -o ionized
quit

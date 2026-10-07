mol new ./1djs_human_model_0.pdb

# Remove waters, ions, altlocs
set sel [atomselect top "protein or resname LIG"]
$sel writepdb clean_complex.pdb

# Remove waters
mol new clean_complex.pdb
set sel [atomselect top "not water"]
$sel writepdb tmp.pdb

mol new tmp.pdb
set sel [atomselect top all]
foreach ch [lsort -unique [$sel get chain]] {
    set s [atomselect top "chain $ch"]
    $s set segid $ch
}

set his [atomselect top "resname HIS"]
$his set resname HSE

set all [atomselect top all]
$all writepdb prefix.pdb

# Save new PDB
set all [atomselect top "all"]
$all writepdb clean.pdb

package require psfgen

# Load topologies
topology top_all36_prot.rtf
topology top_all36_water_ions.rtf

resetpsf

mol new clean.pdb

foreach ch {A B} {
    set sel [atomselect top "chain $ch"]
    $sel writepdb chain_$ch.pdb
}

segment A {
    pdb chain_A.pdb
}
segment B {
    pdb chain_B.pdb
}

coordpdb chain_A.pdb A
coordpdb chain_B.pdb B

guesscoord

# Write output
writepsf complex.psf
writepdb complex.pdb

package require solvate
package require autoionize

solvate complex.psf complex.pdb -t 10 -o solv
autoionize -psf solv.psf -pdb solv.pdb -neutralize -salt 0.15

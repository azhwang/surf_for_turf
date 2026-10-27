
package require psfgen
resetpsf
topology /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/1fq9/goldfish/sim/toppar_charmm/top_all36_prot.rtf
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
segment PROA { first NTER; last CTER; pdb /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/scaffold/receptor_fgfr2.pdb }
coordpdb /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/scaffold/receptor_fgfr2.pdb PROA
segment PROC { first NTER; last CTER; pdb /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/fgfr2/human_stab/fgf_PROC.pdb }
coordpdb /Users/alexanderwang/Desktop/Surf_for_Turf/rebuild_v2/onetoone/fgfr2/human_stab/fgf_PROC.pdb PROC
patch DISU PROA:179 PROA:231
patch DISU PROA:278 PROA:342
regenerate angles dihedrals
guesscoord
writepsf solute_raw.psf
writepdb solute_raw.pdb
quit

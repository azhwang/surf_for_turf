#!/usr/bin/env bash
#
# analyze_exploratory.sh <replicate dir> -- EXPLORATORY analyses of one 50 ns extension
# replicate (NOT part of ANALYSIS_PLAN_extensions.md; added 2026-10-15 after the planned
# analysis was reported). Same window (ext11..extNN) and GB settings as analyze_extension.sh.
#
# Writes <rep>/ext_analysis/:
#   hep_run/    MM/GBSA, FGFR stripped: receptor = heparin(s), ligand = both FGFs
#               (FGF-heparin energy), frames every 0.5 ns
#   fgfr_run/   MM/GBSA, heparin stripped: receptor = both FGFRs, ligand = both FGFs
#               (FGF-FGFR energy without heparin), frames every 0.5 ns
#   decomp_run/ whole-assembly split with per-residue decomposition (idecomp=1), every 1 ns
#
# Needs analyze_extension.sh to have run first (complex_chamber.parm7, mmpbsa_v2/whole_*.parm7).
set -euo pipefail
REP="$1"
cd "$REP"
OUT=ext_analysis
mkdir -p "$OUT/explo"
segs=( $(ls ext[0-9][0-9].dcd | sort) ); nseg=${#segs[@]}
win=( "${segs[@]:10}" )
[ "${#win[@]}" -gt 0 ] || { echo "fewer than 11 segments" >&2; exit 1; }

# component ranges from the PSF: order is always FGFR (PROA,PROB), heparin (CAR*), FGF (PROC,PROD)
eval "$(python3 - <<'EOF'
res=[]; last=None
L=open('step3_input.psf').read().split('\n'); i=next(k for k,l in enumerate(L) if '!NATOM' in l)
for l in L[i+1:i+1+int(L[i].split()[0])]:
    p=l.split(); k=(p[1],p[2])
    if k!=last: res.append(p[1]); last=k
seg={}
for n,s in enumerate(res,1): seg.setdefault(s,[n,n]); seg[s][1]=n
car=sorted(s for s in seg if s.startswith('CAR'))
r0,r1=seg['PROA'][0],seg['PROB'][1]; h0,h1=seg[car[0]][0],seg[car[-1]][1]; f0,f1=seg['PROC'][0],seg['PROD'][1]
assert r1+1==h0 and h1+1==f0, (r0,r1,h0,h1,f0,f1)
nR,nH,nF=r1-r0+1,h1-h0+1,f1-f0+1
print(f'FGFR=":{r0}-{r1}"; HEP=":{h0}-{h1}"; NR={nR}; NH={nH}; NF={nF}')
print(f'HEP_LIG=":{nH+1}-{nH+nF}"; FGFR_LIG=":{nR+1}-{nR+nF}"')
EOF
)"
STRIP=":CLA,POT,WAT"

cat > $OUT/explo/gb.in <<EOF
&general
  startframe=1, endframe=999999, interval=1, verbose=2,
/
&gb
  igb=5, saltcon=0.150,
/
EOF
cat > $OUT/explo/decomp.in <<EOF
&general
  startframe=1, endframe=999999, interval=1, verbose=1,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=1, dec_verbose=0, print_res="all", csv_format=1,
/
EOF

run_split () {   # name other_mask lig_mask(post-strip) stride input
  local s=$1 other=$2 lig=$3 stride=$4 inp=$5
  [ -f $OUT/${s}_run/FINAL_RESULTS_${s}.dat ] && [ "$(cat $OUT/${s}_run/nseg 2>/dev/null)" = "$nseg" ] && return 0
  mkdir -p $OUT/${s}_run
  local strip="$STRIP"; [ -n "$other" ] && strip="$STRIP,${other#:}"
  if [ -n "$other" ]; then
    ( cd $OUT/${s}_run && rm -f ${s}_complex.parm7 ${s}_receptor.parm7 ${s}_ligand.parm7 && \
      ante-MMPBSA.py -p ../../complex_chamber.parm7 -c ${s}_complex.parm7 -r ${s}_receptor.parm7 \
        -l ${s}_ligand.parm7 -s "$strip" -n "$lig" --radii mbondi2 > ante.log 2>&1 )
    cp_=${s}_complex.parm7; rp_=${s}_receptor.parm7; lp_=${s}_ligand.parm7
  else
    cp_=../../mmpbsa_v2/whole_complex.parm7; rp_=../../mmpbsa_v2/whole_receptor.parm7; lp_=../../mmpbsa_v2/whole_ligand.parm7
  fi
  cpptraj -p complex_chamber.parm7 > $OUT/${s}_run/cpptraj.log <<EOF
$(printf "trajin %s $stride\n" "${win[@]}")
autoimage
strip $strip
trajout $OUT/${s}_run/dry.nc netcdf
run
EOF
  local df=FINAL_DECOMP_${s}.dat
  ( cd $OUT/${s}_run && MMPBSA.py -O -i ../explo/$inp -o FINAL_RESULTS_${s}.dat -do $df \
      -cp $cp_ -rp $rp_ -lp $lp_ -y dry.nc -eo FRAMES_${s}.csv > mmpbsa.log 2>&1 && echo "$nseg" > nseg )
  rm -f $OUT/${s}_run/dry.nc
}

echo "FGFR $FGFR HEP $HEP (nR=$NR nH=$NH nF=$NF); hep ligand $HEP_LIG; fgfr ligand $FGFR_LIG" > $OUT/explo/info.txt
run_split hep  "$FGFR" "$HEP_LIG"  "1 last 5"  gb.in
run_split fgfr "$HEP"  "$FGFR_LIG" "1 last 5"  gb.in
run_split decomp ""    ""          "1 last 10" decomp.in
echo "done: $REP"

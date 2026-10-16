#!/usr/bin/env bash
#
# analyze_r90.sh <replicate dir> -- Arg90 heparin engagement for one extension replicate, as fixed in
# ANALYSIS_PLAN_chimeras.md (the same measurement as ext_results/explo_r90_contacts.txt).
#
# Arg90 = FGF2 155-aa numbering = core position 63 (1E0O) / 66 (1FQ9). For each FGF copy, heavy-atom
# contacts (4.5 A) of Arg90 with any heparin residue and with each FGFR chain, every 0.5 ns over
# extension segments 11-50; window means go to <rep>/ext_analysis/r90.txt.
# Requires AmberTools on PATH and <rep>/complex_chamber.parm7 (from the 1 ns MM/GBSA pipeline).
set -euo pipefail
REP="$1"
cd "$REP"
mkdir -p ext_analysis
segs=( $(ls ext[0-9][0-9].dcd | sort) )
[ "${#segs[@]}" -ge 11 ] || { echo "fewer than 11 segments in $REP" >&2; exit 1; }
win=( "${segs[@]:10}" )

eval "$(python3 - <<'EOF'
res=[]; last=None
L=open('step3_input.psf').read().split('\n'); i=next(k for k,l in enumerate(L) if '!NATOM' in l)
for l in L[i+1:i+1+int(L[i].split()[0])]:
    p=l.split(); k=(p[1],p[2])
    if k!=last: res.append((p[1],p[3])); last=k
seg={}
for n,(s,_) in enumerate(res,1): seg.setdefault(s,[n,n]); seg[s][1]=n
core=seg['PROC'][1]-seg['PROC'][0]+1
pos={127:63,129:66}[core]
rC,rD=seg['PROC'][0]+pos-1, seg['PROD'][0]+pos-1
assert res[rC-1][1]=='ARG' and res[rD-1][1]=='ARG', (res[rC-1],res[rD-1])
car=sorted(s for s in seg if s.startswith('CAR'))
print(f'R_C={rC}; R_D={rD}; HEP=":{seg[car[0]][0]}-{seg[car[-1]][1]}"')
print(f'FA=":{seg["PROA"][0]}-{seg["PROA"][1]}"; FB=":{seg["PROB"][0]}-{seg["PROB"][1]}"')
EOF
)"
out=ext_analysis/r90_contacts.dat
cpptraj -p complex_chamber.parm7 > ext_analysis/cpptraj_r90.log <<EOF
$(printf 'trajin %s 1 last 5\n' "${win[@]}")
autoimage
nativecontacts name cH :${R_C}&!@H= ${HEP}&!@H= distance 4.5 out $out
nativecontacts name cA :${R_C}&!@H= ${FA}&!@H= distance 4.5 out $out
nativecontacts name cB :${R_C}&!@H= ${FB}&!@H= distance 4.5 out $out
nativecontacts name dH :${R_D}&!@H= ${HEP}&!@H= distance 4.5 out $out
nativecontacts name dA :${R_D}&!@H= ${FA}&!@H= distance 4.5 out $out
nativecontacts name dB :${R_D}&!@H= ${FB}&!@H= distance 4.5 out $out
run
EOF
# nativecontacts writes native + non-native columns per name; total = sum of the pair
awk -v rep="$REP" -v rc="$R_C" -v rd="$R_D" '!/^#/{for(i=2;i<=NF;i++)s[i]+=$i; n++}
  END{printf "%s  Arg90 residues %s/%s  frames %d\nFGF C: heparin %.1f  FGFR-A %.1f  FGFR-B %.1f\nFGF D: heparin %.1f  FGFR-A %.1f  FGFR-B %.1f\n",
      rep, rc, rd, n, (s[2]+s[3])/n, (s[4]+s[5])/n, (s[6]+s[7])/n, (s[8]+s[9])/n, (s[10]+s[11])/n, (s[12]+s[13])/n}' $out \
  > ext_analysis/r90.txt
cat ext_analysis/r90.txt

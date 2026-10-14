#!/usr/bin/env bash
#
# analyze_extension.sh <replicate dir> -- analysis of one 50 ns extension replicate,
# exactly as fixed in ANALYSIS_PLAN_extensions.md.
#
# Writes <rep>/ext_analysis/:
#   timeseries.dat   every frame (100 ps) of ext01..extNN: heparin-site salt bridges for
#                    FGF C and FGF D (canonical 5 K/R and all K/R), FGFR-FGFR COM distance
#   <split>_run/FINAL_RESULTS_<split>.dat   MM/GBSA (iface1, iface2, whole) on frames every
#                    0.5 ns of ext11..extNN, no decomposition; FRAMES_<split>.csv = per-frame
#                    energies (for the 10 ns block averages in the plan)
#   info.txt         masks used and number of segments
#
# Needs the 1 ns v2 MM/GBSA outputs in <rep> (complex_chamber.parm7, mmpbsa_v2/split_vars.sh,
# mmpbsa_v2/*_{complex,receptor,ligand}.parm7). Requires AmberTools on PATH.
set -euo pipefail
REP="$1"
cd "$REP"
OUT=ext_analysis
mkdir -p "$OUT"

segs=( $(ls ext[0-9][0-9].dcd 2>/dev/null | sort) )
nseg=${#segs[@]}
[ "$nseg" -gt 0 ] || { echo "no ext*.dcd in $REP" >&2; exit 1; }

# masks from this replicate's own PSF: canonical heparin-site K/R (core positions differ by
# template), each FGF copy's heparin, FGFR chains
eval "$(python3 - <<'EOF'
res=[]; last=None
L=open('step3_input.psf').read().split('\n'); i=next(k for k,l in enumerate(L) if '!NATOM' in l)
for l in L[i+1:i+1+int(L[i].split()[0])]:
    p=l.split(); k=(p[1],p[2])
    if k!=last: res.append((p[1],p[3])); last=k
seg={}
for n,(s,_) in enumerate(res,1): seg.setdefault(s,[n,n]); seg[s][1]=n
car=sorted(s for s in seg if s.startswith('CAR'))
core_len=seg['PROC'][1]-seg['PROC'][0]+1
site={129:[104,105,110,114,120], 127:[101,102,107,111,117]}[core_len]   # 1FQ9 / 1E0O
def site_mask(s): return ','.join(str(seg[s][0]+p-1) for p in site)
for s in ('PROC','PROD'):
    names=[res[seg[s][0]+p-2][1] for p in site]
    assert all(n in ('LYS','ARG') for n in names), (s,names)
hepC = car[0]; hepD = car[1] if len(car) == 2 else car[0]     # 2:2:2: own heparin; 2:2:1: shared
def rng(s): return f"{seg[s][0]}-{seg[s][1]}"
print(f'SITE_C=":{site_mask("PROC")}"; SITE_D=":{site_mask("PROD")}"')
print(f'FGF_C=":{rng("PROC")}"; FGF_D=":{rng("PROD")}"; HEP_C=":{rng(hepC)}"; HEP_D=":{rng(hepD)}"')
print(f'FGFR_A=":{rng("PROA")}"; FGFR_B=":{rng("PROB")}"')
EOF
)"
source mmpbsa_v2/split_vars.sh
N='@NZ,NE,NH1,NH2'
O='@OS22,OS23,OS24,OS62,OS63,OS64,O61,O62'
{
  echo "segments: $nseg (${segs[0]} .. ${segs[$((nseg-1))]})"
  echo "site C $SITE_C  site D $SITE_D  heparin C $HEP_C  heparin D $HEP_D"
  echo "FGFR $FGFR_A $FGFR_B ; strip $STRIP_MASK ; iface1 other $IFACE1_OTHER_MASK ; iface2 other $IFACE2_OTHER_MASK"
} > $OUT/info.txt

trajin_all=$(printf 'trajin %s\n' "${segs[@]}")

# --- 1. structural time series, every frame ---
if [ ! -f $OUT/timeseries.dat ]; then
cpptraj -p complex_chamber.parm7 > $OUT/cpptraj_ts.log <<EOF
$trajin_all
autoimage
nativecontacts name siteC ${SITE_C}${N} ${HEP_C}${O} distance 4.0 out $OUT/contacts.dat
nativecontacts name siteD ${SITE_D}${N} ${HEP_D}${O} distance 4.0 out $OUT/contacts.dat
nativecontacts name allC ${FGF_C}&:LYS,ARG${N} ${HEP_C}${O} distance 4.0 out $OUT/contacts.dat
nativecontacts name allD ${FGF_D}&:LYS,ARG${N} ${HEP_D}${O} distance 4.0 out $OUT/contacts.dat
distance rdist ${FGFR_A}@CA ${FGFR_B}@CA out $OUT/rdist.dat
run
EOF
python3 - $OUT <<'EOF'
import sys
d=sys.argv[1]
c=[l.split() for l in open(f'{d}/contacts.dat') if not l.startswith('#')]
r=[l.split() for l in open(f'{d}/rdist.dat') if not l.startswith('#')]
with open(f'{d}/timeseries.dat','w') as f:
    f.write('#frame time_ns siteC siteD allC allD fgfr_dist\n')
    for a,b in zip(c,r):
        v=[float(x) for x in a[1:]]
        f.write(f"{int(a[0])} {1+int(a[0])*0.1:.1f} {v[0]+v[1]:.0f} {v[2]+v[3]:.0f} {v[4]+v[5]:.0f} {v[6]+v[7]:.0f} {float(b[1]):.3f}\n")
EOF
fi

# --- 2. MM/GBSA over the window (ext11 onwards), every 0.5 ns ---
win=( "${segs[@]:10}" )
[ "${#win[@]}" -gt 0 ] || { echo "fewer than 11 segments; MM/GBSA skipped"; exit 0; }
trajin_win=$(printf 'trajin %s 1 last 5\n' "${win[@]}")
cat > $OUT/mmpbsa_ext.in <<EOF
&general
  startframe=1, endframe=999999, interval=1, verbose=2,
/
&gb
  igb=5, saltcon=0.150,
/
EOF
splits=(whole); [ "$IFACE_SPLIT" = 1 ] && splits=(iface1 iface2 whole)
for s in "${splits[@]}"; do
  [ -f $OUT/${s}_run/FINAL_RESULTS_${s}.dat ] && [ "$(cat $OUT/${s}_run/nseg 2>/dev/null)" = "$nseg" ] && continue
  case $s in iface1) other="$IFACE1_OTHER_MASK";; iface2) other="$IFACE2_OTHER_MASK";; *) other="";; esac
  strip="$STRIP_MASK"; [ -n "$other" ] && strip="$STRIP_MASK,${other#:}"
  cpptraj -p complex_chamber.parm7 > $OUT/cpptraj_${s}.log <<EOF
$trajin_win
autoimage
strip $strip
trajout $OUT/${s}_dry.nc netcdf
run
EOF
  mkdir -p $OUT/${s}_run
  ( cd $OUT/${s}_run && MMPBSA.py -O -i ../mmpbsa_ext.in -o FINAL_RESULTS_${s}.dat \
      -cp ../../mmpbsa_v2/${s}_complex.parm7 -rp ../../mmpbsa_v2/${s}_receptor.parm7 \
      -lp ../../mmpbsa_v2/${s}_ligand.parm7 -y ../${s}_dry.nc -eo FRAMES_${s}.csv > mmpbsa.log 2>&1 \
      && echo "$nseg" > nseg )
  rm -f $OUT/${s}_dry.nc
done
echo "done: $REP ($nseg segments)"

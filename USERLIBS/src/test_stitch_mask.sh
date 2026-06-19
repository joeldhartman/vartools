#! /bin/bash
#
# Regression test for the -stitch user library: a point masked out of the
# shift FIT must still receive the fitted per-segment shift when the
# correction is APPLIED.  (Masking affects only the fit, not the correction.)
#
# Setup: two segments overlapping in time.  Segment 0 is flat at mag 10.0.
# Segment 1 is the same baseline plus a +0.3 offset, with a 0.1 mag dip at
# t=4,5,6 that is masked out of the fit (mask=0).  After stitching, the +0.3
# offset must be removed from ALL segment-1 points, including the masked dip
# points, so the dip depth (0.1) is preserved: dip points 10.2 -> 9.9,
# out-of-transit points 10.3 -> 10.0.
#
# Prior to the fix the masked dip points kept their original value (10.2)
# while the unmasked points were shifted, turning the 0.1 dip into a +0.2
# bump and corrupting the relative photometry.  All four shift models
# (median, poly, harmseries, and median+groupbytime) share the same apply
# step, so each is checked.
#
# Requires a vartools binary with the stitch extension available (auto-loaded
# from the install USERLIBS dir, or set VARTOOLS to a specific binary).

set -e

VARTOOLS=${VARTOOLS:-vartools}
TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

LC=$TMP/stitch_masktest.lc
# columns: t mag err lcnum mask
cat > $LC <<'EOF'
0.0 10.0 0.01 0 1
1.0 10.0 0.01 0 1
2.0 10.0 0.01 0 1
3.0 10.0 0.01 0 1
4.0 10.0 0.01 0 1
5.0 10.0 0.01 0 1
6.0 10.0 0.01 0 1
7.0 10.0 0.01 0 1
8.0 10.0 0.01 0 1
9.0 10.0 0.01 0 1
0.0 10.3 0.01 1 1
1.0 10.3 0.01 1 1
2.0 10.3 0.01 1 1
3.0 10.3 0.01 1 1
4.0 10.2 0.01 1 0
5.0 10.2 0.01 1 0
6.0 10.2 0.01 1 0
7.0 10.3 0.01 1 1
8.0 10.3 0.01 1 1
9.0 10.3 0.01 1 1
EOF

# List file carries the harmseries period as a per-star scalar (col 2).
LIST=$TMP/stitch_list.txt
echo "$LC 5.0" > $LIST

FMT='t:1,mag:2,err:3,lcnum:4:int,mask:5:double'
OUT=$TMP/out
mkdir -p $OUT
OUTLC=$OUT/$(basename $LC)

# Check the corrected light curve: segment-1 masked dip points (t=4,5,6,
# output rows 15-17) must be 9.9, and the segment-1 out-of-transit points
# (rows 11-14,18-20) must be 10.0.  Tolerance 1e-4.
check_output () {
    local out=$1 method=$2
    awk -v method="$method" '
      function approx(a,b) { d=a-b; if(d<0)d=-d; return d<1e-4 }
      NR>=15 && NR<=17 {
        if(!approx($2,9.9)) { printf "FAIL [%s]: masked dip row %d mag=%s, expected 9.9\n", method, NR, $2; bad=1 }
      }
      (NR>=11 && NR<=14) || (NR>=18 && NR<=20) {
        if(!approx($2,10.0)) { printf "FAIL [%s]: OOT row %d mag=%s, expected 10.0\n", method, NR, $2; bad=1 }
      }
      END { exit bad }
    ' $out
}

run_method () {
    local method=$1; shift
    rm -f $OUTLC
    $VARTOOLS -l $LIST -inlistvars 'per:2:double' -inputlcformat "$FMT" \
        -stitch mag err mask lcnum "$@" -o $OUT > /dev/null 2>&1
    check_output $OUTLC $method
    echo "PASS [$method]: masked dip points received the per-segment shift"
}

run_method median           median
run_method poly             poly 1
run_method harmseries       harmseries per 1
run_method median_grouptime median groupbytime 100

echo "All -stitch masked-shift regression checks passed."

#! /bin/bash

function ReportVartoolsError () {
    local testexamplenumber=$1
    local testcommand=$2
    local output1=$3
    local expectedoutput=$4
    local exitstatus=$5

        cat > /dev/stderr <<EOF
Unit test exited with error code $exitstatus for test number $testexamplenumber

Testcommand:
------------
EOF
       cat $testcommand > /dev/stderr
       cat > /dev/stderr <<EOF

Expected output:
----------------
EOF
       cat $expectedoutput > /dev/stderr
        cat > /dev/stderr <<EOF

Actual output:
--------------
EOF
       cat $output1 > /dev/stderr
       rm $expectedoutput $output1 $testcommand
       exit 1;
}


function CompareOutput () {
    local testexamplenumber=$1
    local testcommand=$2
    local output1=$3
    local expectedoutput=$4

    diff $output1 $expectedoutput > /dev/null || (
        cat > /dev/stderr <<EOF
Unit test produced unexpected output for test number $testexamplenumber

Testcommand:
------------
EOF
       cat $testcommand > /dev/stderr
       cat > /dev/stderr <<EOF

Expected output:
----------------
EOF
       cat $expectedoutput > /dev/stderr
        cat > /dev/stderr <<EOF

Actual output:
--------------
EOF
       cat $output1 > /dev/stderr
       rm $expectedoutput $output1 $testcommand
       exit 1;
   )
   return 0
}

if (( $# != 1 )) ; then
cat > /dev/stderr <<EOF
Usage: $0 VARTOOLSBINARY

where VARTOOLSBINARY is the path to vartools binary to test
EOF
exit 1
fi

VARTOOLS=$1

goodout=$(mktemp /tmp/vartools_unittest.XXXXXX)
testout=$(mktemp /tmp/vartools_unittest.XXXXXX)
testc=$(mktemp /tmp/vartools_unittest.XXXXXX)

# -addnoise example 1
testnumber=1
echo "$testnumber. Testing -addnoise example 1" > /dev/stderr

cat > $testc <<EOF
gawk '{print \$1, 0., 0.005}' EXAMPLES/1 | \
./vartools -i - -header -randseed 1 \
    -addnoise wavelet gamma fix 0.99 sig_red fix 0.005 sig_white fix 0.005 \
    -o EXAMPLES/OUTDIR1/noisesim.txt
EOF

cat > $goodout <<EOF
#Name
stdin
EOF

gawk '{print $1, 0., 0.005}' EXAMPLES/1 | \
$VARTOOLS -i - -header -randseed 1 \
    -addnoise wavelet gamma fix 0.99 sig_red fix 0.005 sig_white fix 0.005 \
    -o EXAMPLES/OUTDIR1/noisesim.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -addnoise example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -addnoise example 2" > /dev/stderr

cat > $testc <<EOF
gawk '{print \$1, 0., 0.005}' EXAMPLES/1 | \
./vartools -i - -header -randseed 1 \
    -addnoise squareexp rho fix 0.01 sig_red fix 0.005 sig_white fix 0.001 \
    -o EXAMPLES/OUTDIR1/noisesim.txt
EOF

cat > $goodout <<EOF
#Name
stdin
EOF

gawk '{print $1, 0., 0.005}' EXAMPLES/1 | \
$VARTOOLS -i - -header -randseed 1 \
    -addnoise squareexp rho fix 0.01 sig_red fix 0.005 sig_white fix 0.001 \
    -o EXAMPLES/OUTDIR1/noisesim.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -alarm example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -alarm example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -header -alarm
EOF

cat > $goodout <<EOF
#Name Alarm_0
EXAMPLES/2 168.89337
EOF

$VARTOOLS -i EXAMPLES/2 -header -alarm \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -aov example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -aov example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline -ascii \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 1 EXAMPLES/OUTDIR1 \
        whiten clip 5. 1
EOF

cat > $goodout <<EOF
Name               = EXAMPLES/2
Period_1_0         =     1.23583047
AOV_1_0            = 18330.55450
AOV_SNR_1_0        = 7260.34519
AOV_NEG_LN_FAP_1_0 = 7633.23914
Period_2_0         =     0.13011657
AOV_2_0            =  45.37532
AOV_SNR_2_0        =  20.21064
AOV_NEG_LN_FAP_2_0 = 339.21684
Period_3_0         =     0.17657032
AOV_3_0            =  49.35941
AOV_SNR_3_0        =  27.39833
AOV_NEG_LN_FAP_3_0 = 368.38239
Period_4_0         =     0.12675512
AOV_4_0            =  15.24895
AOV_SNR_4_0        =   7.57564
AOV_NEG_LN_FAP_4_0 = 103.14339
Period_5_0         =     0.12056969
AOV_5_0            =  16.69317
AOV_SNR_5_0        =   9.40720
AOV_NEG_LN_FAP_5_0 = 115.00693

EOF

$VARTOOLS -i EXAMPLES/2 -oneline -ascii \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 1 EXAMPLES/OUTDIR1 \
        whiten clip 5. 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -aov_harm example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -aov_harm example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline -ascii \
    -aov_harm 1 0.1 10. 0.1 0.01 2 1 EXAMPLES/OUTDIR1 \
        whiten clip 5. 1
EOF

cat > $goodout <<EOF
Name                     = EXAMPLES/2
Period_1_0               =     1.23533969
AOV_HARM_1_0             =    592360
AOV_HARM_SNR_1_0         =   97480.2
AOV_HARM_NEG_LOG_FAP_1_0 =   9730.81
Mean_AOV_HARM_1_0        =   6.36347
RMS_AOV_HARM_1_0         =   6.07665
Period_2_0               =     0.49981672
AOV_HARM_2_0             =   67.3973
AOV_HARM_SNR_2_0         =   14.1748
AOV_HARM_NEG_LOG_FAP_2_0 =   60.3343
Mean_AOV_HARM_2_0        =   4.65647
RMS_AOV_HARM_2_0         =   4.42623

EOF

$VARTOOLS -i EXAMPLES/2 -oneline -ascii \
    -aov_harm 1 0.1 10. 0.1 0.01 2 1 EXAMPLES/OUTDIR1 \
        whiten clip 5. 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -autcorrelation example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -autocorrelation example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -header \
    -autocorrelation 0.0 10. 0.05 EXAMPLES/OUTDIR1
EOF

cat > $goodout <<EOF
#Name
EXAMPLES/2
EOF

$VARTOOLS -i EXAMPLES/2 -header \
    -autocorrelation 0.0 10. 0.05 EXAMPLES/OUTDIR1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -binlc example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -binlc example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -header \
    -binlc 1 binsize 0.01 0 \
    -o EXAMPLES/OUTDIR1/2.bin.txt
EOF

cat > $goodout <<EOF
#Name
EXAMPLES/2
EOF

$VARTOOLS -i EXAMPLES/2 -header \
    -binlc 1 binsize 0.01 0 \
    -o EXAMPLES/OUTDIR1/2.bin.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -binlc example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -binlc example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -header \
    -LS 0.1 10. 0.1 1 0 \
    -Phase ls \
    -binlc 1 nbins 100 0 \
    -o EXAMPLES/OUTDIR1/2.phasebin.txt
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308
EOF

$VARTOOLS -i EXAMPLES/2 -header \
    -LS 0.1 10. 0.1 1 0 \
    -Phase ls \
    -binlc 1 nbins 100 0 \
    -o EXAMPLES/OUTDIR1/2.phasebin.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -BLS example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -BLS example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -ascii -oneline \
    -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
        1 EXAMPLES/OUTDIR1/ 1 EXAMPLES/OUTDIR1/ 0 fittrap \
        nobinnedrms ophcurve EXAMPLES/OUTDIR1/ -0.1 1.1 0.001
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12319314
BLS_Tc_1_0                   = 53727.297609654204
BLS_SN_1_0                   =   7.21628
BLS_SR_1_0                   =   0.00238
BLS_SDE_1_0                  =   6.31221
BLS_Depth_1_0                =   0.01225
BLS_Qtran_1_0                =   0.04041
BLS_Qingress_1_0             =   0.26351
BLS_OOTmag_1_0               =  10.16692
BLS_i1_1_0                   =   0.98003
BLS_i2_1_0                   =   1.02044
BLS_deltaChi2_1_0            = -24208.38967
BLS_fraconenight_1_0         =   0.42879
BLS_Npointsintransit_1_0     =   180
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   142
BLS_Npointsaftertransit_1_0  =   163
BLS_Rednoise_1_0             =   0.00147
BLS_Whitenoise_1_0           =   0.00489
BLS_SignaltoPinknoise_1_0    =  14.91256
BLS_Period_invtransit_0      =     1.14590298
BLS_deltaChi2_invtransit_0   = -3262.95028
BLS_MeanMag_0                =  10.16740

EOF

$VARTOOLS -i EXAMPLES/3.transit -ascii -oneline \
    -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
        1 EXAMPLES/OUTDIR1/ 1 EXAMPLES/OUTDIR1/ 0 fittrap \
        nobinnedrms ophcurve EXAMPLES/OUTDIR1/ -0.1 1.1 0.001 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -BLSFixPer example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -BLSFixPer example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -ascii -oneline \
    -rms \
    -BLSFixPer fix 2.12345 q 0.01 0.1 200 0 0 1 fittrap \
    -rms
EOF

cat > $goodout <<EOF
Name                             = EXAMPLES/3.transit
Mean_Mag_0                       =  10.16727
RMS_0                            =   0.00542
Expected_RMS_0                   =   0.00104
Npoints_0                        =  3417
BLSFixPer_Period_1               =     2.12345000
BLSFixPer_Tc_1                   = 53727.29676321477
BLSFixPer_SR_1                   =   0.00238
BLSFixPer_Depth_1                =   0.01189
BLSFixPer_Qtran_1                =   0.03626
BLSFixPer_Qingress_1             =   0.20623
BLSFixPer_OOTmag_1               =  10.16687
BLSFixPer_i1_1                   =   0.98158
BLSFixPer_i2_1                   =   1.01785
BLSFixPer_deltaChi2_1            = -24228.56603
BLSFixPer_fraconenight_1         =   0.43055
BLSFixPer_Npointsintransit_1     =   166
BLSFixPer_Ntransits_1            =     4
BLSFixPer_Npointsbeforetransit_1 =   129
BLSFixPer_Npointsaftertransit_1  =   144
BLSFixPer_Rednoise_1             =   0.00151
BLSFixPer_Whitenoise_1           =   0.00489
BLSFixPer_SignaltoPinknoise_1    =  14.08946
BLSFixPer_deltaChi2_invtransit_1 = -2934.30109
BLSFixPer_MeanMag_1              =  10.16740
Mean_Mag_2                       =  10.16678
RMS_2                            =   0.00489
Expected_RMS_2                   =   0.00104
Npoints_2                        =  3417

EOF

$VARTOOLS -i EXAMPLES/3.transit -ascii -oneline \
    -rms \
    -BLSFixPer fix 2.12345 q 0.01 0.1 200 0 0 1 fittrap \
    -rms \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -changeerror example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -changeerror example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/4 -ascii -oneline \
    -chi2 \
    -changeerror \
    -chi2
EOF

cat > $goodout <<EOF
Name                = EXAMPLES/4
Chi2_0              =      5.19874
Weighted_Mean_Mag_0 =  10.35137
Mean_Mag_1          =  10.35142
RMS_1               =   0.00209
Npoints_1           =  3227
Chi2_2              =      1.00031
Weighted_Mean_Mag_2 =  10.35142

EOF

$VARTOOLS -i EXAMPLES/4 -ascii -oneline \
    -chi2 \
    -changeerror \
    -chi2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -changevariable example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -changevariable example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list \
    -LS 0.1 100.0 0.1 1 0 \
    -expr 'phase=t' \
    -changevariable t phase \
    -Phase ls \
    -changevariable t t \
    -o EXAMPLES/OUTDIR1 nameformat "%s.phase.txt" \
        columnformat "t:%17.9f,mag:%9.5f,err:%9.5f,phase:%9.5f" \
    -header
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0
EXAMPLES/1    77.76775250 -5710.91013    0.99392   38.81513
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308
EXAMPLES/3    18.29829471  -26.09202    0.03822   11.87999
EXAMPLES/4     0.99383709  -90.59677    0.12488   11.63276
EXAMPLES/5     7.06979568  -86.82830    0.10042   11.03518
EXAMPLES/6    22.21935786  -59.44675    0.07034    9.96198
EXAMPLES/7     0.14747089  -12.56404    0.01933    4.85468
EXAMPLES/8     0.93696087  -79.26500    0.09115    9.35161
EXAMPLES/9     7.06979568  -48.29326    0.05781   11.13720
EXAMPLES/10     0.96906857  -52.91818    0.06257    9.66926
EOF

$VARTOOLS -l EXAMPLES/lc_list \
    -LS 0.1 100.0 0.1 1 0 \
    -expr 'phase=t' \
    -changevariable t phase \
    -Phase ls \
    -changevariable t t \
    -o EXAMPLES/OUTDIR1 nameformat "%s.phase.txt" \
        columnformat "t:%17.9f,mag:%9.5f,err:%9.5f,phase:%9.5f" \
    -header \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -chi2 example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -chi2 example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -header -l EXAMPLES/lc_list -chi2
EOF

cat > $goodout <<EOF
#Name Chi2_0 Weighted_Mean_Mag_0
EXAMPLES/1  34711.71793  10.24430
EXAMPLES/2   1709.50065  10.11178
EXAMPLES/3     27.06322  10.16684
EXAMPLES/4      5.19874  10.35137
EXAMPLES/5      8.26418  10.43932
EXAMPLES/6      3.94650  10.52748
EXAMPLES/7     10.39941  10.56951
EXAMPLES/8      4.19887  10.61132
EXAMPLES/9      2.67020  10.73129
EXAMPLES/10      3.72218  10.87763
EOF

$VARTOOLS -header -l EXAMPLES/lc_list -chi2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -chi2bin example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -chi2bin example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -header -l EXAMPLES/lc_list \
    -chi2bin 5 5.0 10.0 60.0 1440 14400
EOF

cat > $goodout <<EOF
#Name Chi2Bin_5.0_0 Weight_Mean_Mag_Bin_5.0_0 Chi2Bin_10.0_0 Weight_Mean_Mag_Bin_10.0_0 Chi2Bin_60.0_0 Weight_Mean_Mag_Bin_60.0_0 Chi2Bin_1440.0_0 Weight_Mean_Mag_Bin_1440.0_0 Chi2Bin_14400.0_0 Weight_Mean_Mag_Bin_14400.0_0
EXAMPLES/1 228254.48109  10.23970 405745.83067  10.23903 1619850.72683  10.23607 6446329.92664  10.20336 14569840.53581  10.20009
EXAMPLES/2  10710.26878  10.10721  19525.73651  10.10720  82101.02128  10.10680 172928.13393  10.10747  15967.17229  10.11131
EXAMPLES/3     49.89202  10.16683     65.41749  10.16684    151.04489  10.16681    297.75291  10.16657    395.88712  10.16666
EXAMPLES/4     15.04175  10.35126     23.17347  10.35126     74.10398  10.35130    247.44014  10.35135    164.13887  10.35133
EXAMPLES/5     23.59164  10.43913     37.11690  10.43911    128.16520  10.43912    350.68371  10.43902    402.40959  10.43918
EXAMPLES/6     10.03222  10.52737     15.47193  10.52737     41.62336  10.52735     78.01005  10.52738     75.62758  10.52748
EXAMPLES/7     21.01770  10.56942     28.78196  10.56942     55.11081  10.56944     37.05850  10.56945     10.66248  10.56953
EXAMPLES/8     11.74891  10.61119     19.46837  10.61119     75.72787  10.61118    235.99781  10.61118     89.56209  10.61117
EXAMPLES/9      5.60637  10.73120      7.84047  10.73120     22.83600  10.73119     54.28169  10.73115     42.56913  10.73126
EXAMPLES/10      8.95585  10.87751     13.58515  10.87751     36.27533  10.87748     73.19929  10.87747     82.72547  10.87762
EOF

$VARTOOLS -header -l EXAMPLES/lc_list \
    -chi2bin 5 5.0 10.0 60.0 1440 14400 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -clip example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -clip example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -oneline -i EXAMPLES/5 \
    -rms \
    -clip 3. 1 \
    -rms \
    -o EXAMPLES/OUTDIR1/5.clip.txt
EOF

cat > $goodout <<EOF
Name           = EXAMPLES/5
Mean_Mag_0     =  10.43962
RMS_0          =   0.00288
Expected_RMS_0 =   0.00114
Npoints_0      =  3903
Nclip_1        =    52
Mean_Mag_2     =  10.43960
RMS_2          =   0.00266
Expected_RMS_2 =   0.00114
Npoints_2      =  3851

EOF

$VARTOOLS -oneline -i EXAMPLES/5 \
    -rms \
    -clip 3. 1 \
    -rms \
    -o EXAMPLES/OUTDIR1/5.clip.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -converttime example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -converttime example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 -quiet \
    -converttime input jd inputsubtract 2400000. \
        output hjd outputsubtract 2400000. \
        radec fix 88.079166 32.5533 \
    -o EXAMPLES/OUTDIR1/1.hjdutc
EOF

cat > $goodout <<EOF
EOF

$VARTOOLS -i EXAMPLES/1 -quiet \
    -converttime input jd inputsubtract 2400000. \
        output hjd outputsubtract 2400000. \
        radec fix 88.079166 32.5533 \
    -o EXAMPLES/OUTDIR1/1.hjdutc \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -converttime example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -converttime example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1.UTC -quiet \
    -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 \
    -converttime input jd inputsys-utc \
        output bjd outputsubtract 2400000. outputsys-tdb \
        radec fix 88.079166 32.5533 \
        ephemfile $CSPICE_EPHEM_FILE \
        leapsecfile $CSPICE_LEAPSEC_FILE \
        planetdatafile $CSPICE_PLANETDATA_FILE \
        observatory flwo \
    -o EXAMPLES/OUTDIR1/1.bjdtdb
EOF

cat > $goodout <<EOF
EOF

$VARTOOLS -i EXAMPLES/1.UTC -quiet \
    -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 \
    -converttime input jd inputsys-utc \
        output bjd outputsubtract 2400000. outputsys-tdb \
        radec fix 88.079166 32.5533 \
        ephemfile $CSPICE_EPHEM_FILE \
        leapsecfile $CSPICE_LEAPSEC_FILE \
        planetdatafile $CSPICE_PLANETDATA_FILE \
        observatory flwo \
    -o EXAMPLES/OUTDIR1/1.bjdtdb \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -copylc example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -copylc example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -LS 0.1 10. 0.1 1 0 \
    -copylc 100 \
    -expr 'mag=err*gauss()' \
    -LS 0.1 10. 0.1 1 0 \
    -header
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0 LS_Period_1_3 Log10_LS_Prob_1_3 LS_Periodogram_Value_1_3 LS_SNR_1_3
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308     1.30155234   -0.27586    0.00405    5.65274
EXAMPLES/2_copy1.1     1.23440877 -4000.59209    0.99619   45.98308     0.33922684   -0.11086    0.00364    6.04487
EXAMPLES/2_copy1.2     1.23440877 -4000.59209    0.99619   45.98308     0.29541406   -0.00000    0.00196    4.46575
EXAMPLES/2_copy1.3     1.23440877 -4000.59209    0.99619   45.98308     0.31263418   -0.02265    0.00322    5.20972
EXAMPLES/2_copy1.4     1.23440877 -4000.59209    0.99619   45.98308     0.10555514   -0.72736    0.00483    6.67175
EXAMPLES/2_copy1.5     1.23440877 -4000.59209    0.99619   45.98308     1.16943989   -0.17087    0.00381    5.69068
EXAMPLES/2_copy1.6     1.23440877 -4000.59209    0.99619   45.98308     0.10526938   -0.00119    0.00281    5.85720
EXAMPLES/2_copy1.7     1.23440877 -4000.59209    0.99619   45.98308     0.43628473   -0.02488    0.00324    4.83923
EXAMPLES/2_copy1.8     1.23440877 -4000.59209    0.99619   45.98308     0.22298997   -0.39686    0.00428    5.71502
EXAMPLES/2_copy1.9     1.23440877 -4000.59209    0.99619   45.98308     0.60637624   -0.63588    0.00468    8.30357
EXAMPLES/2_copy1.10     1.23440877 -4000.59209    0.99619   45.98308     0.42906346   -0.00226    0.00288    4.61452
EXAMPLES/2_copy1.11     1.23440877 -4000.59209    0.99619   45.98308     0.79557803   -0.28802    0.00407    5.84913
EXAMPLES/2_copy1.12     1.23440877 -4000.59209    0.99619   45.98308     0.10880413   -0.13341    0.00371    6.01680
EXAMPLES/2_copy1.13     1.23440877 -4000.59209    0.99619   45.98308     0.11377872   -0.00180    0.00286    4.32774
EXAMPLES/2_copy1.14     1.23440877 -4000.59209    0.99619   45.98308     2.33887977   -0.00074    0.00277    4.31805
EXAMPLES/2_copy1.15     1.23440877 -4000.59209    0.99619   45.98308     0.13070211   -0.10909    0.00363    5.23612
EXAMPLES/2_copy1.16     1.23440877 -4000.59209    0.99619   45.98308     0.12467776   -1.49602    0.00594    8.11634
EXAMPLES/2_copy1.17     1.23440877 -4000.59209    0.99619   45.98308     0.31421314   -3.30096    0.00844   12.07729
EXAMPLES/2_copy1.18     1.23440877 -4000.59209    0.99619   45.98308     0.30527086   -1.00245    0.00524    5.33263
EXAMPLES/2_copy1.19     1.23440877 -4000.59209    0.99619   45.98308     1.08766087   -1.27288    0.00563    5.84371
EXAMPLES/2_copy1.20     1.23440877 -4000.59209    0.99619   45.98308     2.65872658   -1.56740    0.00604   12.00997
EXAMPLES/2_copy1.21     1.23440877 -4000.59209    0.99619   45.98308     0.14919473   -0.01698    0.00317    4.77122
EXAMPLES/2_copy1.22     1.23440877 -4000.59209    0.99619   45.98308     0.10473771   -0.56756    0.00457    6.91490
EXAMPLES/2_copy1.23     1.23440877 -4000.59209    0.99619   45.98308     0.15308613   -2.08825    0.00677    8.82486
EXAMPLES/2_copy1.24     1.23440877 -4000.59209    0.99619   45.98308     0.15286045   -0.21157    0.00391    5.29428
EXAMPLES/2_copy1.25     1.23440877 -4000.59209    0.99619   45.98308     0.10431623   -0.51571    0.00449    6.89875
EXAMPLES/2_copy1.26     1.23440877 -4000.59209    0.99619   45.98308     0.23371225   -0.36089    0.00422    4.63263
EXAMPLES/2_copy1.27     1.23440877 -4000.59209    0.99619   45.98308     0.19249444   -0.03989    0.00335    4.82121
EXAMPLES/2_copy1.28     1.23440877 -4000.59209    0.99619   45.98308     0.46016422   -0.02513    0.00324    4.95631
EXAMPLES/2_copy1.29     1.23440877 -4000.59209    0.99619   45.98308     0.10235966   -0.18527    0.00385    5.54652
EXAMPLES/2_copy1.30     1.23440877 -4000.59209    0.99619   45.98308     0.13819236   -0.02953    0.00328    4.84721
EXAMPLES/2_copy1.31     1.23440877 -4000.59209    0.99619   45.98308     0.15631709   -0.00075    0.00277    5.59295
EXAMPLES/2_copy1.32     1.23440877 -4000.59209    0.99619   45.98308     2.00690974   -1.77429    0.00633    9.51435
EXAMPLES/2_copy1.33     1.23440877 -4000.59209    0.99619   45.98308     0.12358801   -0.00396    0.00295    5.60542
EXAMPLES/2_copy1.34     1.23440877 -4000.59209    0.99619   45.98308     0.15742460   -0.01580    0.00316    5.10675
EXAMPLES/2_copy1.35     1.23440877 -4000.59209    0.99619   45.98308     0.19978870   -0.00009    0.00260    4.08781
EXAMPLES/2_copy1.36     1.23440877 -4000.59209    0.99619   45.98308     0.57819890   -0.01744    0.00317    6.18595
EXAMPLES/2_copy1.37     1.23440877 -4000.59209    0.99619   45.98308     1.51741956   -0.00726    0.00303    4.92978
EXAMPLES/2_copy1.38     1.23440877 -4000.59209    0.99619   45.98308     0.15107868   -0.13744    0.00372    5.12093
EXAMPLES/2_copy1.39     1.23440877 -4000.59209    0.99619   45.98308     0.13589821   -0.02081    0.00321    5.73873
EXAMPLES/2_copy1.40     1.23440877 -4000.59209    0.99619   45.98308     2.01994162   -0.40871    0.00430    5.08717
EXAMPLES/2_copy1.41     1.23440877 -4000.59209    0.99619   45.98308     0.34296693   -0.00925    0.00307    6.00696
EXAMPLES/2_copy1.42     1.23440877 -4000.59209    0.99619   45.98308     0.11242176   -0.38220    0.00426    5.66921
EXAMPLES/2_copy1.43     1.23440877 -4000.59209    0.99619   45.98308     0.15663193   -0.03027    0.00328    5.71402
EXAMPLES/2_copy1.44     1.23440877 -4000.59209    0.99619   45.98308     0.29682348   -0.06002    0.00345    7.19634
EXAMPLES/2_copy1.45     1.23440877 -4000.59209    0.99619   45.98308     0.10307191   -0.50396    0.00447    6.66329
EXAMPLES/2_copy1.46     1.23440877 -4000.59209    0.99619   45.98308     0.28749631   -0.01008    0.00308    4.58968
EXAMPLES/2_copy1.47     1.23440877 -4000.59209    0.99619   45.98308     0.47783565   -0.21388    0.00391    4.57496
EXAMPLES/2_copy1.48     1.23440877 -4000.59209    0.99619   45.98308     0.18615859   -0.18509    0.00385    4.24486
EXAMPLES/2_copy1.49     1.23440877 -4000.59209    0.99619   45.98308     1.09147723   -0.03929    0.00334    4.78926
EXAMPLES/2_copy1.50     1.23440877 -4000.59209    0.99619   45.98308     3.07991099   -0.30066    0.00410    7.02316
EXAMPLES/2_copy1.51     1.23440877 -4000.59209    0.99619   45.98308     0.82293918   -1.26071    0.00561    8.58879
EXAMPLES/2_copy1.52     1.23440877 -4000.59209    0.99619   45.98308     0.34679042   -0.51171    0.00448    6.79062
EXAMPLES/2_copy1.53     1.23440877 -4000.59209    0.99619   45.98308     0.13619571   -0.04056    0.00335    6.91699
EXAMPLES/2_copy1.54     1.23440877 -4000.59209    0.99619   45.98308     0.21438388   -0.81053    0.00495    5.47653
EXAMPLES/2_copy1.55     1.23440877 -4000.59209    0.99619   45.98308     0.19588855   -1.67512    0.00619   10.19914
EXAMPLES/2_copy1.56     1.23440877 -4000.59209    0.99619   45.98308     0.42438064   -1.17493    0.00549   10.12036
EXAMPLES/2_copy1.57     1.23440877 -4000.59209    0.99619   45.98308     0.44887592   -0.20760    0.00390    6.38369
EXAMPLES/2_copy1.58     1.23440877 -4000.59209    0.99619   45.98308     0.57605743   -0.97678    0.00520    7.56130
EXAMPLES/2_copy1.59     1.23440877 -4000.59209    0.99619   45.98308     0.25539492   -0.76212    0.00488    8.07138
EXAMPLES/2_copy1.60     1.23440877 -4000.59209    0.99619   45.98308     0.37614391   -0.01746    0.00317    4.81228
EXAMPLES/2_copy1.61     1.23440877 -4000.59209    0.99619   45.98308     0.16226970   -0.04374    0.00337    5.64005
EXAMPLES/2_copy1.62     1.23440877 -4000.59209    0.99619   45.98308     0.19166421   -0.19546    0.00387    5.49266
EXAMPLES/2_copy1.63     1.23440877 -4000.59209    0.99619   45.98308     0.14056530   -0.02250    0.00322    4.62794
EXAMPLES/2_copy1.64     1.23440877 -4000.59209    0.99619   45.98308     0.12089818   -0.14223    0.00373    4.70295
EXAMPLES/2_copy1.65     1.23440877 -4000.59209    0.99619   45.98308     0.10135908   -0.00480    0.00297    5.33695
EXAMPLES/2_copy1.66     1.23440877 -4000.59209    0.99619   45.98308     0.35348978   -0.00124    0.00282    4.20429
EXAMPLES/2_copy1.67     1.23440877 -4000.59209    0.99619   45.98308     0.10537636   -0.45705    0.00439    6.02637
EXAMPLES/2_copy1.68     1.23440877 -4000.59209    0.99619   45.98308     0.22789085   -0.32197    0.00414    5.55893
EXAMPLES/2_copy1.69     1.23440877 -4000.59209    0.99619   45.98308     0.23946960   -0.06220    0.00346    5.71592
EXAMPLES/2_copy1.70     1.23440877 -4000.59209    0.99619   45.98308     0.66044800   -0.02767    0.00326    5.67690
EXAMPLES/2_copy1.71     1.23440877 -4000.59209    0.99619   45.98308     0.10649470   -0.01080    0.00309    4.60399
EXAMPLES/2_copy1.72     1.23440877 -4000.59209    0.99619   45.98308     0.64403936   -0.56662    0.00457    5.77561
EXAMPLES/2_copy1.73     1.23440877 -4000.59209    0.99619   45.98308     0.13186563   -0.20002    0.00388    6.11299
EXAMPLES/2_copy1.74     1.23440877 -4000.59209    0.99619   45.98308     0.13691506   -0.00052    0.00273    5.44508
EXAMPLES/2_copy1.75     1.23440877 -4000.59209    0.99619   45.98308     0.10526938   -0.15978    0.00378    5.32760
EXAMPLES/2_copy1.76     1.23440877 -4000.59209    0.99619   45.98308     0.14031169   -0.48672    0.00444    5.80922
EXAMPLES/2_copy1.77     1.23440877 -4000.59209    0.99619   45.98308     3.38120663   -0.20783    0.00390    7.26696
EXAMPLES/2_copy1.78     1.23440877 -4000.59209    0.99619   45.98308     0.22379209   -0.24263    0.00398    5.27832
EXAMPLES/2_copy1.79     1.23440877 -4000.59209    0.99619   45.98308     0.10362126   -0.16216    0.00379    5.93410
EXAMPLES/2_copy1.80     1.23440877 -4000.59209    0.99619   45.98308     0.10219153   -0.07625    0.00352    6.19129
EXAMPLES/2_copy1.81     1.23440877 -4000.59209    0.99619   45.98308     0.10119421   -0.87677    0.00505    9.38442
EXAMPLES/2_copy1.82     1.23440877 -4000.59209    0.99619   45.98308     0.33127903   -0.32275    0.00414    5.30451
EXAMPLES/2_copy1.83     1.23440877 -4000.59209    0.99619   45.98308     1.40122077   -0.07351    0.00351    5.52852
EXAMPLES/2_copy1.84     1.23440877 -4000.59209    0.99619   45.98308     0.79761797   -0.00000    0.00229    4.78307
EXAMPLES/2_copy1.85     1.23440877 -4000.59209    0.99619   45.98308     0.12198863   -0.08352    0.00355    5.63235
EXAMPLES/2_copy1.86     1.23440877 -4000.59209    0.99619   45.98308     0.45948451   -0.11837    0.00366    4.86883
EXAMPLES/2_copy1.87     1.23440877 -4000.59209    0.99619   45.98308     0.17765335   -1.06828    0.00533    6.87710
EXAMPLES/2_copy1.88     1.23440877 -4000.59209    0.99619   45.98308     0.11978090   -0.04885    0.00339    5.02613
EXAMPLES/2_copy1.89     1.23440877 -4000.59209    0.99619   45.98308     0.11598472   -0.00108    0.00280    4.44741
EXAMPLES/2_copy1.90     1.23440877 -4000.59209    0.99619   45.98308     0.13800843   -0.10738    0.00363    4.64128
EXAMPLES/2_copy1.91     1.23440877 -4000.59209    0.99619   45.98308     0.35591649   -0.00001    0.00247    4.31954
EXAMPLES/2_copy1.92     1.23440877 -4000.59209    0.99619   45.98308     0.21768440   -0.12485    0.00368    6.08624
EXAMPLES/2_copy1.93     1.23440877 -4000.59209    0.99619   45.98308     0.10573454   -0.36367    0.00422    5.94058
EXAMPLES/2_copy1.94     1.23440877 -4000.59209    0.99619   45.98308     0.34912571   -0.15915    0.00378    5.31854
EXAMPLES/2_copy1.95     1.23440877 -4000.59209    0.99619   45.98308     0.57287479   -1.31290    0.00568    7.33831
EXAMPLES/2_copy1.96     1.23440877 -4000.59209    0.99619   45.98308     0.36813137   -0.94701    0.00516    6.39972
EXAMPLES/2_copy1.97     1.23440877 -4000.59209    0.99619   45.98308     0.15919704   -0.33958    0.00418    5.20312
EXAMPLES/2_copy1.98     1.23440877 -4000.59209    0.99619   45.98308     0.16432700   -0.36919    0.00423    6.12480
EXAMPLES/2_copy1.99     1.23440877 -4000.59209    0.99619   45.98308     0.28564831   -0.45396    0.00438    6.51046
EXAMPLES/2_copy1.100     1.23440877 -4000.59209    0.99619   45.98308     0.27407137   -0.00087    0.00278    6.17280
EOF

$VARTOOLS -i EXAMPLES/2 -LS 0.1 10. 0.1 1 0 \
    -copylc 100 \
    -expr 'mag=err*gauss()' \
    -LS 0.1 10. 0.1 1 0 \
    -header \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -decorr example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -decorr example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -rms \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms
EOF

cat > $goodout <<EOF
#Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0 Decorr_constant_term_1 Decorr_constant_term_err_1 LCColumn_1_coeff_1_1 LCColumn_1_coeff_err_1_1 LCColumn_1_coeff_2_1 LCColumn_1_coeff_err_2_1 Decorr_chi2_1 Mean_Mag_2 RMS_2 Expected_RMS_2 Npoints_2
EXAMPLES/1  10.24745   0.15944   0.00101  3122     10.0830375984825      0.0000325849746      0.0097933162509      0.0000059117875      0.0002554062775      0.0000001956124      6.68601  10.24728   0.00211   0.00101  3122
EXAMPLES/2  10.11802   0.03663   0.00102  3313     10.1089753109696      0.0000321843725      0.0003489758977      0.0000058446481     -0.0000050541834      0.0000001932089   1755.64028  10.12417   0.03657   0.00102  3313
EXAMPLES/3  10.16674   0.00490   0.00104  3417     10.1659431955329      0.0000324178381      0.0000836121707      0.0000059111825     -0.0000001270601      0.0000001955122     26.40822  10.16662   0.00485   0.00104  3417
EXAMPLES/4  10.35142   0.00209   0.00114  3227     10.3504285453973      0.0000364204458      0.0001938021136      0.0000065445962     -0.0000058852944      0.0000002184623      4.92236  10.35144   0.00205   0.00114  3227
EXAMPLES/5  10.43962   0.00288   0.00114  3903     10.4381409961476      0.0000333590373      0.0001754549866      0.0000061441544     -0.0000036936281      0.0000002040797      7.82090  10.43987   0.00285   0.00114  3903
EXAMPLES/6  10.52762   0.00209   0.00121  3933     10.5268081102353      0.0000351765599      0.0001431023590      0.0000064831117     -0.0000043442275      0.0000002152375      3.83110  10.52772   0.00206   0.00121  3933
EXAMPLES/7  10.56966   0.00349   0.00116  3626     10.5695886705431      0.0000361386139      0.0000047832341      0.0000067249822     -0.0000006353454      0.0000002256379     10.40848  10.56981   0.00348   0.00116  3626
EXAMPLES/8  10.61152   0.00225   0.00125  3957     10.6116126344903      0.0000363856952     -0.0001259754440      0.0000067138507      0.0000051862416      0.0000002229581      4.06144  10.61174   0.00221   0.00125  3957
EXAMPLES/9  10.73139   0.00187   0.00133  3954     10.7308002296744      0.0000386306851      0.0000885383168      0.0000071460623     -0.0000023248367      0.0000002373138      2.62157  10.73147   0.00186   0.00133  3954
EXAMPLES/10  10.87781   0.00236   0.00143  3974     10.8767761937821      0.0000414614974      0.0001687838294      0.0000076663594     -0.0000048171315      0.0000002542888      3.59696  10.87794   0.00233   0.00143  3974
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -rms \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -dftclean example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -dftclean example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline -ascii \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 1 clip 5. 1
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/2
DFTCLEAN_DSPEC_PEAK_FREQ_0_0 =     0.81189711
DFTCLEAN_DSPEC_PEAK_POW_0_0  = 0.000687634
DFTCLEAN_DSPEC_PEAK_SNR_0_0  =   59.8532

EOF

$VARTOOLS -i EXAMPLES/2 -oneline -ascii \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 1 clip 5. 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -dftclean example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -dftclean example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/4 -oneline -ascii \
    -Injectharm fix 0.697516 0 ampfix 0.1 phaserand 0 0 \
    -Injectharm fix 2.123456 0 ampfix 0.05 phaserand 0 0 \
    -Injectharm fix 0.426515 0 ampfix 0.01 phaserand 0 0 \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 3 clip 5. 1 \
        outwfunc EXAMPLES/OUTDIR1 \
        clean 0.5 5.0 outcbeam EXAMPLES/OUTDIR1 \
        outcspec EXAMPLES/OUTDIR1 \
        findcleanpeaks 3 clip 5. 1 \
        verboseout
EOF

cat > $goodout <<EOF
Name                            = EXAMPLES/4
Injectharm_Period_0             =     0.69751600
Injectharm_Fundamental_Amp_0    =   0.10000
Injectharm_Fundamental_Phase_0  =   0.84019
Injectharm_Period_1             =     2.12345600
Injectharm_Fundamental_Amp_1    =   0.05000
Injectharm_Fundamental_Phase_1  =   0.39438
Injectharm_Period_2             =     0.42651500
Injectharm_Fundamental_Amp_2    =   0.01000
Injectharm_Fundamental_Phase_2  =   0.78310
DFTCLEAN_DSPEC_PEAK_FREQ_0_3    =     1.43890675
DFTCLEAN_DSPEC_PEAK_POW_0_3     = 0.0033349
DFTCLEAN_DSPEC_PEAK_SNR_0_3     =    67.932
DFTCLEAN_DSPEC_PEAK_FREQ_1_3    =     0.43408360
DFTCLEAN_DSPEC_PEAK_POW_1_3     = 0.00294915
DFTCLEAN_DSPEC_PEAK_SNR_1_3     =   59.9866
DFTCLEAN_DSPEC_PEAK_FREQ_2_3    =     2.43569132
DFTCLEAN_DSPEC_PEAK_POW_2_3     = 0.00209661
DFTCLEAN_DSPEC_PEAK_SNR_2_3     =   42.4262
DFTCLEAN_DSPEC_AVESPEC_3        = 3.68371e-05
DFTCLEAN_DSPEC_STDSPEC_3        = 4.85494e-05
DFTCLEAN_DSPEC_AVESPEC_NOCLIP_3 = 8.70341e-05
DFTCLEAN_DSPEC_STDSPEC_NOCLIP_3 = 0.000277645
DFTCLEAN_CSPEC_PEAK_FREQ_0_3    =     1.43086817
DFTCLEAN_CSPEC_PEAK_POW_0_3     = 0.00295622
DFTCLEAN_CSPEC_PEAK_SNR_0_3     =      8650
DFTCLEAN_CSPEC_PEAK_FREQ_1_3    =     0.47427653
DFTCLEAN_CSPEC_PEAK_POW_1_3     = 0.000488863
DFTCLEAN_CSPEC_PEAK_SNR_1_3     =   1429.81
DFTCLEAN_CSPEC_PEAK_FREQ_2_3    =     2.34726688
DFTCLEAN_CSPEC_PEAK_POW_2_3     = 2.55442e-05
DFTCLEAN_CSPEC_PEAK_SNR_2_3     =   74.0001
DFTCLEAN_CSPEC_AVESPEC_3        = 3.41729e-07
DFTCLEAN_CSPEC_STDSPEC_3        = 2.56197e-07
DFTCLEAN_CSPEC_AVESPEC_NOCLIP_3 = 0.000309941
DFTCLEAN_CSPEC_STDSPEC_NOCLIP_3 = 5.52828e-05

EOF

$VARTOOLS -i EXAMPLES/4 -oneline -ascii \
    -Injectharm fix 0.697516 0 ampfix 0.1 phaserand 0 0 \
    -Injectharm fix 2.123456 0 ampfix 0.05 phaserand 0 0 \
    -Injectharm fix 0.426515 0 ampfix 0.01 phaserand 0 0 \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 3 clip 5. 1 \
        outwfunc EXAMPLES/OUTDIR1 \
        clean 0.5 5.0 outcbeam EXAMPLES/OUTDIR1 \
        outcspec EXAMPLES/OUTDIR1 \
        findcleanpeaks 3 clip 5. 1 \
        verboseout \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -ensemblerescalesig example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -ensemblerescalesig example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -chi2 -ensemblerescalesig 3. -chi2
EOF

cat > $goodout <<EOF
#Name Chi2_0 Weighted_Mean_Mag_0 SigmaRescaleFactor_1 Chi2_2 Weighted_Mean_Mag_2
EXAMPLES/1  34711.71793  10.24430   0.32888   3754.43599  10.24786
EXAMPLES/2   1709.50065  10.11178   0.34008    197.71560  10.11842
EXAMPLES/3     27.06322  10.16684   0.36339      3.57371  10.16674
EXAMPLES/4      5.19874  10.35137   0.35266      0.64655  10.35142
EXAMPLES/5      8.26418  10.43932   0.38819      1.24536  10.43964
EXAMPLES/6      3.94650  10.52748   0.40863      0.65899  10.52762
EXAMPLES/7     10.39941  10.56951   0.42085      1.84188  10.56968
EXAMPLES/8      4.19887  10.61132   0.42869      0.77164  10.61154
EXAMPLES/9      2.67020  10.73129   0.44677      0.53297  10.73140
EXAMPLES/10      3.72218  10.87763   0.48353      0.87025  10.87782
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -chi2 -ensemblerescalesig 3. -chi2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -expr example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -expr example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 -expr 'mag=sqrt(mag+5)' -o EXAMPLES/1.add
EOF

cat > $goodout <<EOF
EXAMPLES/1
EOF

$VARTOOLS -i EXAMPLES/1 -expr 'mag=sqrt(mag+5)' -o EXAMPLES/1.add \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -expr example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -expr example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -LS 0.1 10. 0.1 1 0 \
    -rms -chi2 \
    -expr 'mag2=mag' \
    -Killharm ls 0 0 0 \
    -rms -chi2 \
    -expr \
        'mag=(Npoints_5*(Chi2_6-Chi2_2)<-10000)*mag+
            (Npoints_5*(Chi2_6-Chi2_2)>=-10000)*mag2' \
    -o EXAMPLES/OUTDIR1 nameformat '%s.cleanharm'
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0 Mean_Mag_1 RMS_1 Expected_RMS_1 Npoints_1 Chi2_2 Weighted_Mean_Mag_2 Killharm_Mean_Mag_4 Killharm_Period_1_4 Killharm_Per1_Fundamental_Sincoeff_4 Killharm_Per1_Fundamental_Coscoeff_4 Killharm_Per1_Amplitude_4 Mean_Mag_5 RMS_5 Expected_RMS_5 Npoints_5 Chi2_6 Weighted_Mean_Mag_6
EXAMPLES/1     0.97821072 -5530.60640    0.76156   29.55922  10.24745   0.15944   0.00101  3122  34711.71793  10.24430  10.29306     0.97821072  -0.08947  -0.17167   0.38717  10.24061   0.08469   0.00101  3122   8276.59620  10.24430
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308  10.11802   0.03663   0.00102  3313   1709.50065  10.11178  10.12217     1.23440877   0.05008  -0.00222   0.10026  10.11176   0.00231   0.00102  3313      6.51484  10.11178
EXAMPLES/3     1.14786351  -23.48814    0.03471   10.70790  10.16674   0.00490   0.00104  3417     27.06322  10.16684  10.16696     1.14786351  -0.00104   0.00073   0.00253  10.16671   0.00479   0.00104  3417     26.12389  10.16684
EXAMPLES/4     0.99383709  -90.59677    0.12488   11.63276  10.35142   0.00209   0.00114  3227      5.19874  10.35137  10.35322     0.99383709  -0.00080   0.00242   0.00509  10.35138   0.00204   0.00114  3227      4.54950  10.35137
EXAMPLES/5     7.06979568  -86.82830    0.10042   11.03518  10.43962   0.00288   0.00114  3903      8.26418  10.43932  10.43952     7.06979568   0.00020   0.00127   0.00257  10.43961   0.00278   0.00114  3903      7.43429  10.43932
EXAMPLES/6     0.96009571  -58.38646    0.06910    9.76922  10.52762   0.00209   0.00121  3933      3.94650  10.52748  10.52778     0.96009571  -0.00046   0.00070   0.00167  10.52756   0.00202   0.00121  3933      3.67381  10.52748
EXAMPLES/7     0.14747089  -12.56404    0.01933    4.85468  10.56966   0.00349   0.00116  3626     10.39941  10.56951  10.56957     0.14747089   0.00061  -0.00023   0.00130  10.56965   0.00346   0.00116  3626     10.19836  10.56951
EXAMPLES/8     0.93696087  -79.26500    0.09115    9.35161  10.61152   0.00225   0.00125  3957      4.19887  10.61132  10.61128     0.93696087  -0.00029   0.00087   0.00184  10.61152   0.00215   0.00125  3957      3.81615  10.61132
EXAMPLES/9     7.06979568  -48.29326    0.05781   11.13720  10.73139   0.00187   0.00133  3954      2.67020  10.73129  10.73143     7.06979568  -0.00005   0.00066   0.00132  10.73137   0.00182   0.00133  3954      2.51585  10.73129
EXAMPLES/10     0.96906857  -52.91818    0.06257    9.66926  10.87781   0.00236   0.00143  3974      3.72218  10.87763  10.87782     0.96906857  -0.00078  -0.00048   0.00183  10.87775   0.00230   0.00143  3974      3.48930  10.87763
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -LS 0.1 10. 0.1 1 0 \
    -rms -chi2 \
    -expr 'mag2=mag' \
    -Killharm ls 0 0 0 \
    -rms -chi2 \
    -expr \
        'mag=(Npoints_5*(Chi2_6-Chi2_2)<-10000)*mag+
            (Npoints_5*(Chi2_6-Chi2_2)>=-10000)*mag2' \
    -o EXAMPLES/OUTDIR1 nameformat '%s.cleanharm' \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -expr example 3
testnumber=$((testnumber+1))
echo "$testnumber. Testing -expr example 3" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 \
    -expr 'flux=10^(-0.4*(mag-25.0))' \
    -stats flux median \
    -expr 'flux=flux/STATS_flux_MEDIAN_1' \
    -stats flux,mag median,stddev \
    -oneline
EOF

cat > $goodout <<EOF
Name                = EXAMPLES/1
STATS_flux_MEDIAN_1 = 842674.79516438092
STATS_flux_MEDIAN_3 = 1
STATS_flux_STDDEV_3 = 0.12908654281197893
STATS_mag_MEDIAN_3  = 10.18585
STATS_mag_STDDEV_3  = 0.15946976931434603

EOF

$VARTOOLS -i EXAMPLES/1 \
    -expr 'flux=10^(-0.4*(mag-25.0))' \
    -stats flux median \
    -expr 'flux=flux/STATS_flux_MEDIAN_1' \
    -stats flux,mag median,stddev \
    -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -FFT example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -FFT example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/11 \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr \
      'fftreal=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftreal' \
    -expr \
      'fftimag=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftimag' \
    -IFFT fftreal fftimag mag_filter NULL \
    -o EXAMPLES/11.highpassfftfilter columnformat t,mag_filter 
EOF

cat > $goodout <<EOF
EXAMPLES/11   0.00070   0.76113   0.00100 100000
EOF

$VARTOOLS -i EXAMPLES/11 \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr \
      'fftreal=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftreal' \
    -expr \
      'fftimag=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftimag' \
    -IFFT fftreal fftimag mag_filter NULL \
    -o EXAMPLES/11.highpassfftfilter columnformat t,mag_filter \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -FFT example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -FFT example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 \
    -resample splinemonotonic gaps percentile_sep 80 bspline \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr 'fftreal1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal' \
    -expr 'fftimag1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag' \
    -IFFT fftreal1 fftimag1 mag_filter NULL \
    -expr 'fftreal2=fftreal-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal)' \
    -expr 'fftimag2=fftimag-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag)' \
    -IFFT fftreal2 fftimag2 mag_filter2 NULL \
    -resample linear file fix EXAMPLES/2 column 1 \
    -expr 'mag_filter=mag_filter+Mean_Mag_2' \
    -o EXAMPLES/2.fftfilter columnformat t,mag_filter,mag_filter2,mag 
EOF

cat > $goodout <<EOF
EXAMPLES/2  10.11155   0.02901   0.00028 66186
EOF

$VARTOOLS -i EXAMPLES/2 \
    -resample splinemonotonic gaps percentile_sep 80 bspline \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr 'fftreal1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal' \
    -expr 'fftimag1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag' \
    -IFFT fftreal1 fftimag1 mag_filter NULL \
    -expr 'fftreal2=fftreal-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal)' \
    -expr 'fftimag2=fftimag-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag)' \
    -IFFT fftreal2 fftimag2 mag_filter2 NULL \
    -resample linear file fix EXAMPLES/2 column 1 \
    -expr 'mag_filter=mag_filter+Mean_Mag_2' \
    -o EXAMPLES/2.fftfilter columnformat t,mag_filter,mag_filter2,mag \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -findblends example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -findblends example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list_testblend -header \
    -LS 0.1 10. 0.1 1 0 \
    -findblends 2.0 fixcolumn LS_Period_1_0 
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0 FindBlends_Period_1 FindBlends_LCname_1 FindBlends_FluxAmp_1
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308     1.23440877 EXAMPLES/2     82314
EXAMPLES/2.testblend     1.23440877 -1481.99951    0.87328   43.35308     1.23440877 EXAMPLES/2     82314
EOF

$VARTOOLS -l EXAMPLES/lc_list_testblend -header \
    -LS 0.1 10. 0.1 1 0 \
    -findblends 2.0 fixcolumn LS_Period_1_0 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -fluxtomag example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -fluxtomag example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/kplr000757076-2009166043257_llc.fits \
    -readformat 0 1 10 11 \
    -fluxtomag 25.0 0 \
    -o EXAMPLES/OUTDIR1/kplr000757076-2009166043257_llc.asc.txt
EOF

cat > $goodout <<EOF
EXAMPLES/kplr000757076-2009166043257_llc.fits
EOF

$VARTOOLS -i EXAMPLES/kplr000757076-2009166043257_llc.fits \
    -readformat 0 1 10 11 \
    -fluxtomag 25.0 0 \
    -o EXAMPLES/OUTDIR1/kplr000757076-2009166043257_llc.asc.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -GetLSAmpThresh example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -GetLSAmpThresh example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline \
    -LS 0.1 10. 0.1 1 0 \
    -Killharm ls 0 0 0 fitonly \
    -GetLSAmpThresh ls 0.1 -100 harm 0 0
EOF

cat > $goodout <<EOF
Name                                 = EXAMPLES/2
LS_Period_1_0                        =     1.23440877
Log10_LS_Prob_1_0                    = -4000.59209
LS_Periodogram_Value_1_0             =    0.99619
LS_SNR_1_0                           =   45.98308
Killharm_Mean_Mag_1                  =  10.12217
Killharm_Period_1_1                  =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_1 =   0.05008
Killharm_Per1_Fundamental_Coscoeff_1 =  -0.00222
Killharm_Per1_Amplitude_1            =   0.10026
LS_AmplitudeScaleFactor_2            =   0.02425
LS_MinimumAmplitude_2                =   0.00243

EOF

$VARTOOLS -i EXAMPLES/2 -oneline \
    -LS 0.1 10. 0.1 1 0 \
    -Killharm ls 0 0 0 fitonly \
    -GetLSAmpThresh ls 0.1 -100 harm 0 0 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -if example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -if example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -rms \
    -if 'RMS_0>10*Expected_RMS_0' \
        -if 'RMS_0 > 0.1' \
            -stats mag stddev \
        -else \
            -stats mag pct30 \
        -fi \
    -elif 'Npoints_0>3900' \
        -stats mag kurtosis \
    -else \
        -rms \
    -fi \
    -header
EOF

cat > $goodout <<EOF
#Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0 STATS_mag_STDDEV_3 STATS_mag_PCT30.00_5 STATS_mag_KURTOSIS_8 Mean_Mag_10 RMS_10 Expected_RMS_10 Npoints_10
EXAMPLES/1  10.24745   0.15944   0.00101  3122 0.15946976931434603 0 0   0.00000   0.00000   0.00000     0
EXAMPLES/2  10.11802   0.03663   0.00102  3313 0 10.0855 0   0.00000   0.00000   0.00000     0
EXAMPLES/3  10.16674   0.00490   0.00104  3417 0 0 0  10.16674   0.00490   0.00104  3417
EXAMPLES/4  10.35142   0.00209   0.00114  3227 0 0 0  10.35142   0.00209   0.00114  3227
EXAMPLES/5  10.43962   0.00288   0.00114  3903 0 0 4.3605470261901687   0.00000   0.00000   0.00000     0
EXAMPLES/6  10.52762   0.00209   0.00121  3933 0 0 4.0802664679027298   0.00000   0.00000   0.00000     0
EXAMPLES/7  10.56966   0.00349   0.00116  3626 0 0 0  10.56966   0.00349   0.00116  3626
EXAMPLES/8  10.61152   0.00225   0.00125  3957 0 0 3.3159700561447285   0.00000   0.00000   0.00000     0
EXAMPLES/9  10.73139   0.00187   0.00133  3954 0 0 6.745669809252357   0.00000   0.00000   0.00000     0
EXAMPLES/10  10.87781   0.00236   0.00143  3974 0 0 4.0178189593727724   0.00000   0.00000   0.00000     0
EOF

$VARTOOLS -l EXAMPLES/lc_list -rms \
    -if 'RMS_0>10*Expected_RMS_0' \
        -if 'RMS_0 > 0.1' \
            -stats mag stddev \
        -else \
            -stats mag pct30 \
        -fi \
    -elif 'Npoints_0>3900' \
        -stats mag kurtosis \
    -else \
        -rms \
    -fi \
    -header \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Injectharm example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Injectharm example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 -randseed 1 -oneline \
    -Injectharm rand 1.0 5.0 \
        0 amplogrand 0.01 0.1 phaserand \
        0 1 EXAMPLES/OUTDIR1 \
    -LS 0.1 10.0 0.1 1 0
EOF

cat > $goodout <<EOF
Name                           = EXAMPLES/3
Injectharm_Period_0            =     4.36075087
Injectharm_Fundamental_Amp_0   =   0.02480
Injectharm_Fundamental_Phase_0 =   0.78310
LS_Period_1_1                  =     4.38128183
Log10_LS_Prob_1_1              = -1910.50094
LS_Periodogram_Value_1_1       =    0.92429
LS_SNR_1_1                     =   41.57518

EOF

$VARTOOLS -i EXAMPLES/3 -randseed 1 -oneline \
    -Injectharm rand 1.0 5.0 \
        0 amplogrand 0.01 0.1 phaserand \
        0 1 EXAMPLES/OUTDIR1 \
    -LS 0.1 10.0 0.1 1 0 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Injectharm example 2; This version differs from the example on the
#  website by not running with the -parallel option. Using random numbers
#  and -parallel together leads to non-reproducible output which would break
#  this test.
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Injectharm example 2" > /dev/stderr

cat > $testc <<EOF
echo EXAMPLES/4 | \
    gawk '{amp = 0.25; \
        for(i=1; i <= 10; i += 1) { \
            print \$1, amp; amp = amp*0.5; \
        }}' | \
    ./vartools -l - -header -numbercolumns \
    -randseed 2 \
    -Injectharm fix 0.514333 10 \
        amplist phaserand \
        ampfix 0.47077 amprel phasefix 0.60826 phaserel \
        ampfix 0.35916 amprel phasefix 0.26249 phaserel \
        ampfix 0.23631 amprel phasefix -0.06843 phaserel \
        ampfix 0.16353 amprel phasefix 0.60682 phaserel \
        ampfix 0.10621 amprel phasefix 0.28738 phaserel \
        ampfix 0.06203 amprel phasefix 0.95751 phaserel \
        ampfix 0.03602 amprel phasefix 0.58867 phaserel \
        ampfix 0.02900 amprel phasefix 0.22322 phaserel \
        ampfix 0.01750 amprel phasefix 0.94258 phaserel \
        ampfix 0.00768 amprel phasefix 0.66560 phaserel \
        0 0 \
    -o EXAMPLES/OUTDIR1 nameformat InjectRRLyrae.%d.txt \
    -LS 0.1 10.0 0.01 2 0 \
    -aov_harm 2 0.1 10.0 0.1 0.01 2 0 | \
    gawk '{print \$1, \$2, \$3, \$25, \$28, \$33, \$35}'
EOF

cat > $goodout <<EOF
#1_Name 2_Injectharm_Period_0 3_Injectharm_Fundamental_Amp_0 25_LS_Period_1_2 28_LS_SNR_1_2 33_Period_1_3 35_AOV_HARM_SNR_1_3
EXAMPLES/4 0.51433300 0.25000 0.51433699 27.44130 0.51431981 127.663
EXAMPLES/4 0.51433300 0.12500 0.51510351 33.13344 0.51449853 161.823
EXAMPLES/4 0.51433300 0.06250 0.51399704 32.28413 0.51415838 135.203
EXAMPLES/4 0.51433300 0.03125 0.51416696 27.09540 0.51424337 114.268
EXAMPLES/4 0.51433300 0.01562 0.51399704 24.76528 0.51432840 122.375
EXAMPLES/4 0.51433300 0.00781 0.51416696 24.94100 0.51432840 88.0958
EXAMPLES/4 0.51433300 0.00391 0.51382724 24.25570 0.51390355 49.235
EXAMPLES/4 0.51433300 0.00195 1.06567664 22.61650 1.06783462 21.7288
EXAMPLES/4 0.51433300 0.00098 0.99383709 13.90822 0.99317043 14.0428
EXAMPLES/4 0.51433300 0.00049 0.99479057 12.94206 0.99412262 13.758
EOF

echo EXAMPLES/4 | \
    gawk '{amp = 0.25; \
        for(i=1; i <= 10; i += 1) { \
            print $1, amp; amp = amp*0.5; \
        }}' | \
    $VARTOOLS -l - -header -numbercolumns \
    -randseed 2 \
    -Injectharm fix 0.514333 10 \
        amplist phaserand \
        ampfix 0.47077 amprel phasefix 0.60826 phaserel \
        ampfix 0.35916 amprel phasefix 0.26249 phaserel \
        ampfix 0.23631 amprel phasefix -0.06843 phaserel \
        ampfix 0.16353 amprel phasefix 0.60682 phaserel \
        ampfix 0.10621 amprel phasefix 0.28738 phaserel \
        ampfix 0.06203 amprel phasefix 0.95751 phaserel \
        ampfix 0.03602 amprel phasefix 0.58867 phaserel \
        ampfix 0.02900 amprel phasefix 0.22322 phaserel \
        ampfix 0.01750 amprel phasefix 0.94258 phaserel \
        ampfix 0.00768 amprel phasefix 0.66560 phaserel \
        0 0 \
    -o EXAMPLES/OUTDIR1 nameformat InjectRRLyrae.%d.txt \
    -LS 0.1 10.0 0.01 2 0 \
    -aov_harm 2 0.1 10.0 0.1 0.01 2 0 | \
    gawk '{print $1, $2, $3, $25, $28, $33, $35}' \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Injecttransit example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Injecttransit example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 -oneline -randseed 1 \
    -Injecttransit lograndfreq 0.2 2.0 \
        Rpfix 1.0 Mpfix 1.0 \
        phaserand sinirand \
        eomega efix 0. ofix 0. \
        Mstarfix 1.0 Rstarfix 1.0 \
        quad ldfix 0.3471 0.3180 \
        1 EXAMPLES/OUTDIR1 \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 1 EXAMPLES/OUTDIR1 \
        1 fittrap
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3
Injecttransit_Period_0       =     0.72240757
Injecttransit_Rp_0           =   1.00000
Injecttransit_Mp_0           =   1.00000
Injecttransit_phase_0        =   0.39438
Injecttransit_sin_i_0        =   0.96568
Injecttransit_h_0            =   0.00000
Injecttransit_k_0            =   0.00000
Injecttransit_Mstar_0        =   1.00000
Injecttransit_Rstar_0        =   1.00000
Injecttransit_ld_1_0         =   0.34710
Injecttransit_ld_2_0         =   0.31800
BLS_Period_1_1               =     0.72239761
BLS_Tc_1_1                   = 53725.889333375417
BLS_SN_1_1                   =  37.03134
BLS_SR_1_1                   =   0.00159
BLS_SDE_1_1                  =   5.93300
BLS_Depth_1_1                =   0.00935
BLS_Qtran_1_1                =   0.05315
BLS_Qingress_1_1             =   0.19369
BLS_OOTmag_1_1               =  10.16683
BLS_i1_1_1                   =   0.96375
BLS_i2_1_1                   =   1.01691
BLS_deltaChi2_1_1            = -10847.39428
BLS_fraconenight_1_1         =   0.29405
BLS_Npointsintransit_1_1     =   179
BLS_Ntransits_1_1            =    10
BLS_Npointsbeforetransit_1_1 =   227
BLS_Npointsaftertransit_1_1  =   186
BLS_Rednoise_1_1             =   0.00164
BLS_Whitenoise_1_1           =   0.00490
BLS_SignaltoPinknoise_1_1    =  14.75185
BLS_Period_invtransit_1      =     1.25145489
BLS_deltaChi2_invtransit_1   = -2806.09675
BLS_MeanMag_1                =  10.16718

EOF

$VARTOOLS -i EXAMPLES/3 -oneline -randseed 1 \
    -Injecttransit lograndfreq 0.2 2.0 \
        Rpfix 1.0 Mpfix 1.0 \
        phaserand sinirand \
        eomega efix 0. ofix 0. \
        Mstarfix 1.0 Rstarfix 1.0 \
        quad ldfix 0.3471 0.3180 \
        1 EXAMPLES/OUTDIR1 \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 1 EXAMPLES/OUTDIR1 \
        1 fittrap \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Jstet example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Jstet example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -Jstet 0.5 EXAMPLES/dates_tfa
EOF

cat > $goodout <<EOF
#Name Jstet_0 Kurtosis_0 Lstet_0
EXAMPLES/1  98.13279   0.96779  94.97154
EXAMPLES/2  30.19309   0.94719  28.59852
EXAMPLES/3   0.65597   0.92816   0.60885
EXAMPLES/4   0.34402   0.84500   0.29070
EXAMPLES/5   0.58730   0.92120   0.54102
EXAMPLES/6   0.34455   0.93794   0.32317
EXAMPLES/7   0.41754   0.92501   0.38623
EXAMPLES/8   0.46381   0.96124   0.44583
EXAMPLES/9   0.22075   0.80997   0.17880
EXAMPLES/10   0.25784   0.92806   0.23929
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -Jstet 0.5 EXAMPLES/dates_tfa \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Killharm example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Killharm example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline \
    -LS 0.1 10. 0.1 1 0 \
    -rms -chi2 \
    -Killharm ls 0 0 0 \
    -rms -chi2
EOF

cat > $goodout <<EOF
Name                                 = EXAMPLES/2
LS_Period_1_0                        =     1.23440877
Log10_LS_Prob_1_0                    = -4000.59209
LS_Periodogram_Value_1_0             =    0.99619
LS_SNR_1_0                           =   45.98308
Mean_Mag_1                           =  10.11802
RMS_1                                =   0.03663
Expected_RMS_1                       =   0.00102
Npoints_1                            =  3313
Chi2_2                               =   1709.50065
Weighted_Mean_Mag_2                  =  10.11178
Killharm_Mean_Mag_3                  =  10.12217
Killharm_Period_1_3                  =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_3 =   0.05008
Killharm_Per1_Fundamental_Coscoeff_3 =  -0.00222
Killharm_Per1_Amplitude_3            =   0.10026
Mean_Mag_4                           =  10.11176
RMS_4                                =   0.00231
Expected_RMS_4                       =   0.00102
Npoints_4                            =  3313
Chi2_5                               =      6.51484
Weighted_Mean_Mag_5                  =  10.11178

EOF

$VARTOOLS -i EXAMPLES/2 -oneline \
    -LS 0.1 10. 0.1 1 0 \
    -rms -chi2 \
    -Killharm ls 0 0 0 \
    -rms -chi2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Killharm example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Killharm example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/M3.V006.lc -oneline \
    -Killharm fix 1 0.514333 10 0 1 \
        EXAMPLES/OUTDIR1/ fitonly outRphi
EOF

cat > $goodout <<EOF
Name                            = EXAMPLES/M3.V006.lc
Killharm_Mean_Mag_0             =  15.77123
Killharm_Period_1_0             =     0.51433300
Killharm_Per1_Fundamental_Amp_0 =   0.38041
Killharm_Per1_Fundamental_Phi_0 =  -0.07662
Killharm_Per1_Harm_R_2_1_0      =   0.47077
Killharm_Per1_Harm_Phi_2_1_0    =   0.60826
Killharm_Per1_Harm_R_3_1_0      =   0.35917
Killharm_Per1_Harm_Phi_3_1_0    =   0.26249
Killharm_Per1_Harm_R_4_1_0      =   0.23631
Killharm_Per1_Harm_Phi_4_1_0    =  -0.06843
Killharm_Per1_Harm_R_5_1_0      =   0.16353
Killharm_Per1_Harm_Phi_5_1_0    =   0.60682
Killharm_Per1_Harm_R_6_1_0      =   0.10621
Killharm_Per1_Harm_Phi_6_1_0    =   0.28738
Killharm_Per1_Harm_R_7_1_0      =   0.06203
Killharm_Per1_Harm_Phi_7_1_0    =   0.95751
Killharm_Per1_Harm_R_8_1_0      =   0.03602
Killharm_Per1_Harm_Phi_8_1_0    =   0.58867
Killharm_Per1_Harm_R_9_1_0      =   0.02900
Killharm_Per1_Harm_Phi_9_1_0    =   0.22322
Killharm_Per1_Harm_R_10_1_0     =   0.01750
Killharm_Per1_Harm_Phi_10_1_0   =   0.94258
Killharm_Per1_Harm_R_11_1_0     =   0.00768
Killharm_Per1_Harm_Phi_11_1_0   =   0.66560
Killharm_Per1_Amplitude_0       =   1.11128

EOF

$VARTOOLS -i EXAMPLES/M3.V006.lc -oneline \
    -Killharm fix 1 0.514333 10 0 1 \
        EXAMPLES/OUTDIR1/ fitonly outRphi \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -linfit example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -linfit example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 \
    -stats t min \
    -expr t0=STATS_t_MIN_0 \
    -linfit 'a*(t-t0)^2+b*(t-t0)+c' 'a,b,c' \
    -oneline
EOF

cat > $goodout <<EOF
Name          = EXAMPLES/1
STATS_t_MIN_0 = 53725.173920000001
Linfit_a_2    = 0.00025540627746042932
Linfit_erra_2 = 1.9561241332987699e-07
Linfit_b_2    = 0.0097933162509034055
Linfit_errb_2 = 5.9117874714733109e-06
Linfit_c_2    = 10.083037598482507
Linfit_errc_2 = 3.2584974556662493e-05

EOF

$VARTOOLS -i EXAMPLES/1 \
    -stats t min \
    -expr t0=STATS_t_MIN_0 \
    -linfit 'a*(t-t0)^2+b*(t-t0)+c' 'a,b,c' \
    -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -LS example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -LS example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -oneline -ascii \
    -LS 0.1 10. 0.1 5 1 EXAMPLES/OUTDIR1 whiten clip 5. 1
EOF

cat > $goodout <<EOF
Name                     = EXAMPLES/2
LS_Period_1_0            =     1.23440877
Log10_LS_Prob_1_0        = -4000.59209
LS_Periodogram_Value_1_0 =    0.99619
LS_SNR_1_0               =  116.60562
LS_Period_2_0            =     0.55747493
Log10_LS_Prob_2_0        =  -85.74406
LS_Periodogram_Value_2_0 =    0.11152
LS_SNR_2_0               =   14.48656
LS_Period_3_0            =     0.54669773
Log10_LS_Prob_3_0        =  -42.14997
LS_Periodogram_Value_3_0 =    0.05837
LS_SNR_3_0               =    8.23508
LS_Period_4_0            =     0.33412568
Log10_LS_Prob_4_0        =  -33.07258
LS_Periodogram_Value_4_0 =    0.04868
LS_SNR_4_0               =    9.78596
LS_Period_5_0            =     0.25922584
Log10_LS_Prob_5_0        =  -24.23642
LS_Periodogram_Value_5_0 =    0.03682
LS_SNR_5_0               =    9.51180

EOF

$VARTOOLS -i EXAMPLES/2 -oneline -ascii \
    -LS 0.1 10. 0.1 5 1 EXAMPLES/OUTDIR1 whiten clip 5. 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -MandelAgolTransit example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -MandelAgolTransit example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -MandelAgolTransit bls quad 0.3471 0.3180 \
        1 1 1 1 0 0 1 0 0 0 0 1 EXAMPLES/OUTDIR1
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12312625
BLS_Tc_1_0                   = 53727.297046247397
BLS_SN_1_0                   =  38.39425
BLS_SR_1_0                   =   0.00237
BLS_SDE_1_0                  =   4.77204
BLS_Depth_1_0                =   0.01136
BLS_Qtran_1_0                =   0.03000
BLS_i1_1_0                   =   0.98500
BLS_i2_1_0                   =   1.01000
BLS_deltaChi2_1_0            = -24130.93833
BLS_fraconenight_1_0         =   0.42662
BLS_Npointsintransit_1_0     =   146
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   106
BLS_Npointsaftertransit_1_0  =   120
BLS_Rednoise_1_0             =   0.00156
BLS_Whitenoise_1_0           =   0.00490
BLS_SignaltoPinknoise_1_0    =  12.89679
BLS_Period_invtransit_0      =     1.14599569
BLS_deltaChi2_invtransit_0   = -3289.67397
BLS_MeanMag_0                =  10.16740
MandelAgolTransit_Period_1   =     2.12328176
MandelAgolTransit_T0_1       = 53727.29695831
MandelAgolTransit_r_1        =   0.09789
MandelAgolTransit_a_1        =   9.35954
MandelAgolTransit_bimpact_1  =   0.33094
MandelAgolTransit_inc_1      =  87.97368
MandelAgolTransit_e_1        =   0.00000
MandelAgolTransit_omega_1    =   0.00000
MandelAgolTransit_mconst_1   =  10.16687
MandelAgolTransit_ldcoeff1_1 =   0.34710
MandelAgolTransit_ldcoeff2_1 =   0.31800
MandelAgolTransit_chi2_1     =  27.06054

EOF

$VARTOOLS -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -MandelAgolTransit bls quad 0.3471 0.3180 \
        1 1 1 1 0 0 1 0 0 0 0 1 EXAMPLES/OUTDIR1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -medianfilter example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -medianfilter example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 -oneline -chi2 \
    -savelc \
    -medianfilter 0.05 \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianhighpass \
    -restorelc 1 \
    -medianfilter 0.05 replace \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianlowpass
EOF

cat > $goodout <<EOF
Name                = EXAMPLES/1
Chi2_0              =  34711.71793
Weighted_Mean_Mag_0 =  10.24430
Chi2_3              =      5.95454
Weighted_Mean_Mag_3 =  -0.00009
Chi2_7              =  34727.65120
Weighted_Mean_Mag_7 =  10.24440

EOF

$VARTOOLS -i EXAMPLES/1 -oneline -chi2 \
    -savelc \
    -medianfilter 0.05 \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianhighpass \
    -restorelc 1 \
    -medianfilter 0.05 replace \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianlowpass \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -microlens example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -microlens example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/4.microlensinject -oneline \
    -microlens f0 auto f1 auto u0 auto t0 auto tmax auto \
        omodel EXAMPLES/OUTDIR1
EOF

cat > $goodout <<EOF
Name                   = EXAMPLES/4.microlensinject
Microlens_f0_0         = 7.242316197338e-05
Microlens_f1_0         = 7.5541525219661e-05
Microlens_u0_0         = 7.242316197338e-05
Microlens_t0_0         = 3.9109521538222
Microlens_tmax_0       = 53740.494617109
Microlens_chi2perdof_0 = 4.4674961258953

EOF

$VARTOOLS -i EXAMPLES/4.microlensinject -oneline \
    -microlens f0 auto f1 auto u0 auto t0 auto tmax auto \
        omodel EXAMPLES/OUTDIR1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -nonlinfit example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -nonlinfit example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 \
    -stats t min,max \
    -expr t1=STATS_t_MIN_0 \
    -expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \
    -expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \
    -nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \
        'c=(t1+Dt*0.3):(0.1*Dt),d=(Dt*0.1):(0.1*Dt)' \
        linfit a,b amoeba omodel EXAMPLES/OUTDIR1/ \
    -oneline
EOF

cat > $goodout <<EOF
Name                     = EXAMPLES/3
STATS_t_MIN_0            = 53725.173920000001
STATS_t_MAX_0            = 53756.281021000003
Nonlinfit_HasConverged_4 = 1
Nonlinfit_c_BestFit_4    = 53731.405628420733
Nonlinfit_c_Err_4        = 0.00090197861656010614
Nonlinfit_d_BestFit_4    = 1.4997722086319381
Nonlinfit_d_Err_4        = 0.0011353562531547948
Nonlinfit_a_BestFit_4    = 10.167328282857785
Nonlinfit_a_Err_4        = 1.7957444863607521e-05
Nonlinfit_b_BestFit_4    = 0.10065045055297883
Nonlinfit_b_Err_4        = 4.8262250106865319e-05
Nonlinfit_BestFit_Chi2_4 = 90077.830712742856

EOF

$VARTOOLS -i EXAMPLES/3 \
    -stats t min,max \
    -expr t1=STATS_t_MIN_0 \
    -expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \
    -expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \
    -nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \
        'c=(t1+Dt*0.3):(0.1*Dt),d=(Dt*0.1):(0.1*Dt)' \
        linfit a,b amoeba omodel EXAMPLES/OUTDIR1/ \
    -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -nonlinfit example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -nonlinfit example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 \
    -randseed 1 \
    -stats t min,max \
    -expr t1=STATS_t_MIN_0 \
    -expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \
    -expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \
    -nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \
        'a=10.167:0.0002,b=0.1:0.0008,c=(t1+Dt*0.2):(0.005),d=(Dt*0.05):(0.016)' \
        mcmc Nlinkstotal 10000 outchains EXAMPLES/OUTDIR1/ \
    -oneline
EOF

cat > $goodout <<EOF
Name                     = EXAMPLES/3
STATS_t_MIN_0            = 53725.173920000001
STATS_t_MAX_0            = 53756.281021000003
Nonlinfit_a_BestFit_4    = 10.16733595824326
Nonlinfit_b_BestFit_4    = 0.10068199797234947
Nonlinfit_c_BestFit_4    = 53731.405627923843
Nonlinfit_d_BestFit_4    = 1.4997732382607727
Nonlinfit_BestFit_Chi2_4 = 90077.83071503813
Nonlinfit_a_MEDIAN_4     = 10.167339206527549
Nonlinfit_a_STDDEV_4     = 2.1634712351800478e-05
Nonlinfit_b_MEDIAN_4     = 0.10071313795015273
Nonlinfit_b_STDDEV_4     = 7.5501039462404719e-05
Nonlinfit_c_MEDIAN_4     = 53731.40592424116
Nonlinfit_c_STDDEV_4     = 0.0010974050260291169
Nonlinfit_d_MEDIAN_4     = 1.4987003116627102
Nonlinfit_d_STDDEV_4     = 0.001632159684888244

EOF

$VARTOOLS -i EXAMPLES/3 \
    -randseed 1 \
    -stats t min,max \
    -expr t1=STATS_t_MIN_0 \
    -expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \
    -expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \
    -nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \
        'a=10.167:0.0002,b=0.1:0.0008,c=(t1+Dt*0.2):(0.005),d=(Dt*0.05):(0.016)' \
        mcmc Nlinkstotal 10000 outchains EXAMPLES/OUTDIR1/ \
    -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -o example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -o example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -LS 0.1 100.0 0.1 1 0 \
    -expr phase=t \
    -changevariable t phase \
    -Phase ls \
    -o EXAMPLES/OUTDIR1 \
        nameformat "file_%s_%05d_simout.txt" \
        columnformat "t:%11.5f,phase:%8.5f,mag:%7.4f,err:%7.4f"
EOF

cat > $goodout <<EOF
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_Periodogram_Value_1_0 LS_SNR_1_0
EXAMPLES/1    77.76775250 -5710.91013    0.99392   38.81513
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308
EXAMPLES/3    18.29829471  -26.09202    0.03822   11.87999
EXAMPLES/4     0.99383709  -90.59677    0.12488   11.63276
EXAMPLES/5     7.06979568  -86.82830    0.10042   11.03518
EXAMPLES/6    22.21935786  -59.44675    0.07034    9.96198
EXAMPLES/7     0.14747089  -12.56404    0.01933    4.85468
EXAMPLES/8     0.93696087  -79.26500    0.09115    9.35161
EXAMPLES/9     7.06979568  -48.29326    0.05781   11.13720
EXAMPLES/10     0.96906857  -52.91818    0.06257    9.66926
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -LS 0.1 100.0 0.1 1 0 \
    -expr phase=t \
    -changevariable t phase \
    -Phase ls \
    -o EXAMPLES/OUTDIR1 \
        nameformat "file_%s_%05d_simout.txt" \
        columnformat "t:%11.5f,phase:%8.5f,mag:%7.4f,err:%7.4f" \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Phase example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Phase example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -header \
    -Phase fix 1.2354 \
    -o EXAMPLES/OUTDIR1/2.phase.txt
EOF

cat > $goodout <<EOF
#Name
EXAMPLES/2
EOF

$VARTOOLS -i EXAMPLES/2 -header \
    -Phase fix 1.2354 \
    -o EXAMPLES/OUTDIR1/2.phase.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Phase example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Phase example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -Phase bls T0 bls 0.5 \
    -o EXAMPLES/OUTDIR1/3.phase.txt \
    -binlc 1 nbins 200 0 \
    -o EXAMPLES/OUTDIR1/3.phasebin.txt
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12312625
BLS_Tc_1_0                   = 53727.297046247397
BLS_SN_1_0                   =  38.39425
BLS_SR_1_0                   =   0.00237
BLS_SDE_1_0                  =   4.77204
BLS_Depth_1_0                =   0.01136
BLS_Qtran_1_0                =   0.03000
BLS_i1_1_0                   =   0.98500
BLS_i2_1_0                   =   1.01000
BLS_deltaChi2_1_0            = -24130.93833
BLS_fraconenight_1_0         =   0.42662
BLS_Npointsintransit_1_0     =   146
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   106
BLS_Npointsaftertransit_1_0  =   120
BLS_Rednoise_1_0             =   0.00156
BLS_Whitenoise_1_0           =   0.00490
BLS_SignaltoPinknoise_1_0    =  12.89679
BLS_Period_invtransit_0      =     1.14599569
BLS_deltaChi2_invtransit_0   = -3289.67397
BLS_MeanMag_0                =  10.16740

EOF

$VARTOOLS -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -Phase bls T0 bls 0.5 \
    -o EXAMPLES/OUTDIR1/3.phase.txt \
    -binlc 1 nbins 200 0 \
    -o EXAMPLES/OUTDIR1/3.phasebin.txt \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -python example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -python example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list \
             -inputlcformat t:1,mag:2,err:3 \
             -python 'b = numpy.var(mag)' \
                      invars mag outvars b outputcolumns b
EOF

cat > $goodout <<EOF
EXAMPLES/1 0.025422461711037084
EXAMPLES/2 0.0013420988067623005
EXAMPLES/3 2.3966645306408949e-05
EXAMPLES/4 4.3733138204733634e-06
EXAMPLES/5 8.2971716866526236e-06
EXAMPLES/6 4.3664615059428104e-06
EXAMPLES/7 1.216345131566495e-05
EXAMPLES/8 5.0623773543353351e-06
EXAMPLES/9 3.4861868515750583e-06
EXAMPLES/10 5.5813996936871234e-06
EOF

$VARTOOLS -l EXAMPLES/lc_list \
             -inputlcformat t:1,mag:2,err:3 \
             -python 'b = numpy.var(mag)' \
                      invars mag outvars b outputcolumns b \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -python example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -python example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 \
    -LS 0.1 100. 0.1 1 0 \
    -if 'Log10_LS_Prob_1_0<-100' \
        -Phase ls phasevar ph \
        -python 'plotlc(Name,"EXAMPLES/",t,ph,mag,LS_Period_1_0)' \
            init file EXAMPLES/plotlc.py \
    -fi
EOF

cat > $goodout <<EOF
EXAMPLES/1    77.76775250 -5710.91013    0.99392   38.81513
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308
EXAMPLES/3    18.29829471  -26.09202    0.03822   11.87999
EXAMPLES/4     0.99383709  -90.59677    0.12488   11.63276
EXAMPLES/5     7.06979568  -86.82830    0.10042   11.03518
EXAMPLES/6    22.21935786  -59.44675    0.07034    9.96198
EXAMPLES/7     0.14747089  -12.56404    0.01933    4.85468
EXAMPLES/8     0.93696087  -79.26500    0.09115    9.35161
EXAMPLES/9     7.06979568  -48.29326    0.05781   11.13720
EXAMPLES/10     0.96906857  -52.91818    0.06257    9.66926
EOF

$VARTOOLS -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 \
    -LS 0.1 100. 0.1 1 0 \
    -if 'Log10_LS_Prob_1_0<-100' \
        -Phase ls phasevar ph \
        -python 'plotlc(Name,"EXAMPLES/",t,ph,mag,LS_Period_1_0)' \
            init file EXAMPLES/plotlc.py \
    -fi \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -python example 3
testnumber=$((testnumber+1))
echo "$testnumber. Testing -python example 3" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 \
    -LS 0.1 100. 0.1 1 0 \
    -Phase ls phasevar ph \
    -python \
        'for i in range(0,len(mag)):
            plotlc(Name[i],"EXAMPLES/",t[i],ph[i],mag[i],LS_Period_1_0[i])' \
        init file EXAMPLES/plotlc.py \
        process_all_lcs
EOF

cat > $goodout <<EOF
EXAMPLES/1    77.76775250 -5710.91013    0.99392   38.81513
EXAMPLES/2     1.23440877 -4000.59209    0.99619   45.98308
EXAMPLES/3    18.29829471  -26.09202    0.03822   11.87999
EXAMPLES/4     0.99383709  -90.59677    0.12488   11.63276
EXAMPLES/5     7.06979568  -86.82830    0.10042   11.03518
EXAMPLES/6    22.21935786  -59.44675    0.07034    9.96198
EXAMPLES/7     0.14747089  -12.56404    0.01933    4.85468
EXAMPLES/8     0.93696087  -79.26500    0.09115    9.35161
EXAMPLES/9     7.06979568  -48.29326    0.05781   11.13720
EXAMPLES/10     0.96906857  -52.91818    0.06257    9.66926
EOF

$VARTOOLS -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 \
    -LS 0.1 100. 0.1 1 0 \
    -Phase ls phasevar ph \
    -python \
        'for i in range(0,len(mag)):
            plotlc(Name[i],"EXAMPLES/",t[i],ph[i],mag[i],LS_Period_1_0[i])' \
        init file EXAMPLES/plotlc.py \
        process_all_lcs \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -resample example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -resample example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -resample linear \
    -o EXAMPLES/2.resample.example1
EOF

cat > $goodout <<EOF
EXAMPLES/2
EOF

$VARTOOLS -i EXAMPLES/2 -resample linear \
    -o EXAMPLES/2.resample.example1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -resample example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -resample example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/2 -resample splinemonotonic \
    tstart fix 53726 tstop fix 53756 Npoints fix 1000 \
    -o EXAMPLES/2.resample.example2
EOF

cat > $goodout <<EOF
EXAMPLES/2
EOF

$VARTOOLS -i EXAMPLES/2 -resample splinemonotonic \
    tstart fix 53726 tstop fix 53756 Npoints fix 1000 \
    -o EXAMPLES/2.resample.example2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -resample example 3
testnumber=$((testnumber+1))
echo "$testnumber. Testing -resample example 3" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/4 -resample linear \
    file fix EXAMPLES/8 \
    -o EXAMPLES/4.resample.example3
EOF

cat > $goodout <<EOF
EXAMPLES/4
EOF

$VARTOOLS -i EXAMPLES/4 -resample linear \
    file fix EXAMPLES/8 \
    -o EXAMPLES/4.resample.example3 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -resample example 4
testnumber=$((testnumber+1))
echo "$testnumber. Testing -resample example 4" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/1 -resample splinemonotonic \
    tstart fix 53725 tstop fix 53757 delt fix 0.001 \
    gaps percentile_sep 80 bspline nbreaks 15 order 3 \
    extrap nearest \
    -o EXAMPLES/1.resample.example4
EOF

cat > $goodout <<EOF
EXAMPLES/1
EOF

$VARTOOLS -i EXAMPLES/1 -resample splinemonotonic \
    tstart fix 53725 tstop fix 53757 delt fix 0.001 \
    gaps percentile_sep 80 bspline nbreaks 15 order 3 \
    extrap nearest \
    -o EXAMPLES/1.resample.example4 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -rescalesig example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -rescalesig example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/4 -oneline \
    -chi2 -rescalesig -chi2
EOF

cat > $goodout <<EOF
Name                 = EXAMPLES/4
Chi2_0               =      5.19874
Weighted_Mean_Mag_0  =  10.35137
SigmaRescaleFactor_1 =   0.43858
Chi2_2               =      1.00000
Weighted_Mean_Mag_2  =  10.35137

EOF

$VARTOOLS -i EXAMPLES/4 -oneline \
    -chi2 -rescalesig -chi2 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


# -restorelc example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -restorelc example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header -numbercolumns \
    -nobuffer -parallel 4 \
    -savelc \
    -clip 5.0 1 \
    -savelc \
    -LS 0.1 100. 0.1 3 0 clip 5. 1 \
    -aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
    -restorelc 1 \
    -clip 10.0 1 \
    -BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \
    -restorelc 2 \
    -changeerror \
    -autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/ | \
sort -k 1,1
EOF

cat > $goodout <<EOF
#1_Name 2_Nclip_1 3_LS_Period_1_3 4_Log10_LS_Prob_1_3 5_LS_Periodogram_Value_1_3 6_LS_SNR_1_3 7_LS_Period_2_3 8_Log10_LS_Prob_2_3 9_LS_Periodogram_Value_2_3 10_LS_SNR_2_3 11_LS_Period_3_3 12_Log10_LS_Prob_3_3 13_LS_Periodogram_Value_3_3 14_LS_SNR_3_3 15_Period_1_4 16_AOV_1_4 17_AOV_SNR_1_4 18_AOV_NEG_LN_FAP_1_4 19_Nclip_6 20_BLS_Period_1_7 21_BLS_Tc_1_7 22_BLS_SN_1_7 23_BLS_SR_1_7 24_BLS_SDE_1_7 25_BLS_Depth_1_7 26_BLS_Qtran_1_7 27_BLS_i1_1_7 28_BLS_i2_1_7 29_BLS_deltaChi2_1_7 30_BLS_fraconenight_1_7 31_BLS_Npointsintransit_1_7 32_BLS_Ntransits_1_7 33_BLS_Npointsbeforetransit_1_7 34_BLS_Npointsaftertransit_1_7 35_BLS_Rednoise_1_7 36_BLS_Whitenoise_1_7 37_BLS_SignaltoPinknoise_1_7 38_BLS_Period_2_7 39_BLS_Tc_2_7 40_BLS_SN_2_7 41_BLS_SR_2_7 42_BLS_SDE_2_7 43_BLS_Depth_2_7 44_BLS_Qtran_2_7 45_BLS_i1_2_7 46_BLS_i2_2_7 47_BLS_deltaChi2_2_7 48_BLS_fraconenight_2_7 49_BLS_Npointsintransit_2_7 50_BLS_Ntransits_2_7 51_BLS_Npointsbeforetransit_2_7 52_BLS_Npointsaftertransit_2_7 53_BLS_Rednoise_2_7 54_BLS_Whitenoise_2_7 55_BLS_SignaltoPinknoise_2_7 56_BLS_Period_invtransit_7 57_BLS_deltaChi2_invtransit_7 58_BLS_MeanMag_7 59_Mean_Mag_9 60_RMS_9 61_Npoints_9
EXAMPLES/1     0    77.76775250 -5710.91013    0.99392   38.81513     0.97821072 -5530.60640    0.76156   29.55922     1.01657193 -5490.04096    0.71728   27.79505    34.52470117 46218.43332 3912.55946 7221.76774     0     3.05219780 53725.708054615388  26.66040   0.11243   5.37165   0.35851   0.10000   0.12500   0.22000 -53161091.85017   0.35500   285     4   309   448   0.11815   0.12237   6.02371     1.01277747 53725.829693408894  23.40282   0.09809   4.50317   0.36748   0.10500   0.59500   0.69500 -40461972.19159   0.35470   227     4    63   268   0.10708   0.12692   6.78011    20.00000000 -23119393.58400  10.24430  10.24745   0.15944  3122
EXAMPLES/10     0     0.96906857  -52.91818    0.06257    9.66926    22.21935786  -50.53617    0.05981    9.19953     6.48064604  -30.84071    0.03729    5.36490    22.24937863  37.60833   8.30015 110.84854     0     0.37312764 53725.443504723371  19.58111   0.00040   3.77173   0.00133   0.10500   0.67000   0.77000 -412.30826   0.13440   404    18   428   397   0.00082   0.00233   5.88202     4.77975095 53727.922276796293  18.26272   0.00041   4.03736   0.00145   0.09000   0.53000   0.61500 -443.36990   0.40628   413     3   567   424   0.00056   0.00233   4.21955     0.99080932 -758.00480  10.87763  10.87781   0.00236  3974
EXAMPLES/2     0     1.23440877 -4000.59209    0.99619   45.98308     0.55252400 -3750.68873    0.70251   32.20447     0.35510389 -2989.13043    0.24101   10.55229     1.23583047 9274.25316 654.12320 4979.79170     0     4.94437027 53725.45822129061  29.36340   0.02676   6.83057   0.06573   0.10500   0.00500   0.10500 -3122684.45271   0.24901   871     6    26   107   0.01976   0.02549   8.10016     0.54996672 53725.455777945892  17.41964   0.01592   3.48461   0.05536   0.10500   0.46000   0.56000 -1105077.99512   0.20882   340    13   336   354   0.03069   0.03369   6.35792     5.22727867 -2033854.28019  10.11178  10.11802   0.03663  3313
EXAMPLES/3     2    18.29829471  -27.22364    0.03970   12.17097     1.14786351  -23.40687    0.03457   10.48313     1.05806466  -20.65933    0.03088    9.27350    17.56461059  35.43692  10.52015 103.28794     0    16.40794224 53737.602936245486  34.21975   0.00085   5.49747   0.00246   0.10500   0.70500   0.80500 -3094.85163   0.43976   398     2   165   211   0.00079   0.00483   4.02468     0.94160522 53726.019010685653  30.41691   0.00074   4.31994   0.00264   0.10500   0.84500   0.94500 -2364.56373   0.33407   310     8   360   221   0.00122   0.00485   5.14786    15.40203327 -2404.93187  10.16684  10.16674   0.00486  3415
EXAMPLES/4     6     0.99383709  -84.73842    0.11774   12.31519     1.16943989  -77.73898    0.10777   11.20711     0.93414718  -73.54166    0.10184   10.54794    34.70790113  92.61914  17.51651 276.74730     0     2.32488927 53726.812966932586  30.95354   0.00093   4.33228   0.00569   0.05000   0.68000   0.72500 -2977.65425   1.00000    61     1    98     8   0.00071   0.00195   7.53922     1.98505107 53725.580855469459  30.25312   0.00092   4.22373   0.00651   0.04000   0.18500   0.22000 -2870.35671   0.99585    50     2   253     8   0.00076   0.00197  10.72951    13.53502538 -934.10617  10.35137  10.35141   0.00204  3221
EXAMPLES/5     1     7.06979568  -86.77716    0.10039   11.00843     0.87625637  -73.24175    0.08453    9.12747     1.01990495  -51.62791    0.05972    6.18553     8.37612050  71.69357  10.25483 217.96768     0    20.00000000 53745.073920000003  18.12676   0.00064   4.90228   0.00173   0.10000   0.94500   1.04000 -1632.92132   0.43363   698     2   298   583   0.00073   0.00282   3.29460     0.87900978 53725.929868414547  16.45821   0.00056   3.95272   0.00180   0.10000   0.81000   0.90500 -1254.84860   0.20057   432     9   233   439   0.00112   0.00284   4.50717    20.00000000 -2993.85612  10.43932  10.43962   0.00287  3902
EXAMPLES/6     1    22.21935786  -59.61922    0.07054    9.98686     0.96009571  -58.05122    0.06871    9.70196     0.48987561  -30.86340    0.03740    4.84451    22.85589452  46.00325   9.33718 137.68537     0     0.37298915 53725.448999501263  21.62511   0.00039   3.90357   0.00131   0.10500   0.68500   0.78500 -549.87354   0.15460   394    18   423   393   0.00081   0.00206   6.02355     2.86102607 53727.984878110394  20.09105   0.00040   4.04948   0.00142   0.10500   0.93000   1.03000 -573.18316   0.47471   412     5   392   608   0.00054   0.00205   5.45946     1.92423528 -802.58344  10.52748  10.52761   0.00208  3932
EXAMPLES/7     1     0.14747089  -12.46647    0.01922    4.80537     0.12327234  -10.31697    0.01649    3.96906     0.68505751   -9.94923    0.01602    3.82622     1.80654114  19.70992   5.85501  52.63608     0     3.63778582 53726.552674147853  34.15570   0.00052   4.92465   0.00409   0.03500   0.36000   0.39000 -916.31144   0.46930   107     3   229    87   0.00082   0.00342   7.08279     0.37591921 53725.362670617018  30.73327   0.00046   3.80264   0.00337   0.01500   0.48000   0.49000 -710.54027   0.35387    66    13    78    64   0.00135   0.00347   5.93376     3.45042962 -1325.03039  10.56951  10.56966   0.00348  3625
EXAMPLES/8     0     0.93696087  -79.26500    0.09115    9.35161    16.37215842  -77.78479    0.08943    9.16090     1.11896047  -61.96259    0.07120    7.14271    12.29039040 129.90291  17.81185 389.85231     0     3.73236282 53725.425854490481  18.73923   0.00057   4.91213   0.00208   0.10500   0.01500   0.11500 -1101.88000   0.32552   519     5   195   397   0.00095   0.00216   4.79798     0.53079447 53725.602536531653  12.08756   0.00037   2.39737   0.00115   0.10500   0.75500   0.85500 -455.24278   0.24208   425    11   371   452   0.00120   0.00223   3.05435    13.71981339 -2076.09183  10.61132  10.61152   0.00225  3957
EXAMPLES/9     0     7.06979568  -48.29326    0.05781   11.13720     0.87625637  -42.47599    0.05105    9.73682     0.76430224  -25.04394    0.03107    5.59685     6.97310721  37.95931  10.42873 111.95866     0     7.48568220 53730.095756047165  28.05210   0.00039   5.83543   0.00122   0.08500   0.61500   0.69500 -449.11329   0.37314   481     3   239   397   0.00041   0.00184   4.85572     0.87594502 53725.920663129713  23.60989   0.00033   4.36861   0.00142   0.09500   0.80500   0.89500 -316.76326   0.44009   240     7   241   418   0.00047   0.00184   6.60947     3.40368315 -493.75117  10.73129  10.73139   0.00187  3954
EOF

$VARTOOLS -l EXAMPLES/lc_list -header -numbercolumns \
    -nobuffer -parallel 4 \
    -savelc \
    -clip 5.0 1 \
    -savelc \
    -LS 0.1 100. 0.1 3 0 clip 5. 1 \
    -aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
    -restorelc 1 \
    -clip 10.0 1 \
    -BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \
    -restorelc 2 \
    -changeerror \
    -autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/ | \
sort -k 1,1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -restricttimes example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -restricttimes example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 -stats t min,max \
             -restricttimes JDrange 53740 53750 \
             -stats t min,max -oneline
EOF

cat > $goodout <<EOF
Name                  = EXAMPLES/3
STATS_t_MIN_0         = 53725.173920000001
STATS_t_MAX_0         = 53756.281021000003
RestrictTimes_MinJD_1 = 53740
RestrictTimes_MaxJD_1 = 53750
STATS_t_MIN_2         = 53740.336210000001
STATS_t_MAX_2         = 53745.478681000001

EOF

$VARTOOLS -i EXAMPLES/3 -stats t min,max \
             -restricttimes JDrange 53740 53750 \
             -stats t min,max -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -restricttimes example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -restricttimes example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 -stats mag min,max \
             -restricttimes expr '(mag>10.16311)&&(mag<10.17027)' \
             -stats mag min,max -oneline
EOF

cat > $goodout <<EOF
Name            = EXAMPLES/3
STATS_mag_MIN_0 = 10.141400000000001
STATS_mag_MAX_0 = 10.1921
STATS_mag_MIN_2 = 10.163119999999999
STATS_mag_MAX_2 = 10.170260000000001

EOF

$VARTOOLS -i EXAMPLES/3 -stats mag min,max \
             -restricttimes expr '(mag>10.16311)&&(mag<10.17027)' \
             -stats mag min,max -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -restricttimes example 3
testnumber=$((testnumber+1))
echo "$testnumber. Testing -restricttimes example 3" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 -stats mag pct20.0,pct80.0 \
             -restricttimes expr \
               '(mag>STATS_mag_PCT20_00_0)&&(mag<STATS_mag_PCT80_00_0)' \
             -stats mag min,max -oneline
EOF

cat > $goodout <<EOF
Name                 = EXAMPLES/3
STATS_mag_PCT20.00_0 = 10.16311
STATS_mag_PCT80.00_0 = 10.17027
STATS_mag_MIN_2      = 10.163119999999999
STATS_mag_MAX_2      = 10.170260000000001

EOF

$VARTOOLS -i EXAMPLES/3 -stats mag pct20.0,pct80.0 \
             -restricttimes expr \
               '(mag>STATS_mag_PCT20_00_0)&&(mag<STATS_mag_PCT80_00_0)' \
             -stats mag min,max -oneline \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -restricttimes example 4
testnumber=$((testnumber+1))
echo "$testnumber. Testing -restricttimes example 4" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -oneline \
             -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
                  0 0 0 fittrap nobinnedrms \
             -expr 'minph=0.5-BLS_Qtran_1_0/2.0' \
             -expr 'maxph=0.5+BLS_Qtran_1_0/2.0' \
             -expr 'ph=t' \
             -changevariable t ph \
             -Phase bls T0 bls 0.5 \
             -restricttimes exclude JDrangebylc expr minph expr maxph \
             -o EXAMPLES/3.cliptransit columnformat t,mag,err,ph
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12319314
BLS_Tc_1_0                   = 53727.297609654204
BLS_SN_1_0                   =   7.21628
BLS_SR_1_0                   =   0.00238
BLS_SDE_1_0                  =   6.31221
BLS_Depth_1_0                =   0.01225
BLS_Qtran_1_0                =   0.04041
BLS_Qingress_1_0             =   0.26351
BLS_OOTmag_1_0               =  10.16692
BLS_i1_1_0                   =   0.98003
BLS_i2_1_0                   =   1.02044
BLS_deltaChi2_1_0            = -24208.38967
BLS_fraconenight_1_0         =   0.42879
BLS_Npointsintransit_1_0     =   180
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   142
BLS_Npointsaftertransit_1_0  =   163
BLS_Rednoise_1_0             =   0.00147
BLS_Whitenoise_1_0           =   0.00489
BLS_SignaltoPinknoise_1_0    =  14.91256
BLS_Period_invtransit_0      =     1.14590298
BLS_deltaChi2_invtransit_0   = -3262.95028
BLS_MeanMag_0                =  10.16740
RestrictTimes_MinJD_6        = 0.47979532860751506
RestrictTimes_MaxJD_6        = 0.52020467139248494

EOF

$VARTOOLS -i EXAMPLES/3.transit -oneline \
             -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
                  0 0 0 fittrap nobinnedrms \
             -expr 'minph=0.5-BLS_Qtran_1_0/2.0' \
             -expr 'maxph=0.5+BLS_Qtran_1_0/2.0' \
             -expr 'ph=t' \
             -changevariable t ph \
             -Phase bls T0 bls 0.5 \
             -restricttimes exclude JDrangebylc expr minph expr maxph \
             -o EXAMPLES/3.cliptransit columnformat t,mag,err,ph \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -rms example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -rms example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list -header \
    -rms
EOF

cat > $goodout <<EOF
#Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0
EXAMPLES/1  10.24745   0.15944   0.00101  3122
EXAMPLES/2  10.11802   0.03663   0.00102  3313
EXAMPLES/3  10.16674   0.00490   0.00104  3417
EXAMPLES/4  10.35142   0.00209   0.00114  3227
EXAMPLES/5  10.43962   0.00288   0.00114  3903
EXAMPLES/6  10.52762   0.00209   0.00121  3933
EXAMPLES/7  10.56966   0.00349   0.00116  3626
EXAMPLES/8  10.61152   0.00225   0.00125  3957
EXAMPLES/9  10.73139   0.00187   0.00133  3954
EXAMPLES/10  10.87781   0.00236   0.00143  3974
EOF

$VARTOOLS -l EXAMPLES/lc_list -header \
    -rms \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -rmsbin example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -rmsbin example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -header -l EXAMPLES/lc_list \
    -rmsbin 5 5.0 10.0 60.0 1440 14400
EOF

cat > $goodout <<EOF
#Name RMSBin_5.0_0 Expected_RMS_Bin_5.0_0 RMSBin_10.0_0 Expected_RMS_Bin_10.0_0 RMSBin_60.0_0 Expected_RMS_Bin_60.0_0 RMSBin_1440.0_0 Expected_RMS_Bin_1440.0_0 RMSBin_14400.0_0 Expected_RMS_Bin_14400.0_0
EXAMPLES/1   0.15943   0.00048   0.15944   0.00035   0.15947   0.00016   0.15879   0.00007   0.13152   0.00003
EXAMPLES/2   0.03659   0.00047   0.03659   0.00035   0.03649   0.00016   0.02723   0.00006   0.00378   0.00003
EXAMPLES/3   0.00292   0.00048   0.00255   0.00035   0.00180   0.00016   0.00125   0.00006   0.00062   0.00003
EXAMPLES/4   0.00146   0.00053   0.00135   0.00040   0.00113   0.00018   0.00092   0.00007   0.00040   0.00003
EXAMPLES/5   0.00180   0.00051   0.00162   0.00038   0.00135   0.00017   0.00093   0.00007   0.00057   0.00003
EXAMPLES/6   0.00130   0.00054   0.00117   0.00040   0.00090   0.00018   0.00059   0.00007   0.00025   0.00003
EXAMPLES/7   0.00208   0.00051   0.00180   0.00038   0.00124   0.00018   0.00047   0.00007   0.00011   0.00003
EXAMPLES/8   0.00151   0.00056   0.00142   0.00041   0.00130   0.00019   0.00096   0.00007   0.00033   0.00003
EXAMPLES/9   0.00105   0.00060   0.00090   0.00044   0.00068   0.00020   0.00047   0.00008   0.00021   0.00003
EXAMPLES/10   0.00139   0.00064   0.00123   0.00047   0.00092   0.00021   0.00065   0.00008   0.00030   0.00003
EOF

$VARTOOLS -header -l EXAMPLES/lc_list \
    -rmsbin 5 5.0 10.0 60.0 1440 14400 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -SoftenedTransit example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -SoftenedTransit example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -SoftenedTransit bls 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1 0
EOF

cat > $goodout <<EOF
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12312625
BLS_Tc_1_0                   = 53727.297046247397
BLS_SN_1_0                   =  38.39425
BLS_SR_1_0                   =   0.00237
BLS_SDE_1_0                  =   4.77204
BLS_Depth_1_0                =   0.01136
BLS_Qtran_1_0                =   0.03000
BLS_i1_1_0                   =   0.98500
BLS_i2_1_0                   =   1.01000
BLS_deltaChi2_1_0            = -24130.93833
BLS_fraconenight_1_0         =   0.42662
BLS_Npointsintransit_1_0     =   146
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   106
BLS_Npointsaftertransit_1_0  =   120
BLS_Rednoise_1_0             =   0.00156
BLS_Whitenoise_1_0           =   0.00490
BLS_SignaltoPinknoise_1_0    =  12.89679
BLS_Period_invtransit_0      =     1.14599569
BLS_deltaChi2_invtransit_0   = -3289.67397
BLS_MeanMag_0                =  10.16740
SoftenedTransit_Period_1     =     2.12322112
SoftenedTransit_T0_1         = 53727.29783160
SoftenedTransit_eta_1        =     0.06171206
SoftenedTransit_cval_1       =   -10.87159958
SoftenedTransit_delta_1      =    -0.01206461
SoftenedTransit_mconst_1     =    10.16686817
SoftenedTransit_chi2perdof_1 =    27.04335183

EOF

$VARTOOLS -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -SoftenedTransit bls 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1 0 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -Starspot example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -Starspot example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3.starspot -oneline \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 0 \
    -Starspot aov 0.0298 0.08745 20. 85. 30. 0. -1 \
        1 0 0 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1/
EOF

cat > $goodout <<EOF
Name                   = EXAMPLES/3.starspot
Period_1_0             =     3.07960303
AOV_1_0                = 2861.35783
AOV_SNR_1_0            = 605.83431
AOV_NEG_LN_FAP_1_0     = 4755.85353
Period_2_0             =     3.08383230
AOV_2_0                = 2559.76138
AOV_SNR_2_0            = 541.88234
AOV_NEG_LN_FAP_2_0     = 4578.45788
Period_3_0             =     3.10356260
AOV_3_0                = 2165.11113
AOV_SNR_3_0            = 458.19878
AOV_NEG_LN_FAP_3_0     = 4314.25066
Period_4_0             =     6.36005672
AOV_4_0                = 1514.99364
AOV_SNR_4_0            = 320.34471
AOV_NEG_LN_FAP_4_0     = 3762.76328
Period_5_0             =     3.19393583
AOV_5_0                = 1026.91829
AOV_SNR_5_0            = 216.85084
AOV_NEG_LN_FAP_5_0     = 3185.86267
Starspot_Period_1      =     3.12218969
Starspot_a_1           =   0.02980
Starspot_b_1           =   0.08745
Starspot_alpha_1       =  22.51312
Starspot_inclination_1 =  69.03963
Starspot_chi_1         =  30.00411
Starspot_psi0_1        =   0.00000
Starspot_mconst_1      =  10.16641
Starspot_chi2perdof_1  =  26.58796

EOF

$VARTOOLS -i EXAMPLES/3.starspot -oneline \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 0 \
    -Starspot aov 0.0298 0.08745 20. 85. 30. 0. -1 \
        1 0 0 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1/ \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -stats example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -stats example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/3 \
    -oneline \
    -expr 'mag2=mag+0.01*gauss()' \
    -stats mag,mag2 \
        mean,weightedmean,median,stddev,meddev,medmeddev,MAD,kurtosis,skewness,pct10,pct20,pct80,pct90,max,min,sum
EOF

cat > $goodout <<EOF
Name                      = EXAMPLES/3
STATS_mag_MEAN_1          = 10.166743412350014
STATS_mag_WEIGHTEDMEAN_1  = 10.166840511593556
STATS_mag_MEDIAN_1        = 10.16667
STATS_mag_STDDEV_1        = 0.0048962905656505422
STATS_mag_MEDDEV_1        = 0.0048968410484820264
STATS_mag_MEDMEDDEV_1     = 0.0028399999999990655
STATS_mag_MAD_1           = 0.0042117199999986143
STATS_mag_KURTOSIS_1      = 4.9831502545166657
STATS_mag_SKEWNESS_1      = 0.14038264034605824
STATS_mag_PCT10.00_1      = 10.160911997073457
STATS_mag_PCT20.00_1      = 10.16311
STATS_mag_PCT80.00_1      = 10.17027
STATS_mag_PCT90.00_1      = 10.172283998536729
STATS_mag_MAX_1           = 10.1921
STATS_mag_MIN_1           = 10.141400000000001
STATS_mag_SUM_1           = 34739.762240000091
STATS_mag2_MEAN_1         = 10.166882723132543
STATS_mag2_WEIGHTEDMEAN_1 = 10.166911774719013
STATS_mag2_MEDIAN_1       = 10.166930994622186
STATS_mag2_STDDEV_1       = 0.011231063371529723
STATS_mag2_MEDDEV_1       = 0.01123116713766443
STATS_mag2_MEDMEDDEV_1    = 0.0076167519498877567
STATS_mag2_MAD_1          = 0.011295643141683544
STATS_mag2_KURTOSIS_1     = 2.9418281705218856
STATS_mag2_SKEWNESS_1     = 0.001194929476051031
STATS_mag2_PCT10.00_1     = 10.152328371334622
STATS_mag2_PCT20.00_1     = 10.157470019145888
STATS_mag2_PCT80.00_1     = 10.176434258298739
STATS_mag2_PCT90.00_1     = 10.18130759497164
STATS_mag2_MAX_1          = 10.208880822707791
STATS_mag2_MIN_1          = 10.129160734154709
STATS_mag2_SUM_1          = 34740.238264943939

EOF

$VARTOOLS -i EXAMPLES/3 \
    -oneline \
    -expr 'mag2=mag+0.01*gauss()' \
    -stats mag,mag2 \
        mean,weightedmean,median,stddev,meddev,medmeddev,MAD,kurtosis,skewness,pct10,pct20,pct80,pct90,max,min,sum \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -SYSREM example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -SYSREM example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/trendlist_tfa -header \
    -rms \
    -SYSREM 2 1 EXAMPLES/3 5. 5. 8. 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1/sysrem.trends 1 \
    -rms
EOF

cat > $goodout <<EOF
#Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0 SYSREM_MeanMag_1 SYSREM_Trend_0_Coeff_1 SYSREM_Trend_1_Coeff_1 SYSREM_Trend_2_Coeff_1 SYSREM_RMS_1 Mean_Mag_2 RMS_2 Expected_RMS_2 Npoints_2
EXAMPLES/4  10.35142   0.00209   0.00114  3227  10.35130  -0.83218  -0.86967  -0.03662   0.00087  10.35137   0.00116   0.00114  3227
EXAMPLES/5  10.43962   0.00288   0.00114  3903  10.43976  -0.03154  -0.05758  -0.02459   0.00177  10.43962   0.00197   0.00114  3903
EXAMPLES/6  10.52762   0.00209   0.00121  3933  10.52774   0.28190   0.30516   0.02490   0.00104  10.52762   0.00115   0.00121  3933
EXAMPLES/7  10.56966   0.00349   0.00116  3626  10.56988   0.14423   0.11560  -0.02852   0.00239  10.56969   0.00254   0.00116  3626
EXAMPLES/8  10.61152   0.00225   0.00125  3957  10.61160   0.01083   0.00955   0.00017   0.00178  10.61152   0.00181   0.00125  3957
EXAMPLES/9  10.73139   0.00187   0.00133  3954  10.73148   0.14601   0.13900  -0.00614   0.00124  10.73139   0.00142   0.00133  3954
EXAMPLES/10  10.87781   0.00236   0.00143  3974  10.87792   0.23422   0.25903   0.02672   0.00105  10.87781   0.00120   0.00143  3974
EOF

$VARTOOLS -l EXAMPLES/trendlist_tfa -header \
    -rms \
    -SYSREM 2 1 EXAMPLES/3 5. 5. 8. 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1/sysrem.trends 1 \
    -rms \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -TFA example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -TFA example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list_tfa -oneline -rms \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 1 0 0
EOF

cat > $goodout <<EOF
Name           = EXAMPLES/3.transit
Mean_Mag_0     =  10.16727
RMS_0          =   0.00542
Expected_RMS_0 =   0.00104
Npoints_0      =  3417
TFA_MeanMag_1  =  10.16714
TFA_RMS_1      =   0.00471

EOF

$VARTOOLS -l EXAMPLES/lc_list_tfa -oneline -rms \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 1 0 0 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -TFA_SR example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -TFA_SR example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list_tfa_sr_harm -oneline -rms \
    -LS 0.1 10. 0.1 1 0 \
    -savelc \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 harm 0 0 period ls \
    -o EXAMPLES/OUTDIR1 nameformat 2.test_tfa_sr_harm \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 
EOF

cat > $goodout <<EOF
Name                                  = EXAMPLES/2
Mean_Mag_0                            =  10.11802
RMS_0                                 =   0.03663
Expected_RMS_0                        =   0.00102
Npoints_0                             =  3313
LS_Period_1_1                         =     1.23440877
Log10_LS_Prob_1_1                     = -4000.59209
LS_Periodogram_Value_1_1              =    0.99619
LS_SNR_1_1                            =   45.98308
Killharm_Mean_Mag_3                   =  10.12217
Killharm_Period_1_3                   =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_3  =   0.05008
Killharm_Per1_Fundamental_Coscoeff_3  =  -0.00222
Killharm_Per1_Amplitude_3             =   0.10026
Mean_Mag_4                            =  10.11176
RMS_4                                 =   0.00231
Expected_RMS_4                        =   0.00102
Npoints_4                             =  3313
TFA_MeanMag_6                         =  10.11766
TFA_RMS_6                             =   0.03555
Killharm_Mean_Mag_7                   =  10.12211
Killharm_Period_1_7                   =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_7  =   0.04802
Killharm_Per1_Fundamental_Coscoeff_7  =  -0.00268
Killharm_Per1_Amplitude_7             =   0.09620
Mean_Mag_8                            =  10.11169
RMS_8                                 =   0.00725
Expected_RMS_8                        =   0.00102
Npoints_8                             =  3313
TFA_SR_MeanMag_10                     =  10.11788
TFA_SR_RMS_10                         =   0.03642
Killharm_Mean_Mag_12                  =  10.12210
Killharm_Period_1_12                  =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_12 =   0.04986
Killharm_Per1_Fundamental_Coscoeff_12 =  -0.00237
Killharm_Per1_Amplitude_12            =   0.09984
Mean_Mag_13                           =  10.11166
RMS_13                                =   0.00210
Expected_RMS_13                       =   0.00102
Npoints_13                            =  3313

EOF

$VARTOOLS -l EXAMPLES/lc_list_tfa_sr_harm -oneline -rms \
    -LS 0.1 10. 0.1 1 0 \
    -savelc \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 harm 0 0 period ls \
    -o EXAMPLES/OUTDIR1 nameformat 2.test_tfa_sr_harm \
    -Killharm ls 0 0 0 \
    -rms -restorelc 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -TFA_SR example 2
testnumber=$((testnumber+1))
echo "$testnumber. Testing -TFA_SR example 2" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list_tfa_sr_bin -oneline -rms \
    -aov Nbin 20 0.1 10. 0.1 0.01 1 0 \
    -savelc \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 bin 100 period aov \
    -o EXAMPLES/OUTDIR1 nameformat %s.test_tfa_sr_bin \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1 
EOF

cat > $goodout <<EOF
Name                                  = EXAMPLES/2
Mean_Mag_0                            =  10.11802
RMS_0                                 =   0.03663
Expected_RMS_0                        =   0.00102
Npoints_0                             =  3313
Period_1_1                            =     1.23583047
AOV_1_1                               = 18330.55450
AOV_SNR_1_1                           = 3093.61938
AOV_NEG_LN_FAP_1_1                    = 7633.23914
Killharm_Mean_Mag_3                   =  10.12175
Killharm_Period_1_3                   =     1.23583047
Killharm_Per1_Fundamental_Sincoeff_3  =   0.04514
Killharm_Per1_Fundamental_Coscoeff_3  =   0.02155
Killharm_Per1_Harm_2_Sincoeff_3       =   0.00033
Killharm_Per1_Harm_2_Coscoeff_3       =   0.00006
Killharm_Per1_Harm_3_Sincoeff_3       =   0.00017
Killharm_Per1_Harm_3_Coscoeff_3       =   0.00007
Killharm_Per1_Harm_4_Sincoeff_3       =   0.00015
Killharm_Per1_Harm_4_Coscoeff_3       =  -0.00015
Killharm_Per1_Harm_5_Sincoeff_3       =   0.00012
Killharm_Per1_Harm_5_Coscoeff_3       =   0.00023
Killharm_Per1_Harm_6_Sincoeff_3       =  -0.00011
Killharm_Per1_Harm_6_Coscoeff_3       =  -0.00006
Killharm_Per1_Amplitude_3             =   0.10012
Mean_Mag_4                            =  10.11187
RMS_4                                 =   0.00200
Expected_RMS_4                        =   0.00102
Npoints_4                             =  3313
TFA_MeanMag_6                         =  10.11766
TFA_RMS_6                             =   0.03555
Killharm_Mean_Mag_7                   =  10.12162
Killharm_Period_1_7                   =     1.23583047
Killharm_Per1_Fundamental_Sincoeff_7  =   0.04303
Killharm_Per1_Fundamental_Coscoeff_7  =   0.02062
Killharm_Per1_Harm_2_Sincoeff_7       =   0.00119
Killharm_Per1_Harm_2_Coscoeff_7       =   0.00122
Killharm_Per1_Harm_3_Sincoeff_7       =   0.00029
Killharm_Per1_Harm_3_Coscoeff_7       =  -0.00134
Killharm_Per1_Harm_4_Sincoeff_7       =  -0.00144
Killharm_Per1_Harm_4_Coscoeff_7       =  -0.00052
Killharm_Per1_Harm_5_Sincoeff_7       =   0.00032
Killharm_Per1_Harm_5_Coscoeff_7       =  -0.00036
Killharm_Per1_Harm_6_Sincoeff_7       =  -0.00006
Killharm_Per1_Harm_6_Coscoeff_7       =  -0.00079
Killharm_Per1_Amplitude_7             =   0.09712
Mean_Mag_8                            =  10.11188
RMS_8                                 =   0.00723
Expected_RMS_8                        =   0.00102
Npoints_8                             =  3313
TFA_SR_MeanMag_10                     =  10.11796
TFA_SR_RMS_10                         =   0.03655
Killharm_Mean_Mag_12                  =  10.12173
Killharm_Period_1_12                  =     1.23583047
Killharm_Per1_Fundamental_Sincoeff_12 =   0.04509
Killharm_Per1_Fundamental_Coscoeff_12 =   0.02147
Killharm_Per1_Harm_2_Sincoeff_12      =   0.00034
Killharm_Per1_Harm_2_Coscoeff_12      =   0.00010
Killharm_Per1_Harm_3_Sincoeff_12      =   0.00017
Killharm_Per1_Harm_3_Coscoeff_12      =  -0.00001
Killharm_Per1_Harm_4_Sincoeff_12      =   0.00005
Killharm_Per1_Harm_4_Coscoeff_12      =  -0.00017
Killharm_Per1_Harm_5_Sincoeff_12      =   0.00009
Killharm_Per1_Harm_5_Coscoeff_12      =   0.00014
Killharm_Per1_Harm_6_Sincoeff_12      =  -0.00012
Killharm_Per1_Harm_6_Coscoeff_12      =  -0.00007
Killharm_Per1_Amplitude_12            =   0.09998
Mean_Mag_13                           =  10.11183
RMS_13                                =   0.00183
Expected_RMS_13                       =   0.00102
Npoints_13                            =  3313

Name                                  = EXAMPLES/3.starspot
Mean_Mag_0                            =  10.18727
RMS_0                                 =   0.02402
Expected_RMS_0                        =   0.00104
Npoints_0                             =  3417
Period_1_1                            =     3.07960303
AOV_1_1                               = 2861.35783
AOV_SNR_1_1                           = 605.83431
AOV_NEG_LN_FAP_1_1                    = 4755.85353
Killharm_Mean_Mag_3                   =  10.18395
Killharm_Period_1_3                   =     3.07960303
Killharm_Per1_Fundamental_Sincoeff_3  =  -0.01015
Killharm_Per1_Fundamental_Coscoeff_3  =  -0.02571
Killharm_Per1_Harm_2_Sincoeff_3       =   0.00906
Killharm_Per1_Harm_2_Coscoeff_3       =   0.01009
Killharm_Per1_Harm_3_Sincoeff_3       =  -0.00185
Killharm_Per1_Harm_3_Coscoeff_3       =  -0.00049
Killharm_Per1_Harm_4_Sincoeff_3       =  -0.00401
Killharm_Per1_Harm_4_Coscoeff_3       =   0.00008
Killharm_Per1_Harm_5_Sincoeff_3       =   0.00318
Killharm_Per1_Harm_5_Coscoeff_3       =  -0.00045
Killharm_Per1_Harm_6_Sincoeff_3       =  -0.00149
Killharm_Per1_Harm_6_Coscoeff_3       =  -0.00022
Killharm_Per1_Amplitude_3             =   0.05709
Mean_Mag_4                            =  10.18396
RMS_4                                 =   0.00526
Expected_RMS_4                        =   0.00104
Npoints_4                             =  3417
TFA_MeanMag_6                         =  10.18678
TFA_RMS_6                             =   0.02225
Killharm_Mean_Mag_7                   =  10.18344
Killharm_Period_1_7                   =     3.07960303
Killharm_Per1_Fundamental_Sincoeff_7  =  -0.01043
Killharm_Per1_Fundamental_Coscoeff_7  =  -0.02290
Killharm_Per1_Harm_2_Sincoeff_7       =   0.00982
Killharm_Per1_Harm_2_Coscoeff_7       =   0.00950
Killharm_Per1_Harm_3_Sincoeff_7       =  -0.00044
Killharm_Per1_Harm_3_Coscoeff_7       =  -0.00123
Killharm_Per1_Harm_4_Sincoeff_7       =  -0.00238
Killharm_Per1_Harm_4_Coscoeff_7       =   0.00067
Killharm_Per1_Harm_5_Sincoeff_7       =   0.00269
Killharm_Per1_Harm_5_Coscoeff_7       =  -0.00121
Killharm_Per1_Harm_6_Sincoeff_7       =  -0.00152
Killharm_Per1_Harm_6_Coscoeff_7       =  -0.00013
Killharm_Per1_Amplitude_7             =   0.05327
Mean_Mag_8                            =  10.18380
RMS_8                                 =   0.00805
Expected_RMS_8                        =   0.00104
Npoints_8                             =  3417
TFA_SR_MeanMag_10                     =  10.18717
TFA_SR_RMS_10                         =   0.02379
Killharm_Mean_Mag_12                  =  10.18425
Killharm_Period_1_12                  =     3.07960303
Killharm_Per1_Fundamental_Sincoeff_12 =  -0.00973
Killharm_Per1_Fundamental_Coscoeff_12 =  -0.02591
Killharm_Per1_Harm_2_Sincoeff_12      =   0.00864
Killharm_Per1_Harm_2_Coscoeff_12      =   0.00976
Killharm_Per1_Harm_3_Sincoeff_12      =  -0.00198
Killharm_Per1_Harm_3_Coscoeff_12      =   0.00004
Killharm_Per1_Harm_4_Sincoeff_12      =  -0.00381
Killharm_Per1_Harm_4_Coscoeff_12      =  -0.00022
Killharm_Per1_Harm_5_Sincoeff_12      =   0.00297
Killharm_Per1_Harm_5_Coscoeff_12      =  -0.00065
Killharm_Per1_Harm_6_Sincoeff_12      =  -0.00139
Killharm_Per1_Harm_6_Coscoeff_12      =  -0.00004
Killharm_Per1_Amplitude_12            =   0.05705
Mean_Mag_13                           =  10.18388
RMS_13                                =   0.00449
Expected_RMS_13                       =   0.00104
Npoints_13                            =  3417

EOF

$VARTOOLS -l EXAMPLES/lc_list_tfa_sr_bin -oneline -rms \
    -aov Nbin 20 0.1 10. 0.1 0.01 1 0 \
    -savelc \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 bin 100 period aov \
    -o EXAMPLES/OUTDIR1 nameformat %s.test_tfa_sr_bin \
    -Killharm aov 5 0 0 \
    -rms -restorelc 1  \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -TFA_SR example 3
testnumber=$((testnumber+1))
echo "$testnumber. Testing -TFA_SR example 3" > /dev/stderr

cat > $testc <<EOF
./vartools -l EXAMPLES/lc_list_tfa_sr_decorr -oneline -rms \
    -savelc \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \
        decorr 0 1 1 2 \
        25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 bin 100 \
    -o EXAMPLES/OUTDIR1 nameformat %s.test_tfa_sr_decorr \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1
EOF

cat > $goodout <<EOF
Name                        = EXAMPLES/1
Mean_Mag_0                  =  10.24745
RMS_0                       =   0.15944
Expected_RMS_0              =   0.00101
Npoints_0                   =  3122
Decorr_constant_term_2      =     10.0830375984825
Decorr_constant_term_err_2  =      0.0000325849746
LCColumn_1_coeff_1_2        =      0.0097933162509
LCColumn_1_coeff_err_1_2    =      0.0000059117875
LCColumn_1_coeff_2_2        =      0.0002554062775
LCColumn_1_coeff_err_2_2    =      0.0000001956124
Decorr_chi2_2               =      6.68601
Mean_Mag_3                  =  10.24728
RMS_3                       =   0.00211
Expected_RMS_3              =   0.00101
Npoints_3                   =  3122
TFA_MeanMag_5               =  10.24577
TFA_RMS_5                   =   0.15642
Decorr_constant_term_6      =     10.0894399398392
Decorr_constant_term_err_6  =      0.0000325849746
LCColumn_1_coeff_1_6        =      0.0092745178509
LCColumn_1_coeff_err_1_6    =      0.0000059117875
LCColumn_1_coeff_2_6        =      0.0002510444516
LCColumn_1_coeff_err_2_6    =      0.0000001956124
Decorr_chi2_6               =    945.72347
Mean_Mag_7                  =  10.24429
RMS_7                       =   0.02623
Expected_RMS_7              =   0.00101
Npoints_7                   =  3122
TFA_SR_MeanMag_9            =  10.24745
TFA_SR_RMS_9                =   0.15944
Decorr_constant_term_11     =     10.0830881231895
Decorr_constant_term_err_11 =      0.0000325849746
LCColumn_1_coeff_1_11       =      0.0097704933097
LCColumn_1_coeff_err_1_11   =      0.0000059117875
LCColumn_1_coeff_2_11       =      0.0002562034043
LCColumn_1_coeff_err_2_11   =      0.0000001956124
Decorr_chi2_11              =      5.15430
Mean_Mag_12                 =  10.24732
RMS_12                      =   0.00188
Expected_RMS_12             =   0.00101
Npoints_12                  =  3122

EOF

$VARTOOLS -l EXAMPLES/lc_list_tfa_sr_decorr -oneline -rms \
    -savelc \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1 \
    -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 xycol 2 3 1 0 0 \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1 \
    -TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \
        decorr 0 1 1 2 \
        25.0 xycol 2 3 1 \
        1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \
        0 0.001 100 bin 100 \
    -o EXAMPLES/OUTDIR1 nameformat %s.test_tfa_sr_decorr \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms -restorelc 1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout

# -wwz example 1
testnumber=$((testnumber+1))
echo "$testnumber. Testing -wwz example 1" > /dev/stderr

cat > $testc <<EOF
./vartools -i EXAMPLES/8 -oneline \
    -wwz maxfreq 2.0 freqsamp 0.25 tau0 auto tau1 auto dtau 0.1 \
        outfulltransform EXAMPLES/OUTDIR1/ pm3d \
        outmaxtransform EXAMPLES/OUTDIR1
EOF

cat > $goodout <<EOF
Name                = EXAMPLES/8
MaxWWZ_0            = 345.87167006071132
MaxWWZ_Freq_0       = 0.30645161290322581
MaxWWZ_TShift_0     = 53735.173920000001
MaxWWZ_Power_0      = 243.91899737448088
MaxWWZ_Amplitude_0  = 0.0019655476933700114
MaxWWZ_Neffective_0 = 1651.192187564548
MaxWWZ_AverageMag_0 = 10.611015587079519
Med_WWZ_0           = 135.53839042847693
Med_Freq_0          = 0.20967741935483872
Med_Power_0         = 116.72893440099938
Med_Amplitude_0     = 0.00096335163669205139
Med_Neffective_0    = 2188.0750923512305
Med_AverageMag_0    = 10.611486630808924

EOF

$VARTOOLS -i EXAMPLES/8 -oneline \
    -wwz maxfreq 2.0 freqsamp 0.25 tau0 auto tau1 auto dtau 0.1 \
        outfulltransform EXAMPLES/OUTDIR1/ pm3d \
        outmaxtransform EXAMPLES/OUTDIR1 \
> $testout

lastcode=$?

if (( $lastcode != 0 )) ; then
    ReportVartoolsError $testnumber $testc $testout $goodout $lastcode
fi

CompareOutput $testnumber $testc $testout $goodout


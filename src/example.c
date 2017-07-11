/*     This file is part of VARTOOLS version 1.31                      */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

/* This function outputs examples for various commands */

void example(char *c, ProgramData *p)
{
  int commandfound = 0;
  int l = 0;
  OutText s;
  int i;
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  if(!strncmp(c,"-addnoise",9) && strlen(c) == 9)
    {
#ifdef _HAVE_GSL
      printtostring(&s,
		    "\nExample 1.\n\n");
      printtostring(&s,
		    "gawk '{print $1, 0., 0.005}' EXAMPLES/1 | \\\n");
      printtostring(&s,
		    "vartools -i - -header -randseed 1 \\\n");
      printtostring(&s,
		    "\t-addnoise wavelet gamma fix 0.99 sig_red fix 0.005 sig_white fix 0.005 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/noisesim.txt\n\n");
      printtostring(&s,
		    "Simulate a light curve with time-correlated noise and the same time sampling as the light curve EXAMPLES/1. The red-noise component has power spectral density proportional 1/f^0.99 and has standard deviation 0.005. The white-noise component has standard deviation 0.005. The simulated light curve with the noise added is output to the file EXAMPLES/OUTDIR/noisesim.txt. Use different values for -randseed to simulate different light curves.\n\n");
      printtostring(&s,
		    "Example 2.\n");
#endif
      printtostring(&s,
		    "\ngawk '{print $1, 0., 0.005}' EXAMPLES/1 | \\\n");
      printtostring(&s,
		    "vartools -i - -header -randseed 1 \\\n");
      printtostring(&s,
		    "\t-addnoise squareexp rho fix 0.01 sig_red fix 0.005 sig_white fix 0.001 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/noisesim.txt\n\n");
      printtostring(&s,
		    "Simulate a light curve with time-correlated noise and the same time sampling as the light curve EXAMPLES/1. We use a squared-exponential model for the red-noise component with a correlation time-scale of 0.01 days and standard deviation 0.005 mag. An additional white-noise component is included with a standard deviation of 0.001. The simulated light curve with the noise added is output to the file EXAMPLES/OUTDIR/noisesim.txt. Use different values for -randseed to simulate different light curves.\n\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-alarm",6) && strlen(c) == 6)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header -alarm\n\n");
      printtostring(&s,
		    "Computes the alarm variability statistic for the light curve EXAMPLES/1.\n\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-aov",4) && strlen(c) == 4)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-aov Nbin 20 0.1 10. 0.1 0.01 5 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\twhiten clip 5. 1\n\n");
      printtostring(&s,
		    "Runs the phase-binning AoV period-finding algorithm on the light curve EXAMPLES/2. 20 phase-bins are used. Periods between 0.1 and 10.0 days are searched. The coarse search is done at a frequency resolution of 0.1/T (T is the time-span of the lc, 31.1d in this case). The fine search around the peaks is done at a frequency resolution of 0.01/T. The top 5 peaks are identified, between each cycle the best-fit signal is removed and the periodogram is regenerated. The periodogram is output to the directory EXAMPLES/OUTDIR1. The filename will be 2.aov. An iterative 5-sigma clipping is applied when identifying peaks in the periodogram.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-aov_harm",9) && strlen(c) == 9)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-aov_harm 1 0.1 10. 0.1 0.01 2 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\twhiten clip 5. 1\n\n");
      printtostring(&s,
		    "Runs the harmonic-fitting AoV period-finding algorithm on the light curve EXAMPLES/2. 1 harmonic is used (i.e. the model is a simple sine-curve). Periods between 0.1 and 10.0 days are searched. The coarse search is done at a frequency resolution of 0.1/T (T is the time-span of the lc, 31.1d in this case). The fine search around the peaks is done at a frequency resolution of 0.01/T. The top 2 peaks are identified, between each cycle the best-fit signal is removed and the periodogram is regenerated. The periodogram is output to the directory EXAMPLES/OUTDIR1. The filename will be 2.aov_harm. An iterative 5-sigma clipping is applied when identifying peaks in the periodogram.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-autocorrelation",16) && strlen(c) == 16)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-autocorrelation 0.0 10. 0.05 EXAMPLES/OUTDIR1\n\n");
      printtostring(&s,
		    "Compute the discrete auto-correlation function (DACF) of the light curve EXAMPLES/2. The DACF is calculated between time-lags of 0 and 10.0 days with a time-step of 0.05 days. It is output to the directory EXAMPLES/OUTDIR1 with the filename 2.autocorr\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-binlc",6) && strlen(c) == 6)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-binlc median binsize 0.01 tcenter \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/2.bin.txt\n\n");
      printtostring(&s,
		    "Median-bin the light curve EXAMPLES/2 in time. A binsize of 0.01 days is used. The output time for each bin is taken to be the center of the bin. After binning the light curve, the -o command outputs the binned light curve to the file EXAMPLES/OUTDIR1/2.bin.txt\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-Phase ls \\\n");
      printtostring(&s,
		    "\t-binlc median nbins 100 tcenter \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/2.phasebin.txt\n\n");
      printtostring(&s,
		    "A slightly more involved example. First the -LS command is used to apply the Lomb-Scargle period-finding algorithm to the light curve EXAMPLES/2 (the period is searched between 0.1 and 10.0 days; the scan is done at a frequency resolution of 0.1/T where T=31.1d is the time-span of the LC; only 1 peak in the periodogram is found, and the periodogram is not output). The light curve is then phased using the period identified by the -LS command. We then bin the phased light curve as in Example 1, though here we specify the number of bins to use rather than the size of the bins. The phase-binned light curve is then output to the file EXAMPLES/OUTDIR1/2.phasebin.txt\n\n");
      printtostring(&s,
		    "Example 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-expr 'rmsbin=mag' -expr 'npoints=1' \\\n");
      printtostring(&s,
		    "\t-binlc median binsize 0.01 bincolumns rmsbin:stddev,npoints:sum \\\n");
      printtostring(&s,
		    "\t\ttcenter \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/2.bin.txt columnformat t,mag,err,rmsbin,npoints\n\n");
      printtostring(&s,
		    "In addition to median-binning the light curves as in Example 1, also record the standard deviation of the points in each bin (to be stored in the rmsbin variable) and the number of points that were in each bin (to be stored in the npoints variable). The output time for each bin is taken to be the center of the bin. After binning the light curve, the -o command outputs the binned light curve to the file EXAMPLES/OUTDIR1/2.bin.txt and includes the rmsbin and npoints variables.\n\n");
      commandfound=1;
    }
  if(!strncmp(c,"-BLS",4) && strlen(c) == 4)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.transit -oneline \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1/ 1 EXAMPLES/OUTDIR1/ 0 fittrap \\\n");
      printtostring(&s,
		    "\t\tnobinnedrms ophcurve EXAMPLES/OUTDIR1/ -0.1 1.1 0.001\n\n");
      printtostring(&s,
		    "Apply the Box-Least Squares (BLS) transit-search algorithm to the light curve EXAMPLES/3.transit (an LC with a transit injected). At each trial frequency we scan fractional transit durations between 0.01 and 0.1 (the option q fixes the range, if we had used r then the range of transit durations to search would depend on the trial frequency). Periods between 0.1 and 20.0 days are searched. A total of 100000 frequencies are searched, and at each frequency we use 200 phase-bins. The time-zone is set to 0 (this only affects how the BLS_fraconenight output statistic is calculated). We only search for one peak in the BLS spectrum. The BLS spectrum is output to the directory EXAMPLES/OUTDIR1, with the filename 3.transit.bls. The best-fit box-transit model is written to the directory EXAMPLES/OUTDIR1 with the filename 3.transit.bls.model. We do not subtract the best-fit transit before passing to the next command (in this example there are no further commands). By specifying fittrap we have the routine fit a trapezoid-shaped transit model after finding the top peak in the BLS spectrum. By specifying nobinnedrms we speed up the routine compared to the default behavior. We output a model transit phase curve to EXAMPLES/OUTDIR1/ (it will have the filename 3.transit.bls.phcurve), this can be used for overplotting the model on the data in the file bls.model (which is only sampled at the observed phases). The phase curve will go from -0.1 to 1.1 with a step-size of 0.001 in phase. By giving the -oneline option before the -BLS command we have the statistics written with one line per statistic.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-BLSFixPer",10) && strlen(c) == 10)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.transit -oneline \\\n");
      printtostring(&s,
		    "\t-rms \\\n");
      printtostring(&s,
		    "\t-BLSFixPer fix 2.12345 q 0.01 0.1 200 0 0 1 fittrap \\\n");
      printtostring(&s,
		    "\t-rms\n\n");
      printtostring(&s,
		    "Fit a Box-shape transit to the light curve EXAMPLES/3.transit at the period of 2.12345 days. Allow for fractional transit durations between 0.01 and 0.1, and use 200 phase-bins. Set the time-zone to 0 (this only affects the BLSFixPer_fraconenight statistic. Do not output the model, but subtract the best-fit model from the light curve before passing it on to the next command (in this example, that is the final -rms command). The fittrap option causes the procedure to fit a trapezoid-shaped transit. The calls to -rms before and after -BLSFixPer show how subtracting the BLS model reduces the scatter of the light curve.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-changeerror",12) && strlen(c) == 12)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/4 -oneline \\\n");
      printtostring(&s,
		    "\t-chi2 \\\n");
      printtostring(&s,
		    "\t-changeerror \\\n");
      printtostring(&s,
		    "\t-chi2\n\n");
      printtostring(&s,
		    "Calculate the chi^2 per degree of freedom of the light curve EXAMPLES/4. In between the two calls to -chi2 were replace the formal errors in the light curve with the RMS of the light curve. As a result, the second call to -chi2 yields chi^2 per degree of freedom = 1.0\n");
      commandfound=1;
    }
  if(!strcmp(c,"-changevariable"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list \\\n");
      printtostring(&s,
		    "\t-LS 0.1 100.0 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-expr \'phase=t\' \\\n");
      printtostring(&s,
		    "\t-changevariable t phase \\\n");
      printtostring(&s,
		    "\t-Phase ls \\\n");
      printtostring(&s,
		    "\t-changevariable t t \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 nameformat \"%%s.phase.txt\" \\\n");
      printtostring(&s,
		    "\t\tcolumnformat \"t:%%17.9f,mag:%%9.5f,err:%%9.5f,phase:%%9.5f\" \\\n");
      printtostring(&s,
		    "\t-header\n\n");
      printtostring(&s,
		    "Use -LS to find periodic signals in the light curves from EXAMPLES/lc_list. Phase the light curves with this period, and output the light curves including time, magnitude, magnitude uncertainty, and phase, to the directory EXAMPLES/OUTDIR1, appending phase.txt to the end of each filename. Note here we first use the -expr command to set the new vector variable \"phase\" equal to \"t\" (this is done on a per-point basis). We then change the time variable to \"phase\", so that the subsequent -Phase command stores the phase in the variable \"phase\" rather than in the variable \"t\". Note the -Phase command expects the time variable to store the times on input, this is why we need the preceding -expr command. We then change the time variable back to \"t\" so that the output light curves will be sorted by time rather than Phase. To have vartools write more than just the default t, mag, and err columns to the output light curves we give the columnformat keyword.\n\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      commandfound=1;
    }
  if(!strncmp(c,"-chi2",5) && strlen(c) == 5)
    {
      printtostring(&s,
		    "\nvartools -header -l EXAMPLES/lc_list -chi2\n\n");
      printtostring(&s,
		    "Calculate chi2 per degree of freedom and the weighted mean magnitude for all light curves in the list EXAMPLES/lc_list.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-chi2bin",8) && strlen(c) == 8)
    {
      printtostring(&s,
		    "\nvartools -header -l EXAMPLES/lc_list \\\n");
      printtostring(&s,
		    "\t-chi2bin 5 5.0 10.0 60.0 1440 14400\n\n");
      printtostring(&s,
		    "Apply a set of moving-mean filters to the light curves in the list EXAMPLES/lc_list and calculate chi^2 per degree of freedom and the weighted mean magnitude for each filter. We use 5 filters of 5.0, 10.0, 60.0, 1440.0, and 14400.0 minutes. Note that the \"filter\" here refers to replacing each point in the light curve with the mean of all points that are within the specified number of minutes of that point. The formal errors are decreased by the amount expected for white noise, as a result the chi2 values for the light curves in this example (which have non-negligible red-noise) increase as the filter-size is increased.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-clip",5) && strlen(c) == 5)
    {
      printtostring(&s,
		    "\nvartools -header -i EXAMPLES/5 \\\n");
      printtostring(&s,
		    "\t-rms \\\n");
      printtostring(&s,
		    "\t-clip 3. 1 \\\n");
      printtostring(&s,
		    "\t-rms \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/5.clip.txt\n\n");
      printtostring(&s,
		    "Calculate the RMS of the light curve EXAMPLES/5, apply iterative 3 sigma clipping, calculate the RMS of the clipped light curve, and output the clipped light curve to EXAMPLES/OUTDIR1/5.clip.txt\n");
      commandfound=1;
    }
  if(!strncmp(c,"-converttime",12) && strlen(c) == 12)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1 -quiet \\\n");
      printtostring(&s,
		    "\t-converttime input jd inputsubtract 2400000. \\\n");
      printtostring(&s,
		    "\t\toutput hjd outputsubtract 2400000. \\\n");
      printtostring(&s,
		    "\t\tradec fix 88.079166 32.5533 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/1.hjdutc\n\n");
      printtostring(&s,
		    "Convert the light curve EXAMPLES/1 from JD-2400000 to HJD-2400000 (Heliocentric Julian Date; assuming an elliptical orbit for the Earth-Moon Barycenter with linear perturbations to the orbital elements) and output to EXAMPLES/OUTDIR1/1.hjdutc.\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1.UTC -quiet \\\n");
      printtostring(&s,
		    "\t-readformat 0 inpututc '%%Y-%%M-%%DT%%h:%%m:%%s' 1 2 3 \\\n");
      printtostring(&s,
		    "\t-converttime input jd inputsys-utc \\\n");
      printtostring(&s,
		    "\t\toutput bjd outputsubtract 2400000. outputsys-tdb \\\n");
      printtostring(&s,
		    "\t\tradec fix 88.079166 32.5533 \\\n");
      printtostring(&s,
		    "\t\tephemfile CSPICEKERNELS/de421.bsp \\\n");
      printtostring(&s,
		    "\t\tleapsecfile CSPICEKERNELS/naif0009.tls \\\n");
      printtostring(&s,
		    "\t\tplanetdatafile CSPICEKERNELS/pck00009.tpc \\\n");
      printtostring(&s,
		    "\t\tobservatory flwo \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/1.bjdtdb\n\n");
      printtostring(&s,
		    "Convert the time in the light curve EXAMPLES/1.UTC from UTC to Barycentric Julian Date (BJD) in the Barycentric Dynamical Time (TDB) reference system. We use the inpututc keyword to the -readformat command to specify that the input is in UTC and to provide its format; the UTC is automatically converted to JD on input. The -converttime command converts from the input JD on the UTC system to the output BJD on the TDB system. We subtract 2400000 from the output time (providing an easier to read number). For conversion to BJD we must additionally specify the RA and DEC of the source, we also need to provide CSPICE ephemeris, leap-second and planetary-data files (not included with this distribution). For completeness we also give the observatory where the observations were made, however this is a very minor correction. The time-converted light curve is output to EXAMPLES/OUTDIR1/1.bjdtdb.\n\n");
      commandfound = 1;
    }
  if(!strcmp(c,"-copylc"))
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-copylc 100 \\\n");
      printtostring(&s,
		    "\t-expr 'mag=err*gauss()' \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-header\n\n");
      printtostring(&s,
		    "Calculate the -LS periodogram for the light curve EXAMPLES/2, then make 100 copies of the light curve replacing the magnitudes with Gaussian random noise, and run the -LS periodogram for each simulation. This is an example of how one might carry out Monte Carlo simulations with VARTOOLS to determine the bandwidth correction to the false alarm probability for a given light curve sampling.  In practice one would want to run substantially more than 100 simulations.\n\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-decorr",7) && strlen(c) == 7)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-rms \\\n");
      printtostring(&s,
		    "\t-decorr 1 1 1 0 1 1 2 0 \\\n");
      printtostring(&s,
		    "\t-rms\n\n");
      printtostring(&s,
		    "Fit a quadratic polynomial to the light curves given in the file EXAMPLES/lc_list. To do this we use no global terms and 1 lc term. For the lc term we use the first column in the light curve (the JD) and fit a second order polynomial in this term to each light curve. We fit for the zero-point term and correct the light curve (so that commands after -decorr will receive light curves with the best-fit quadratic polynomial removed, note that when the light curve is corrected the mean is kept constant), we also subtract the first term in the signal that we are decorrelating against (use JD~-~JD0 rather than JD, since JD*JD runs into round-off problems whereas (JD~-~JD0)*(JD~-~JD0) does not), but we do not output the corrected light curves to the disk. The rms is determined before and after the fit. To interpret the output, note that for light curve 1 we find that the best-fit quadratic has the form: 0.0002554062775*(JD~-~53725.173920)*(JD~-~53725.173920) +~0.0097933162509*(JD~-~53725.173920) +~10.0830375984825, and that fitting this equation reduces the RMS from 0.15944 mag to 0.00211 mag (a quadratic signal was injected into this particular light curve).\n");
      commandfound=1;
    }
  if(!strncmp(c,"-dftclean",9) && strlen(c) == 9)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -oneline \\\n");
      printtostring(&s,
		    "\t-dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\tfinddirtypeaks 1 clip 5. 1\n\n");
      printtostring(&s,
		    "Compute the Discrete Fourier Transform (DFT) of the light curve EXAMPLES/2. We will over-sample the DFT by a factor of 4, and compute the DFT up to a maximum frequency of 10 cycles per day. The resulting power spectrum will be output to the directory EXAMPLES/OUTDIR1 with the filename 4.dftclean.dspec. We search for the highest peak in the power spectrum using a 5-sigma iterative clipping to determine the spectroscopic signal-to-noise ratio of the peak. In this example the highest peak is at 0.81189711 cycles/day corresponding to a period of 1.23168 days, which is close to the injected period of 1.2354 days.\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/4 -oneline \\\n");
      printtostring(&s,
		    "\t-Injectharm fix 0.697516 0 ampfix 0.1 phaserand 0 0 \\\n");
      printtostring(&s,
		    "\t-Injectharm fix 2.123456 0 ampfix 0.05 phaserand 0 0 \\\n");
      printtostring(&s,
		    "\t-Injectharm fix 0.426515 0 ampfix 0.01 phaserand 0 0 \\\n");
      printtostring(&s,
		    "\t-dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\tfinddirtypeaks 3 clip 5. 1 \\\n");
      printtostring(&s,
		    "\t\toutwfunc EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\tclean 0.5 5.0 outcbeam EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\toutcspec EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\tfindcleanpeaks 3 clip 5. 1 \\\n");
      printtostring(&s,
		    "\t\tverboseout\n\n");
      printtostring(&s,
		    "A slightly more involved example to illustrate the use of the clean algorithm. Here 3 harmonic signals are added to the light curve EXAMPLES/4 with periods of 0.697516, 2.123456 and 0.426515 days, and amplitudes of 0.1, 0.05 and 0.01 mag. We then calculate the DFT as in Example 1, but this time we determine the 3 highest peaks in the DFT and we also output the window function to EXAMPLES/OUTDIR1 (the filename will be 4.dftclean.wfunc). We then apply the CLEAN deconvolution algorithm to the DFT using a gain of 0.5 and a S/N threshold of 5.0. The clean-beam, and clean power spectrum are written to EXAMPLES/OUTDIR1 (with filenames 4.dftclean.cbeam and 4.dftclean.cspec respectively). We search the clean power spectrum for 3 peaks, and include the average and standard deviation of the dirty and clean power spectra in the output table of statistics. The frequencies found in the clean spectrum are a little closer to the injected frequencies.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-ensemblerescalesig",19) && strlen(c) == 19)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-chi2 -ensemblerescalesig 3. -chi2\n\n");
      printtostring(&s,
		    "Transform the magnitude uncertainties for the light curves in the list EXAMPLES/lc_list. The -chi2 commands before and after -ensemblerescalesig are used to demonstrate that the chi^2 per dof values change, however, for this small set of light curves, which span a limited range of magnitudes and include several artificial variables, the -ensemblerescalesig command is not recommended. This command should only be used in cases where the light curves span a large range of magnitudes, and one can easily see that the \"floor\" in a magnitude-RMS plot falls well above the expected floor. The example here just illustrates the method of calling this procedure.\n");
      commandfound = 1;
    }
  if(!strcmp(c,"-expr"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1 -expr 'mag=sqrt(mag+5)' -o EXAMPLES/1.add\n\n");
      printtostring(&s,
		    "Use the vartools analytic expression evaluation command -expr to add a constant (in this case 5) to all the magnitude values in the light curve EXAMPLES/1, and then take the square root. Write out the result to EXAMPLES/1.add.\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-rms -chi2 \\\n");
      printtostring(&s,
		    "\t-expr 'mag2=mag' \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -chi2 \\\n");
      printtostring(&s,
		    "\t-expr \\\n");
      printtostring(&s,
		    "\t\t\'mag=(Npoints_5*(Chi2_6-Chi2_2)<-10000)*mag+\n");
      printtostring(&s,"\t\t    (Npoints_5*(Chi2_6-Chi2_2)>=-10000)*mag2\' \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 nameformat '%%s.cleanharm'\n\n");
      printtostring(&s,
		    "An example of using the vartools analytic expression evaluation command, together with other commands, to fit a sinusoid signal to various light curves and subtract the signal only for cases where the sinusoid shows a significant delta chi2 improvement over the flat model. Here the -LS command is used to find a periodic signal, the first call to -expr copies the light curve magnitudes to a new vector variable mag2, -Killharm fits and subtracts a sinusoid model from mag, -chi2 and -rms are calculated on the residuals (currently stored in mag), and the second call to -expr sets the magnitude to the residuals if there is a significant chi2 improvement (Npoints_5*(Chi2_6-Chi2_2)<-10000), or back to the original magnitudes if the improvement is not significant (>=-10000). Note that the '<' and '>=' will evaluate to 1 (0) when true (false). One could alternatively use an -if, -else, and -fi construct to achieve a similar result.\n");
      printtostring(&s,
		    "\nExample 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1 \\\n");
      printtostring(&s,
		    "\t-expr 'flux=10^(-0.4*(mag-25.0))' \\\n");
      printtostring(&s,
		    "\t-stats flux median \\\n");
      printtostring(&s,
		    "\t-expr 'flux=flux/STATS_flux_MEDIAN_1' \\\n");
      printtostring(&s,
		    "\t-stats flux,mag median,stddev \\\n");
      printtostring(&s,
		    "\t-oneline\n\n");
      printtostring(&s,
		    "In this example we use the -expr command to convert the magnitudes in a light curve to fluxes. We then use -stats to compute the median flux, and the second -expr command to normalize the flux by its median value. The final call to stats then calculates the median and standard deviations of the magnitudes and fluxes.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-findblends",11) && strlen(c) == 11)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list_testblend -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-findblends 2.0 fixcolumn LS_Period_1_0 \n\n");
      printtostring(&s,
		    "This example illustrates a basic usage of the -findblends command. The list EXAMPLES/lc_list_testblend contains two light curves as well as the x and y coordinates of each light curve. These light curves are searched for a sinusoidal variation using the -LS command. We then check for variability blends using the -findblends command. Here we consider any stars within 2 pixels of each other as potential blends. We use the period determined by -LS as the variability period for the -findblends command. We do this by giving the \"fixcolumn\" keyword, and then giving the name of the output column which stores the period. Running this example gives \"EXAMPLES/2\" as the source of the variability for both light curves since that light curve has a higher variability amplitude (in flux) than EXAMPLES/2.testblend.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-fluxtomag",10) && strlen(c) == 10)
    {
      printtostring(&s,
		    "\nvartools -i Q1_public/kplr000757076-2009166043257_llc.fits \\\n");
      printtostring(&s,
		    "\t-readformat 0 1 10 11 \\\n");
      printtostring(&s,
		    "\t-fluxtomag 25.0 0 \\\n");
      printtostring(&s,
		    "\t-o kplr000757076-2009166043257_llc.asc.txt\n\n");
      printtostring(&s,
		    "Read in the binary fits Q1 Kepler public light curve for KIC 757076 (not included), convert it to flux using a zero-point magnitude of 25.0, and output the result to an ascii text file.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-GetLSAmpThresh",15) && strlen(c) == 15)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -oneline \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 fitonly \\\n");
      printtostring(&s,
		    "\t-GetLSAmpThresh ls 0.1 -100 harm 0 0\n\n");
      printtostring(&s,
		    "Example illustrating a use of the -GetLSAmpThresh command. Read in the light curve EXAMPLES/2 and search for a period with -LS, Use the -Killharm command as above to get the peak-to-peak amplitude of the signal at the period found by the -LS command. Then run -GetLSAmpThresh to get the minimum amplitude that the signal in EXAMPLES/2 could have had, and still have resulted in an LS detection with Log10_LS_Prob~<~-100.0. The \"ls\" keyword tells the command to take the input period from the last -LS command, and 0.1 is the minimum period searched for by -LS (this sets the FAP scale). We give \"harm~0~0\" to fit a sinusoid for the signal (more complicated signal forms are possible). Running this command shows that the signal could have had an amplitude as low as 0.00248 and still have had Log10_LS_Prob~<~-100.0. For comparison, the real signal has an amplitude of 0.1~mag, and easily passes the threshold.\n");
      commandfound=1;
    }
  if(!strcmp(c,"-if") ||
     !strcmp(c,"-elif") ||
     !strcmp(c,"-else") ||
     !strcmp(c,"-fi"))
    {
      printtostring(&s,
		    "./vartools -l EXAMPLES/lc_list -rms \\\n");
      printtostring(&s,
		    "\t-if 'RMS_0>10*Expected_RMS_0' \\\n");
      printtostring(&s,
		    "\t\t-if 'RMS_0 > 0.1' \\\n");
      printtostring(&s,
		    "\t\t\t-stats mag stddev \\\n");
      printtostring(&s,
		    "\t\t-else \\\n");
      printtostring(&s,
		    "\t\t\t-stats mag pct30 \\\n");
      printtostring(&s,
		    "\t\t-fi \\\n");
      printtostring(&s,
		    "\t-elif 'Npoints_0>3900' \\\n");
      printtostring(&s,
		    "\t\t-stats mag kurtosis \\\n");
      printtostring(&s,
		    "\t-else \\\n");
      printtostring(&s,
		    "\t\t-rms \\\n");
      printtostring(&s,
		    "\t-fi \\\n");
      printtostring(&s,
		    "\t-header\n\n");
      printtostring(&s,
		    "Example illustrating the use of an -if, -elif, -else, -fi construct. Here we compute the -rms of light curves given in the file EXAMPLES/lc_list. Subsequent processing depends on the values returned by -rms. Note the use of nested -if statements.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-Injectharm",11) && strlen(c) == 11)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 -randseed 1 -oneline \\\n");
      printtostring(&s,
		    "\t-Injectharm rand 1.0 5.0 \\\n");
      printtostring(&s,
		    "\t\t0 amplogrand 0.01 0.1 phaserand \\\n");
      printtostring(&s,
		    "\t\t0 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10.0 0.1 1 0\n\n");
      printtostring(&s,
		    "Inject a sinusoid into the light curve EXAMPLES/3 and then conduct a search for a periodic sinusoidal signal using -LS. In this example we adopt a random period between 1.0 and 5.0 days, we include 0 harmonic overtones (Nharm=0 implies only the fundamental is used) and take the amplitude of the fundamental signal from a uniform-log distribution, and we also adopt a random phase. We use no sub-harmonics and we output the model to the directory EXAMPLES/OUTDIR1 (it will have the filename \"EXAMPLES/OUTDIR1/3.injectharm.model\", the model values are in the third column). By giving the option \"-randseed~1\" we explicitly seed the random number generator with the value 1. Use \"-randseed~time.\" to get a non-repeatable test.\n\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\necho EXAMPLES/4 | \\\n");
      printtostring(&s,
		    "\tgawk '{amp = 0.25; \\\n");
      printtostring(&s,
		    "\t\tfor(i=1; i <= 10; i += 1) { \\\n");
      printtostring(&s,
		    "\t\t\tprint $1, amp; amp = amp*0.5; \\\n");
      printtostring(&s,
		    "\t\t}}' | \\\n");
      printtostring(&s,
		    "\tvartools -l - -header -numbercolumns -parallel 4 \\\n");
      printtostring(&s,
		    "\t-Injectharm fix 0.514333 10 \\\n");
      printtostring(&s,
		    "\t\tamplist column 2 phaserand \\\n");
      printtostring(&s,
		    "\t\tampfix 0.47077 amprel phasefix 0.60826 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.35916 amprel phasefix 0.26249 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.23631 amprel phasefix -0.06843 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.16353 amprel phasefix 0.60682 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.10621 amprel phasefix 0.28738 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.06203 amprel phasefix 0.95751 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.03602 amprel phasefix 0.58867 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.02900 amprel phasefix 0.22322 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.01750 amprel phasefix 0.94258 phaserel \\\n");
      printtostring(&s,
		    "\t\tampfix 0.00768 amprel phasefix 0.66560 phaserel \\\n");
      printtostring(&s,
		    "\t\t0 0 \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10.0 0.01 2 0 \\\n");
      printtostring(&s,
		    "\t-aov_harm 2 0.1 10.0 0.1 0.01 2 0\n\n");
      printtostring(&s,
		    "Inject an RR Lyrae signal into a light curve and recover it. The initial gawk command creates a list with 10 rows, each row gives the name of the light curve to inject the signal into (EXAMPLES/4) and the amplitude of the injected signal. The call to vartools uses -Injectharm to inject a signal with a fixed period of 0.514333 days, and with 10 harmonics. For the fundamental mode we take the amplitude from the input list, and we allow the phase to be random. We then fix the relative amplitudes and phases of the 10 harmonics (these values give an RR-Lyrae shaped signal, see \"vartools -example -Killharm\" for an example of how these coefficients can be determined). After injecting the signal it is then recovered using the -LS and -aov_harm commands. To speed up the process we use \"-parallel 4\" which processes up to four light curves simultaneously (use a higher or lower number depending on the number of CPU cores on your machine). Note that the order of rows in the output is arbitrary when using parallel.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-Injecttransit",14) && strlen(c) == 14)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 -oneline -randseed 1 \\\n");
      printtostring(&s,
		    "\t-Injecttransit lograndfreq 0.2 2.0 \\\n");
      printtostring(&s,
		    "\t\tRpfix 1.0 Mpfix 1.0 \\\n");
      printtostring(&s,
		    "\t\tphaserand sinirand \\\n");
      printtostring(&s,
		    "\t\teomega efix 0. ofix 0. \\\n");
      printtostring(&s,
		    "\t\tMstarfix 1.0 Rstarfix 1.0 \\\n");
      printtostring(&s,
		    "\t\tquad ldfix 0.3471 0.3180 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\t1 fittrap\n\n");
      printtostring(&s,
		    "Inject a transit into the light curve \"EXAMPLES/3\" and recover \\ it using -BLS. We draw the period from a uniform-log random distribution in frequency with limits of 0.2 c/d and 2.0 c/d. We fix the planet radius to 1.0 R_J and the planet mass to 1.0 M_J. We adopt a random phase, and draw sin(i) from a random orientation distribution (subject to the constraint that a transit must occur). We fix the eccentricity and argument of periastron to 0, we fix the star mass and radius to 1.0 M_sun and 1.0 R_sun, respectively, and we adopt a quadratic limb darkening law, with coefficients 0.3471 and 0.3180. We output the model to the directory EXAMPLES/OUTDIR1 (the filename will be EXAMPLES/OUTDIR1/3.injecttransit.model). After injecting the transit we use the -BLS command to recover the transit.\n");
      commandfound=1;
    }
  if(!strncmp(c, "-Jstet", 6) && strlen(c) == 6)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-Jstet 0.5 EXAMPLES/dates_tfa\n\n");
      printtostring(&s,
		    "Calculate Stetson's J statistic for the light curves in the list EXAMPLES/lc_list. Use 0.5 days to distinguish between \"near\" and \"far\" observations.\n");
      commandfound=1;
    }
  if(!strncmp(c, "-Killharm", 9) && strlen(c) == 9)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -oneline \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-rms -chi2 \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -chi2\n\n");
      printtostring(&s,
		    "Search for a periodic signal in the light curve EXAMPLES/2 using the Lomb-Scargle algorithm, and then fit and remove a sinusoid using -Killharm. We include calls to -rms and -chi2 before and after calling -Killharm to show how these statistics change after subtracting the best-fit sinusoid. For the -Killharm command we take the period from the last ls command, we only fit the fundamental (no harmonics or sub-harmonics), and we do not output the best-fit model.\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/M3.V006.lc -oneline \\\n");
      printtostring(&s,
		    "\t-Killharm fix 1 0.514333 10 0 1 \\\n");
      printtostring(&s,
		    "\t\tEXAMPLES/OUTDIR1/ fitonly outRphi\n\n");
      printtostring(&s,
		    "Fit a harmonic series to the RR Lyrae light curve EXAMPLES/M3.V006.lc. We fix the period to 0.514333 days, and we fit 10 harmonics plus the fundamental. We do not fit sub-harmonics. The best-fit model is output to EXAMPLES/OUTDIR1 (the filename will be EXAMPLES/OUTDIR1/M3.V006.lc.killharm.model). We do not subtract the model (the fitonly keyword) and we give relative amplitudes and phases (the amplitudes and phases in this format can be used in the -Injectharm command to inject a harmonic series with the fixed signal shape, but random overall amplitude and phase. See \"vartools -example -Injectharm\").\n");
      commandfound=1;
    }
  if(!strcmp(c,"-linfit"))
    {
      printtostring(&s,
		    "./vartools -i EXAMPLES/1 \\\n");
      printtostring(&s,
		    "\t-stats t min \\\n");
      printtostring(&s,
		    "\t-expr t0=STATS_t_MIN_0 \\\n");
      printtostring(&s,
		    "\t-linfit 'a*(t-t0)^2+b*(t-t0)+c' 'a,b,c' \\\n");
      printtostring(&s,
		    "\t-oneline\n\n");
      printtostring(&s,
		    "Fit a quadratic function in time to the light curve EXAMPLES/1 using the -linfit command. We first determine the minimum time value using the -stats command, and store it in the variable t0 with the -expr command. For the -linfit command we fit the function 'a*(t-t0)^2+b*(t-t0)+c' where we subtract the t0 to avoid round-off errors accumulating and producing a bad fit. We then give the list of free parameters as a,b,c. The output table gives the best-fit value and standard error for each free parameter.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-LS",3) && strlen(c) == 3)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -oneline \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 5 1 EXAMPLES/OUTDIR1 whiten clip 5. 1\n\n");
      printtostring(&s,
		    "Run the Lomb-Scargle period-finding algorithm on the light curve EXAMPLES/2. Search for periods between 0.1 and 10.0 days at a frequency resolution of 0.1/T (T is the time-span of the lc, 31.1d in this case). Report the top 5 peaks, and output the periodogram to EXAMPLES/OUTDIR1 (the filename will be EXAMPLES/OUTDIR1/2.ls). Pre-whiten the light curve and re-apply L-S before finding the next peak, and use a 5 sigma iterative clipping in determining the spectroscopic S/N.\n");
      commandfound =1;
    }
  if(!strncmp(c,"-MandelAgolTransit",18) && strlen(c) == 18)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.transit -oneline \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \\\n");
      printtostring(&s,
		    "\t-MandelAgolTransit bls quad 0.3471 0.3180 \\\n");
      printtostring(&s,
		    "\t\t1 1 1 1 0 0 1 0 0 0 0 1 EXAMPLES/OUTDIR1\n\n");
      printtostring(&s,
		    "Use -BLS to identify a transit signal in the light curve EXAMPLES/3.transit and fit a Mandel-Agol transit model to it. For the -MandelAgolTransit command we indicate that the initial parameter values should be determined based on the results from the -BLS command. We use a quadratic limb-darkening law, with parameters 0.3471 and 0.3180. We set flags to vary the ephemeris (P and T0), Rp/R*, a/R* and the impact parameter. We do not vary either the eccentricity or argument of periastron. We vary the mean out-of-transit magnitude, and we do not vary either of the limb-darkening coefficients. We do not fit and RV curve or subtract the best-fit model from the light curve. We output the best-fit model to the directory EXAMPLES/OUTDIR1 (the filename will be EXAMPLES/OUTDIR1/3.transit.mandelagoltransit.model.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-medianfilter",13) && strlen(c) == 13)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1 -header -chi2 \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-medianfilter 0.05 \\\n");
      printtostring(&s,
		    "\t-chi2 -o EXAMPLES/OUTDIR1/1.medianhighpass \\\n");
      printtostring(&s,
		    "\t-restorelc 1 \\\n");
      printtostring(&s,
		    "\t-medianfilter 0.05 replace \\\n");
      printtostring(&s,
		    "\t-chi2 -o EXAMPLES/OUTDIR1/1.medianlowpass\n\n");
      printtostring(&s,
		    "Apply a high-pass, and a low-pass median filter to the quadratically varying light curve EXAMPLES/1. Before applying the high-pass filter we save the light curve using the -savelc command. For the high-pass filter we use -medianfilter with a time-scale of 0.05 days. We output the filtered light curve to the file EXAMPLES/OUTDIR1/1.medianhighpass. We then restore the light curve to its state at the first -savelc command and apply the low-pass filter (by including the keyword \"replace\" in the -medianfilter command). The low-pass filtered light curve is output to the file EXAMPLES/OUTDIR1/1.medianlowpass. We include calls to -chi2 before the filtering and after each of the filters to show how the chi2 per degree of freedom is affected.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-microlens",10) && strlen(c) == 10)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/4.microlensinject -oneline \\\n");
      printtostring(&s,
		    "\t-microlens f0 auto f1 auto u0 auto t0 auto tmax auto \\\n");
      printtostring(&s,
		    "\t\tomodel EXAMPLES/OUTDIR1\n\n");
      printtostring(&s,
		    "Fit a simple microlensing model to the simulated light curve EXAMPLES/4.microlensinject. The initial values for all parameters are set automatically. The best-fit model is output to the directory EXAMPLES/OUTDIR1 (the filename is 4.microlensinject.microlens).\n");
      commandfound=1;
    }
  if(!strcmp(c,"-nonlinfit"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring_nowrap(&s,
		    "./vartools -i EXAMPLES/3 \\\n");
      printtostring_nowrap(&s,
		    "\t-stats t min,max \\\n");
      printtostring_nowrap(&s,
		    "\t-expr t1=STATS_t_MIN_0 \\\n");
      printtostring_nowrap(&s,
		    "\t-expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \\\n");
      printtostring_nowrap(&s,
		    "\t-expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \\\n");
      printtostring_nowrap(&s,
		    "\t-nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \\\n");
      printtostring_nowrap(&s,
		    "\t\t'c=(t1+Dt*0.3):(0.1*Dt),d=(Dt*0.1):(0.1*Dt)' \\\n");
      printtostring_nowrap(&s,
		    "\t\tlinfit a,b amoeba omodel EXAMPLES/OUTDIR1/ \\\n");
      printtostring_nowrap(&s,
		    "\t-oneline\n\n");
      printtostring(&s,
		    "Add a Gaussian function to a light curve (done with the third '-expr' call), and then fit a Gaussian to the light curve using the -nonlinfit command. We first determine the minimum and maximum time values using the stats commands. The parameters of the injected Gaussian (the peak time, and the standard deviation) are then given in terms of these parameters.  For the -nonlinfit command we fit a function of the form 'a+b*exp(-(t-c)^2/(2*d^2))' with the free parameters being a, b, c and d. Two of these parameters, c and d, enter in a nonlinear way into the function and thus must be varied using a nonlinear fitting algorithm. After the function we list these two nonlinear parameters, together with initial guesses for the optimal values and initial step-sizes to use in the search. For 'c' the initial value is '(t1+Dt*0.3)' while the step-size is '(0.1*Dt)'. For 'd' the initial value is '(Dt*0.1)' with a step-size of '(0.1*Dt)'. Obviously for this artificial example we could have simply set the initial values to the values used in creating the artificial signal, but we choose slightly different values to illustrate the search procedure. Because the parameters 'a' and 'b' enter linearly into the model, we can optimize these parameters using linear least squares. We do this by giving the 'linfit' keyword and then listing the two parameters. Note that we could also have varied all four parameters a,b,c,d through the non-linear algorithm (in which case we would need to provide initial guesses and step-sizes for a and b on the non-linear parameter string, and we would not provide the 'linfit' keyword or the parameter string following it). We use the 'amoeba' algorithm to perform the non-linear search. This is an implementation of the downhill simplex method, a greedy optimization algorithm which will find a local chi^2 minimum that is close to the initial starting values. We output the model to the directory 'EXAMPLES/OUTDIR1'. In this case the model will be in the file 'EXAMPLES/OUTDIR1/3.nonlinfit.model'.\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring_nowrap(&s,
		    "./vartools -i EXAMPLES/3 \\\n");
      printtostring_nowrap(&s,
		    "\t-stats t min,max \\\n");
      printtostring_nowrap(&s,
		    "\t-expr t1=STATS_t_MIN_0 \\\n");
      printtostring_nowrap(&s,
		    "\t-expr 'Dt=(STATS_t_MAX_0-STATS_t_MIN_0)' \\\n");
      printtostring_nowrap(&s,
		    "\t-expr 'mag=mag+0.1*exp(-0.5*((t-(t1+Dt*0.2))/(Dt*0.05))^2)' \\\n");
      printtostring_nowrap(&s,
		    "\t-nonlinfit 'a+b*exp(-(t-c)^2/(2*d^2))' \\\n");
      printtostring_nowrap(&s,
		    "\t\t'a=10.167:0.0002,b=0.1:0.0008,c=(t1+Dt*0.2):(0.005),d=(Dt*0.05):(0.016)' \\\n");
      printtostring_nowrap(&s,
		    "\t\tmcmc Nlinkstotal 10000 outchains EXAMPLES/OUTDIR1/ \\\n");
      printtostring_nowrap(&s,
		    "\t-oneline\n\n");
      printtostring(&s,
		    "Similar to Example 1, in this case we use an MCMC procedure to explore the chi^2 landscape. We will run a total of 10000 links in the MCMC chain, and output the chain to the directory EXAMPLES/OUTDIR1. The file will be EXAMPLES/OUTDIR1/3.mcmc in this case.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-o",2) && strlen(c) == 2)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 100.0 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-expr phase=t \\\n");
      printtostring(&s,
		    "\t-changevariable t phase \\\n");
      printtostring(&s,
		    "\t-Phase ls \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\tnameformat \"file_%%s_%%05d_simout.txt\" \\\n");
      printtostring(&s,
		    "\t\tcolumnformat \"t:%%11.5f,phase:%%8.5f,mag:%%7.4f,err:%%7.4f\"\n\n");
      printtostring(&s,
		    "Example illustrating the use of the \"nameformat\" and \"columnformat\" keywords for the -o command. Light curves are read-in from the list, the -LS command is used to find the periods. The -expr command then defines a new vector \"phase\" which is initialized to the times in the light curves. The -changevariable command causes subsequent commands to use phase in cases where the time would normally be used. This, together with the following -Phase command, causes the vector \"phase\" to store the light curve phase for the period found with -LS. The light curves are then output to the directory EXAMPLES/OUTDIR1. The nameformat keyword gives the rule for naming the output files. The first light curve (\"EXAMPLES/1\") will yield output to the file \"EXAMPLES/OUTDIR1/file_1_00001_simout.txt\", the second (\"EXAMPLES/2\") to the file \"EXAMPLES/OUTDIR1/file_2_00002_simout.txt\", and so on. If the nameformat had not been given, the first file would have been output to \"EXAMPLES/OUTDIR1/1\" and so on. The columnformat keyword specifies how the data will be formatted in the output light curve. Here we indicate that four quantities, the time, phase, magnitude, and error will be included in the output. We also give printf like formatting rules for each of these to make the output easier to read. If columnformat had not been given, then only t, mag and err would have been output, and they would have all been output using the formats %%17.9f, %%9.5f, and %%9.5f respectively.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-Phase",6) && strlen(c) == 6)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -header \\\n");
      printtostring(&s,
		    "\t-Phase fix 1.2354 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/2.phase.txt\n\n");
      printtostring(&s,
		    "Phase the light curve EXAMPLES/2 with a period 1.2354 d and output the result to EXAMPLES/OUTDIR1/2.phase.txt.\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.transit -oneline \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \\\n");
      printtostring(&s,
		    "\t-Phase bls T0 bls 0.0 phasevar ph startphase -0.5 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/3.phase.txt \\\n");
      printtostring(&s,
		    "\t\tcolumnformat t,mag,err,ph \\\n");
      printtostring(&s,
		    "\t-changevariable t ph \\\n");
      printtostring(&s,
		    "\t-binlc 1 nbins 200 0 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1/3.phasebin.txt\n\n");
      printtostring(&s,
		    "Example illustrating the use of the -Phase command. We apply -BLS to identify a transit signal in the light curve EXAMPLES/3.transit. We then phase the light curve taking the period from bls, and the time of zero phase from bls. We set the phase of mid-transit to 0.0, and store the phases to the variable ph, rather than overwriting the times. We use the \"startphase -0.5\" term to have the phases run from -0.5 to 0.5, rather than from 0 to 1. We output the result to EXAMPLES/OUTDIR1/3.phase.txt using the \"columnformat\" keyword to include the phases in the fourth column of the output. We then use the \"-changevariable\" command to switch the time variable to \"ph\", and then median-bin the phased light curve using 200 phase bins, and output the result to EXAMPLES/OUTDIR1/3.phasebin.txt.\n");
      commandfound=1;
    }
#ifdef _HAVE_PYTHON
  if(!strcmp(c,"-python"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 -header \\\n");
      printtostring(&s,
		    "\t-python 'b = numpy.var(mag)' invars mag outvars b outputcolumns b\n\n");
      printtostring(&s,
		    "Example of using python to calculate the variance in the magnitudes for each light curve read from the list file EXAMPLES/lc_list. The python expression \"b = numpy.var(mag)\" will be evaluated for each light curve. The variable \"mag\" will be passed as an input to the python command (it will be treated as a numpy array within python), and the variable \"b\" will be read out from the command. The value of this variable, which stores the variance, will be included as a column in the output ASCII table. The full output from this command is:\n\n#Name PYTHON_b_0\n"
		    "EXAMPLES/1 0.025422461711037084\n"
		    "EXAMPLES/2 0.0013420988067623005\n"
		    "EXAMPLES/3 2.3966645306408949e-05\n"
		    "EXAMPLES/4 4.3733138204733634e-06\n"
		    "EXAMPLES/5 8.2971716866526236e-06\n"
		    "EXAMPLES/6 4.3664615059428104e-06\n"
		    "EXAMPLES/7 1.216345131566495e-05\n"
		    "EXAMPLES/8 5.0623773543353351e-06\n"
		    "EXAMPLES/9 3.4861868515750583e-06\n"
		    "EXAMPLES/10 5.5813996936871234e-06\n");
      printtostring(&s,
		    "\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\n> cat EXAMPLES/plotlc.py\n\n");
      printtostring(&s,
		    "import matplotlib.pyplot as plt\n"
		    "\n"
		    "def plotlc(lcname,outdir,t,ph,mag,P):\n"
		    "\tlcbasename = lcname.split('/')[-1]\n"
		    "\tplt.figure(1)\n"
		    "\tplt.subplot(211)\n"
		    "\tplt.gca().invert_yaxis()\n"
		    "\ttcorr = t - t[0]\n"
		    "\tplt.plot(tcorr, mag, 'bo', markersize=0.5)\n"
		    "\tplt.ylabel('magnitude')\n"
		    "\tplt.title(lcbasename+' P='+str(P))\n"
		    "\tplt.xlabel('time - '+str(t[0]))\n"
		    "\tplt.subplot(212)\n"
		    "\tplt.gca().invert_yaxis()\n"
		    "\tplt.plot(ph, mag, 'bo', markersize=0.5)\n"
		    "\tplt.ylabel('magnitude')\n"
		    "\tplt.xlabel('phase')\n"
		    "\tplt.savefig(outdir+'/'+lcbasename+'.png',format=\"png\")\n"
		    "\tplt.close()\n");

      printtostring(&s,
		    "\n> vartools -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 100. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-if 'Log10_LS_Prob_1_0<-100' \\\n");
      printtostring(&s,
		    "\t\t-Phase ls phasevar ph \\\n");
      printtostring(&s,
		    "\t\t-python 'plotlc(Name,\"EXAMPLES/\",t,ph,mag,LS_Period_1_0)' \\\n");
      printtostring(&s,
		    "\t\t\tinit file EXAMPLES/plotlc.py \\\n");
      printtostring(&s,
		    "\t-fi\n");
      printtostring(&s,
		    "\nExample of running the Lomb-Scargle periodogram on light curves in the list file EXAMPLES/lc_list, and then using matplotlib.pyplot in python to create .png plots of those light curves which have a log10 false alarm probability from LS of less than -100.  The routine for plotting a light curve is stored in the file EXAMPLES/plotlc.py, where we import the matplotlib.pyplot module and then define a python function for making a plot. We use the \"init file\" option to the -python command in VARTOOLS to load this code on initialization (technically the code is incorporated into a module together with a function that VARTOOLS constructs for calling the python commands to be run on each light curve). The expression 'plotlc(Name,\"EXAMPLES/\",t,ph,mag,LS_Period_1_0)' will then be executed on each light curve, where the variable \"Name\" stores the full light curve name, \"ph\" stores the phases which are calculated with the -Phase command, and \"LS_Period_1_0\" is the period returned by the -LS command.\n\nWhen run on the files in EXAMPLES/lc_list this will produce the plots EXAMPLES/1.png and EXAMPLES/2.png showing the light curves vs time and phase.\n\nNote that to run this example vartools will need to be compiled against a version of python which is able to import matplotlib.\n");
      printtostring(&s,
		    "\nExample 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\n> vartools -l EXAMPLES/lc_list -inputlcformat t:1,mag:2,err:3 -header \\\n");
      printtostring(&s,
		    "\t-LS 0.1 100. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-Phase ls phasevar ph \\\n");
      printtostring(&s,
		    "\t-python \\\n");
      printtostring(&s,
		    "\t\t'for i in range(0,len(mag)):\n\t\t\tplotlc(Name[i],\"EXAMPLES/\",t[i],ph[i],mag[i],LS_Period_1_0[i])' \\\n");
      printtostring(&s,
		    "\t\tinit file EXAMPLES/plotlc.py \\\n");
      printtostring(&s,
		    "\t\tprocess_all_lcs\n");
      printtostring(&s,
		    "\nSame as in example 2, except here we plot all of the light curves (no \"-if\" command to VARTOOLS), and we use the \"process_all_lcs\" keyword to send all of the light curves to python at once. In this case VARTOOLS supplies these as lists of numpy arrays, so we use the for loop to cycle through all of the light curves calling plotlc on each light curve in turn.\n");
      commandfound = 1;
    }
#endif
  if(!strcmp(c,"-resample"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -resample linear \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/2.resample.example1\n\n");
      printtostring(&s,
		    "Resample the light curve EXAMPLES/2 using linear interpolation and extrapolation, with default options for determining the resample times. Output the resampled light curve to the file EXAMPLES/2.resample.example1.\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/2 -resample splinemonotonic \\\n");
      printtostring(&s,
		    "\ttstart fix 53726 tstop fix 53756 Npoints fix 1000 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/2.resample.example2\n\n");
      printtostring(&s,
		    "Resample the light curve EXAMPLES/2 using monotonic spline interpolation. Specify the start and stop times for the resampled points, and also specify the number of resampled points. Output the resampled light curve to the file EXAMPLES/2.resample.example2.\n");
      printtostring(&s,
		    "\nExample 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/4 -resample linear \\\n");
      printtostring(&s,
		    "\tfile fix EXAMPLES/8 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/4.resample.example3\n\n");
      printtostring(&s,
		    "Resample the light curve EXAMPLES/4 onto the same time-base as the light curve EXAMPLES/8 using linear interpolation. Output the result to EXAMPLES/4.resample.example3.\n");
      printtostring(&s,
		    "\nExample 4:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/1 -resample splinemonotonic \\\n");
      printtostring(&s,
		    "\ttstart fix 53725 tstop fix 53757 delt fix 0.001 \\\n");
      printtostring(&s,
		    "\tgaps percentile_sep 80 bspline nbreaks 15 order 3 \\\n");
      printtostring(&s,
		    "\textrap nearest \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/1.resample.example4\n\n");
      printtostring(&s,
		    "Resample the light curve EXAMPLES/1 onto a uniform time-base between t=53725 and t=53757 with a step-size of 0.001 days. Use different methods of interpolation for resampled points that are close to observed times, for those that are far from observed times, and for those that are extrapolations beyond the observed time-range. For the close points we will use monotonic spline interpolation, while for the far points we use B-spline interpolation (using 15 breaks and a third order spline). The cutoff between near and far separations is taken to be the 80th percentile of the time separations in the input light curve. For extrapolated points we will use the nearest neighbor resampling. The result is output to the file EXAMPLES/1.resample.example4.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-rescalesig",11) && strlen(c) == 11)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/4 -oneline \\\n");
      printtostring(&s,
		    "\t-chi2 -rescalesig -chi2\n\n");
      printtostring(&s,
		    "Rescale the formal errors in the light curve EXAMPLES/4 such that chi2 per degree of freedom equals 1. We call -chi2 before and after -rescalesig to demonstrate that this rescaling has been done.\n");
      commandfound=1;
    }
  if(!strcmp(c,"-restricttimes"))
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 -stats t min,max \\\n");
      printtostring(&s,
		    "\t-restricttimes JDrange 53740 53750 \\\n");
      printtostring(&s,
		    "\t-stats t min,max -oneline\n\n");
      printtostring(&s,
		    "Filter the light curve to remove any points outside the range 53740 < t < 5370. The calls to -stats before and after the call to -restricttimes show how the command affects the timespan of the light curve. The output from this command is:\n\n");
      printtostring(&s,
"Name                  = EXAMPLES/3\n"
"STATS_t_MIN_0         = 53725.173920000001\n"
"STATS_t_MAX_0         = 53756.281021000003\n"
"RestrictTimes_MinJD_1 = 53740\n"
"RestrictTimes_MaxJD_1 = 53750\n"
"STATS_t_MIN_2         = 53740.336210000001\n"
		    "STATS_t_MAX_2         = 53745.478681000001\n");
      printtostring(&s,
		    "\nExample 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 -stats mag min,max \\\n");
      printtostring(&s,
		    "\t-restricttimes expr '(mag>10.16311)&&(mag<10.17027)' \\\n");
      printtostring(&s,
		    "\t-stats mag min,max -oneline\n\n");
      printtostring(&s,
		    "Filter the light curve to remove any points outside the range 10.16311 < mag < 10.17027. The calls to -stats before and after the call to -restricttimes show how the commands affects the light curve range. The output from this command is:\n\n");
      printtostring(&s,
		    "Name            = EXAMPLES/3\n"
		    "STATS_mag_MIN_0 = 10.141400000000001\n"
		    "STATS_mag_MAX_0 = 10.1921\n"
		    "STATS_mag_MIN_2 = 10.163119999999999\n"
		    "STATS_mag_MAX_2 = 10.170260000000001\n\n");
      printtostring(&s,
		    "\nExample 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 -stats mag pct20.0,pct80.0 \\\n");
      printtostring(&s,
		    "\t-restricttimes expr \\\n");
      printtostring(&s,
		    "\t\t'(mag>STATS_mag_PCT20_00_0)&&(mag<STATS_mag_PCT80_00_0)' \\\n");
      printtostring(&s,
		    "\t-stats mag min,max -oneline\n\n");
      printtostring(&s,
		    "Same as Example 2, but here we make the cut in terms of the magnitude percentiles, which are calculated with the -stats command, rather than giving the magnitude range explicitly.\n\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-restorelc",10) && strlen(c) == 10)
    {
      printtostring(&s,
		    "\n\nSee \"vartools -example -savelc\" for an example of how to use the -restorelc command.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-rms",4) && strlen(c) == 4)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header \\\n");
      printtostring(&s,
		    "\t-rms\n\n");
      printtostring(&s,
		    "Calculate the mean magnitude, RMS and expected RMS based on the formal magnitude uncertainties for the light curves in EXAMPLES/lc_list.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-rmsbin",7) && strlen(c) == 7)
    {
      printtostring(&s,
		    "\nvartools -header -l EXAMPLES/lc_list \\\n");
      printtostring(&s,
		    "\t-rmsbin 5 5.0 10.0 60.0 1440 14400\n\n");
      printtostring(&s,
		    "Apply a set of moving-mean filters to the light curves in the list EXAMPLES/lc_list and calculate mean, RMS, and expected RMS assuming white noise for each filter. We use 5 filters of 5.0, 10.0, 60.0, 1440.0, and 14400.0 minutes. Note that the \"filter\" here refers to replacing each point in the light curve with the mean of all points that are within the specified number of minutes of that point.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-savelc",7) && strlen(c) == 7)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list -header -numbercolumns \\\n");
      printtostring(&s,
		    "\t-nobuffer -parallel 4 \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-clip 5.0 1 \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-LS 0.1 100. 0.1 3 0 clip 5. 1 \\\n");
      printtostring(&s,
		    "\t-aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \\\n");
      printtostring(&s,
		    "\t-restorelc 1 \\\n");
      printtostring(&s,
		    "\t-clip 10.0 1 \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \\\n");
      printtostring(&s,
		    "\t-restorelc 2 \\\n");
      printtostring(&s,
		    "\t-changeerror \\\n");
      printtostring(&s,
		    "\t-autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/\n\n");
      printtostring(&s,
		    "An example of running a battery of variability selection algorithms on a number of light curves in parallel. We first save the initial state of the light curve using the -savelc command, then apply iterative 5-sigma clipping to the light curve. We save the 5-sigma clipped light curve. We then run the -LS and -aov period finding algorithms. We then restore the light curve to its state before the 5-sigma clipping and apply 10-sigma clipping and run BLS (BLS would be sensitive to eclipses. To search for eclipses we would want to use a less aggressive sigma-clipping). After BLS we restore the light curve to its state after the 5-sigma clipping was applied and replace the errors in the light curve with the RMS. Finally we run the -autocorrelation command on the light curve which will output the autocorrelation function to the EXAMPLES/OUTDIR1 directory (see for example EXAMPLES/OUTDIR1/2.autocorr which is periodic and has a first peak at 1.23 days).\n");
      commandfound=1;
    }
  if(!strncmp(c,"-SoftenedTransit",16) && strlen(c) == 16)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.transit -oneline \\\n");
      printtostring(&s,
		    "\t-BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \\\n");
      printtostring(&s,
		    "\t-SoftenedTransit bls 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1 0\n\n");
      printtostring(&s,
		    "Fit a Protopapas et al. 2005 softened transit model to the light curve EXAMPLES/3.transit. We first identify the transit signal using -BLS and then use the results from BLS to initialize the parameters for -SoftenedTransit. We vary the ephemeris, eta, cval, delta and mconst. We do not subtract the best-fit model from the light curve. The model is output to the directory EXAMPLES/OUTDIR1 (the filename will be 3.transit.softenedtransit.model). We do not fit a harmonic series to the light curve simultaneously with the transit.\n");
      commandfound=1;
    }
  if(!strncmp(c,"-Starspot",9) && strlen(c) == 9)
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3.starspot -oneline \\\n");
      printtostring(&s,
		    "\t-aov Nbin 20 0.1 10. 0.1 0.01 5 0 \\\n");
      printtostring(&s,
		    "\t-Starspot aov 0.0298 0.08745 20. 85. 30. 0. -1 \\\n");
      printtostring(&s,
		    "\t\t1 0 0 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1/\n\n");
      printtostring(&s,
		    "Find the rotation period of the simulated light curve EXAMPLES/3.starspot using the -aov command, and then fit a Dorren 1987 single-starspot model to the light curve. Take the period from the -aov command, take 0.0298 and 0.08745 for the parameters a and b, and set the initial values of the spot radius, stellar inclination, spot latitude, and spot longitude to 20., 85., 30., and 0. degrees, respectively. Have the routine automatically determine the initial value for the unspotted photosphere magnitude (by giving -1 for the initial value). We set the flags to vary the period, spot radius, inclination, latitude, longitude, and unspotted magnitude. We keep a and b fixed. We do not subtract the model from the light curve. We output the best-fit model to the directory EXAMPLES/OUTDIR1 (the filename will be 3.starspot.starspot.model). The light curve was simulated with parameters P = 3.12345, a=0.0298, b=0.08745, alpha=30., i=90., chi=45. Due to degeneracies in the model, the recovered values are substantially different from the input values.\n");
      commandfound=1;
    }
  if(!strcmp(c,"-stats"))
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/3 \\\n");
      printtostring(&s,
		    "\t-oneline \\\n");
      printtostring(&s,
		    "\t-expr 'mag2=mag+0.01*gauss()' \\\n");
      printtostring(&s,
		    "\t-stats mag,mag2 \\\n");
      printtostring_nowrap(&s,
		    "\t\tmean,weightedmean,median,stddev,meddev,medmeddev,MAD,kurtosis,skewness,pct10,pct20,pct80,pct90,max,min,sum\n\n");
      printtostring(&s,
		    "Calculate a variety of statistics for the magnitudes in a light curve, and for the magnitudes after adding gaussian noise to them. The call to -expr defines a new vector mag2 which is equal to mag with some gaussian noise added. In the call to stats we first tell it which variables to compute the statistics for (mag and mag2), we then give the statistics to compute. Note that the pct## statistics are percentiles (in this case the 10th, 20th, 80th and 90th percentiles).\n");
      commandfound=1;
    }
  if(!strncmp(c,"-SYSREM",7) && strlen(c) == 7)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/trendlist_tfa -header \\\n");
      printtostring(&s,
		    "\t-rms \\\n");
      printtostring(&s,
		    "\t-SYSREM 2 1 EXAMPLES/3 5. 5. 8. 1 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1/sysrem.trends 1 \\\n");
      printtostring(&s,
		    "\t-rms\n\n");
      printtostring(&s,
		    "Example illustrating the syntax of the -SYSREM command. We apply SYSREM to the light curves in the list EXAMPLES/trendlist_tfa, using 2 initial \"color\"-like terms (corresponding to the second and third columns in EXAMPLES/trendlist_tfa) and one initial \"airmass\"-like term. For the \"airmass\"-like term was use the file EXAMPLES/3. We apply 5 sigma-clipping for determining the mean magnitudes of the light curves and for determining whether or not points contribute to the iterative SYSREM fit. We set the saturation magnitude to 8.0 and we correct the light curves. We output the SYSREM model light curves to the directory EXAMPLES/OUTDIR1 and we output the trends to the file EXAMPLES/OUTDIR1/sysrem.trends. We calculate the RMS before and after applying -SYSREM to illustrate the effect of SYSREM on the light curves. In practice one would apply SYSREM to a much larger set of light curves, and would use the airmass, seeing, temperature and other known image-dependent parameters (i.e. parameters which differ from image to image, but are fixed for all stars in a given image) for the initial \"airmass\"-like terms, and parameters such as the ra/dec, star colors, and other known star-dependent parameters (i.e. parameters which differ from star to star, but are fixed for a given star over all images) for the initial \"color\"-like terms.\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-TFA",4) && strlen(c) == 4)
    {
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list_tfa -oneline -rms \\\n");
      printtostring(&s,
		    "\t-TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 0 0\n\n");
      printtostring(&s,
		    "Apply the trend-filtering algorithm to the light curves in the list EXAMPLES/lc_list_tfa (the only light curve in that list is EXAMPLES/3.transit). We use the list of template trend light-curves given in EXAMPLES/trendlist_tfa (the second and third columns are the x and y positions for each star, we use the xycol keyword to explicitly state this, though for this particular set of vartools commands this would have been the default). We have TFA only include trend stars that are more than 25 pixels from the source (the coordinates for EXAMPLES/3.transit are given in EXAMPLES/lc_list_tfa). We pass the corrected light curve to the subsequent commands and do not output the TFA coefficients or the TFA model.\n\n");
      commandfound = 1;
    }
  if(!strncmp(c,"-TFA_SR",7) && strlen(c) == 7)
    {
      printtostring(&s,
		    "\nExample 1:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list_tfa_sr_harm -oneline -rms \\\n");
      printtostring(&s,
		    "\t-LS 0.1 10. 0.1 1 0 \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 0 0 \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\t0 0.001 100 harm 0 0 period ls \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 nameformat 2.test_tfa_sr_harm \\\n");
      printtostring(&s,
		    "\t-Killharm ls 0 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \n\n");
      printtostring(&s,
		    "Example illustrating the use of signal-reconstruction TFA with a harmonic series for the model signal. We process light curves in the list EXAMPLES/lc_list_tfa_sr_harm (EXAMPLES/2 is the only light curve in this list). We first use -LS to find the period of the signal. We save a copy of the light curve and then use the -Killharm command followed by the -rms command to determine the amplitude of the signal and the RMS of the residual. We restore the light curve to its initial state (this is needed because -Killharm removes the signal from the light curve when the \"fitonly\" keyword is not used). We then filter the light curve using the -TFA command (no signal-reconstruction; see \"vartools -TFA\" or \"vartools -example -TFA\" for an explanation of the syntax). After applying -TFA we run -Killharm and -rms to determine the amplitude of the signal in the trend-filtered light curve and the RMS of the residual. Note that the amplitude is lower, and the residual RMS is higher after applying -TFA than before applying -TFA. This is because non-reconstructive TFA filters the signal in addition to filtering the noise from the light curve. We restore the light curve to its initial (pre-filtered) state using the -restorelc command, and then apply signal-reconstruction TFA. The syntax for the -TFA_SR command is similar to that for the -TFA command, however this time we choose to output the TFA coefficients and the TFA model to the directory EXAMPLES/OUTDIR1. For the signal-reconstruction specific parameters we set dotfafirst to 0, tfathresh to 0.001 and maxiter to 100 (these parameters only matter if we are using \"decorr\", \"bin\" or \"signal\" which can be used to iteratively fit the signal and filter the light curve, for the \"harm\" mode there is no iteration). We then give the \"harm\" keyword to indicate that we will be using a harmonic series for the signal. We use 0 harmonics and 0 sub-harmonics (only the fundamental is used, i.e. we fit a simple sine-curve). We specify \"period\" to indicate that we wish to fit a period other than the time-span of the light curve, and we specify \"ls\" to have the period taken from the previous -LS command. We output the filtered light curve to the file EXAMPLES/OUTDIR1/2.test_tfa_sr_harm, and we then use the -Killharm and -rms commands to determine the amplitude of the signal after -TFA_SR and the residual RMS. In this case the amplitude is comparable to the value before filtering, but the residual RMS is lower than without filtering.\n\n");
      printtostring(&s,
		    "Example 2:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list_tfa_sr_bin -oneline -rms \\\n");
      printtostring(&s,
		    "\t-aov Nbin 20 0.1 10. 0.1 0.01 1 0 \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-Killharm aov 5 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 0 0 \\\n");
      printtostring(&s,
		    "\t-Killharm aov 5 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\t0 0.001 100 bin 100 period aov \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 nameformat %%s.test_tfa_sr_bin \\\n");
      printtostring(&s,
		    "\t-Killharm aov 5 0 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \n\n");
      printtostring(&s,
		    "This example is the same as Example 1, except that this time we use phase-binning to determine the signal to preserve while filtering the noise, rather than fitting a harmonic series. In addition to filtering EXAMPLES/2, we also filter EXAMPLES/3.starspot. In this example we use the -aov command to search for the period rather than -LS, and we include 5 higher-order harmonics in the -Killharm commands (to give a somewhat better fit to the non-sinusoidal starspot signal). In this case the values for dotfafirst tfathresh and maxiter which are supplied to the -TFA_SR command matter. In this case -TFA_SR will iterate between binning the light curve to define the signal, and applying -TFA to the residual light curve to filter the noise. Setting dotfafirst=0 has the signal measured first before applying TFA (this is generally the best approach), we set tfathresh=0.001 and maxiter=100 to set the conditions for when the iteration will terminate. We use 100 bins and do the binning in phase rather than time (by giving the \"period\" keyword), and we use \"aov\" is the source for the period.\n\n");
      printtostring(&s,
		    "Example 3:\n");
      printtostring(&s,
		    "----------\n");
      printtostring(&s,
		    "\nvartools -l EXAMPLES/lc_list_tfa_sr_decorr -oneline -rms \\\n");
      printtostring(&s,
		    "\t-savelc \\\n");
      printtostring(&s,
		    "\t-decorr 1 1 1 0 1 1 2 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 0 0 \\\n");
      printtostring(&s,
		    "\t-decorr 1 1 1 0 1 1 2 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1 \\\n");
      printtostring(&s,
		    "\t-TFA_SR EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa \\\n");
      printtostring(&s,
		    "\t\tdecorr 0 1 1 2 \\\n");
      printtostring(&s,
		    "\t\t25.0 xycol 2 3 1 \\\n");
      printtostring(&s,
		    "\t\t1 EXAMPLES/OUTDIR1 1 EXAMPLES/OUTDIR1 \\\n");
      printtostring(&s,
		    "\t\t0 0.001 100 bin 100 \\\n");
      printtostring(&s,
		    "\t-o EXAMPLES/OUTDIR1 nameformat %%s.test_tfa_sr_decorr \\\n");
      printtostring(&s,
		    "\t-decorr 1 1 1 0 1 1 2 0 \\\n");
      printtostring(&s,
		    "\t-rms -restorelc 1\n\n");
      printtostring(&s,
		    "This example illustrates the use of simultaneous decorrelation against light-curve specific trends. In this case we process the light curve EXAMPLES/1 and give the \"decorr\" keyword to the -TFA_SR command. We set iterativeflag=0 after the \"decorr\" keyword to have the routine simultaneously fit the light-curve specific trends and the TFA templates given in EXAMPLES/trendlist_tfa (this is more correct than if we iterated, but would run significantly slower than the iterative procedure if we were processing several light curves). We read in one column to decorrelate from the light curve (set Nlcterms=1), we take that to be the first column in the light curve (set lccolumn1=1; the JD in this case) and we fit a 2nd order polynomial in that term (set lcorder1=2). We still have to provide a source for the signal in addition to the decorrelation, in this case we use binning with 100 bins (as in Example 2), but we do not provide a period (so the binning is done in time rather than phase). In this example we use the -decorr command rather than -aov and -Killharm to illustrate how -TFA vs. -TFA_SR affects the light curve. Note that -TFA_SR does not reduce the signal, the way -TFA does, but does reduce the residual RMS compared with not applying trend-filtering at all.\n\n");
      commandfound=1;
    }
  if(!strcmp(c,"-wwz"))
    {
      printtostring(&s,
		    "\nvartools -i EXAMPLES/8 -oneline \\\n");
      printtostring(&s,
		    "\t-wwz maxfreq 2.0 freqsamp 0.25 tau0 auto tau1 auto dtau 0.1 \\\n");
      printtostring(&s,
		    "\t\toutfulltransform EXAMPLES/OUTDIR1/ pm3d \\\n");
      printtostring(&s,
		    "\t\toutmaxtransform EXAMPLES/OUTDIR1\n\n");
      printtostring(&s,
		    "Calculate the Weighted-Wavelet Z-Transform for the light curve stored in EXAMPLES/8. The transform will be calculated for frequencies between 0 and 2.0 cycles per day, using a frequency step-size of 0.25/T where T is the difference between the last and first observation times in the light curve. We will consider time-shifts running from the first observation to the last (done by setting both tau0 and tau1 to auto) in steps of 0.1 days. The full 2-d transform will be output to the file EXAMPLES/OUTDIR1/8.wwz, while the maximum Z frequency for each time-shift (i.e. the projection of the wavelet onto the time-shift axis) will be output to the file EXAMPLES/OUTDIR1/8.mwwz. For the full 2-d transform we will output the file in the gnuplot \"pm3d\" format. A plot of the full transform may then be generated in gnuplot using the commands:\n\n");
      printtostring(&s,
	"gnuplot> set pm3d map\n");
      printtostring(&s,
	"gnuplot> unset key\n");
      printtostring(&s,
	"gnuplot> splot \"EXAMPLES/OUTDIR1/8.wwz\" u 1:2:3\n\n");
      printtostring(&s,
	"A significant periodic signal is seen in this light curve around the time 53735.17392 with a frequency of 0.30645 cycles per day. This signal is not present near the beginning or end of the time series.\n");
      commandfound=1;
    }

#ifdef DYNAMICLIB
  if(!commandfound) {
    for(i=0; i < p->NUserLib; i++) {
      if(!strcmp(c,p->UserLib[i].commandname)) {
	if(p->UserLib[i].ShowExample_function != NULL) {
	  p->UserLib[i].ShowExample_function(stdout);
	  commandfound = 1;
	  break;
	}
	commandfound = 0;
	break;
      }
    }
  }
#endif

  if(!commandfound)
    error2(ERR_EXAMPLE_BADCOMMAND,c);
  printtostring(&s, "\n");
  fprintf(stderr,s.s);
  if(s.s != NULL)
    free(s.s);
  exit(ERR_USAGE);
}

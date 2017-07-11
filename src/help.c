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
/*     This file is part of VARTOOLS version 1.152                      */
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

/* Functions for displaying syntax and help for commands - these are part of the vartools program by J. Hartman */

int listcommands_noexit(char *c, ProgramData *p, OutText *s)
{
  int commandfound = 0;
  int i;

  if(c == NULL || (!strncmp(c,"-addnoise",9) && strlen(c) == 9))
    {
      printtostring(s,
		    "-addnoise\n");
      printtostring(s,
		    "\t<   \"white\"\n");
      printtostring(s,
		    "\t\t\t<\"sig_white\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t  | \"squareexp\"\n");
      printtostring(s,
		    "\t\t\t<\"rho\" <\"fix\" val | \"list\" [\"column\" col]>>\n"); 
      printtostring(s,
		    "\t\t\t<\"sig_red\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t<\"sig_white\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t[\"bintime\" <\"fix\" val | \"list\" [\"column\" col]>]\n");
      printtostring(s,
		    "\t  | \"exp\"\n");
      printtostring(s,
		    "\t\t\t<\"rho\" <\"fix\" val | \"list\" [\"column\" col]>>\n"); 
      printtostring(s,
		    "\t\t\t<\"sig_red\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t<\"sig_white\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t[\"bintime\" <\"fix\" val | \"list\" [\"column\" col]>]\n");
      printtostring(s,
		    "\t  | \"matern\"\n");
      printtostring(s,
		    "\t\t\t<\"nu\" <\"fix\" val | \"list\" [\"column\" col]>>\n"); 
      printtostring(s,
		    "\t\t\t<\"rho\" <\"fix\" val | \"list\" [\"column\" col]>>\n"); 
      printtostring(s,
		    "\t\t\t<\"sig_red\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t<\"sig_white\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t[\"bintime\" <\"fix\" val | \"list\" [\"column\" col]>]\n");
#ifdef _HAVE_GSL
      printtostring(s,
		    "\t  | \"wavelet\n");
      printtostring(s,
		    "\t\t\t<\"gamma\" <\"fix\" val | \"list\" [\"column\" col]>>\n"); 
      printtostring(s,
		    "\t\t\t<\"sig_red\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
      printtostring(s,
		    "\t\t\t<\"sig_white\" <\"fix\" val | \"list\" [\"column\" col]>>\n");
#endif
      printtostring(s,
	"\t>\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-alarm",6) && strlen(c) == 6))
    {
      printtostring(s, "-alarm\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-aov",4) && strlen(c) == 4))
    {
      printtostring(s,
		    "-aov [\"Nbin\" Nbin] minp maxp subsample finetune Npeaks operiodogram\n");
      printtostring(s,
		    "\t[outdir] [\"whiten\"] [\"clip\" clip clipiter] [\"uselog\"]\n");
      printtostring(s,
		    "\t[\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\" | \"fix\" period\n");
      printtostring(s,
		    "\t\t| \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t| \"fixcolumn\" <colname | colnum>>]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-aov_harm",9) && strlen(c) == 9))
    {
      printtostring(s,"-aov_harm Nharm minp maxp subsample finetune Npeaks operiodogram [outdir]\n");
      printtostring(s,
		    "\t[\"whiten\"] [\"clip\" clip clipiter]\n");
      printtostring(s,
		    "\t[\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\" | \"fix\" period\n");
      printtostring(s,
		    "\t\t| \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t| \"fixcolumn\" <colname | colnum>>]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-autocorrelation",16) && strlen(c) == 16))
    {
      printtostring(s,"-autocorrelation start stop step outdir\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-binlc",6) && strlen(c) == 6))
    {
      printtostring(s,"-binlc <\"average\" | \"median\" | \"weightedaverage\">\n");
      printtostring(s,"\t<\"binsize\" binsize | \"nbins\" nbins>\n");
      printtostring(s,"\t[\"bincolumns\" var1[:stats1][,var2[:stats2],...]]\n");
      printtostring(s,"\t[\"firstbinshift\" firstbinshift]\n");
      printtostring(s,"\t<\"tcenter\" | \"taverage\" | \"tmedian\" | \"tnoshrink\" [\"bincolumnsonly\"]>\n");
      
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-BLS",4) && strlen(c) == 4))
    {
      printtostring(s,"-BLS < \"r\" rmin rmax | \"q\" qmin qmax\n");
      printtostring(s,"\t\t| \"density\" rho minexpdurfrac maxexpdurfrac >\n");
      printtostring(s,"\tminper maxper nfreq nbins\n");
      printtostring(s,"\ttimezone Npeak outperiodogram [outdir] omodel [model_outdir]\n");
      printtostring(s,"\tcorrectlc [\"fittrap\"] [\"nobinnedrms\"]\n");
      printtostring(s,"\t[\"ophcurve\" outdir phmin phmax phstep]\n");
      printtostring(s,"\t[\"ojdcurve\" outdir jdstep]\n");
      printtostring(s,"\t[\"stepP\" | \"steplogP\"]\n");
      printtostring(s,"\t[\"adjust-qmin-by-mindt\" [\"reduce-nbins\"]]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-BLSFixPer",10) && strlen(c) == 10))
    {
      printtostring(s,"-BLSFixPer <\"aov\" | \"ls\" | \"list\" [\"column\" col]\n");
      printtostring(s,"\t\t| \"fix\" period |  \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,"\t\t| \"expr\" expr>\n");
      printtostring(s,"\t<\"r\" rmin rmax | \"q\" qmin qmax >\n");
      printtostring(s,"\tnbins timezone omodel [model_outdir] correctlc [\"fittrap\"]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-BLSFixDurTc")))
    {
      printtostring(s,"-BLSFixDurTc\n");
      printtostring(s,"\t<\"duration\" <\"fix\" dur | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,"\t\t| \"list\" [\"column\" col]>>\n");
      printtostring(s,"\t<\"Tc\" <\"fix\" Tc | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,"\t\t| \"list\" [\"column\" col]>>\n");
      printtostring(s,"\t[\"fixdepth\" <\"fix\" depth | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,"\t\t| \"list\" [\"column\" col]>\n");
      printtostring(s,"\t\t[\"qgress\" <\"fix\" qgress | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,"\t\t\t| \"list\" [\column\" col]>]]\n");
      printtostring(s,"\tminper maxper nfreq timezone\n");
      printtostring(s,"\tNpeak outperiodogram [outdir] omodel [model_outdir]\n");
      printtostring(s,"\tcorrectlc [\"fittrap\"]\n");
      printtostring(s,"\t[\"ophcurve\" outdir phmin phmax phstep]\n");
      printtostring(s,"\t[\"ojdcurve\" outdir jdstep]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-changeerror",12) && strlen(c) == 12))
    {
      printtostring(s,"-changeerror\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-changevariable")))
    {
      printtostring(s,"-changevariable <\"t\" | \"mag\" | \"err\" | \"id\"> var\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-chi2",5) && strlen(c) == 5))
    {
      printtostring(s,"-chi2\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-chi2bin",8) && strlen(c) == 8))
    {
      printtostring(s,"-chi2bin Nbin bintime1...bintimeN\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-clip",5) && strlen(c) == 5))
    {
      printtostring(s,"-clip sigclip iter [\"niter\" n] [\"median\"]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-converttime",12) && strlen(c) == 12))
    {
      printtostring(s,"-converttime\n");
      printtostring(s,"\t<\"input\" <\"mjd\" | \"jd\" | \"hjd\" | \"bjd\" >>\n");
      printtostring(s,"\t[\"inputsubtract\" value] [\"inputsys-tdb\" | \"inputsys-utc\"]\n");
      printtostring(s,"\t<\"output\" <\"mjd\" | \"jd\" | \"hjd\" | \"bjd\" >>\n");
      printtostring(s,"\t[\"outputsubtract\" value] [\"outputsys-tdb\" | \"outputsys-utc\"]\n");
      printtostring(s,"\t[\"radec\" <\"list\" [\"column\" col] | \"fix\" raval decval>\n");
      printtostring(s,"\t[\"epoch\" epoch]]\n");
      printtostring(s,"\t[\"ppm\" <\"list\" [\"column\" col] | \"fix\" mu_ra mu_dec>]\n");
      printtostring(s,"\t[\"input-radec\" <\"list\" [\"column\" col] | \"fix\" raval decval>\n");
      printtostring(s,"\t[\"epoch\" epoch]]\n");
      printtostring(s,"\t[\"input-ppm\" <\"list\" [\"column\" col] | \"fix\" mu_ra mu_dec>]\n");      
      printtostring(s,"\t[\"ephemfile\" file] [\"leapsecfile\" file] [\"planetdatafile\" file]\n");
      printtostring(s,"\t[\"observatory\" < code | \"show-codes\">\n");
      printtostring(s,"\t\t| \"coords\"\n");
      printtostring(s,"\t\t\t<\"fix\" latitude[deg] longitude[deg_E] altitude[m]\n");
      printtostring(s,"\t\t\t| \"list\" [\"column\" collat collong colalt]\n");
      printtostring(s,"\t\t\t| \"fromlc\" collat collong colalt>]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-copylc")))
    {
      printtostring(s,"-copylc Ncopies\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-decorr",7) && strlen(c) == 7))
    {
      printtostring(s,"-decorr correctlc zeropointterm subtractfirstterm Nglobalterms globalfile1\n");
      printtostring(s,"\torder1 ... Nlcterms lccolumn1 lcorder1 ... omodel\n");
      printtostring(s,"\t[model_outdir]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-dftclean",9) && strlen(c) == 9))
    {
      printtostring(s,"-dftclean nbeam [\"maxfreq maxf\"] [\"outdspec\" dspec_outdir]\n");
      printtostring(s,"\t[\"finddirtypeaks\" Npeaks [\"clip\" clip clipiter]]\n");
      printtostring(s,"\t[\"outwfunc\" wfunc_outdir]\n");
      printtostring(s,"\t[\"clean\" gain SNlimit [\"outcbeam\" cbeam_outdir]\n");
      printtostring(s,"\t[\"outcspec\" cspec_outdir]\n");
      printtostring(s,"\t[\"findcleanpeaks\" Npeaks [\"clip\" clip clipiter]]]\n");
      printtostring(s,"\t[\"useampspec\"] [\"verboseout\"]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-difffluxtomag",14) && strlen(c) == 14))
    {
      printtostring(s,"-difffluxtomag mag_constant offset [\"magcolumn\" col]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-ensemblerescalesig",19) && strlen(c) == 19))
    {
      printtostring(s,"-ensemblerescalesig sigclip\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-expr")))
    {
      printtostring(s,"-expr var\"=\"expression\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-findblends",11) && strlen(c) == 11))
    {
      printtostring(s,"-findblends matchrad [\"radec\"]\n");
      printtostring(s,"\t[\"xycol\" xcol ycol]\n");
      printtostring(s,"\t<\"fix\" period | \"list\" [\"column\" col]\n");
      printtostring(s,"\t\t| \"fixcolumn\" <colname | colnum>>\n");
      printtostring(s,"\t[\"starlist\" starlistfile] [\"zeromag\" zeromagval] [\"nofluxconvert\"]\n");
      printtostring(s,"\t[\"Nharm\" Nharm] [\"omatches\" outputmatchfile]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-fluxtomag",10) && strlen(c) == 10))
    {
      printtostring(s,"-fluxtomag mag_constant offset\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-GetLSAmpThresh",15) && strlen(c) == 15))
    {
      printtostring(s,"-GetLSAmpThresh <\"ls\" | \"list\" [\"column\" col]> minp thresh\n");
      printtostring(s,"\t<\"harm\" Nharm Nsubharm | \"file\" listfile> [\"noGLS\"]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-if")) || (!strcmp(c,"-elif")) ||
     (!strcmp(c,"-else")) || (!strcmp(c,"-fi")))
    {
      printtostring(s,"-if <expression> [-command1 ... -commandN]\n");
      printtostring(s,"\t[-elif <expression> [-command1 ... -commandN]]\n");
      printtostring(s,"\t...\n");
      printtostring(s,"\t[-elif <expression> [-command1 ... -commandN]]\n");
      printtostring(s,"\t[-else [-command1 ... -commandN]]\n");
      printtostring(s,"\t-fi\n");
      commandfound=1;
    }
  if(c == NULL || (!strncmp(c,"-Injectharm",11) && strlen(c) == 11))
    {
      printtostring(s,"-Injectharm <\"list\" [\"column\" col] | \"fix\" per\n");
      printtostring(s,"\t| \"rand\" minp maxp\n");
      printtostring(s,"\t| \"logrand\" minp maxp | \"randfreq\" minf maxf\n");
      printtostring(s,"\t| \"lograndfreq\" minf maxf>\n");
      printtostring(s,"\tNharm (<\"amplist\" [\"column\" col]\n");
      printtostring(s,"\t| \"ampfix\" amp | \"amprand\" minamp maxamp \n");
      printtostring(s,"\t| \"amplogrand\" minamp maxamp> [\"amprel\"]\n");
      printtostring(s,"\t<\"phaselist\" [\"column\" col]\n");
      printtostring(s,"\t| \"phasefix\" phase | \"phaserand\"> [\"phaserel\"])0...Nharm Nsubharm\n");
      printtostring(s,"\t(<\"amplist\" [\"column\" col] | \"ampfix\" amp\n");
      printtostring(s,"\t| \"amprand\" minamp maxamp \n");
      printtostring(s,"\t| \"amplogrand\" minamp maxamp> [\"amprel\"]\n");
      printtostring(s,"\t<\"phaselist\" [\"column\" col]\n");
      printtostring(s,"\t| \"phasefix\" phase | \"phaserand\"> [\"phaserel\"])1...Nsubharm\n\tomodel [modeloutdir]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-Injecttransit",14) && strlen(c) == 14))
    {
      printtostring(s,"-Injecttransit <\"Plist\" [\"column\" col] | \"Pfix\" per\n");
      printtostring(s,"\t\t| \"Pexpr\" expr | \"Prand\" minp maxp\n");
      printtostring(s,"\t\t| \"Plogrand\" minp maxp | \"randfreq\" minf maxf\n");
      printtostring(s,"\t\t| \"lograndfreq\" minf maxf>\n");
      printtostring(s,"\t<\"Rplist\" [\"column\" col] | \"Rpfix\" Rp | \"Rpexpr\" expr\n");
      printtostring(s,"\t\t| \"Rprand\" minRp maxRp | \"Rplogrand\" minRp maxRp>\n");
      printtostring(s,"\t<\"Mplist\" [\"column\" col] | \"Mpfix\" Mp | \"Mpexpr\" expr\n");
      printtostring(s,"\t\t| \"Mprand\" minMp maxMp | \"Mplogrand\" minMp maxMp>\n");
      printtostring(s,"\t<\"phaselist\" [\"column\" col] | \"phasefix\" phase\n");
      printtostring(s,"\t\t| \"phasexpr\" expr | \"phaserand>\n");
      printtostring(s,"\t<\"sinilist\" [\"column\" col] | \"sinifix\" sin_i\n");
      printtostring(s,"\t\t| \"siniexpr\" expr | \"sinirand\">\n");
      printtostring(s,"\t<\"eomega\" <\"elist\" [\"column\" col] | \"efix\" e \"expr\" expr | \"erand\">\n");
      printtostring(s,"\t\t<\"olist\" [\"column\" col] | \"ofix\" omega | \"oexpr\" expr | \"orand\">\n");
      printtostring(s,"\t| \"hk\" <\"hlist\" [\"column\" col] | \"hfix\" h | \"hexpr\" expr | \"hrand\">\n");
      printtostring(s,"\t\t<\"klist\" [\"column\" col] | \"kfix\" k | \"kexpr\" expr | \"krand\">>\n");      
      printtostring(s,"\t<\"Mstarlist\" [\"column\" col] | \"Mstarfix\" Mstar | \"Mstarexpr\" expr>\n");
      printtostring(s,"\t<\"Rstarlist\" [\"column\" col] | \"Rstarfix\" Rstar | \"Rstarexpr\" expr>\n");
      printtostring(s,"\t<\"quad\" | \"nonlin\"> <\"ldlist\" [\"column\" col]\n");
      printtostring(s,"\t| \"ldfix\" ld1 ... ldn | \"ldexpr\" ld1 ... ldn>\n");
      printtostring(s,"\t[\"dilute\" <\"list\" [\"column\" col] | \"fix\" dilute | \"expr\" dilutexpr>]\n");
      printtostring(s,"\t omodel [modeloutdir]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-Jstet",6) && strlen(c) == 6))
    {
      printtostring(s,"-Jstet timescale dates\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-Killharm",9) && strlen(c) == 9))
    {
      printtostring(s,"-Killharm <\"aov\" | \"ls\" | \"both\" | \"injectharm\" \n");
      printtostring(s,"\t| \"fix\" Nper per1 ... perN\n");
      printtostring(s,"\t| \"list\" Nper [\"column\" col1]> Nharm Nsubharm\n");
      printtostring(s,"\tomodel [model_outdir] [\"fitonly\"]\n");
      printtostring(s,"\t[\"outampphase\" | \"outampradphase\" | \"outRphi\" | \"outRradphi\"]\n");
      printtostring(s,"\t[\"clip\" val]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-linfit")))
    {
      printtostring(s,
		    "-linfit function paramlist [\"modelvar\" varname]\n");
      printtostring(s,
		    "\t[\"correctlc\"]\n");
      printtostring(s,
		    "\t[\"omodel\" model_outdir [\"format\" nameformat]]\n");
      commandfound=1;
    }
  if(c == NULL || (!strncmp(c,"-LS",3) && strlen(c) == 3))
    {
      printtostring(s,
		    "-LS minp maxp subsample Npeaks operiodogram [outdir] [\"noGLS\"] [\"whiten\"]\n");
      printtostring(s,
		    "\t[\"clip\" clip clipiter] [\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\"\n");
      printtostring(s,
		    "\t\t| \"fix\" period | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t| \"fixcolumn\" <colname | colnum>>]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-MandelAgolTransit",18) && strlen(c) == 18))
    {
      printtostring(s,"-MandelAgolTransit <bls | blsfixper\n");
      printtostring(s,"\t\t| P0 T00 r0 a0 <\"i\" inclination | \"b\" bimpact> e0 omega0 mconst0>\n");
      printtostring(s,"\t<\"quad\" | \"nonlin\"> ldcoeff1_0 ... ldcoeffn_0 fitephem\n");
      printtostring(s,"\tfitr fita fitinclterm fite fitomega fitmconst fitldcoeff1 ... fitldcoeffn\n");
      printtostring(s,"\tfitRV [RVinputfile RVmodeloutfile K0 gamma0 fitK fitgamma]\n");
      printtostring(s,"\tcorrectlc omodel [model_outdir] [\"modelvar\" var]\n");
      printtostring(s,"\t[\"ophcurve\" curve_outdir phmin phmax phstep]\n");
      printtostring(s,"\t[\"ojdcurve\" curve_outdir jdstep]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-medianfilter",13) && strlen(c) == 13))
    {
      printtostring(s,"-medianfilter time [\"average\" | \"weightedaverage\"] [\"replace\"]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-microlens",10) && strlen(c) == 10))
    {
      printtostring(s,
		    "-microlens\n");
      printtostring(s,
		    "\t[\"f0\"\n");
      printtostring(s,
		    "\t\t[\"fix\" fixval | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t\t| \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t\t| \"auto\"]\n");
      printtostring(s,
		    "\t\t[\"step\" initialstepsize] [\"novary\"]]\n");
      printtostring(s,
		    "\t[\"f1\"\n");
      printtostring(s,
		    "\t\t[\"fix\" fixval | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t\t| \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t\t| \"auto\"]\n");
      printtostring(s,
		    "\t\t[\"step\" initialstepsize] [\"novary\"]]\n");
      printtostring(s,
		    "\t[\"u0\"\n");
      printtostring(s,
		    "\t\t[\"fix\" fixval | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t\t| \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t\t| \"auto\"]\n");
      printtostring(s,
		    "\t\t[\"step\" initialstepsize] [\"novary\"]]\n");
      printtostring(s,
		    "\t[\"t0\"\n");
      printtostring(s,
		    "\t\t[\"fix\" fixval | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t\t| \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t\t| \"auto\"]\n");
      printtostring(s,
		    "\t\t[\"step\" initialstepsize] [\"novary\"]]\n");
      printtostring(s,
		    "\t[\"tmax\"\n");
      printtostring(s,
		    "\t\t[\"fix\" fixval | \"list\" [\"column\" col]\n");
      printtostring(s,
		    "\t\t\t| \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t\t| \"auto\"]\n");
      printtostring(s,
		    "\t\t[\"step\" initialstepsize] [\"novary\"]]\n");
      printtostring(s,
		    "\t[\"correctlc\"] [\"omodel\" outdir]\n\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-nonlinfit")))
    {
      printtostring(s,
		    "-nonlinfit function paramlist [\"linfit\" linfitparams]\n");
      printtostring(s,
		    "\t[\"errors\" error_expr]\n");
      printtostring(s,
		    "\t[\"covariance\"\n");
      printtostring(s,
		    "\t\t<\"squareexp\" amp_var rho_var\n");
      printtostring(s,
		    "\t\t | \"exp\" amp_var rho_var\n");
      printtostring(s,
		    "\t\t | \"matern\" amp_var rho_var nu_var>]\n");
      printtostring(s,
		    "\t[\"priors\" priorlist] [\"constraints\" constraintlist]\n");
      printtostring(s,
		    "\t<\"amoeba\" [\"tolerance\" tol] [\"maxsteps\" steps]\n");
      printtostring(s,
		    "\t | \"mcmc\" [\"Naccept\" N | \"Nlinkstotal\" N]\n");
      printtostring(s,
		    "\t\t\t[\"fracburnin\" frac] [\"eps\" eps] [\"skipamoeba\"]\n");
      printtostring(s,
		    "\t\t\t[\"chainstats\" exprlist statslist]\n");
      printtostring(s,
		    "\t\t\t[\"maxmemstore\" maxmem]\n");
      printtostring(s,
		    "\t\t\t[\"outchains\" outdir [\"format\" format] [\"printevery\" N]] >\n");
      printtostring(s,
		    "\t[\"modelvar\" varname] [\"correctlc\"]\n");
      printtostring(s,
		    "\t[\"omodel\" model_outdir [\"format\" nameformat]]\n");
      commandfound=1;
    }
  if(c == NULL || (!strncmp(c,"-o",2) && strlen(c) == 2))
    {
      printtostring(s,"-o <outdir | outname> [\"nameformat\" formatstring]\n");
      printtostring(s,"\t[\"columnformat\" formatstring] [\"delimiter\" delimchar]\n");
      printtostring(s,"\t[\"fits\"] [\"noclobber\"]\n");
      commandfound = 1;
    }
  if(c == NULL || ((!strncmp(c,"-Phase",6) || !strncmp(c,"-phase",6)) && strlen(c) == 6))
    {
      printtostring(s,
		    "-Phase <\"aov\" | \"ls\" | \"bls\" | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t| \"list\" [\"column\" col] | \"fix\" period>\n"); 
      printtostring(s,
		    "\t[\"T0\" <\"bls\" phaseTc | \"fixcolumn\" <colname | colnum>\n");
      printtostring(s,
		    "\t\t| \"list\" [\"column\" col] | \"fix\" T0>]\n");
      printtostring(s,
		    "\t[\"phasevar\" var] [\"startphase\" startphase]\n");
      commandfound = 1;
    }
#ifdef _HAVE_PYTHON
  if(c == NULL || (!strcmp(c,"-python")))
    {
      printtostring(s,
		    "-python <\"fromfile\" commandfile | commandstring>\n");
      printtostring(s,
		    "\t[\"init\" <\"file\" initializationfile | initializationstring>\n");
      printtostring(s,
		    "\t\t| \"continueprocess\" prior_python_command_number]\n");
      printtostring(s,
		    "\t[\"vars\" variablelist\n");
      printtostring(s,
		    "\t\t| [\"invars\" inputvariablelist] [\"outvars\" outputvariablelist]]\n");
      printtostring(s,
		    "\t[\"outputcolumns\" variablelist] [\"process_all_lcs\"]\n");
      commandfound = 1;
    }
#endif
  if(c == NULL || (!strcmp(c,"-resample")))
    {
      printtostring(s,
		    "-resample\n");
      printtostring(s,
		    "\t<\"nearest\" |\n");
      printtostring(s,
                    "\t  \"linear\"  |\n");
      printtostring(s,
                    "\t  \"spline\"  [\"left\" yp1] [\"right\" ypn] |\n");
      printtostring(s,
                    "\t  \"splinemonotonic\" |\n");
      printtostring(s,
                    "\t  \"bspline\" [\"nbreaks\" nbreaks] [\"order\" order] >\n");
      printtostring(s,
		    "\t[\"file\" <\"fix\" times_file [\"column\" time_column] |\n");
      printtostring(s,
                    "\t\t\"list\" [\"listcolumn\" col] [\"tcolumn\" time_column] > |\n");
      printtostring(s,
                    "\t[\"tstart\" <\"fix\" tstart | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,
		    "\t\t\"list\" [\"column\" col] | \"expr\" expression > ]\n");
      printtostring(s,
                    "\t[\"tstop\" <\"fix\" tstop | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,
                    "\t\t\"list\" [\"column\" col] | \"expr\" expression > ]\n");
      printtostring(s,
                    "\t[[\"delt\" <\"fix\" delt | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,
                    "\t\t\"list\" [\"column\" col] | \"expr\" expression > ]\n");
      printtostring(s,
                    "\t | [\"Npoints\" <\"fix\" Np | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,
                    "\t\t\"list\" [\"column\" col] | \"expr\" expression > ]]]\n");
      printtostring(s,
		    "\t[\"gaps\" \n");
      printtostring(s,
		    "\t\t<\"fix\" time_sep | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,
		    "\t\t\t\"list\" [\"column\" col] | \"expr\" expression |\n");
      printtostring(s,
		    "\t\t\t\"frac_min_sep\" val | \"frac_med_sep\" val | \"percentile_sep\" val>\n");
      printtostring(s,
		    "\t\t<\"nearest\" |\n");
      printtostring(s,
                    "\t\t  \"linear\"  |\n");
      printtostring(s,
                    "\t\t  \"spline\"  [\"left\" yp1] [\"right\" ypn] |\n");
      printtostring(s,
                    "\t\t  \"splinemonotonic\" |\n");
      printtostring(s,
                    "\t\t  \"bspline\" [\"nbreaks\" nbreaks] [\"order\" order] >]\n");
      printtostring(s,
		    "\t[\"extrap\" \n");
      printtostring(s,
		    "\t\t<\"nearest\" |\n");
      printtostring(s,
                    "\t\t  \"linear\"  |\n");
      printtostring(s,
                    "\t\t  \"spline\"  [\"left\" yp1] [\"right\" ypn] |\n");
      printtostring(s,
                    "\t\t  \"splinemonotonic\" |\n");
      printtostring(s,
                    "\t\t  \"bspline\" [\"nbreaks\" nbreaks] [\"order\" order] >]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-rescalesig",11) && strlen(c) == 11))
    {
      printtostring(s,"-rescalesig\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-restorelc",10) && strlen(c) == 10))
    {
      printtostring(s,"-restorelc savenumber\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-restricttimes")))
    {
      printtostring(s,"-restricttimes [\"exclude\"]\n");
      printtostring(s,"\t< \"JDrange\" minJD maxJD |\n");
      printtostring(s,"\t  \"JDrangebylc\"\n");
      printtostring(s,"\t\t<\"fix\" minJD | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,"\t\t \"expr\" expression>\n");
      printtostring(s,"\t\t<\"fix\" maxJD | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> |\n");
      printtostring(s,"\t\t \"expr\" expression> |\n");
      printtostring(s,"\t  \"JDlist\" JDfilename | \n");
      printtostring(s,"\t  \"imagelist\" imagefilename | \n");
      printtostring(s,"\t  \"expr\" expression>\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-restoretimes")))
    {
      printtostring(s,"-restoretimes prior_restricttimes_command\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-rms",4) && strlen(c) == 4))
    {
      printtostring(s,"-rms\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-rmsbin",7) && strlen(c) == 7))
    {
      printtostring(s,"-rmsbin Nbin bintime1...bintimeN\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-savelc",7) && strlen(c) == 7))
    {
      printtostring(s,"-savelc\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-SoftenedTransit",16) && strlen(c) == 16))
    {
      printtostring(s,"-SoftenedTransit <bls | blsfixper | P0 T00 eta0 delta0 mconst0 cval0>\n");
      printtostring(s,"\tfitephem fiteta fitcval fitdelta fitmconst correctlc\n");
      printtostring(s,"\tomodel [model_outdir] fit_harm [<\"aov\" | \"ls\" | \"bls\" \n");
      printtostring(s,"\t| \"list\" [\"column\" col] | \"fix\" Pharm> nharm nsubharm]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-Starspot",9) && strlen(c) == 9))
    {
      printtostring(s,"-Starspot\n");
      printtostring(s,
		    "\t<aov | ls | list [\"column\" col] | \"fix\" period |\n");
      printtostring(s,"\t\t\"fixcolumn\" <colname | colnum>>\n");
      printtostring(s,"\ta0 b0 alpha0 i0 chi0 psi00 mconst0 fitP fita fitb\n");
      printtostring(s,"\tfitalpha fiti fitchi fitpsi fitmconst correctlc omodel [model_outdir]\n");
      commandfound = 1;
    }
  if(c == NULL || !strcmp(c,"-stats"))
    {
      printtostring(s,"-stats var1,var2,... stats1,stats2,...\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-SYSREM",7) && strlen(c) == 7))
    {
      printtostring(s,"-SYSREM Ninput_color [\"column\" col1] Ninput_airmass initial_airmass_file\n");
      printtostring(s,"\tssima_clip1 sigma_clip2 saturation correctlc omodel [model_outdir]\n");
      printtostring(s,"\totrends [trend_outfile] useweights\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-TFA",4) && strlen(c) == 4))
    {
      printtostring(s,"-TFA trendlist [\"readformat\" Nskip jdcol magcol]\n");
      printtostring(s,"\tdates_file pixelsep [\"xycol\" xcol ycol]\n");
      printtostring(s,"\tcorrectlc ocoeff [coeff_outdir] omodel [model_outdir]\n");
      commandfound = 1;
    }
  if(c == NULL || (!strncmp(c,"-TFA_SR",7) && strlen(c) == 7))
    {
      printtostring(s,"-TFA_SR trendlist [\"readformat\" Nskip jdcol magcol] dates_file\n");
      printtostring(s,"\t[\"decorr\" iterativeflag Nlcterms lccolumn1 lcorder1 ...] pixelsep\n");
      printtostring(s,"\t[\"xycol\" colx coly]\n");
      printtostring(s,"\tcorrectlc ocoeff [coeff_outdir] omodel [model_outdir] dotfafirst\n");
      printtostring(s,"\ttfathresh maxiter <\"bin\" nbins [\"period\" <\"aov\" | \"ls\" \n");
      printtostring(s,"\t| \"bls\" | \"list\" [\"column\" col] | \"fix\" period>]\n");
      printtostring(s,"\t| \"signal\" filename\n");
      printtostring(s,"\t| \"harm\" Nharm Nsubharm [\"period\" <\"aov\" | \"ls\" \n");
      printtostring(s,"\t| \"bls\" | \"list\" [\"column\" col] | \"fix\" period>]>\n");
      commandfound = 1;
    }
  if(c == NULL || (!strcmp(c,"-wwz")))
    {
      printtostring(s,"-wwz <\"maxfreq\" <\"auto\" | maxfreq>> <\"freqsamp\" freqsamp>\n");
      printtostring(s,"\t<\"tau0\" <\"auto\" | tau0>> <\"tau1\" <\"auto\" | tau1>>\n");
      printtostring(s,"\t<\"dtau\" <\"auto\" | dtau>> [\"c\" cval]\n");
      printtostring(s,"\t[\"outfulltransform\" outdir");
#ifdef USECFITSIO
      printtostring(s," [\"fits\" | \"pm3d\"]");
#else
      printtostring(s," [\"pm3d\"]");
#endif
      printtostring(s," [\"format\" format]]\n");
      printtostring(s,"\t[\"outmaxtransform\" outdir [\"format\" format]]\n");
      commandfound = 1;
    }
  printtostring(s,"\n");
  return(commandfound);
}

void listcommands(char *c, ProgramData *p)
{
  OutText s;
  int val;
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;
  val = listcommands_noexit(c, p, &s);
  if(!val)
    error2(ERR_HELP_BADCOMMAND,c);
  if(c != NULL)
    fprintf(stderr,"Command syntax:\n\n");
  if(s.s != NULL)
    {
      fprintf(stderr,s.s);
      free(s.s);
    }
    exit(ERR_USAGE);
}

void usage_common(OutText *s)
{
  char s2[MAXLEN];

  printtostring(s, "Usage: vartools <-i lcname ");
#ifdef _USEBINARY_LC
  printtostring(s, " [\"binary\"] ");
#endif
  printtostring(s, "| -l lclist [\"column\" col]");
#ifdef _USEBINARY_LC
  printtostring(s, " [\"binary\"]");
#endif
  printtostring(s, " [\"opencommand\" command]>\n");
  printtostring(s, "\t[-header] [-headeronly] [-tab]\n");
  printtostring(s, "\t[-readformat Nskip [\"stringid\" colstringid] [\"inpututc\" format]\n");
  printtostring(s, "\t\tcol_time col_mag col_sig]\n");
  printtostring(s, "\t[-inputlcformat var1:col1[:type1[:fmt1][,var2:col2[:vtype2[:fmt2]],...]]\n");
  printtostring(s, "\t\t[\"skipnum\" Nskip] [\"skipchar\" chars]>]]\n");
  printtostring(s, "\t[-inlistvars var1:col1[:vtype1[:fmt1][,var2:col2[:vtyp2[:fmt2]],....]]]\n");
  printtostring(s, "\t[-basename] [-readall] [-binaryperiodogram] [-jdtol jdtol] [-matchstringid]\n");
  printtostring(s, "\t[-nobuffer] [-help [\"all\" | commandname]] [-quiet] [-randseed seed]\n");
  printtostring(s, "\t[-bufferlines nlines] [-numbercolumns] [-listcommands [commandname]]\n");
  printtostring(s, "\t[-redirectstats statsfile [\"append\"]] [-oneline]\n");
  printtostring(s, "\t[-example commandname] [-showinputlistformat]\n");
  printtostring(s, "\t");
    printtostring(s, "[-showinputlcformat] ");
    printtostring(s, "[-log-command-line] [-skipmissing] [-noskipempty]");
    printtostring(s, "\n");
#ifdef PARALLEL
  printtostring(s, "\t[-parallel Nproc] [-functionlist]\n");
#endif
#ifdef DYNAMICLIB
  printtostring(s, "\t[-L libraryfile] [-F libraryfile] [-f functodefine]\n");
#endif
#ifdef VARTOOLS_VERSION
    printtostring(s, "\t[-version]\n");
#endif
  printtostring(s, "\tCommand1... CommandN\n\n");
}

void usage(char *argv)
{
  OutText s;
  char s2[MAXLEN];
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  usage_common(&s);
  printtostring(&s, "\tFor a list of commands type:\n");
  printtostring(&s, "\t\tvartools -listcommands\n\n");
  printtostring(&s, "\tTo view the syntax for a specific command type:\n");
  printtostring(&s, "\t\tvartools -listcommands $COMMANDNAME\n");
  printtostring(&s, "\tTo see the syntax for the -LS command, for example, type:\n");
  printtostring(&s, "\t\tvartools -listcommands -LS\n\n");
  printtostring(&s, "\tTo get detailed help type:\n");
  printtostring(&s, "\t\tvartools -help all\n\n");
  printtostring(&s, "\tFor help on a particular command or option type:\n");
  printtostring(&s, "\t\tvartools -help $COMMANDNAME\n");
  printtostring(&s, "\tTo see help on the -readformat option, for example, type:\n");
  printtostring(&s, "\t\tvartools -help -readformat\n\n");
  printtostring(&s, "\tTo see an example use for a particular command type:\n");
  printtostring(&s, "\t\tvartools -example $COMMANDNAME\n\n");
#ifdef VARTOOLS_VERSION
  sprintf(s2,"\tVARTOOLS version %s\n\n",VARTOOLS_VERSION);
  printtostring(&s, s2);
#endif
  fprintf(stderr,s.s);
  if(s.s != NULL)
    free(s.s);
  exit(ERR_USAGE);
}

void help(char *c, ProgramData *p)
{
  int commandfound = 0;
  int i, all;
  OutText s;
  char s2[MAXLEN];
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  if(c == NULL)
    usage(NULL);

  printtostring(&s,"\n");

  if(!strcmp(c,"all"))
    all = 1;
  else
    all = 0;

  if(all) {
    usage_common(&s);
#ifdef VARTOOLS_VERSION
    sprintf(s2,"VARTOOLS version %s\n\n", VARTOOLS_VERSION);
    printtostring(&s, s2);
#endif
    printtostring(&s, "This program provides tools for calculating variability/periodicity statistics of light curves as well as tools for modifying light curves. The program is run by issuing a sequence of commands to perform actions on light curves, each command is executed in turn with the resulting light curves passed to the next command. The results of each command are sent to stdout as an ascii table, or can be output to a file if the -redirectstats option is given. The program can execute on one light curve at a time, or on all light curves at once (for some commands all light curves must be operated on together).\n\n");
    commandfound = 1;
  }
  if(all == 1 || (!strncmp(c,"-i",2) && strlen(c) == 2) || (!strncmp(c,"-l",2) && strlen(c) == 2))
    {
      printtostring(&s,"<-i lcname ");
#ifdef _USEBINARY_LC
      printtostring(&s,"[\"binary\"] ");
#endif
      printtostring(&s,"| -l lclist [\"column\" col]");
#ifdef _USEBINARY_LC
      printtostring(&s," [\"binary\"]");
#endif
      printtostring(&s," [\"opencommand\" command]");
      printtostring(&s,">\n\n");
      printtostring(&s,"Either provide an individual light curve, or a list of light curves. If lcname or lclist is \"-\" then input will be taken from stdin. Some commands can only be used with lists. The list should contain a single light curve filename per line, by default in the first column. You can change the column of the file-list by giving the \"column\" keyword and specifying the column number. Additional columns may be necessary to specify, for example, the periods to use for pre-whitening or any other values required by some commands.");
#ifdef _USEBINARY_LC
      printtostring(&s," Give the optional \"binary\" keyword to either command after the filename if the input light curves are in Penev's binary format. Note that light curves in binary fits format will be identified by the suffix .fits appearing in the filename.");
#endif
      printtostring(&s," Use the \"opencommand\" option to apply a shell command to each light curve file before reading it in. In this case 'command' is a string which is passed to the shell with instances of %%s replaced by the filename as read in from the input list, and instances of %%%% replaced by %%. The command will be executed by the shell (after substituting %%s and %%%%) and the light curve will be read from the stdout of that command.");
      printtostring(&s,"\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-example",8) && strlen(c) == 8))
    {
      printtostring(&s,"-example <\"command\">\n\n");
      printtostring(&s,"Show an example usage of the specified command.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-header",7) && strlen(c) == 7))
    {
      printtostring(&s,"-header\n\n");
      printtostring(&s,"Use this option to provide a one line header at the start of the output table.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-headeronly",11) && strlen(c) == 11))
    {
      printtostring(&s,"-headeronly\n\n");
      printtostring(&s,"Use this option to output the one line header and then quit without processing any light curves.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-showinputlistformat"))
    {
      printtostring(&s,"-showinputlistformat\n\n");
      printtostring(&s,"Print the expected format of the input light curve list and exit. This command was called \"inputlistformat\" in older versions of VARTOOLS.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-showinputlcformat"))
    {
      printtostring(&s,"-showinputlcformat\n\n");
      printtostring(&s,"Print the expected format of the input light curve(s) and exit.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-tab",4) && strlen(c) == 4))
    {
      printtostring(&s,"-tab\n\n");
      printtostring(&s,"Use this option to use a tab-delimited starbase format for the output table rather than the default space-delimited format.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-basename",9) && strlen(c) == 9))
    {
      printtostring(&s,"-basename\n\n");
      printtostring(&s,"Use this option to print only the basename of each light curve in the output table.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-readformat",11) && strlen(c) == 11))
    {
      printtostring(&s,"-readformat Nskip [\"stringid\" colstringid] [\"inpututc\" format] col_time col_mag col_sig\n\n");
      printtostring(&s,"This option is now deprecated, it is suggested to use the -inputlcformat option instead.\n\nUse this option to specify the format of the input light curves. Nskip is the number of lines to skip at the beginning of each file (this includes any lines that begin with '#' which are otherwise automatically ignored), the default value is 0. If you need to read in a column of strings from each light curve specify \"stringid\" and then the column, by default no column of strings is read in. The stringid column must be specified in this fashion, however, if the matchstringid option is set. If the keyword \"inpututc\" is given, then the input time is taken to be a string giving the UTC of the observation, this will be converted to JD on input. The user must specify the format of the UTC string, this is taken as a format string for a scanf command, with %%Y parsed as the year, %%M as the month, %%D as the day, %%h as the hour, %%m as the minute, and %%s as the second. Note that the year, month, day, hour, and minute will all be converted to integers while seconds are treated as floating point. For example, if the UTC in the light curve has the format 2011-05-18T04:08:03 one would give the format \"%%Y-%%M-%%DT%%h:%%m:%%s\". col_time, col_mag and col_sig are the columns that contain the time, magnitude and magnitude uncertainties (or differential flux and differential flux uncertainty if the light curve will be passed through the -difffluxtomag command), the default values are 1, 2 and 3. Use 0 to not read in a column, if col_time is set to 0 then input times will be set to the line number in the file (starting from 0), if col_mag is set to 0 all magnitudes are set to 0.0, if col_sig is set to 0 all uncertainties are set to 1.0 (use the \"-changeerror\" command to set the uncertainties to the RMS of each light curve). The time is assumed to be in days, though for most commands this is not important.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-functionlist")))
    {
      printtostring(&s,"-functionlist\n\n");
      printtostring(&s,"Show the list of supported functions, operators, constants, and special variables understood by the VARTOOLS analytic expression evaluator.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-inputlcformat")))
    {
      printtostring(&s, "-inputlcformat var1:col1[:type1:[fmt1]][,var2:col2[:type2:[fmt2]],...]\n");
      printtostring(&s, "\t[\"skipnum\" Nskip] [\"skipchar\" skipchar1,skipchar2,...]>\n\n");
      printtostring(&s, "Use this option to specify the format of the input light curves. Provide a comma-separated list of variables to read-in. Here \"var1\", \"var2\" etc are symbolic names for each variable, these can be any alphanumeric strings, the first character should not be a number. The special variable names \"t\", \"mag\", \"err\" and \"id\" are used for the time, magnitude, magnitude uncertainty, and string image identifier (useful for matching observations from different light curves). Following the variable name, \"col1\", \"col2\", etc are the column numbers in the light curve file to read into the associated variables");
#ifdef USECFITSIO
      printtostring(&s, ", or they may be the column header names if the input light curves are in binary fits table format");
#endif
#ifdef _USEBINARY_LC
#ifndef USECFITSIO
      printtostring(&s, ", or they may be the column header names if the input light curves are in");
#else
      printtostring(&s, " or");
#endif
      printtostring(&s, " Penev's binary format");
#endif
#ifdef USECFITSIO
      printtostring(&s, ". If column header names are given, the column names must correspond to the same column numbers for all light curves processed in a single vartools call. The program will only determine the column numbers from the header of the first light curve processed, and will not check to make sure that subsequent light curves use the same columns");
#elif _USEBINARY_LC
      printtostring(&s, ". If column header names are given, the column names must correspond to the same column numbers for all light curves processed in a single vartools call. The program will only determine the column numbers from the header of the first light curve processed, and will not check to make sure that subsequent light curves use the same columns");
#endif
      printtostring(&s, ". If the column number is zero, then the variable will be created, but it will not be read-in from the light curve. Two additional optional control parameters are allowed. This includes the variable type (options are \"double\", \"float\", \"int\", \"long\", \"short\", \"string\", \"char\", or \"utc\"), and a scanf-type format string. The default type is \"double\" except for the variable \"id\" which has a default type of \"string\". Note that for the variables \"t\", \"mag\" and \"err\" the type must be \"double\" or \"utc\", while for \"id\" it must be \"string\". The variable type \"utc\" is a special type which is used to read-in a UTC date and immediately convert it into JD (stored as a double). In this case the user must provide the format of the UTC string, this is taken as a format string for a scanf command, with %%Y parsed as the year, %%M as the month, %%D as the day, %%h as the hour, %%m as the minute, and %%s as the second. Note that the year, month, day, hour, and minute will all be converted to integers while seconds are treated as floating point. For example, if the UTC in the light curve has the format 2011-05-18T04:08:03 one would give the format \"%%Y-%%M-%%DT%%h:%%m:%%s\". If the column number is set to 0, then one must specify the type, while the format will be used to indicate how to initialize the variable. In this case fmt should be an analytic expression which can include any previously defined variables as well as the special variable \"NR\" which is the integer record number of a point in the light curve (0 for the first observation, 1 for the second, etc.). See \"vartools -functionlist\" for a list of supported functions, constants and operators. For the special variables \"t\", \"mag\", and \"err\" if the column number is 0 it is not necessary to specify the type or format; the defaults are t=NR, mag=0, and err=1. The optional keywords \"skipnum\" and \"skipchar\" can be used to specify a number of lines to skip at the start of each light curve, and to change the comment character(s). Lines for which the first non-white-space character matches one of the comment characters will be skipped (the default is \"#\"). Multiple comment characters can be specified using a comma-separated list.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-inlistvars")))
    {
      printtostring(&s, "-inlistvars var1:col1[:vtype1[:fmt1][,var2:col2[:vtyp2[:fmt2]],....]]\n");
      printtostring(&s,"Use this option to read one or more variables from the input light curve list. Provide a comma-separated list of variables to read-in. Here \"var1\", \"var2\" etc are symbolic names for each variable, these can be any alphanumeric strings, the first character should not be a number. The special variable names \"t\", \"mag\", \"err\" and \"id\" are reserved for light curves themselves, and cannot be used here. Following the variable name, \"col1\", \"col2\", etc are the column numbers in the list file to read into the associated variables. If the column number is zero, then the variable will be created, but it will not be read-in from the list. Two additional optional control parameters are allowed. This includes the variable type (options are \"double\", \"float\", \"int\", \"long\", \"short\", \"string\", \"char\", or \"utc\"), and a scanf-type format string. The default type is \"double\". If the column number is set to 0, then one must specify the type, while the format will be used to indicate how to initialize the variable. In this case fmt should be an analytic expression which can include any previously defined variables as well as the special variable \"NF\" which is the integer record number of the line in the list file (0 for the first light curve, 1 for the second, etc.). See \"vartools -functionlist\" for a list of supported functions, constants and operators.\n\n");
      commandfound=1;
    }
  if(all == 1 || (!strncmp(c,"-readall",8) && strlen(c) == 8))
    {
      printtostring(&s,"-readall\n\n");
      printtostring(&s,"Use this option to force the program to read in all the light curves at once. If not specified, this mode will only be used if a command that requires storing all the light curves in memory is issued.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-ascii",6) && strlen(c) == 6))
    {
      printtostring(&s,"-ascii\n\n");
      printtostring(&s,"Deprecated. By default all periodograms are now output in ascii format. This option does nothing, and is only included so to prevent scripts made with older versions of vartools from failing. Use the \"-binaryperiodogram\" option to output periodograms in binary.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-binaryperiodogram")))
    {
      printtostring(&s,"-binaryperiodogram\n\n");
      printtostring(&s,"This option causes periodograms to be output in binary format rather than the default ascii format. The format will be an integer giving the size of the periodogram, followed by double arrays of this size, one for each column included in the ascii version of the periodogram.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-bufferlines")))
    {
      printtostring(&s,"-bufferlines Nlines\n\n");
      printtostring(&s,"If this is specified then Nlines will be buffered within VARTOOLS before being written to standard out. This buffering is separate from the system buffering done on stdout itself, which may be turned off with the -nobuffer option. The -bufferlines option only has an effect if the -parallel option is used, and will help speed up the processing by preventing threads from waiting on each other to output results. The tradeoff with using a larger buffer is that more memory is used, and results will be output less frequently.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-jdtol",6) && strlen(c) == 6))
    {
      printtostring(&s,"-jdtol jdtol\n\n");
      sprintf(s2,"Time measurements within jdtol of each other are considered equal. The default value is %f days.\n\n", DEFAULT_JDTOL);
      printtostring(&s,s2);
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-matchstringid",14) && strlen(c) == 14))
    {
      printtostring(&s,"-matchstringid\n\n");
      printtostring(&s,"If this option is specified then commands that require image-point matching will use a string-id for each image rather than by comparing the jd values. If this option is set then the \"stringid\" column for each light curve must be specified with the -readformat option.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-nobuffer",9) && strlen(c) == 9))
    {
      printtostring(&s,"-nobuffer\n\n");
      printtostring(&s,"If this is specified then stdout will not be buffered, by default it is buffered.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-quiet",6) && strlen(c) == 6))
    {
      printtostring(&s,"-quiet\n\n");
      printtostring(&s,"If this is specified the output table will not be generated, however any calls to output light curves, periodograms etc. will still be executed.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-randseed",9) && strlen(c) == 9))
    {
      printtostring(&s,"-randseed seed\n\n");
      printtostring(&s,"Use this to set the seed for the random number generator. The seed should be an integer, or you may use the word \"time\" to seed the generator with the system clock. If this is not specified the seed value will be 1.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-numbercolumns",14) && strlen(c) == 14))
    {
      printtostring(&s,"-numbercolumns\n\n");
      printtostring(&s,"Prefix each column name in the header with its column number.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-redirectstats",14) && strlen(c) == 14))
    {
      printtostring(&s,"-redirectstats statsfile [\"append\"]\n\n");
      printtostring(&s,"Output the statistics to the file statsfile rather than to stdout. If the \"append\" keyword is given, then the statistics will be appended to that file. This is useful if you wish to output processed light curves to stdout as part of a pipeline.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-oneline",8) && strlen(c) == 8))
    {
      printtostring(&s,"-oneline\n\n");
      printtostring(&s,"Output each statistic on a separate line rather than using the default of outputing a table. This option can provide more readable output when processing a single light curve. It is not suggested when processing a list of light curves.\n\n");
      commandfound = 1;
    }
#ifdef PARALLEL
  if(all == 1 || (!strncmp(c,"-parallel",9) && strlen(c) == 9))
    {
      printtostring(&s,"-parallel Nproc\n\n");
      printtostring(&s,"Process up to Nproc light curves in parallel. Note that if this option is used light curve results will be output in the order that they are finished processing, which is not necessarily the order of the input list.\n\n");
      commandfound = 1;
    }
#endif
#ifdef DYNAMICLIB
  if(all == 1 || (!strcmp(c,"-L")))
    {
      printtostring(&s,"-L libraryfile\n\n");
      printtostring(&s,"Dynamically load a user compiled library defining a new light curve processing command. Here libraryfile is the name of the library. This option allows users to develop their own commands to be incorporated into vartools. See the ReadME file in the USERLIBS directory included with this distribution for examples of how to write, compile, and use your own vartools command libraries. Explicitly loading the library with the \"-L\" option is needed only if the library is not installed in the userlibs install. See also the \"-F\" option which can be used to load libraries defining new functions for use with the analytic expression evaluator.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-F")))
    {
      printtostring(&s,"-F libraryfile\n\n");
      printtostring(&s,"Dynamically load a user compiled library defining new functions which may be used in the evaluation of analytic expressions. Here libraryfile is the name of the library. This option allows users to define their own functions that may be used with the \"-expr\", \"-linfit\", \"-if\", and other commands that handle analytic expressions. See the ReadME file in the USERFUNCS directory included with this distribution for examples of how to write, compile, and use your own analytic function libraries. See also the \"-L\" option which can be used to load libraries defining new vartools processing commands.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-f")))
    {
      printtostring(&s,"-f functodefine\n\n");
      printtostring(&s,"Define an analytic function. Here functodefine should have the format \"funcname(arg1,arg2,....,argN)=function_expression\". Once defined, the function may then be used in subsequent analytic expressions given on the command line.\n\n");
      commandfound = 1;
    }
#endif
  if(all == 1 || !strcmp(c,"-log-command-line"))
    {
      printtostring(&s,"-log-command-line\n\n");
      printtostring(&s,"Print the command-line syntax given to vartools as a comment to the output ascii table, before giving the header.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-skipmissing"))
    {
      printtostring(&s,"-skipmissing\n\n");
      printtostring(&s,"Do not abort if a missing or unreadable light curve file is encountered. Instead skip the light curve and proceed with others in the list.\n\n");
      commandfound=1;
    }
  if(all == 1 || !strcmp(c,"-noskipempty"))
    {
      printtostring(&s,"-noskipempty\n\n");
      printtostring(&s,"By default empty light curves are skipped and not included in the output table. To not skip these, and include them in the output, give this option. Note that this option does not have an effect if the -readall option is used.\n\n");
      commandfound=1;
    }
#ifdef VARTOOLS_VERSION
  if(all == 1 || !strcmp(c,"-version"))
    {
      printtostring(&s,"-version\n\n");
      printtostring(&s,"Output the VARTOOLS version number at the start of the output ascii table.\n\n");
      commandfound = 1;
    }
#endif
  if(all == 1)
    {
      printtostring(&s," Commands:\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-addnoise",9) && strlen(c) == 9))
    {
      listcommands_noexit("-addnoise",p,&s);
      printtostring(&s,"Add time-correlated noise to the light-curve. The noise is assumed to be Gaussian. The user must specify the model to be assumed for time correlations. Options include:\n\n");
      printtostring(&s,"\"white\" - pure white noise without a time-correlated component. This noise is parameterized by the standard deviation, sig_white.\n\n");
      printtostring(&s,"\"squareexp\" - a Gaussian-process with a squared-exponential covariance matrix. The covariance between times t_i and t_j is given by sig_white*sig_white*delta_ij + sig_red*sig_red*exp(-(t_i-t_j)^2/2/rho^2) where delta_ij is the Kronecker delta function, and sig_white, sig_red and rho are parameters of the model. Note that rho and sig_red must both be greater than zero. The sig_white term allows for an additional noise component that is uncorrelated in time. The optional \"bintime\" keyword can be used to chunk the light curve into timebins of a specified duration and independently simulate the correlated noise over each bin. This can speed up the simulation substantially in cases where the full light curve duration is much longer than the correlation timescale.\n\n");
      printtostring(&s,"\"exp\" - a Gaussian-process with an exponentially decaying covariance matrix. The covariance between times t_i and t_j is given by sig_white*sig_white*delta_ij + sig_red*sig_red*exp(-(t_i-t_j)/rho) where delta_ij is the Kronecker delta function, and sig_white, sig_red and rho are parameters of the model. Note that rho and sig_red must both be greater than zero. The sig_white term allows for an additional noise component that is uncorrelated in time. The optional \"bintime\" keyword can be used to chunk the light curve into timebins of a specified duration and independently simulate the correlated noise over each bin. This can speed up the simulation substantially in cases where the full light curve duration is much longer than the correlation timescale.\n\n");
      printtostring(&s,"\"matern\" - a Gaussian-process with a Matern covariance matrix. The covariance between times t_i and t_j is given by sig_white*sig_white*delta_ij + sig_red*sig_red*C(nu,x)*K_nu(x) with x=sqrt(2*nu)*(|t_i-t_j|)/rho and C(x,y)=(2^(1-x)/Gamma(x))*(y)^x with K_nu being the modified Bessel function of the second kind, and Gamma the usual gamma function. In this case the parameters are nu, rho, sig_white, and sig_red.  All but sig_white must be greater than zero. The sig_white term allows for an additional noise component that is uncorrelated in time. The optional \"bintime\" keyword can be used to chunk the light curve into timebins of a specified duration and independently simulate the correlated noise over each bin. This can speed up the simulation substantially in cases where the full light curve duration is much longer than the correlation timescale.\n\n");
#ifdef _HAVE_GSL
      printtostring(&s,"\"wavelet\" - the sum of a red-noise component with a power-spectral-density proportional to 1/f^gamma (gamma must satisfy -1~<~gamma~<~1) with standard deviation sig_red, and a white noise component with standard deviation sig_white. The red-noise is generated using the method of \"McCoy and Walden, 1996, Journal of Computational and Graphical Statistics, Vol 5, No. 1, pp. 26-56 (section 4.2).\"\n\n");
#endif
      printtostring(&s,"For each of the parameters nu, rho,");
#ifdef _HAVE_GSL
      printtostring(&s," gamma,");
#endif
      printtostring(&s," sig_red, and sig_white, one may either use the \"fix\" keyword and specify the value on the command-line, or use the \"list\" keyword in which case the value is read-in from the input list (optionally give the column in the list with the \"column\" keyword, if not specified the next free column is assumed).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-alarm",6) && strlen(c) == 6))
    {
      listcommands_noexit("-alarm",p,&s);
      printtostring(&s,"Calculate the alarm variability statistic for each light curve. Cite Tamuz, Mazeh, and North 2006, MNRAS, 367, 1521 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-aov",4) && strlen(c) == 4))
    {
      listcommands_noexit("-aov",p,&s);
      printtostring(&s,"Perform an AoV period search on the light curves using phase binning. Specify \"Nbin\" and a number to change the number of bins used from the default value of 8. minp and maxp are the minimum and maximum periods to search. The intial search will use a frequency resolution of subsample/T where T is the baseline of the light cruve. The peak periods will be refined using a resolution of finetune/T. The program will output the Npeaks highest peaks in the periodogram of the light curve. By default the program will output, for each peak, the period, the theta_aov statistic, the signal to noise ratio (theta_aov~-~<theta_aov>)/RMS(theta_aov), and the negative natural logarithm of the formal false alarm probability (this is calculated from the value of theta_aov using the Horne and Baliunas, 1986, ApJ, 302, 757 estimate for the bandwidth penality). operiodogram should be either 0 or 1. If it is set to 1 then the periodogram for each light curve will also be output to the directory outdir, with the suffix \".aov\". The first column in the output is the period, the second column is theta_aov. If the \"whiten\" keyword is given, then the light curve will be whitened at each peak period and the periodogram will be recomputed before searching for the next peak period. The average and RMS theta_aov used for each peak are computed on the whitened periodogram. The output spectrum will contain the Npeaks computed periodograms. Use the \"clip\" keyword to change the clipping parameters for calculating the average and RMS of the power spectrum when computing the SNR value of a peak. clip is the sigma-clipping factor, and clipiter is a flag that is 1 or 0 to toggle iterative clipping. By default iterative 5-sigma clipping is used. Use the keyword \"uselog\" to output the log(theta_aov) SNR instead of the standard statistics (this is the original behavior for this command), i.e. the statistic is (<-ln(theta_aov)>~-~ln(theta_aov))/RMS(-ln(theta_aov)) where the average and RMS are taken over the entire periodogram for a light curve, in this case the output will also include the average and RMS of -ln(theta_aov). Use the \"fixperiodSNR\" option to output the AoV statistic, SNR etc at a specified period. See the help for the \"-LS\" command for an explanation of the syntax. Cite Schwarzenberg-Czerny, A., 1989, MNRAS, 241, 153 and Devor, J., 2005, ApJ, 628, 411 if you use this tool. (The Devor citation is needed because this code is based on his implementation of AoV).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-aov_harm",9) && strlen(c) == 9))
    {
      listcommands_noexit("-aov_harm",p,&s);
      printtostring(&s,"Perform an AoV period search on the light curves using a multi-harmonic model. Nharm is the number of harmonics to use. minp and maxp are the minimum and maximum periods to search. The initial search will use a frequency resolution of subsample/T where T is the baseline of the light cruve. The peak periods will be refined using a resolution of finetune/T. The program will output the Npeaks highest peaks in the periodogram of the light curve. The output for each peak includes the period, the theta_AoV statistic, the signal to noise ratio (theta_aov~-~<theta_aov>)/RMS(theta_aov), the false alarm probability (see the help for the \"-aov command\"), and the average and RMS of theta_aov. operiodogram should be either 0 or 1. If it is set to 1 then the periodogram for each light curve will also be output to the directory outdir, with the suffix \".aov_harm\". The first column in the output is the period, the second column is theta_aov. If the \"whiten\" keyword is given, then the light curve will be whitened at each peak period and the periodogram will be recomputed before searching for the next peak period. The average and RMS theta_aov output for each peak are computed on the whitened periodogram. The output spectrum will contain the Npeaks computed periodograms. Use the \"clip\" keyword to change the clipping parameters for calculating the average and RMS of the power spectrum when computing the SNR value of a peak. clip is the sigma-clipping factor, and clipiter is a flag that is 1 or 0 to toggle iterative clipping. By default iterative 5-sigma clipping is used. Use the \"fixperiodSNR\" option to output the AoV statistic, SNR etc at a specified period. See the help for the \"-LS\" command for an explanation of the syntax. Cite Schwarzenberg-Czerny, A., 1996, ApJ, 460, L107 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-autocorrelation",16) && strlen(c) == 16))
    {
      listcommands_noexit("-autocorrelation",p,&s);
      printtostring(&s,"Calculate the discrete auto-correlation function (Edelson and Krolik 1988, ApJ, 333, 646), these are written out to outdir/basename.autocorr, start, stop and step are the times in days for sampling the auto-correlation. Note that rather than using the variance of the light curve in the denominator, we use the formal uncertainty (we do not subtract the \"measurement error\" from the variance as done in the Edelson and Krolik formula since this could lead to imaginary numbers in the case where \"measurement errors\" are overestimated). If you wish to use the variance in the denominator rather than the formal uncertainty you should issue the -changeerror command before calling this routine. Note that due to binning, when the variance is used in the denominator the autocorrelation function may be smaller than 1 unless the time step used is less than the time difference between any consecutive measurements in the light curve.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-binlc",6) && strlen(c) == 6))
    {
      listcommands_noexit("-binlc",p,&s);
      printtostring(&s,"Bin the light curves in time (or phase if a -Phase command has already been issued). Use the \"average\" keyword to take the average of points in a bin, \"median\" to take the median of the points or \"weightedaverage\" to take the weighted average (for backward compatibility, you can also use the numbers 0, 1 or 2). One should specify either the binsize in units of the time coordinate (e.g. in phase if the light curves have been phased), or the number of bins to split the light curve time-span into. By default all other columns in the light curve will be binned in the same way, you can optionally change the binning type for some of the columns using the \"bincolumns\" keyword. If you do this, then the following string of the form var1[:stats][,var2[:stats2],...] is a comma separated list of variable names, with a statistics name appended to each variable following a \":\". The choices for statistics are the same as for the \"-stats\" command. By default the first bin begins at the initial time in the light curve (t0), if firstbinshift is specified, then this will be shifted by t0 - firstbinshift/binsize. If the \"tcenter\" keyword is given then the output time for each bin is the time at the center of the bin, if \"taverage\" is given then the time will be the average of the times of points that fall within the bin, if \"tmedian\" is given then the time will be the median of the times of points in the bin, and if \"tnoshrink\" is given then the size of the light curve will not be reduced by the binning, and instead all points in the light curve will be replaced by the binned value (for backward compatability, you can also use the numbers 0, 1, 2, or 3). If you only want to apply the binning to some of the columns in the light curve, but not all of them, you can give the \"tnoshrink\" keyword followed by \"bincolumnsonly\". In this case only the variables explicitely listed after \"bincolumns\" will be binned (this includes t, mag and err).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-BLS",4) && strlen(c) == 4))
    {
      listcommands_noexit("-BLS",p,&s);
      printtostring(&s,"This command runs BLS on the light curves. You may set a fixed minimum and maximum q (fraction of orbit in transit) for the search, you can specify a minimum and maximum stellar radius to consider (in solar radii) in which case the qmin and qmax values are calculated for each trial period P as q~=~0.076~*~R**(2/3)~/~P**(2/3) (with P in days - this assumes R~=~M as an approximation for the lower main sequence), or you can specify a stellar density (in grams per cubic centimeter) and a minimum and maximum fraction of the expected transit duration (assuming a circular orbit) to consider. The minimum and maximum periods to search are given in days, nfreq is the number of trial frequencies (10000 is a decent number). nbins is the number of phase bins to break the light curve into (200 is a decent number). timezone is the number to add to the julian date to get the local date (-7 for arizona), this is used to determine which nights different observations come from and thus to determine the fraction of delta-chi2 that comes from a single night. Npeak is the number of peaks in the BLS spectrum to find and report. outperiodogram is a flag that is set to 1 to output the BLS period vs. SN spectrum. outdir is the output directory for the BLS spectrum if outperiodogram is set to 1, \".bls\" will be appended to the filename. omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".bls.model\" will be appended to the filename. correctlc is either 0 or 1, set to 1 it will subtract the transit model from the light curve before passing it to the next command. If the optional \"fittrap\" keyword is specified, the routine will fit a trapezoidal transit to each BLS peak. This allows a refined estimate of the transit time, duration and depth. In this case the output table will also include \"qingress\" which is the fraction of the transit duration covered by ingress, a value of 0 corresponds to a perfectly box-shaped transit, a value of 0.5 corresponds to a V-shaped transit. The optional keyword \"nobinnedrms\" adjusts the way in which the BLS_SN statistic is calculated. When \"nobinnedrms\" is given the procedure runs faster, but the SN will tend to be suppressed for high significance detections. If this option is given use another parameter (delta chi2 or signal to pink noise) for selecting transits rather than BLS_SN. Note that in some cases the default behavior will supress the signal at very high peaks (characterized by the periodogram going to zero in the center of a peak), to avoid this behavior use the \"nobinnedrms\" keyword. If the optional keyword \"ophcurve\" is given then a model phase curve will be output to a file in the directory outdir with suffix \".bls.phcurve\" with phases between phmin and phmax and a uniform step size of phstep. If the optional keyword \"ojdcurve\" is given then a model light curve will be output to a file in the directory outdir with suffix \".bls.jdcurve\" with times between the first and last times in the light curve with a uniform step size of jdstep. By default the BLS spectrum is sampled at uniform frequency steps. To sample it at uniform steps in period or log(period) using the \"stepP\" or \"steplogP\" keywords. Use the optional \"adjust-qmin-by-mindt\" keyword to adaptively set the minimum q value to the maximum of qmin or mindt*frequency where mindt is the minimum time difference between consecutive points in the light curve. If this keyword is given you may also use the \"reduce-nbins\" keyword to adaptively reduce the number of phase bins at each frequency such that there are no more than 2 bins to sample a transit of duration qmin. Cite Kovacs, G., Zucker, S., Mazeh, T. 2002, A&A, 391, 369 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-BLSFixPer",10) && strlen(c) == 10))
    {
      listcommands_noexit("-BLSFixPer",p,&s);
      printtostring(&s,"This command runs BLS at a fixed period on the light curves - it searches merely for the most transit-like signal at the specified period. The period comes either from the last aov command, from the last ls command, is specified in the input-list (by default it will be taken from the next available column, use the \"column\" keyword to specify the column), is fixed to the value given on the command-line for all light curves, it can be set to the value in an output column from a previously executed command, or it can be set to the value of an analytic expression. Like BLS you can either specify a minimum and maximum stellar radius to consider, or a minimum and maximum q to consider. Similarly you must specify the number of bins to use and the timezone.  omodel is 1 to output the model light curve to model_outdir (with \".blsfixper.model\" appended to the end of the file name. correctlc is 1 to subtract the transit model from the light curve before passing it to the next command. \"fittrap\" is an optional keyword, which if specified, fits a trapezoidal transit after running the BLS scan. Cite Kovacs, G., Zucker, S., Mazeh, T. 2002, A&A, 391, 369 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-changeerror",12) && strlen(c) == 12))  
    {
      listcommands_noexit("-changeerror",p,&s);
      printtostring(&s,"Replace the formal errors in a light curve with the RMS of the light curve.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-changevariable")))
    {
      listcommands_noexit("-changevariable",p,&s);
      printtostring(&s,"Change the variable used for time, magnitude, magnitude uncertainty, or image-id in subsequent commands. First specify what is being changed, and then give the name of the variable which will now fulfill this role. Suppose for example that you have read in the variables mag and mag2 using the -inputlcformat option. After \"-changevariable mag mag2\" is issued, subsequent commands will use mag2 for the magnitudes rather than mag. The variable mag will still exist, if you wish to switch back to using mag for the magnitudes, you should issue the command \"-changevariable mag mag\".\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-chi2",5) && strlen(c) == 5))
    {
      listcommands_noexit("-chi2",p,&s);
      printtostring(&s,"Calculate chi2 per dof for the light curves. The output will include chi2 and the error weighted mean magnitude.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-chi2bin",8) && strlen(c) == 8))
    {
      listcommands_noexit("-chi2bin",p,&s);
      printtostring(&s,"Calculate chi2 per dof after applying a moving mean filter (each point in the light curve is replaced by the average of all points within a range centered on the point) to the light curves. (Note that the light curves passed to the next command are unchanged). Nbin filters are used (with Nbin resulting estimates of chi2 and the error weighted mean). The width of each filter is given by 2.0*bintime (bintime is given in minutes)\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-clip",5) && strlen(c) == 5))
    {
      listcommands_noexit("-clip",p,&s);
      printtostring(&s,"Sigma-clip the light curves using a clipping factor of sigclip. iter is a flag that is either 1 for iterative clipping (performed continuously until no further points are removed), or 0 to not do continuous iterative clipping. To clip for a fixed number of iterations give the \"niter\" keyword followed by the number of iterations to use. By default clipping is done with respect to the mean. To use the median instead provide the \"median\" keyword. The output table will include the number of points that were clipped. Note that points with errors <= 0 or NaN magnitude values will be clipped. If sigclip is <= 0, then sigma clipping is not performed, but points with errors <= 0 or NaN magnitude values will be clipped.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-converttime",12) && strlen(c) == 12))
    {
      listcommands_noexit("-converttime",p,&s);
      printtostring(&s,"Convert the time system of the light curves. The user must specify \"input\" followed by the input time system (options are \"mjd\" = modified Julian date, \"jd\" = Julian date, \"hjd\" = Helio-centric Julian date, or \"bjd\" = Bary-centric Julian date). The \"bjd\" option is only available if vartools has been linked to the JPL NAIF cspice library (available from http://naif.jpl.nasa.gov/naif/toolkit_C.html). This program uses the convention MJD~=~JD~-~2400000.5 (Note the 0.5!). Conversions between JD and HJD assume that the observations are made at the Earth-Moon Barycenter, and that this position follows an elliptical orbit about the center of the solar system with linear perturbations to the orbital elements. Differences between BJD and HJD are due mostly to the orbits of Jupiter and Saturn (which are ignored for HJD), and can be as large as 4.2 seconds for observations made between 1980 and 2020. If the input time value has a constant subtracted from it, use the \"inputsubtract\" keyword to specify the constant. For example, if the input time is HJD-2400000, one would give \"input hjd inputsubtract 2400000\". One can optionally specify whether the input JD values have been calculated directly from UTC without accounting for leap-seconds (\"inputsys-utc\", which is the default, and is typically the case if one is using the JD values given in the image headers from most ground-based observatories), or if the time was converted to TDB (barycentric dynamical time, i.e. corrections have been made for leap-seconds) before converting to JD (\"inputsys-tdb\"). In the year 2011 there is approximately a one minute difference between the two systems. The cspice library is necessary to handle TDB. The user must then specify \"output\" followed by the output time system to use. If a constant value should be subtracted from the output times use the \"outputsubtract\" keyword followed by the constant. One may also specify whether the output system should be UTC or TDB based (by default the input system is assumed). If converting to/from HJD or BJD the user must also indicate the source for the RA/DEC coordinates of the target. This can either be \"list\" (in which case the values are taken from the input-list, either from the next available column, or from the user-specified column) or \"fix\" followed by the ra and dec values (in degrees) given on the command line. The epoch may optionally be specified as well, the default is 2000.0. The user may also specify a proper-motion for the target by giving the \"ppm\" keyword and the source for those values. The values are assumed to be in units of mas per year. If for some reason a different radec or ppm was used in determining the input time values (e.g. if the input HJD was calculated for a given image assuming the coordinates of the center of the field, but one wishes to output BJD for the coordinates of a specific star in the field) one may use the \"input-radec\" and \"input-ppm\" keywords. Conversions to/from HJD or BJD assume that the source is infinitely far away and that the observations are made from the surface of the Earth. We do not account for the Shapiro or Einstein time-delays in converting to BJD. See \"Eastman, Siverd, and Gaudi, 2010, PASP, 122, 935\" for an illuminating discussion of time-conversions. Conversions to/from TDB or BJD require JPL NAIF \"kernel\" files which store ancillary time/solar-system data. These include an ephemeris file to determine the position of the Earth with respect to the Solar System Barycenter (e.g. ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de421.bsp), a leap-second file (e.g. ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0009.tls or a later one if available), and a planetary physical data file used to determine the location of the observer with respect to the center of the Earth in the J2000.0 inertial frame (e.g. ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00009.tpc). The files can either be specified on the command-line (with the \"ephemfile\" etc. keywords) or one can set environment variables giving the path to these files (the variables are CSPICE_EPHEM_FILE, CSPICE_LEAPSEC_FILE and CSPICE_PLANETDATA_FILE respectively). Finally you may specify the observatory or the latitude/longitude/altitude where the observations were made (if not specified, then the center of the Earth is assumed; this can introduce up to a 21 millisecond error in the BJD correction). To specify an observatory provide the \"observatory\" keyword followed by the code for the observatory. To see a list of the codes use \"observatory show-codes\". To give latitude/longitude/altitude values use the \"coords\" keyword, and then either use \"fix\" to specify the values on the command-line, \"list\" to read the values from the input-list, or \"fromlc\" to read the values from the indicated columns in the light curves (necessary if one is processing light curves containing data from multiple sites). This command has an internal precision for conversions near J2000.0 of approximately 0.1 milliseconds.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-copylc")))
    {
      listcommands_noexit("-copylc",p,&s);
      printtostring(&s,"Make Ncopies copies of the current light curve. These will each be processed by VARTOOLS starting with the command following this -copylc command. Data in the output ascii table for commands preceeding the -copylc command will also be replicated for the light curve copies. Any files output by preceeding VARTOOLS commands will not be replicated. The suffix \"_copy$copycommandnum.$copynum\" will be appended to each name, where $copycommandnum is the -copylc command which created the copy (in case there are multiple -copylc commands given) and $copynum indicates which copy this is (from 0 to Ncopies-1). The -copylc command cannot be used with the -readall option.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-decorr",7) && strlen(c) == 7))
    {
      listcommands_noexit("-decorr",p,&s);
      printtostring(&s,"Decorrelate the light curves against specified signals. The signals are either global signals (e.g. the airmass) that are input in separate files (with the format: JD signal_value or the format: stringid signal_value if the -matchstringid option was set), or are light curve specific signals (e.g. the sub-pixel position) that are input as additional columns in the light curves. correctlc is a flag that is 0 or 1, if the value is 1 then the light curves passed to the next command will be decorrelated, if the value is 0 then resulting chi2 from decorrelation and the decorrelation coefficients will be output to the table, but the light curves themselves will not be decorrelated. zeropointterm is a flag that is 0 or 1, if the value is 1 then a zeropoint term is included in the decorrelation, if it is 0 then no such term is included. subtractfirstterm is a flag that is 0 or 1, if the value is 1 then the light curves are decorrelated against the signal minus the first signal value (this is useful if one is decorrelating against the JD, for example, in which case one would decorrelate against JD - JD0 where JD0 is the first JD in the light curve), if it is 0 then the light curves are decorrelated against the signal without subtracting the first term. Nglobalterms is the number of global files. globalfile1...globalfileN or the names of those files, order1...orderN are the orders of the polynomials used to decorrelate each signal (each must be greater than 0). Nlcterms is the number of light curve specific signals. The columns of these signals are given by lccolumn1...lccolumnN. The orders of the polynomials are given by lcorder1...lcorderN. omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".decorr.model\" will be appended to the filename.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-dftclean",9) && strlen(c) == 9))
    {
      listcommands_noexit("-dftclean",p,&s);
      printtostring(&s,"Compute the DFT power spectrum of the light curve and optionally apply the clean algorithm (Roberts et al. 1987) to it. nbeam is the number of points per 1/T frequency element to include in the power spectrum. Use \"maxfreq\" to set the maximum frequency to maxf. The default value for the maximum frequency is 1/(2*min_time_separation). Use \"outdspec\" to output the dirty power spectrum, the suffix \".dftclean.dspec\" will be appended to the name of the output file. Use \"finddirtypeaks\" to find the top Npeaks peaks in the dirty power spectrum. By default 5sigma iterative clipping is used in calculating the SNR of a peak. To change this use \"clip\" followed by the new sigma-clipping value and 1 to do iterative clipping or 0 to not do iterative clipping. Use \"outwfunc\" to output the window function, with the suffix \".dftclean.wfunc\". Use \"clean\" to apply the clean algorithm to the power spectrum. Please cite Roberts, D.H., Lehar, J., and Dreher, J.W. 1987, AJ, 93, 4 if you use this algorithm. To use clean you must specify the gain which is a value between 0.1 and 1 (with lower values giving slower convergence). The procedure will continue to clean the spectrum until the last peak is less than SNlimit times the noise. Use \"outcbeam\" to output the clean beam which is convolved with the deconvolved spectrum to produce the final clean power spectrum. The suffix for the clean beam file is \".dftclean.cbeam\". Use \"outcspec\" to output the clean power spectrum, the suffix is \".dftclean.cspec\". Finally use \"findcleanpeaks\" to find the top Npeaks peaks in the clean power spectrum. By default the SNR values output by one of the find commands will be determined on the power spectrum, and will be defined as (peak~-~ave)/std. To use the amplitude spectrum rather than the power spectrum give the option \"useampspec\". To output the average and standard deviation of the spectrum before and after clipping, in addition to the final SNR value, use the option \"verboseout\". This routine uses the FDFT algorithm by Kurtz 1985, MNRAS, 213, 773 to compute the DFT.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-difffluxtomag",14) && strlen(c) == 14))
    {
      listcommands_noexit("-difffluxtomag",p,&s);
      printtostring(&s,"Convert light curves from ISIS differential flux into magnitudes. In this case a light curve list must be used, and the reference magnitudes of the stars (after aperture correction) should be given as an additional column in the light curve list. Use the \"magcolumn\" keyword to specify the column, by default it will be the next available column. mag_constant is the magnitude of a source with a flux of 1 ADU, offset is an additive constant to apply to the output light curves. This command does not yield any output to stdout, you must explicitey call -rms or -chi2 if you want statistics.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-ensemblerescalesig",19) && strlen(c) == 19))
    {
      listcommands_noexit("-ensemblerescalesig",p,&s);
      printtostring(&s,"Transform the magnitude uncertainties by e_mag~-->~sqrt(a*e_mag*e_mag~+~b). The parameters 'a' and 'b' are determined by fitting a linear relation between (expected~rms)^2 and (chi2/dof)*(expected~rms)^2 for all light curves. The result of this transformation is that chi2 per dof is distributed about unity. Outliers in the distribution are clipped during the fitting using a clipping factor of sigclip. This routine requires a list of light curves. If this command is invoked, then all light curves will be read into memory. The output will include the average rescale factor for each light curve (taken to be sqrt(chi2_after/chi2_before). See also the -rescalesig command which does a sigma-rescaling on a light curve-by-light curve basis.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-expr")))
    {
      listcommands_noexit("-expr",p,&s);
      printtostring(&s,"Evaluate an analytic expression. The argument to this command is an equality, with the left-hand-side being the name of a variable to update, and the right-hand-side being an analytic expression. See \"vartools -functionlist\" for the allowed syntax. If the variable on the left-hand-side has not previously been defined, it will be created as a light-curve vector. Note that variables which appear on the right-hand-side can be the name of a vector which is read-in from the light curve and set using the -inputlcformat option, a scalar or vector created by another command (e.g. the fitting parameters, or the output model vector, created by the -linfit command), or any output parameter from a previously executed command (these latter variables have the names of the output columns which can be found by running vartools with the -headeronly command. Note that any leading numbers or '_' characters are removed; also any character that is not a digit, a letter, or a '_' in the header name will be replaced with a '_' in its variable name equivalent. For example if the header name is '2_STATS_mag_PCT20.00_0' its variable name would be 'STATS_mag_PCT20_00_0', where the leading '2_' is removed, and the '.' is replaced with a '_'). The left-hand-side of the equality cannot be a variable associated with an output column.\n\nUse square brackets '[]' after the variable name on the left-hand-side to update only a subset of the components of the vector. Within the square brackets you can either specify a single index value (indexed starting from 0), an expression that evaluates to a single index value, a range of indices (using a ':' to distinguish between the first index and the last index to operate on) or an expression that evaluates to a vector, in which case only components where the resulting vector has a value > 0 will be updated. Examples of each of these are as follows:\n\n1: '-expr mag[0]=5.0' would set the magnitude value for the first observation in the light curve to 5.0, and would not affect any of the other components in the light curve.\n\n2: '-expr mag[i+2]=5.0' would set the magnitude value for the observation indexed by 'i+2' (if 'i' is a scalar variable, if it is a light curve vector it would be treated as in example 8 below) to 5.0.\n\n3: '-expr mag[len(mag)-2]=5.0' would set the magnitude value for the second to last observation to 5.0. Note that len(mag) returns the length of the magnitude vector, and due to indexing from zero, len(mag)-1 is the index of the last point in the vector.\n\n4: '-expr mag[0:2]=5.0' would set the magnitude values for the first and second observations in the light curve to 5 (note that here we follow the python convention where the upper index in the range is not itself included, so mag[0:2] corresponds to mag[0] and mag[1]).\n\n5: '-expr mag[:5]=5.0' would set the magnitude values for the first five observations in the light curve to 5.\n\n6: '-expr mag[5:]=5.0' would set the magnitude values for all but the first five observations to 5.0.\n\n7: '-expr mag[(t>25.0)&&(t<30.0)]=5.0' would set the magnitude values for any points with 25 < t < 30.0 to 5.0.\n\n8: '-expr mag[t-25.0]=5.0' would set the magnitude values for any points with t > 25 to 5.0.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-findblends",11) && strlen(c) == 11))
    {
      listcommands_noexit("-findblends",p,&s);
      printtostring(&s,"Use this routine to find variability blends. The routine operates by matching a list of potential variables to a list of light curves. For each potential variable, the flux amplitude of all matching light curves will be measured and the one with the highest amplitude will be chosen. The name and flux amplitude of the highest amplitude light curve are included in the output statistics table. If this routine is used, light curves must be read from an input list with x and y coordinates given as additional columns in the list file (i.e. you must use \"-l\" for input rather than \"-i\"), this will be used as the list of potential variables. By default the x and y values will be taken from the next available columns in the input list, use the \"xycol\" keyword to specify the columns. matchrad is the matching radius used to match stars that are possibly blended. If \"radec\" is specified, then matchrad is given in arcseconds and the x and y coordinates are ra/dec in degrees, if \"radec\" is not specified then rectangular matching will be used. Use the \"fix\", \"list\" or \"fixcolumn\" keywords to specify the source for the periods to use for each potential variable. If \"fix\" is given then the periods for all potential variables are set to the listed value. If \"list\" is given then the periods should be included as an additional column in the input list file (by default it will be taken from the next available column, use the \"column\" keyword to specify the column). If \"fixcolumn\" is given then the period for each potential variable in the list will be set to an arbitrary previously computed statistic by giving the name or number of the output column as seen with the -header or -numbercolumns options. By default the list of potential variables will be matched to itself. If, however, \"starlist\" is specified, then the list of potential variables will be matched to the list contained in the file starlistfile. The first three columns in this file are the light curve name, and the x and y coordinates. Use \"zeromag\" to specify the zero-point magnitude for converting magnitudes into fluxes, the default value is 25.0. If you do not wish to convert from magnitudes into fluxes (e.g. if the input light curves are already in flux units), use \"nofluxconvert\", make sure you know what you're doing then when you use this option. The amplitude of the light curves are taken to be the peak-to-peak amplitude of a fourier series fit to the light curve with period P. If the period is less than or equal to zero then the amplitude will be measured by fitting the fourier series to the unphased light curve (i.e. the period will be set equal to the total duration of the light curve). Use \"Nharm\" to specify the number of harmonics to include in the Fourier series (if it is 0 then only a sinusoid will be fit, the default value is 2). You can output the names and flux amplitudes of all stars matching to each potential variable by using \"omatches\" and giving the name of the file to output this information to.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-fluxtomag",10) && strlen(c) == 10))
    {
      listcommands_noexit("-fluxtomag",p,&s);
      printtostring(&s,"Convert light curves from flux into magnitudes. mag_constant is the magnitude of a source with a flux of 1 ADU, offset is an additive constant to apply to the output light curves. This command does not yield any output to stdout.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-GetLSAmpThresh",15) && strlen(c) == 15))
    {
      listcommands_noexit("-GetLSAmpThresh",p,&s);
      printtostring(&s,"This routine determines the minimum amplitude that a light curve could have and still be detected at a given period (which is either specified in the input list, or taken from the last LS command) with a LS~log(formal~false~alarm~probability)~<~thresh. The light curve signal is calculated by fitting a fourier series at the period of the light curve with Nharm harmonics and Nsubharm subharmonics if \"harm\" is specified, or it is read in from a file if \"file\" is specified. If reading the signal from a file, listfile is a file with two columns of the form: signal_file signal_amp, with one line for each light curve being operated on. Each signal file should contain the signal magnitude in the third column, signal_amp is the amplitude in magnitudes of the signal (this is to allow for cases where the signal amplitude is greater than just the max observed signal magnitude minus the min observed signal magnitude, e.g. if the signal comes from fitting a fourier series to the light curve at a different period from the one used to do the ls selection). To calculate the minimum amplitude, the signal is subtracted it from the light curve and re-added after scaling it by a number. minP is the minimum period that would be searched (this is needed since it sets the false alarm probability scale). The output is the minimum factor by which the signal could be scaled and still be detectable together with the minimum peak-to-peak amplitude. Note that by default the Generalized Lomb-Scargle periodogram is assumed to be used for this command. If you wish the threshhold to be computed for the traditional Lomb-Scargle periodogram, then you should give the \"noGLS\" keyword.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-if")) || (!strcmp(c,"-elif")) || (!strcmp(c,"-else")) || (!strcmp(c,"-fi")))
    {
      listcommands_noexit("-if",p,&s);
      printtostring(&s,"Make execution of commands conditional upon the evaluation of an expression. If the expression evaluates as 0 the commands following \"-if <expression>\" and preceding the next \"-elif\", \"-else\" or \"-fi\" statement will not be executed, if it evaluates as a number different from 0 when cast to an integer the commands will be executed. The \"-elif <expression>\" and \"-else\" constructs may also be used to provide a set of commands to be executed in case the expressions associated with previous \"-if\" and \"-elif\" statements have not yet evaluated as true. The construct may be terminated with the \"-fi\" statement. Nested \"-if\", \"-elif\", \"-else\", \"-fi\" constructs are allowed. Any conditional constructs that are not explicitely terminated with \"-fi\" are assumed to be terminated after the last command given on the command-line. CAUTION: conditional constructs are ignored by commands which process all light curves simultaneously (e.g. -SYSREM or -findblends) as well as by the -savelc and -restorelc commands.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-Injectharm",11) && strlen(c) == 11))
    {
      listcommands_noexit("-Injectharm",p,&s);
      printtostring(&s,"Inject a harmonic series into the light curves. The harmonic series has the form:\n\n\tA_1*cos(2*pi*(t/P~+~phi_1))\n\t\t+~sum_{k=2}^{Nharm+1}~A_k*cos(2*pi*(t*k/P~+~phi_k))\n\t\t+~sum_{k=2}^{Nsubharm+1}~Asub_k*cos(2*pi*(t/k/P~+~phisub_k))\n\nThe period for the series is either specified for each light curve in the list file (use \"list\", by default it will be taken from the next available column, though the user may optionally specify the column with the \"column\" keyword), it is fixed to a particular value for all light curves (use \"fix\"), it is generated randomly either from a uniform distribution (use \"rand\") or from a uniform log distribution (use \"logrand\"), or the frequency is generated from a uniform distribution (use \"randfreq\") or from a uniform log distribution (use \"lograndfreq\"). Specify the number of harmonics, and then for each harmonic specify how the amplitude and phase of that harmonic is to be generated. Note that the fundamental period in this case corresponds to harmonic 1, so for Nharm=0 there must be one set of amplitude/phase commands, for Nharm=1 there are two sets, etc. amplist, ampfix, amprand and amplogrand specify whether the amplitude for the harmonic is specified as a column in the light curve list file, or if it is fixed to a particular value for all light curves, or if it is randomly generated from a uniform or uniform log distribution. You may optionally specify \"amprel\" after the amplitude command to make the specified amplitude of this harmonic to be relative to the amplitude of the fundamental mode (i.e. the input is R_i1~=~A_i/A_1 rather than A_i). This would be used, for example, if one wishes to simulate a specific signal form (say a Cepheid) but with a random over-all amplitude. phaselist, phasefix and phaserand specify whether the phase of the harmonic is specified as a column in the light curve list file, or if it is fixed to a particular value for all light curves, or if it randomly generated from a uniform distribution between 0 and 1. The phase is the phase at time t=0. If the optional \"phaserel\" term is given then the phase for the harmonic will be relative to the phase of the fundamental (i.e. the input is phi_k1~=~phi_k~-~k*phi_1 rather than phi_k). Then specify the number of sub-harmonics and list how the amplitudes and phases for each sub-harmonic are to be generated (periods k*P with k=2...Nsubharm+1). omodel is a flag that is either 1 or 0 depending on whether or not the model light curve will be output, if it is 1 then the output directory is given in modeloutdir, the suffix \".injectharm.model\" will be appended to the filename.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-Injecttransit",14) && strlen(c) == 14))
    {
      listcommands_noexit("-Injecttransit",p,&s);
      printtostring(&s,"Inject a Mandel-Agol transit model into the light curve. The source for the following parameters must be given: period in days (or freq in 1/day), radius of the planet in jupiter radii, mass of the planet in jupiter masses, phase of the transit at time T=0 (phase = 0 corresponds to transit center), sin_i of the transit, e/omega (omega in degrees) or h/k=esinomega/ecosomega, mass of the star in solar masses, radius of the star in solar radii and the limb darkening coefficients. Each parameter can either be specified as a column in the input list (use the *list keywords, by default it will be taken from the next available column, use the \"column\" keyword to specify the column), fixed to a specific value given on the command line (use the *fix keywords), or set to the value of an analytic expression (use the *expr keywords). In some cases the parameter can be generated from a uniform random or uniform log random distribution. If the \"sinirand\" keyword is used then sin_i is drawn from a uniform orientation distribution with the constraint that there must be a transit (for a circular orbit this corresponds to cos(i) being uniform between 0 and (Rstar~+~Rp)/a). The \"dilute\" keyword can optionally be specified. If it is then the transit signal is scaled by the value dilute. omodel is a flag that is either 1 or 0 depending on whether or not the model light curve will be output, if it is 1 then the output directory is given in modeloutdir, the suffix \".injecttransit.model\" will be appended to the filename.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-Jstet",6) && strlen(c) == 6))
    {
      listcommands_noexit("-Jstet",p,&s);
      printtostring(&s,"Calculate Stetson's J statistic, L statistic and the Kurtosis of each light curve. The timescale is the time in minutes that distinguishes between \"near\" and \"far\" observations. The dates file should contain JDs for all possible observations in the first columns - this is used to calculate the maximum possible weight. Note that the J statistic calculated here differs from Stetson's definition by including an extra factor of (sum(weights)/weight_max). Cite Stetson, P.B. 1996, PASP, 108, 851 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-Killharm",9) && strlen(c) == 9))
    {
      listcommands_noexit("-Killharm",p,&s);
      printtostring(&s,"This command whitens light curves against one or more periods. The mean value of the light curve, the period of the light curve and the cos and sin coefficients are output. The light curves passed to the next command are whitened light curves (unless the keyword \"fitonly\" is given). The origin of the period(s) is either from the most recent previous aov command (either aov, or aov_harm), from the most recent previous LS command, two periods one from aov, the other from LS, the period from the most recent injectharm command, Nper periods per1 through perN that are fixed for all light curves, or Nper periods specified in the input list (the periods are read off in order as additional columns in the input light curve list - a list must be used for this option; use the optional \"column\" keyword to specify the column for the first period, subsequent periods are read in order following that column). The light curves will be whitened using Nharm higher-harmonics (frequencies of 2*f0, 3*f0, ... (Nharm~+~1)*f0) and Nsubharm sub-harmonics (frequencies of f0/2, f0/3, ... f0/(Nsubharm~+~1)). omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".killharm.model\" will be appended to the filename. If \"fitonly\" is specified, then the model is not subtracted from the light curve (a keyword is used rather than a flag to maintain compatability with scripts written before this option was added). By default the a_k and b_k cos and sin coefficients are output. If the keyword \"outampphase\" or \"outampradphase\" is given, then the output will be the amplitudes A_k~=~sqrt(a_k^2~+~b_k^2) and the phases (phi_k~=~atan2(-b_k,~a_k)/2pi for \"outampphase\" or phi_k~=~atan2(-b_k,~a_k) for \"outampradphase\"). If the keyword \"outRphi\" or \"outRradphi\" is given then the output will be the relative amplitudes R_k1~=~A_k/A_1 and phases phi_k1~=~phi_k~-~k*phi_1 (in units of 0 to 1, or in radians for the two keywords respectively). Note that for sub-harmonics, k~=~1/2,~1/3, etc. For the fundamental mode the amplitude A_1 and the phase phi_1 will be given. Finally one can also fit the model, applying a clipping to the residuals, and refit the model to the points which passed the clipping. To do this give the \"clip\" keyword, followed by the number of sigma to use for the clipping.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-linfit")))
    {
      listcommands_noexit("-linfit",p,&s);
      printtostring(&s,"Fit a function that is linear in its free parameters to each light curve. \"function\" is the analytic function to fit (e.g. \'a*t^2+b*t+c\' for a quadratic function of time), paramlist is a comma-separated list of free parameters (e.g. 'a1,a2,a3' for the quadratic function). Note that the free parameters should have names that are not used by any vector variables (e.g. t, mag, err, or other variables defined by the -expr command or -inputlcformat option). Note that these variables may be used by other commands as well (e.g. on the right-hand-side of the -expr command). To store the best-fit model in a vector variable for use by later commands, give the \"modelvar\" keyword followed by the variable name. To subtract the best-fit model from the light curve give the \"correctlc\" keyword. To output the model to a file, give the \"omodel\" keyword followed by the directory. By default output files will have names of the form model_outdir/basefilename.linfit.model. Optionally give the \"format\" keyword and then a rule for specifying the name of the output model (see the \"nameformat\" option to the -o command).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-LS",3) && strlen(c) == 3))
    {
      listcommands_noexit("-LS",p,&s);
      printtostring(&s,
"Perform a Generalized Lomb-Scargle (L-S) search of the light curves for periodic sinusoidal signals. The search is done over frequencies between fmin = 1/maxp to fmax = 1/minp, with a uniform frequency step-size of Delta f = subsample/T, where T is the time-span of the observations.\n\n");
      printtostring(&s,
"The program will find the Npeaks strongest peaks in the L-S periodogram. For each peak it will output the period found, log10(FAP) (the logarithm of the formal false alarm probability), and the spectroscopic signal to noise ratio (SNR) for the peaks (SNR = (LS - <LS>)/RMS(LS); where for the Generalized L-S we use the statistic LS=(chi_0^2 - chi(freq)^2)/chi_0^2 where chi_0^2 is chi^2 about the weighted mean and chi(freq)^2 is chi^2 about the best-fit sinusoidal signal with frequency freq; for the non-Generalized case, LS is the traditional normalized L-S statistic).\n\n");
      printtostring(&s,
"operiodogram is a flag that is either 1 to output the periodogram for each light curve to a separate file, or 0 not to. If it is set to 1, then the periodogram will be output to a file with the name $outdir/$basename\".ls\", where $basename is the base filename of the light curve. The output periodogram has the format: freq LS log10(FAP). Here LS is the same as described above. For the Generalized Lomb-Scargle periodogram it is restricted to lie between 0 and 1 (0 means no signal at all, 1 is a perfect sinusoidal variation), while for the traditional periodogram it is unbounded.\n\n");
      printtostring(&s,
		    "By default the Generalized Lomb-Scargle periodogram due to Zechmeister and K\\\"urster 2009, A&A, 496, 577 is calculated. This differs from the tranditional Lomb-Scargle periodogram in allowing for a floating mean and non-uniform errors. If you wish to compute the traditional periodogram give the \"noGLS\" keyword.\n\n");
      printtostring(&s,
"If the \"whiten\" keyword is given, then the light curve will be whitened at each peak period and the periodogram will be recomputed before searching for the next peak period. If \"whiten\" is used then the RMS used in finding the signal-to-noise ratio is computed on the whitened periodogram. The output spectrum will also contain the Npeaks successively whitened periodograms.\n\n");
      printtostring(&s,
"Use the \"clip\" keyword to change the clipping parameters for calculating the average and RMS of the power spectrum when computing the SNR value of a peak. clip is the sigma-clipping factor, and clipiter is a flag that is 1 or 0 to toggle iterative clipping. By default iterative 5-sigma clipping is used.\n\n");
      printtostring(&s,
"If the \"fixperiodSNR\" keyword is given, then log10(FAP) and the SNR will be given for a specific period as well as for the peaks. This period can either be from the most recent issued -aov or -aov_harm command (keyword \"aov\"), from the most recent -LS command (keyword \"ls\"), from the most recent -Injectharm command (keyword \"injectharm\"), it can be fixed to a particular value for all light curves (keyword \"fix\"), it can be specified as an additional column in the input list (keyword \"list\"; by default the next available column is assumed, use the \"column\" keyword to specify the column), or it can be set to an arbitrary previously computed statistic (using the \"fixcolumn\" option together with the name or number of the output column as seen with the -header or -numbercolumns options).\n\n");
      printtostring(&s,
"Cite Zechmeister and K\\\"urster 2009, A&A, 496, 577, and Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. 1992, Numerical Recipes in C, 2nd ed. (New York: Cambridge University Press) if you use this tool (we use an appropriately modified version of the fasper algorithm from Numerical Recipes to calculate the periodogram). If you compute the traditional L-S periodogram rather than the default Generalized one, then the references are Lomb, N.R. 1976, A&SS, 39, 447, Scargle, J.D. 1982, ApJ, 263, 835, Press, W.H. & Rybicki, G.B. 1989, ApJ, 338, 277, and Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. 1992, Numerical Recipes in C, 2nd ed. (New York: Cambridge University Press).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-MandelAgolTransit",18) && strlen(c) == 18))
    {
      listcommands_noexit("-MandelAgolTransit",p,&s);
      printtostring(&s,"Fit a Mandel and Agol (2002) transit model to the light curve. The initial parameters either are automatically estimated using the results from bls if \"bls\" is specified, a blsfixper command, or they are entered directly in the command line as the values for P0, T00 etc. The parameters are P0 - period, T00 - time of transit, r0 - ratio of planet radius to star radius = Rp/R*, a0 - ratio of semi-major axis to star radius = a/R*, either \"i\" - the inclination angle in degrees, or \"b\" - the impact parameter, e0 - eccentricity, omega0 - argument of periastron in degrees, and mconst0 - the out of transit magnitude. If mconst0 is negative then it will be estimated directly from the light curve. \"quad\" or \"nonlin\" are used to specify the limb darkening law, either quadratic or non-linear, following Claret. The number of initial limb darkening coefficients must be 2 for quadratic limb-darkening or 4 for non-linear limb-darkening. fitephem fitr fita ... should be 1 or 0 depending on whether or not the corresponding parameter is allowed to vary. Here fitephem corresponds to varying P0 and T00. If fitRV is 1 then the program will simultaneously fit an RV curve stored in the file RVinputfile (first column JD, second column RV, third column RV error). The model RV will be output to the file RVmodeloutfile (this will be a model evaluated at a number of evenly spaced phase points rather than at the observed RV phases). If fitting RV specify the initial K0 and gamma0 together with flags fitK and fitgamma to specify whether or not these parameters are allowed to vary in the fit. Set correctlc to 1 to subtract the model from the light curve, set it to 0 to only perform the fit. omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".mandelagoltransit.model\" will be appended to the filename. Use the \"modelvar\" keyword to save the best-fit model light curve to the variable var.  If the optional keyword \"ophcurve\" is given then a model phase curve will be output to a file in the directory outdir with suffix \".mandelagoltransit.phcurve\" with phases between phmin and phmax and a uniform step size of phstep. If the optional keyword \"ojdcurve\" is given then a model light curve will be output to a file in the directory outdir with suffix \".mandelagoltransit.jdcurve\" with times between the first and last times in the light curve and with a uniform step size of jdstep. Cite Mandel, K., & Agol, E. 2002, ApJ, 580, L171 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-medianfilter",13) && strlen(c) == 13))
    {
      listcommands_noexit("-medianfilter",p,&s);
      printtostring(&s,"Apply either a high-pass or a low-pass median filter to the light curve. By default, a high-pass filter is applied, i.e. the median magnitude of all points within 'time' of a given observation is subtracted from that observation. If \"average\" or \"weightedaverage\" is specified then it will be the average magnitude or the weighted average magnitude that is subtracted. If the \"replace\" keyword is given, then each point in the light curve will be replaced with the running median (or the average or the weighted average). In this case the command acts as a low-pass filter.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-microlens",10) && strlen(c) == 10))
    {
      listcommands_noexit("-microlens",p,&s);
      printtostring(&s,
	      "Fit a simple microlensing model to the light curve. This program uses the functional form for the model given by Wozniak, P. R., et al. 2001, AcA, 51, 175. For each parameter you can optionally specify the source for the initial parameter value: either \"fix\" the initial value for all light curves to the value fixval, use \"list\" to specify the value as an additional column in the input list (by default it will be the next available column, use the \"column\" keyword to specify the column), use \"fixcolumn\" to take the initial value from a previously computed statistic, or use \"auto\" to automatically determine the initial parameter value. The default for each parameter is \"auto\". You can then optionally specify an initial step size for each parameter for use in the downhill simplex fit. You may also optionally use \"novary\" to not vary this parameter in the fit. Use the \"correctlc\" keyword to subtract the model fit from the light curve before passing to the next command, and the \"omodel\" keyword to output the model. The suffix will be \"microlens\".\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strcmp(c,"-nonlinfit")))
    {
      listcommands_noexit("-nonlinfit",p,&s);
      printtostring(&s,"Fit a function that is nonlinear in its free parameters to each light curve. \"function\" is the analytic function to fit (e.g. 'a*exp(-(t-t0)^2/2/sigma^2)+b' for a Gaussian function in time).\n\n");
      printtostring(&s,"\"paramlist\" is a comma-separated list of parameters to vary. For each parameter you must specify the initial guess and step-size using the format var=init:step. For example, if the initial value for t0 is 5.0, and its step-size is 2.0, and if the initial value for sigma is 10.0 and its step-size is 8.0 you would give the expression 't0=5.0:2.0,sigma=10.0:8.0'. You may also use more complicated analytic expressions involving statistics computed by prior VARTOOLS commands to initialize the variables. Note that the free parameters should have names that are not used by any vector variables (e.g. t, mag, err, or other variables defined by the -expr command or -inputlcformat option). Note that these variables may be used by other commands as well (e.g. on the right-hand-side of the -expr command).\n\n");
      printtostring(&s,"If there are any parameters in the function which enter linearly, you may optimize these using linear least squares by giving the keyword \"linear\" and then providing the list of parameters. For example, this could be 'linear a,b' in the Gaussian case. Doing this speeds up the optimization for these parameters, however they will be excluded from the posterior distribution used to calculate uncertainties or other statistics if using the mcmc fitting method.\n\n");
      printtostring(&s,"You may use an analytic expression for the magnitude uncertainties used in calculating the likelihood function by giving the keyword \"errors\" followed by the expression. Parameters which are varied in the fit may be included in this expression.\n\n");
      printtostring(&s,"By default this command assumes uncorrelated errors. You may optionally allow for correlated uncertainties by using the \"covariance\" keyword. This command support three different noise models: a squared exponential model where the covariance between times t_i and t_j is proportional to amp_var*exp(-(t_i-t_j)^2/(2*rho_var)^2), an exponential model where the covariance is proportional to amp_var*exp(-|t_i-t_j|/rho_var), and a so-called Matern model where the covariance is given by amp_var*C(nu_var,x)*K_nu(x) with x=sqrt(2*nu_var)*(|t_i-t_j|)/rho_var and C(x,y)=(2^(1-x)/Gamma(x))*(y)^x with K_nu being the modified Bessel function of the second kind, and Gamma the usual gamma function. In this case nu_var>0, rho_var>0 and amp_var>0 are the three parameters characterizing the model, with rho_var being the correlation length scale, amp_var being the amplitude, and nu_var characterizing the shape of the correlation tale (when nu_var=0.5 the correlation is exponential in |t_i-t_j| whereas when nu_var -> infinity it converges to the squared exponential model). To use the squared exponential model or the exponential model provide the \"squareexp\" keyword, or the \"exp\" keyword, following \"covariance\" and then provide the names of the variables that will be used for amp_var and rho_var, respectively. For the Matern model provide the \"matern\" keyword and the names of the variables to be used for amp_var, rho_var, and nu_var. If these variables do not appear in the parameter list, then they should be specified as 'amp_var=expr' where amp_var is the name of the variable to use, and expr is an expression used to determine the fixed value to be used for this parameter. If they are provided in the free parameter list, then just the name of the variable should be given here. The covariance parameters will be forced by the program to be positive and non-zero. Note that linear fitting of a subset of the parameters is not permitted when the covariance option is used.\n\n");
      printtostring(&s,"You may place priors on the variables using the keyword 'priors' and then providing a comma-separated list. Each entry in the list should evaluate to -2*ln(P) where P is the prior probability for a variable given its value. For example, to place a Gaussian prior on t0 with mean 4.0 and standard deviation 3.0 you would use 'prior (t0-4.0)^2/3.0^2'.\n\n");
      printtostring(&s,"To place constraints on the parameters use the keyword 'constraints' and provide a comma-separated list. For example, to require sigma > 0, use 'constraints sigma>0'.\n\n");
      printtostring(&s,"You must specify a method for fitting the function. The allowed values are 'amoeba' to run a down-hill simplex optimization, or 'mcmc' to run a Differential Evolution Markov Chain Monte Carlo routine.\n\n");
      printtostring(&s,"If you use amoeba, you may also optionally specify the convergence tolerance which is minimum fractional change in chi^2 between iterations, or maxsteps which is a maximum number of iterations to try before giving up. A column will be included in the output ascii table indicating whether or not the fit converged.\n\n");
      printtostring(&s,"If you use mcmc, then you may also optionally specify the number of accepted links to run in a given fit ('Naccept') or the total number of links to run ('Nlinkstotal'). By default the Nlinkstotal option is assumed with a value of 100000. You may also specify an initial fraction of the chain to ignored in computing statistics from the posterior distribution ('fracburnin'; the default is 0.1). The optional 'eps' parameter controls the scale of random variations allowed in generating proposed parameter sets using the differential evolution model (the default is 0.001). By default the downhill simplex algorithm is run initially before beginning the mcmc procedure, this can be skipped by giving the \"skipamoeba\" keyword. By default the median and standard deviation for each of the parameters varied in the mcmc will be included in the output ascii table. You may change the statistics, and or the quantities that are used with the 'chainstats' keyword. You must then provide a comma-separated list of analytic expressions to calculate from the chain, and a comma-separated list of statistics to report for each of the expressions. The available statistics are the same as for the '-stats' command. You may set a maximum limit to the total amount of memory to be used by this mcmc run with the 'maxmemstore' keyword followed by the limit in GB. This will limit the length of the chain used for calculating the statistics. The default is 4.0.  If you wish to output the chains themselves, provide the 'outchains' keyword followed by the directory to output the chains to. Each chain will be written to a file of the form outdir/lcname.mcmc, where lcname is the basename of the input light curve file. You can optionally modify the naming convention using the \"format\" keyword (see the \"nameformat\" option to the -o command for syntax). The chain file will have one column per quantity calculated (by default these are just the varied parameters) followed by -2*ln(L) where L is the likelihood for that link. You can reduce the number of links output to the chain file using the 'printevery' keyword.\n\n");
      printtostring(&s,"To store the best-fit model in a vector variable for use by later commands, give the \"modelvar\" keyword followed by the variable name. To subtract the best-fit model from the light curve give the \"correctlc\" keyword. To output the model to a file, give the \"omodel\" keyword followed by the directory. Output files will have names of the form model_outdir/basefilename.nonlinfit.model. Optionally give the \"format\" keyword and then a rule for specifying the name of the output model (see the \"nameformat\" option to the -o command).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-o",2) && strlen(c) == 2))
    {
      listcommands_noexit("-o",p,&s);

      printtostring(&s,"Output the light curves to directory outdir or to the file outname. If a light curve list is used, the directory form will be used, if a single light curve is read in, then the outname form will be used.\n\nThe default output filename for the outdir form is: $outdir/$inname where inname is the base filename of the input light curve. You can optionally specify a format rule for the output name by giving the \"nameformat\" keyword followed by the formatstring. In that case the output filename will be $outdir/$formatstring with instances of %%s replaced with $inname, instances of %%d replaced with the light curve number (starting with 1), instances of %%0nd where n is an integer replaced with the formatted light curve number, and instances of %%%% will be replaced with %%. For example, if the second line in the file \"inlist\" is \"tmp/file2.lc\", the command \"vartools -l inlist -rms -o ./directory nameformat file%%s%%05d.txt\" would result in copying the file \"tmp/file2.lc\" to \"./directory/file2.lc00002.txt\". If a single light curve is read in, then the parameter given to -o is the name of the output light curve. The \"nameformat\" option and formatstring will be ignored if they are given. If \"-\" is given for outname, then the light curve will be output to stdout. In that case you should also use the -quiet command to avoid mixing the output light curve with the output statistics.\n\nBy default the output light curves will have three columns: time, mag, and err. You can use the \"columnformat\" keyword to change this format. The formatstring is a comma-separated list of variable names to output, optionally using a colon after each variable name to specify the printf format to use for that variable. For example, \"columnformat t:%%.17g,mag:%%.5f,err:%%.5f,xpos:%%.3f\" would output the variables t, mag, err, and xpos using formats %%.17g, %%.5f, %%.5f, and %%.3f respectively. Here xpos is a non-default variable that one would have read-in with the -inputlcformat command. If the light curves are output in fits format, then terms after the colon will be used to specify the units of the column in the light curve header.\n\nBy default a single space character is used to delimit columns when outputing ascii data. You can change the character used for delimiting columns by giving the keyword \"delimiter\" followed by the character to use.\n\nBy default light curves are output in ascii format. Give the keyword \"fits\" to output the light curves in binary fits table format. The output light curve will have the extension \".fits\" appended if it is not already present. The keyword \"noclobber\" may be used to prevent overwritting any existing files. VARTOOLS will terminate if it encounters an existing file with noclobber set.\n\n");
      commandfound = 1;
    }
  if(all == 1 || ((!strncmp(c,"-Phase",6) || !strncmp(c,"-phase",6)) && strlen(c) == 6))
    {
      listcommands_noexit("-Phase",p,&s);
      printtostring(&s,"Replace the time variable of a light curve with its phase and sort the light curve by the phase. The period for phasing is either taken from the last aov command, the last ls command, the last BLS command, is set to an arbitrary previously computed statistic by giving the name or number of the output column as seen with the -header or -numbercolumns options (the \"fixcolumn\" option), is specified in the list (in this case the list must be used), or is fixed to the value P specified on the command line. The user may optionally specify a reference time for phase = 0 (T0) by giving \"T0\" followed either by \"bls\" and the phase to assign the time of mid-transit found from the last bls command, by \"fixcolumn\" and the name or number of a previously computed statistic, by \"list\" (in which case the value to use for each light curve is given in the input list), or by \"fix\" and the T0 value to use for all processed light curves. By default t0 is the initial time in the light curve. If the \"phasevar\" keyword is used, then the phases will be stored in the variable var, rather than overwriting the times. The optional \"startphase\" keyword may be used to change the range of phases (by default they run from 0 to 1, using this keyword will cause them to run from startphase to startphase+1).\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-python"))
    {
      listcommands_noexit("-python",p,&s);
      printtostring(&s,"Run arbitrary python routines on a light curve. Either provide a string on the command line with python code to apply to the light curve, or give the \"fromfile\" keyword followed by the name of a file with python code to run. VARTOOLS will take the commands supplied by the user and embed them within a python function which receives variables from VARTOOLS, executes the code supplied by the user, and then returns the variables back to VARTOOLS. VARTOOLS will use the python C library to compile this function and then execute it on each light curve to be processed. Note that in order to access the variable data VARTOOLS will automatically \"import numpy\", so you may use routines within this module without importing it explicitly.\n\nAny code that you would like to execute through python before running the processing function on each light curve (e.g., code that defines your own functions that are then called through the commandstring) can be supplied by providing the \"init\" keyword and then either giving the code as a single string on the command line, or giving the \"file\" keyword followed by the name of a file from which the code will be imported.\n\nBy default VARTOOLS will run every -python command as a separate sub-process, meaning, for example, that any changes that you make to global variables in one -python command will not be seen by other -python commands. Similarly each command will only use its own initialization code. You can optionally combine multiple calls to -python on the VARTOOLS command-line into a single sub-process by using the \"continueprocess\" keyword. Supply the keyword and then the number of the prior -python command whose sub-process a new -python command should use. The first call to \"-python\" on the command line will have prior_python_command_nnumber=1, the second will have prior_python_command_number=2 and so on. If you use the \"continueprocess\" keyword then you will not supply any initialization code for this command. The command will instead inherit the initialization code that was executed for the prior -python command.\n\nBy default all active variables will be passed from VARTOOLS into the python function and will then be returned to VARTOOLS. The input variables are provided as numpy arrays, except in the case of variables storing string data which are supplied as lists. To specify the variables which are input and output from the command you can use the \"vars\" keyword followed by a comma-separated list of variable names, or you can use the \"invars\" and/or \"outvars\" keywords to separately specify the variables that are input and output. If you would like any variables to be output to the ASCII table (e.g., statistics that you compute from the light curves) use the \"outputcolumns\" keyword followed by a comma-separated list of variable names. Note that light curve vectors or variables storing string data cannot be included in this list.\n\nBy default VARTOOLS will only send one light curve at a time to the user function. If, however, you wish to pass all of the light curves to PYTHON at once, you can supply the \"process_all_lcs\" keyword. In that case light curve vectors will be supplied as a list of numpy arrays, while other variables, such as those storing values computed for each light curve by a command, will be provided as numpy arrays.\n\nNote that if VARTOOLS is run with the -parallel option, then a separate python sub-process will be run for each thread. This allows for parallel processing of the light curves through python without using the python global interpreter lock. Note that because we are using separate processes, the initialization code will be executed separately for each thread, meaning extra over-head for each thread used, and any changes that you make to global variables in one thread will not be visible to other threads.\n\nSee \"vartools -example -python\" for examples of how to use this command.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-resample"))
    {
      listcommands_noexit("-resample",p,&s);
      printtostring(&s,"Resample the light curve onto a new time-base. All light curve vectors will be resampled. First specify the method for resampling. Options include \"nearest\" to set resampled values to the value of the observation that is closed in time, \"linear\" to do linear interpolation between points, \"spline\" to do cubic spline interpolation, \"spline_monotonic\" to do cubic spline interpolation with the interpolating function constrained to be monotonic between input observations, or \"bspline\" to interpolate with a Basis-spline function. For \"spline\" interpolation you may optionally specify the left and right boundary conditions for the spline. For \"bspline\" interpolation you may optionally specify the number of breaks (the default is 15) and the order of the spline function (the default is 3). If the number of breaks is < 2, then the routine will increase the breaks until chi squared per degree of freedom is less than one. Caution, this can be quite slow. Light curve vectors containing text data rather than numeric data (e.g., image ids) will be resampled using the \"nearest\" method irrespective of what method is specified for the numeric data.\n\n");
      printtostring(&s,"By default the code will resample the light curve onto a uniform time base running from the first observed time in the light curve to the last time, using a time step equal to the minimum time separation between consecutive points in the input light curve. You may change the start time, the stop time, the time separation (or, instead, the number of points in the resampled light curve) using the \"tstart\", \"tstop\", and \"delt\" or \"Npoints\" keywords. For each quantity you may either \"fix\" it to a specified value for all light curves, you may have it set to a column output from the a previous command (\"fixcolumn\"), you may read the value from the input light curve \"list\", or you can provide an analytic \"expr\"ession which is evaluated to determine the value of the quantity for each light curve. If you wish to use a general time-base (not necessarily uniform) give the \"file\" keyword and the indicate a file to read the new times from. This may either be a \"fix\"ed file that is used for all light curve (use the \"column\" keyword after the file name to give the column in the file containing the times to use, by default the first column is assumed), or you can use a different file for each light curve with the name of the file read-in from the input light curve \"list\". If you read the filenames from the light curve list, by default the next unused column in the list will be assumed. Use the \"listcolumn\" keyword to change which column in the list gives the names of the files. You can also change the column from which the times are read within these files using the \"tcolumn\" keyword.\n\n");
      printtostring(&s,"By default the same resampling method will be used for all points in the light curve. By giving the \"gaps\" keyword you can make the method depend on how far away a resampled time is from the closest observed time. You first need to indicate how the time separation used to distinguish between the near and far points will be determined. The \"fix\" keyword sets it to a particular value given on the command line. The \"fixcolumn\" keyword sets it to the result of a previously executed command. The \"list\" keyword will cause the values to be tken from the input light curve list file (use the \"column\" keyword to specify the column, otherwise the next available column will be assumed). Use the \"expr\" keyword to provide an analytic expression to be evaluated for each light curve. Use the \"frac_min_sep\" keyword to set the time separation to a fixed factor times the minimum separation between subsequent points in the input light curve (e.g. if you give \"frac_min_sep 5.0\" and the minimum separation between points in the light curve is 1 day, then the separation timescale will be set to 5 days). Use the \"frac_med_sep\" keyword to set the time separation to a fixed factor times the median separation between subsequent points in the input light curve. Or use the \"percentile_sep\" keyword to use a percentile of the input separations (e.g., if you give \"percentile_sep 70\", then the N separations in the input light curve will be ordered, and the separation to be used for distinguishing between near and far points will be set equal to the 0.7*N longest separation in the input light curve.).  After specifying how the separation is to be determined, you must next indicate the resampling method to be used for the far resampled points. The method specified before the \"gaps\" keyword will then apply to the near resampled points. The resampling options here are the same as discussed above.\n\n");
      printtostring(&s,"You may also use a different method for resampled points which are extrapolations rather than interpolations. Give the \"extrap\" keyword, followed by the method to use (options are as already discussed).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-rescalesig",11) && strlen(c) == 11))
    {
      listcommands_noexit("-rescalesig",p,&s);
      printtostring(&s,"Rescale the magnitude uncertainties of each light curve so that Chi2 per dof is equal to 1 for every light curve. The rescale factor for each light curve will be included in the output table.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-restorelc",10) && strlen(c) == 10))
    {
      listcommands_noexit("-restorelc",p,&s);
      printtostring(&s,"This command can be used in conjunction with the -savelc command to restore the light curve to a previous state. savenumber is 1, 2, 3, ... to restore the light curve to its state at the first, second, third ... call to -savelc. For example, suppose you want to try running TFA followed by LS on the light curve using 3 different template lists. A command string of the form: \"... -savelc -TFA trendlist1 ... -LS ... -restorelc 1 -TFA trendlist2 ... -LS ... -restorelc 1 -TFA trendlist3 ... -LS ...\" would accomplish this as each time \"-restorelc 1\" is called the light curve is restored to its form at the first -savelc command.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-restricttimes"))
    {
      listcommands_noexit("-restricttimes",p,&s);
      printtostring(&s,"This command is used to filter observations from the light curves based on the time values. If the optional keyword \"exclude\" is given, then specified times will be removed from the light curve(s), otherwise non-specified times will be removed. The times to include or exclude can be done either by specifying a range of times (keyword \"JDrange\" or \"JDrangebylc\"), by providing a list of times (\"JDlist\"), or by providing a list of string IDs (\"imagelist\"). If \"JDrange\" is used then the same range of times is set for all light curves. If \"JDrangebylc\" is used, then the user may specify a different minJD and/or maxJD for each light curve. The options here are \"fix\" to fix the parameter to a value given on the command-line, \"list\" to have the parameter read-in from the input light curve list (optionally giving the column number with the keyword \"column\", by default the next column in the list is used), \"fixcolumn\" to use a previously computed column from the output statistics file, or \"expr\" to set the value to an analytic expression which is evaluated for each light curve at the time this command is executed. If \"JDlist\" is used then the user should supply a file which contains a list of JDs in the first column. By default points with times not given in this file will be removed from all light curves, but if the \"exclude\" option is given, then the times in the file will be removed from the light curves. If \"imagelist\" is used then the user should supply a file which contains a list of string-IDs (image names) in the first column to use in selecting points from the light curves. If \"expr\" is used then the user should supply an analytic expression to be evaluated for each point in the light curve. Only points evaluating to a value greater than zero will be included (or excluded if the \"exclude\" keyword is used).  For example giving the command \"-restricttimes expr '(mag>9.0)&&(mag<9.5)'\" would restrict the light curve to only points in the range 9.0 < mag < 9.5, filtering out all other points.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-restoretimes"))
    {
      listcommands_noexit("-restoretimes",p,&s);
      printtostring(&s,"Restore the observations in the light curve that were filtered out through a prior -restricttimes command. The value for prior_restricttimes_command should be set to the number of the -restricttimes command that you want to restore the times from, where prior_restricttimes_command=1 for the first -restricttimes command given on the command line, =2 for the second -restricttimes command, etc. The restored points are appended to the current light curve, and the light curve is then sorted by time. You can use the -restricttimes and -restoretimes commands to apply modifications to isolated portions of the light curve.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-rms",4) && strlen(c) == 4))
    {
      listcommands_noexit("-rms",p,&s);
      printtostring(&s,"Calculate the rms of the light curves. The output will include RMS, the mean magnitude, the expected RMS and the number of points in the light curve\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-rmsbin",7) && strlen(c) == 7))
    {
      listcommands_noexit("-rmsbin",p,&s);
      printtostring(&s,"Similar to chi2bin, this calculates RMS after applying a moving mean filter.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-savelc",7) && strlen(c) == 7))
    {
      listcommands_noexit("-savelc",p,&s);
      printtostring(&s,"This command can be used in conjunction with the -restorelc command to save a light curve and later restore it to this state. For example, suppose you want to try running TFA followed by LS on the light curve using 3 different template lists. A command string of the form: \"... -savelc -TFA trendlist1 ... -LS ... -restorelc 1 -TFA trendlist2 ... -LS ... -restorelc 1 -TFA trendlist3 ... -LS ...\" would accomplish this as each time \"-restorelc 1\" is called the light curve is restored to its form at the first -savelc command.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-SoftenedTransit",16) && strlen(c) == 16))
    {
      listcommands_noexit("-SoftenedTransit",p,&s);
      printtostring(&s,"Fit a Protopapas, Jimenez and Alcock transit model to the light curve. The initial parameters either come from bls if \"bls\" is specified, a blsfixper command if \"blsfixper\" is specified, otherwise they are entered directly in the command line as the values for P0 T0 etc. If mconst0 is < 0 then it will be estimated directly from the light curve. fitephem fiteta fitcval ... should be 1 or 0 depending on whether or not the corresponding parameter is allowed to vary. Here fitephem is used to fit both P and T0. Set correctlc to 1 to subtract the model from the light curve, set it to 0 to only perform the fit. omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".softenedtransit.model\" will be appended to the filename. If fit_harm is 1 then simultaneously fit a harmonic series to the light curve. The period comes from either the last aov, the last ls, the last bls model or is specified in the command line as Pharm (you must type \"list\" and then the value), finally the number of harmonics and sub-harmonics must also be specified. Cite Protopapas, P., Jimenez, R., & Alcock, C. 2005, MNRAS, 362, 460 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-Starspot",9) && strlen(c) == 9))
    {
      listcommands_noexit("-Starspot",p,&s);
      printtostring(&s,"This command fits a single, circular, uniform temperature spot model to the light curve using the Dorren 1987 model. Use aov or ls to choose the initial period to be from AoV or LS, use \"list\" to specify the initial period in the input list (by default it will be the next available column, use the \"column\" keyword to specify the column), use \"fix\" to fix the period to a value given on the command line, or set the value to any previously computed statistic by giving \"fixcolumn\" and the name or number of the previously computed statistic. The other parameters are initial guesses as in Dorren 1987, mconst is the constant magnitude term. alpha0~=~spot radius is degrees, i0~=~stellar rotation axis inclination in degrees (90~=~rotation axis perpendicular to line-of-sight), chi0~=~spot latitude in degrees (0~=~equator, positive if spot center is between equator and pole), psi00~=~longitude of the spot center in degrees. Set mconst0 to a negative value to have the program automatically estimate the initial mconst. Note that here we take a~=~a_d/(pi~*~(1~-~mu_*~/~3)) and b~=~b_d/(pi~*~(1~-~mu_*~/~3)) where a_d and b_d are the a and b terms from Dorren 1987 (for the Sun at 5000 angstroms one could use a=0.0298 and b=0.08745). fitP... are flags denoting whether or not a parameter is to be varied. Set correctlc to 1 to subtract the model from the light curve, set it to 0 to only perform the fit. omodel is a flag set to 1 or zero that can be used to output the model for the light curve, the output directory is then given in modeloutdir, the suffix \".starspot.model\" will be appended to the filename. Cite Dorren 1987, ApJ, 320, 756 if you use this tool.\n\n");
      commandfound = 1;
    }
  if(all == 1 || !strcmp(c,"-stats"))
    {
      listcommands_noexit("-stats",p,&s);
      printtostring(&s,"Compute statistics on one or more light-curve vectors (e.g. t, mag, err). var1,var2,... is a comma-separated list of variable names to compute the statistics on. stats1,stats2,... is a comma-separated list of one or more statistics to compute for each variable (every statistic is computed for every variable). Available statistics include:\n");
      printtostring(&s,"\tmean\n");
      printtostring(&s,"\tweightedmean - mean of the vector weighted by the light curve uncertainties.\n");
      printtostring(&s,"\tmedian\n");
      printtostring(&s,"\twmedian - median vector weighted by light curve uncertainties\n");
      printtostring(&s,"\tstddev - standard deviation calculated with respect to the mean.\n");
      printtostring(&s,"\tmeddev - standard deviation calculated with respect to the median.\n");
      printtostring(&s,"\tmedmeddev - Median of the absolute deviations from the median.\n");
      printtostring(&s,"\tMAD - 1.483*medmeddev. For a gaussian distribution this equals stddev in the limit of large N.\n");
      printtostring(&s,"\tkurtosis\n");
      printtostring(&s,"\tskewness\n");
      printtostring(&s,"\tpct%%f - %%f percentile, where %%f is a floating point number between 0 and 100. Here 0 corresponds to the minimum value and 100 to the maximum value in the vector.\n");
      printtostring(&s,"\twpct%%f - percentile including light curve uncertainties as weights.\n");
      printtostring(&s,"\tmax - maximum value, equivalent to pct100.\n");
      printtostring(&s,"\tmin - minimum value, equivalent to pct0.\n");
      printtostring(&s,"\tsum - sum of all elements in the vector.\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-SYSREM",7) && strlen(c) == 7))
    {
      listcommands_noexit("-SYSREM",p,&s);
      printtostring(&s,"Run the SYSREM algorithm to identify and remove trends from an ensemble of light curves. If this command is used, the -readall option will automatically be set. Also, if this command is used then one must have as input a light curve list. A total of Ninput_color + Ninput_airmass trends will be searched for. Ninput_color of these will have the colors set initially, Ninput_airmass will have the airmass terms set initially. The initial color terms will be read in as additional columns in the input light curve list, there must be Ninput_color of these columns (by default they will be read starting from the next available column, use the \"column\" keyword to specify the column for the first color term, subsequent terms will be read in order from the following columns). The initial_airmass_file gives a list of the initial airmass trends to use, the first column is JD (or the string-ids for the -matchstringid option) and the subsequent Ninput_airmass columns are the initial airmass trends. As implemented this routine first filters the trends with the airmass terms specified initially and then filters the trends with the color terms specified initially. sigma_clip1 is the sigma clipping used in calculating the mean magnitudes for the light curves. sigma_clip2 is the sigma clipping used in determining whether or not points contribute to the airmass or color terms when doing the fit. Any points with magnitude less than saturation will not contribute. omodel is a flag that is 1 or zero for outputing the model light curves, the output directory is given in model_outdir, the suffix \".sysrem.model\" will be appended to the filename. The format for the model light curve is: JD mag mag_model sig clip, where clip is either 1 if the point was included in determining the trends or 0 if it wasn't. otrends is a flag that is 1 or zero for outputing the final trends. These will be output to the file trend_outfile, the first column is the JD and subsequent columns are for each trend signal. If you use this command, be sure to cite Tamuz, Mazeh and Zucker (2005, MNRAS, 356, 1466).\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-TFA",4) && strlen(c) == 4))
    {
      listcommands_noexit("-TFA",p,&s);
      printtostring(&s,"Run the Trend Filtering Algorithm on the light curves. The trendlist is a file containing a list of light curves to be used as trends, it has the format: trendname trendx trendy, where trendname is the file name and trendx and trendy are the coordinates of the trend star. You can optionally specify the readformat for the light curves in the trendlist, by default Nskip is 0, jdcol is 1 and magcol is 2. Note that if -matchstringid is set then the JDs of the trend light curves are not read-in and instead jdcol corresponds to the column storing the string-ids for point matching between light curves. dates_file is a file giving all the dates in the second column, the unique string-ids are read-in from the first column if the -matchstringid option is used (it is the same format as the dates file for the ISIS image subtraction program where the first column is the filename and the second is JD). Trend stars within pixelsep of the light curve in question will not be used in detrending the light curve. To you this routine you need an input light curve list, the x and y positions of each light curve must be given as columns in the list. By default these are the next available columns, use the \"xycol\" keyword to specify the columns. correctlc is a flag that is 0 or 1 indicating whether the light curves passed to the next command should have the tfa applied. ocoeff is a flag to denote whether or not to output the list of coefficients for each trend for a given light curve, if set to 1 then they will be output to coeff_outdir with \".tfa.coeff\" appended to the end of the filename. omodel is a flag set to 1 or zero that can be used to output the tfa model for the light curve, the output directory is then given in model_outdir, the suffix \".tfa.model\" will be appended to the filename. If you use this routine be sure to cite Kovacs, Bakos and Noyes, 2005, MNRAS, 356, 557\n\n");
      commandfound = 1;
    }
  if(all == 1 || (!strncmp(c,"-TFA_SR",7) && strlen(c) == 7))
    {
      listcommands_noexit("-TFA_SR",p,&s);
      printtostring(&s,"Run the Trend Filtering Algorithm in signal reconstruction mode on the light curves (i.e. iteratively filter the light curve and fit a simple signal to the light curve). In addition to the parameters used by the -TFA command, this command allows for simulataneously decorrelating the light curve against additional light curve specific signals. This is useful, for example, if one wishes to do EPD (external parameter decorrelation) and TFA on a high-amplitude variable star. To do this, specify \"decorr\", the iterativeflag is 1 if the decorrelation and the TFA will be done iteratively (this is faster) or 0 if they will be done simultaneously (more correct, but slower). The Nlcterms, lccolumn and lcorder terms are the same as for the -decorr command. Any global signals to decorrelate against should just be included in the TFA trendlist. Other parameters that are different from the -TFA command include: \"dotfafirst\" is a flag that is 0 or 1, if the flag is 1 then TFA will be applied to the input light curve first and the signal will then be determined on the residual in each iteration, if it is 0 then the signal is determined and subtracted from the light curve and TFA is applied to the residual in each iteration. The iterations will stop once the fractional change in the RMS is less than \"tfathresh\" or if the number of iterations reaches \"maxiter\". The model signal is either taken to be the mean value of the binned light curve (specify \"bins\" and then the number of bins to use), it is a fixed signal form read in from a file (specify \"signal\" and then a filename), or it is a Fourier series that is simultaneously fit to the light curve with TFA (specify \"harm\"). If binning is used, then you can specify \"period\" to make the binning phase-binning, and then specify where the period should be read from. If \"list\" is used then the period is read from the input light curve list (by default it is taken from the next available column, use the \"column\" keyword to specify the column), if \"fix\" is used then the period is fixed to the specified period value for all light curves. If \"signal\" is used, then the filename will be a file containing a list of signal files, one file for each light curve, with each signal file containing the signal in the second column. The quantity a*signal + b is fit to the light curve, where a and b are free parameters. If \"harm\" is used, then the signal will be a Fourier series that is fit simultaneously to the light curve with TFA. In this case there is no TFA iteration. Nharm is the number of harmonics in addition to the fundamental to include in the fourier series, Nsubharm is the number of subharmonics. If \"period\" is specified then period for the fourier series will be taken from the specified source. If \"period\" is not specified then the period for the Fourier series will be set to the time-span of the observations. If you use this routine be sure to cite Kovacs, Bakos and Noyes, 2005, MNRAS, 356, 557\n\n");
      commandfound = 1;
    }
  if(all || (!strcmp(c,"-wwz")))
    {
      listcommands_noexit("-wwz",p,&s);
      printtostring(&s,"Run the Weighted Wavelet Z-Transform as defined by Foster, 1996, AJ, 112, 1709 using an abbreviated Morlet wavelet given by f(z)=exp(i*2*pi*f*(t-tau)-c*(2*pi*f)^2*(t-tau)^2). The wavelet transform will be calculated for frequencies from maxfreq. Give the \"maxfreq\" keyword and then either the maximum frequency (in cycles per day) to consider, or use \"auto\" to 1.0/(2.0*delmin) where delmin is the minimum time seperation between consecutive points in the light curve. The frequency sampling will be given by freqsamp/T where T is the time base-line spanned by the light curve. Use the \"freqsamp\" keyword to give the value of freqsamp. The wavelet transform will also be calculated for time-shifts between tau0 and tau1 in step-sizes of dtau. In each case you may either specify the value to use for that parameter or give \"auto\". If you use \"auto\" for tau0, it will use the minimum time in the light curve. If you use \"auto\" for tau1, it will use the maximum time in the light curve. If you use \"auto\" for dtau, it will use delmin.  By default this routine assumes c=(1/(8*pi^2)). You may change this value by giving the \"c\" keyword and then providing the value.  The routine will output on the ascii table the maximum value of the Z-transform, together with the corresponding frequency, power, amplitude, effective number of points, time-shifts and local mean light curve value. It will also output the median value over all time-shifts of the maximum Z at each time-shift, the median frequency, power, amplitude, effective number of points, and average magnitude. To output the full wavelet transform (at each trial time-shift and frequency) give the \"outfulltransform\" keyword followed by the directory to output these to.");
#ifdef USECFITSIO
      printtostring(&s," You may output these as multi-extension fits image files by giving the \"fits\" keyword.");
#endif
      printtostring(&s," If you give the \"pm3d\" keyword then a blank line will be included after each Time-step scan. The output data file is then in the format expected by the gnuplot pm3d plotting style.");
      printtostring(&s," By default these will be output to files named outdir/lcname.wwz, where lcname is the basename of the input light curve file. You can optionally modify the naming convention using the \"format\" keyword (see the \"nameformat\" option to the -o command for syntax). To output the transform that maximizes Z over frequencies as a function of time-shift, give the \"outmaxtransform\" keyword followed by the directory to output these to. By default these will be output to files named outdir/lcname.mwwz. You can optionally modify the naming convention using the \"format\" keyword. If you use this routine cite Foster, 1996, AJ, 112, 1709.\n\n");
      commandfound = 1;
    }
  if(s.s != NULL)
    {
      if(all) {
	printf(s.s);
      }
      else
	fprintf(stderr,s.s);
    }
#ifdef DYNAMICLIB
    for(i=0; i < p->NUserLib; i++) {
      if(all == 1 || (!strcmp(c,p->UserLib[i].commandname)))
	{
	  if(all) {
	    p->UserLib[i].ShowSyntax_function(stdout);
	    fprintf(stdout,"\n");
	    p->UserLib[i].ShowHelp_function(stdout);
	    fprintf(stdout,"\n");
	  }
	  else {
	    p->UserLib[i].ShowSyntax_function(stderr);
	    fprintf(stderr,"\n");
	    p->UserLib[i].ShowHelp_function(stderr);
	    fprintf(stderr,"\n");
	  }
	  commandfound = 1;
	}
    }
#endif
  if(s.s != NULL)
    free(s.s);
  if(!commandfound)
    error2(ERR_HELP_BADCOMMAND,c);
  exit(ERR_USAGE);
}

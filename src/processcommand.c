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

/* This file contains functions that launch the light curve processing commands for the program vartools by J. Hartman */

void SkipCommand(ProgramData *p, Command *c, int thisindex, int lc, int lc2)
{
  int i;
  if(thisindex > p->max_colcommand)
    return;
  if(p->col_commandstart[thisindex] < 0)
    return;
  for(i = p->col_commandstart[thisindex]; i <= p->col_commandstop[thisindex]; i++) {
    zerooutcolumnvalue(&(p->outcolumns[i]), lc2, lc);
  }
}

void ProcessCommandSingle(ProgramData *p, Command *c, int lc, int thisindex, int lc2)
{
  double d0, d1, d2, d3, d4;
  double *d1ptr, *d2ptr, *d3ptr, *d4ptr;
  int i1, i2, i3, i4, i;
  char *s1;
  char outname[MAXLEN], outname2[MAXLEN], outname3[MAXLEN], outname4[MAXLEN],
    tmpstring[256];

  if(p->isifcommands) {
    if(!TestIf(p->IfStack[lc2], p, c, lc, lc2)) {
      if(c->cnum != CNUM_SAVELC && c->cnum != CNUM_RESTORELC) {
	SkipCommand(p, c, thisindex, lc, lc2);
	return;
      }
    }
  }

  /* Sort the light curve in time, and merge unequal values if needed */
  if(c->require_sort || c->require_distinct) {
    i1 = sortlcbytime(p->NJD[lc2], p->t[lc2], lc2, p);
    if(c->require_distinct && i1) {
      mergeequallctimes(p, lc2);
    }
  }

  /* This routine runs a command on a single light curve */
  switch(c->cnum)
    {

    case CNUM_SORTLC:
      /* Sort the light curve by a vector */
      if(!c->SortLC->issortvar && !c->SortLC->isreverse) {
	i1 = sortlcbytime(p->NJD[lc2], p->t[lc2], lc2, p);
      }
      else if(!c->SortLC->issortvar && c->SortLC->isreverse) {
	i1 = sortlcbytime_rev(p->NJD[lc2], p->t[lc2], lc2, p);
      }
      else if(!c->SortLC->isreverse) {
	switch(c->SortLC->sortdtype) {
	case VARTOOLS_TYPE_DOUBLE:
	  i1 = sortlcbyvardbl(p->NJD[lc2], (*((double ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_INT:
	  i1 = sortlcbyvarint(p->NJD[lc2], (*((int ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_STRING:
	  i1 = sortlcbyvarstring(p->NJD[lc2], (*((char ****) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  i1 = sortlcbyvarfloat(p->NJD[lc2], (*((float ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  i1 = sortlcbyvarchar(p->NJD[lc2], (*((char ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_LONG:
	  i1 = sortlcbyvarlong(p->NJD[lc2], (*((long ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	default:
	  error(ERR_CODEERROR);
	  break;
	}
      } else {
	switch(c->SortLC->sortdtype) {
	case VARTOOLS_TYPE_DOUBLE:
	  i1 = sortlcbyvardbl_rev(p->NJD[lc2], (*((double ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_INT:
	  i1 = sortlcbyvarint_rev(p->NJD[lc2], (*((int ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_STRING:
	  i1 = sortlcbyvarstring_rev(p->NJD[lc2], (*((char ****) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  i1 = sortlcbyvarfloat_rev(p->NJD[lc2], (*((float ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  i1 = sortlcbyvarchar_rev(p->NJD[lc2], (*((char ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	case VARTOOLS_TYPE_LONG:
	  i1 = sortlcbyvarlong_rev(p->NJD[lc2], (*((long ***) (c->SortLC->sortvar->dataptr)))[lc2], lc2, p);
	  break;
	default:
	  error(ERR_CODEERROR);
	  break;
	}
      }
      break;

    case CNUM_DIFFFLUXTOMAG:
      /* Convert from isis differential flux to magnitudes */
      difffluxtomag(p->t[lc2],p->mag[lc2],p->sig[lc2],p->NJD[lc2],c->DiffFluxtomag->magstar[lc][lc2],c->DiffFluxtomag->mag_constant1, c->DiffFluxtomag->offset);
      break;

    case CNUM_FLUXTOMAG:
      /* Convert from isis differential flux to magnitudes */
      fluxtomag(p->t[lc2],p->mag[lc2],p->sig[lc2],p->NJD[lc2],c->Fluxtomag->mag_constant1, c->Fluxtomag->offset);
      break;

    case CNUM_EXPRESSION:
      /* Evaluate an analytic expression */
      RunExpressionCommand(lc, lc2, p, c->ExpressionCommand);
      break;

    case CNUM_LINFIT:
      /* Fit a model that is linear in its free parameters to the light curve */
      DoLinfit(p, c->Linfit, lc2, lc);
      break;

    case CNUM_NONLINFIT:
      /* Fit a model that is nonlinear in its free parameters to the light curve */
      DoNonlinfit(p, c->Nonlinfit, lc2, lc);
      break;

    case CNUM_WWZ:
      /* Run the WWZ transform */
      DoWWZ(p, c->WWZ, lc2, lc);
      break;

    case CNUM_FINDBLENDS:
      /* Find variability blends */
      /* Set the period and potential variable x, y coordinates*/
      if(c->FindBlends->pertype == PERTYPE_FIXCOLUMN)
	{
	  getoutcolumnvalue(c->FindBlends->linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &(c->FindBlends->periods[lc2][0]));
	}
      else if(c->FindBlends->pertype == PERTYPE_FIX)
	{
	  c->FindBlends->periods[lc2][0] = c->FindBlends->fixperiod;
	}
      getoutcolumnvalue(c->FindBlends->linkedcolumn_varname, lc2, lc, VARTOOLS_TYPE_STRING, &(c->FindBlends->varnames[lc2][0]), MAXLEN);
      c->FindBlends->varx[lc2] = c->FindBlends->varxyin[lc2][0];
      c->FindBlends->vary[lc2] = c->FindBlends->varxyin[lc2][1];
      findblends(1,p->NJD,p->t,p->mag,p->sig,c->FindBlends);
      break;

    case CNUM_CLIP:
      /* Clip the light curve */
      c->Clip->Nclip[lc2] = sigclip(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &d1, &d2, &d3, &i1, c->Clip->sigclip, c->Clip->iter, lc2, p, c->Clip->niter, c->Clip->usemedian, c->Clip->markclip, c->Clip->clipvar, c->Clip->noinitmark);
      break;

    case CNUM_CONVERTTIME:
      /* Perform a time conversion */
      converttime(p->NJD[lc2], p->t[lc2], lc2, lc, c->ConvertTime, p);
      break;

    case CNUM_ENSEMBLERESCALESIG:
      /* There must be a bug in the program if we're calling ensemblerescalesig in single light curve process mode! */
      error(ERR_CODEERROR);
      break;

    case CNUM_HARMONICFILTER:
      /* Apply a fourier filter to the light curve */
      doHarmonicFilter(p, c->HarmonicFilter, lc2, lc);
      break;

    case CNUM_PRINT:
      /* Run the -print command */
      RunPrintCommand(p, c->PrintCommand, lc2, lc);
      break;

    case CNUM_RESAMPLE:
      /* Resample the times of observation in the light curve */
      DoResample(p, c->Resample, lc2, lc);
      break;

    case CNUM_RESCALESIG:
      /* Rescale sigma for a light curve so that it has chi2 = 1 */
      /* First get the old chi2 value */
      c->Rescalesig->chi2_old[lc2] = chi2(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &d1, &i1, c->Rescalesig->usemask, c->Rescalesig->maskvar, lc, lc2);
      if(i1 > 1)
	{
	  c->Rescalesig->chi2_old[lc2] /= (double) (i1 - 1);
	}
      else
	c->Rescalesig->chi2_old[lc2] = -1.;
      /* Rescale sigma */
      rescalesigma_chi2(p->NJD[lc2], p->sig[lc2], c->Rescalesig->chi2_old[lc2]);
      /* Get the new chi2 value */
      c->Rescalesig->chi2_new[lc2] = chi2(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &d1, &i2, c->Rescalesig->usemask, c->Rescalesig->maskvar, lc, lc2);
      if(i2 > 1)
	{
	  c->Rescalesig->chi2_new[lc2] /= (double) (i2 - 1);
	}
      else
	c->Rescalesig->chi2_new[lc2] = -1.;
      /* Calculate the rescale factor */
      c->Rescalesig->rescalefactor[lc2] = sqrt(c->Rescalesig->chi2_new[lc2] / c->Rescalesig->chi2_old[lc2]);
      break;

    case CNUM_SAVELC:
      /* Save the light curve */
      dosavelc(p, c->Savelc, lc2, lc);
      break;

    case CNUM_COPYLC:
      /* Setup the light copies */
      docopylccommand(p, c->CopyLC, lc2, lc);
      break;

    case CNUM_RESTORELC:
      /* Restore the saved light curve */
      dorestorelc(p, c[c->Restorelc->saveindex - thisindex].Savelc, c->Restorelc, lc2, lc2, lc);
      break;

    case CNUM_RESTRICTTIMES:
      /* Filter points from the light curves based on the times of
	 observation */
      if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDRANGE) {
	if(c->RestrictTimes->minJDtype == PERTYPE_FIXCOLUMN) {
	  getoutcolumnvalue(c->RestrictTimes->minJD_linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &(c->RestrictTimes->minJD[lc2]));
	  d1 = c->RestrictTimes->minJD[lc2];
	}
	else if(c->RestrictTimes->minJDtype == PERTYPE_FIX) {
	  c->RestrictTimes->minJD[lc2] = c->RestrictTimes->minJDfixval;
	  d1 = c->RestrictTimes->minJD[lc2];
	}
	else if(c->RestrictTimes->minJDtype == PERTYPE_SPECIFIED) {
	  d1 = c->RestrictTimes->minJD[lc];
	}
	else if(c->RestrictTimes->minJDtype == PERTYPE_EXPR) {
	  d1 = EvaluateExpression(lc, lc2, 0, c->RestrictTimes->minJDexpr);
	  c->RestrictTimes->minJD[lc2] = d1;
	}
	if(c->RestrictTimes->maxJDtype == PERTYPE_FIXCOLUMN) {
	  getoutcolumnvalue(c->RestrictTimes->maxJD_linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &(c->RestrictTimes->maxJD[lc2]));
	  d2 = c->RestrictTimes->maxJD[lc2];
	}
	else if(c->RestrictTimes->maxJDtype == PERTYPE_FIX) {
	  c->RestrictTimes->maxJD[lc2] = c->RestrictTimes->maxJDfixval;
	  d2 = c->RestrictTimes->maxJD[lc2];
	}
	else if(c->RestrictTimes->maxJDtype == PERTYPE_SPECIFIED) {
	  d2 = c->RestrictTimes->maxJD[lc];
	}
	else if(c->RestrictTimes->maxJDtype == PERTYPE_EXPR) {
	  d2 = EvaluateExpression(lc, lc2, 0, c->RestrictTimes->maxJDexpr);
	  c->RestrictTimes->maxJD[lc2] = d2;
	}
	RestrictTimes_JDrange_apply(p->NJD[lc2], p->t[lc2], lc2, p, 
				    c->RestrictTimes, d1, d2,
				    c->RestrictTimes->exclude,
				    c->RestrictTimes->markrestrict,
				    c->RestrictTimes->markvar, 
				    c->RestrictTimes->noinitmark);
      }
      else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDLIST) {
	RestrictTimes_JDlist_apply(p->NJD[lc2], p->t[lc2], lc2, p,
				   c->RestrictTimes,
				   c->RestrictTimes->JD_restrictlist,
				   c->RestrictTimes->N_restrictlist,
				   c->RestrictTimes->exclude,
				   c->RestrictTimes->markrestrict,
				   c->RestrictTimes->markvar, 
				   c->RestrictTimes->noinitmark);
      }
      else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_IMAGELIST) {
	RestrictTimes_imagelist_apply(p->NJD[lc2], p->stringid[lc2], 
				      p->stringid_idx[lc2], lc2, p,
				      c->RestrictTimes,
				      c->RestrictTimes->image_restrictlist,
				      c->RestrictTimes->image_restrictlist_indx,
				      c->RestrictTimes->N_restrictlist,
				      c->RestrictTimes->exclude,
				      c->RestrictTimes->markrestrict,
				      c->RestrictTimes->markvar, 
				      c->RestrictTimes->noinitmark);
      }
      else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_EXPR) {
	RestrictTimes_expr_apply(p, c->RestrictTimes, lc2, lc,
				 c->RestrictTimes->markrestrict,
				 c->RestrictTimes->markvar, 
				 c->RestrictTimes->noinitmark);
      }
      break;

    case CNUM_RESTORETIMES:
      RestoreTimes(p, c->RestoreTimes, lc2, lc2);
      break;

    case CNUM_CHI2_NOBIN:
      /* calculate chi2 without binning */
      c->Chi2_NoBin->chi2val[lc2] = chi2(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &c->Chi2_NoBin->wtave[lc2], &i1, c->Chi2_NoBin->usemask, c->Chi2_NoBin->maskvar, lc, lc2);
      if(i1 > 1)
	c->Chi2_NoBin->chi2val[lc2] /= (double) (i1 - 1);
      else
	c->Chi2_NoBin->chi2val[lc2] = -1.;
      break;

    case CNUM_CHI2_BIN:
      /* calculate chi2 with binning */
      for(i=0; i < c->Chi2_Bin->Nbin ; i++)
	{
	  c->Chi2_Bin->chi2binvals[lc2][i] = binnedchi2(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->Chi2_Bin->bintimes[i], &c->Chi2_Bin->wtavebin[lc2][i], &i1, c->Chi2_Bin->usemask, c->Chi2_Bin->maskvar, lc, lc2);
	  if(i1 > 1)
	    c->Chi2_Bin->chi2binvals[lc2][i] /= (double) (i1 - 1);
	  else
	    c->Chi2_Bin->chi2binvals[lc2][i] = -1.;
	}
      break;

    case CNUM_CHANGEERROR:
      /* Replace the formal errors in a light curve with the RMS */
      c->Changeerror->rmsval[lc2] = changeerror(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &c->Changeerror->ave[lc2], &c->Changeerror->ngood[lc2], c->Changeerror->usemask, c->Changeerror->maskvar, lc, lc2);
      break;

    case CNUM_CHANGEVARIABLE:
      /* Switch the time, mag, sig, or string-id variable */
      DoChangeVariable(p, c->Changevariable, lc2);
      break;

    case CNUM_RMS_NOBIN:
      /* calculate RMS without binning */
      c->RMS_NoBin->rmsval[lc2] = rms(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &c->RMS_NoBin->ave[lc2], &c->RMS_NoBin->rmsthy[lc2], &c->RMS_NoBin->ngood[lc2], c->RMS_NoBin->usemask, c->RMS_NoBin->maskvar, lc, lc2);
      break;

    case CNUM_RMS_BIN:
      /* Calculate RMS with binning */
      for(i=0; i < c->RMS_Bin->Nbin ; i++)
	{
	  c->RMS_Bin->rmsbinvals[lc2][i] = binnedrms(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->RMS_Bin->bintimes[i], &d1, &c->RMS_Bin->rmsthybin[lc2][i], &i1, c->RMS_Bin->usemask, c->RMS_Bin->maskvar, lc, lc2);
	}
      break;

    case CNUM_JSTET:
      /* Calculate JSTET */
      getJstet(p->NJD[lc2], c->Jstet->Jstet_time, c->Jstet->wkmax, p->t[lc2], p->mag[lc2], p->sig[lc2], &d1, &c->Jstet->jst[lc2], &c->Jstet->kur[lc2], &c->Jstet->lst[lc2], lc2, lc, c->Jstet->usemask, c->Jstet->maskvar);
      break;

    case CNUM_ADDFITSKEYWORD:
      Run_AddFitsKeyword_Command(p, c->AddFitsKeyword, lc2, lc);
      break;

    case CNUM_ADDNOISE:
      /* Add time-correlated noise to the light curve */
      addnoise(p, c->AddNoise, lc2, lc);
      break;

    case CNUM_ALARM:
      /* Calculate the Alarm */
      c->Alarm->alarmvals[lc2] = doalarm(p->NJD[lc2], p->mag[lc2], p->sig[lc2],
					 lc2, lc, c->Alarm->usemask,
					 c->Alarm->maskvar);
      break;

    case CNUM_AUTOCORR:
      /* Calculate the auto-correlation */
      i1 = 0;
      i2 = 0;
      while(p->lcnames[lc][i1] != '\0')
	{
	  if(p->lcnames[lc][i1] == '/')
	    i2 = i1 + 1;
	  i1++;
	}
      sprintf(outname,"%s/%s%s",c->Autocorr->outdir,&p->lcnames[lc][i2],c->Autocorr->suffix);
      autocorrelation(p->t[lc2], p->mag[lc2], p->sig[lc2], p->NJD[lc2], c->Autocorr->start, c->Autocorr->stop, c->Autocorr->step, outname, lc2, lc, c->Autocorr->usemask, c->Autocorr->maskvar);
      break;


    case CNUM_AOV:
      /* Calculate the AoV with phase binning */
      RunAOVCommand(p, c, c->Aov, lc2, lc, thisindex);
      break;

    case CNUM_HARMAOV:
      /* Calculate the AoV with Harmonics */
      RunAOVHarmCommand(p, c, c->AovHarm, lc2, lc, thisindex);
      break;

    case CNUM_LS:
      /* Calculate the Lomb-Scargle Periodogram */
      RunLombScargleCommand(p, c->Ls, c, lc2, lc, thisindex);
      break;

#ifdef _HAVE_GSL
    case CNUM_FFT:
      /* Calculate light curve statistics */
      RunFFTCommand(p, lc, lc2, c->FFT);
      break;
#endif

    case CNUM_GETLSAMPTHRESH:
      /* Get the amplitude scale-factor for which the signal just passes */
      /* The LS threshhold */
      if(p->NJD[lc2] > 1)
	{
	  if(c->GetLSAmpThresh->pertype == PERTYPE_LS)
	    {
	      i1 = c->GetLSAmpThresh->lastlsindex;
	      c->GetLSAmpThresh->period[lc2][0] = c[i1-thisindex].Ls->peakperiods[lc2][0];
	    }
	  if(c->GetLSAmpThresh->harm_specsigflag)
	    {
	      if(gnu_getline(&(c->GetLSAmpThresh->line),&(c->GetLSAmpThresh->line_size),c->GetLSAmpThresh->listfile) < 0)
		{
		  error2(ERR_GETLSAMPTHRESH_FILETOSHORT,c->GetLSAmpThresh->listfilename);
		}
	      sscanf(c->GetLSAmpThresh->line,"%s %lf",c->GetLSAmpThresh->filename,c->GetLSAmpThresh->amp);
	      if((c->GetLSAmpThresh->infile = fopen(c->GetLSAmpThresh->filename,"r")) == NULL)
		error2(ERR_FILENOTFOUND,c->GetLSAmpThresh->filename);
	      c->GetLSAmpThresh->sizesigfile = 0;
	      while(gnu_getline(&(c->GetLSAmpThresh->line),&(c->GetLSAmpThresh->line_size),c->GetLSAmpThresh->infile) >= 0)
		c->GetLSAmpThresh->sizesigfile = c->GetLSAmpThresh->sizesigfile + 1;
	      if(c->GetLSAmpThresh->sizesigfile != p->NJD[lc2])
		error2(ERR_SIGFILEWRONGLENGTH,c->GetLSAmpThresh->filename);
	    }
	  if(c->GetLSAmpThresh->pertype != PERTYPE_SPECIFIED)
	    getlsampthresh(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->GetLSAmpThresh->period[lc2][0],c->GetLSAmpThresh->harm_specsigflag,c->GetLSAmpThresh->infile,c->GetLSAmpThresh->Nsubharm,c->GetLSAmpThresh->Nharm,c->GetLSAmpThresh->minPer,c->GetLSAmpThresh->thresh,&c->GetLSAmpThresh->ampthresh_scale[lc2],&c->GetLSAmpThresh->amp[lc2],c->GetLSAmpThresh->use_orig_ls);
	  else
	    getlsampthresh(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->GetLSAmpThresh->period[lc][0],c->GetLSAmpThresh->harm_specsigflag,c->GetLSAmpThresh->infile,c->GetLSAmpThresh->Nsubharm,c->GetLSAmpThresh->Nharm,c->GetLSAmpThresh->minPer,c->GetLSAmpThresh->thresh,&c->GetLSAmpThresh->ampthresh_scale[lc2],&c->GetLSAmpThresh->amp[lc2],c->GetLSAmpThresh->use_orig_ls);
	}
      break;


    case CNUM_DECORR:
      /* Decorrelate the light curves */
      /* First get the lc average */
      d1 = rms(p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], &d2, &d3, &i1, 0, NULL, lc, lc2);
      /* Get the output name if we're outputing the model light curve */
      if(c->Decorr->omodel)
	{
	  i3 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i3] != '\0')
	    {
	      if(p->lcnames[lc][i3] == '/')
		i2 = i3 + 1;
	      i3++;
	    }
	  sprintf(outname,"%s/%s%s",c->Decorr->modeloutdir,&p->lcnames[lc][i2],c->Decorr->modelsuffix);
	}
      /* Do the decorrelation only if there is at least 1 degree of freedom left over */
      if(i1 >= c->Decorr->N_decorrterms_total + 1)
	{
	  docorr(p->mag[lc2], p->sig[lc2], p->NJD[lc2], c->Decorr->N_decorrterms, c->Decorr->decorr_terms[lc2], c->Decorr->order, c->Decorr->b[lc2], c->Decorr->b_err[lc2], d2, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc2);

	  if(c->Decorr->correctlc)
	    {
	      magcorr((void *) (p->t[lc2]),VARTOOLS_TYPE_DOUBLE,p->mag[lc2], p->sig[lc2], p->NJD[lc2], c->Decorr->N_decorrterms, c->Decorr->decorr_terms[lc2], c->Decorr->order, c->Decorr->b[lc2], &c->Decorr->chi2val[lc2], &d3, d2, c->Decorr->omodel, outname, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc2);
	    }
	  else
	    magcorr_chi2only(p->t[lc2],p->mag[lc2],p->sig[lc2], p->NJD[lc2], c->Decorr->N_decorrterms,c->Decorr->decorr_terms[lc2],c->Decorr->order,c->Decorr->b[lc2],&c->Decorr->chi2val[lc2], &d3, d2, c->Decorr->omodel, outname, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc2);
	  c->Decorr->chi2val[lc2] /= (i1 - c->Decorr->N_decorrterms_total);
	}
      else
	{
	  for(i2=0;i2<c->Decorr->N_decorrterms_total;i2++)
	    {
	      c->Decorr->b[lc2][i2] = -1.;
	      c->Decorr->b_err[lc2][i2] = -1.;
	    }
	  c->Decorr->chi2val[lc2] = -1.;
	}
      break;

    case CNUM_KILLHARM:
      /* Remove harmonics */

      if(c->Killharm->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->Killharm->modeloutdir,&p->lcnames[lc][i2],c->Killharm->modelsuffix);
	}

      if(c->Killharm->pertype == PERTYPE_AOV)
	{
	  c->Killharm->Nper = 1;
	  i1=c->Killharm->lastaovindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    c->Killharm->periods[lc2][0] = c[i1-thisindex].Aov->peakperiods[lc2][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    c->Killharm->periods[lc2][0] = c[i1-thisindex].AovHarm->peakperiods[lc2][0];

	}
      else if(c->Killharm->pertype == PERTYPE_LS)
	{
	  c->Killharm->Nper = 1;
	  i1 = c->Killharm->lastlsindex;
	  c->Killharm->periods[lc2][0] = c[i1-thisindex].Ls->peakperiods[lc2][0];
	}
      else if(c->Killharm->pertype == PERTYPE_BOTH)
	{
	  c->Killharm->Nper = 2;
	  i1=c->Killharm->lastaovindex;
	  i2=c->Killharm->lastlsindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    c->Killharm->periods[lc2][0] = c[i1-thisindex].Aov->peakperiods[lc2][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    c->Killharm->periods[lc2][0] = c[i1-thisindex].AovHarm->peakperiods[lc2][0];
	  c->Killharm->periods[lc2][1] = c[i2-thisindex].Ls->peakperiods[lc2][0];
	}
      else if(c->Killharm->pertype == PERTYPE_INJECTHARM)
	{
	  c->Killharm->Nper = 1;
	  i1 = c->Killharm->lastaovindex;
	  c->Killharm->periods[lc2][0] = c[i1-thisindex].Injectharm->periodinject[lc2];
	}
      else if(c->Killharm->pertype == PERTYPE_FIX)
	{
	  for(i1 = 0; i1 < c->Killharm->Nper; i1++)
	    {
	      c->Killharm->periods[lc2][i1] = c->Killharm->fixedperiods[i1];
	    }
	}
      if(p->NJD[lc2] > 1)
	{
	  if(c->Killharm->pertype != PERTYPE_SPECIFIED)
	    dokillharms(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->Killharm->Nper,c->Killharm->periods[lc2],c->Killharm->Nsubharm,c->Killharm->Nharm,c->Killharm->subharmA[lc2],c->Killharm->subharmB[lc2],c->Killharm->harmA[lc2],c->Killharm->harmB[lc2],c->Killharm->fundA[lc2],c->Killharm->fundB[lc2],&c->Killharm->mean[lc2], c->Killharm->omodel, outname, c->Killharm->amp[lc2], c->Killharm->fitonly, c->Killharm->outtype, c->Killharm->clip);
	  else
	    dokillharms(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->Killharm->Nper,c->Killharm->periods[lc],c->Killharm->Nsubharm,c->Killharm->Nharm,c->Killharm->subharmA[lc2],c->Killharm->subharmB[lc2],c->Killharm->harmA[lc2],c->Killharm->harmB[lc2],c->Killharm->fundA[lc2],c->Killharm->fundB[lc2],&c->Killharm->mean[lc2], c->Killharm->omodel, outname,c->Killharm->amp[lc2], c->Killharm->fitonly, c->Killharm->outtype, c->Killharm->clip);
	}
      break;

    case CNUM_INJECTHARM:
      /* Inject a harmonic series */

      if(c->Injectharm->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->Injectharm->modeloutdir,&p->lcnames[lc][i2],c->Injectharm->modelsuffix);
	}
      doinjectharm(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],lc2,lc,c->Injectharm,outname);
      break;

    case CNUM_INJECTTRANSIT:
      /* Inject a transit model */
      if(c->Injecttransit->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->Injecttransit->modeloutdir,&p->lcnames[lc][i2],c->Injecttransit->modelsuffix);
	}
      doinjecttransit(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],lc2,lc,c->Injecttransit,outname);
      break;


    case CNUM_STARSPOT:
      /* Fit a starspot model to the light curve and remove it if we're doing that */
      if(c->Starspot->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->Starspot->modeloutdir,&p->lcnames[lc][i2],c->Starspot->modelsuffix);
	}
      if(c->Starspot->pertype == PERTYPE_AOV)
	{
	  i1 = c->Starspot->lastaovindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    c->Starspot->period[lc2][0] = c[i1-thisindex].Aov->peakperiods[lc2][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    c->Starspot->period[lc2][0] = c[i1-thisindex].AovHarm->peakperiods[lc2][0];
	}
      else if(c->Starspot->pertype == PERTYPE_LS)
	{
	  i1 = c->Starspot->lastlsindex;
	  c->Starspot->period[lc2][0] = c[i1-thisindex].Ls->peakperiods[lc2][0];
	}
      else if(c->Starspot->pertype == PERTYPE_FIXCOLUMN)
	{
	  getoutcolumnvalue(c->Starspot->linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &(c->Starspot->period[lc2][0]));
	}
      else if(c->Starspot->pertype == PERTYPE_FIX)
	{
	  c->Starspot->period[lc2][0] = c->Starspot->fixedperiod;
	}
      c->Starspot->a[lc2] = c->Starspot->a0;
      c->Starspot->b[lc2] = c->Starspot->b0;
      c->Starspot->alpha[lc2] = c->Starspot->alpha0;
      c->Starspot->inclination[lc2] = c->Starspot->inclination0;
      c->Starspot->chi[lc2] = c->Starspot->chi0;
      c->Starspot->psi0[lc2] = c->Starspot->psi00;
      c->Starspot->mconst[lc2] = c->Starspot->mconst0;
      if(p->NJD[lc2] > 1)
	{
	  if(c->Starspot->pertype != PERTYPE_SPECIFIED)
	    fitstarspot_amoeba(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],&c->Starspot->period[lc2][0],&c->Starspot->a[lc2],&c->Starspot->b[lc2],&c->Starspot->alpha[lc2],&c->Starspot->inclination[lc2],&c->Starspot->chi[lc2],&c->Starspot->psi0[lc2],&c->Starspot->mconst[lc2],c->Starspot->fitP,c->Starspot->fita,c->Starspot->fitb,c->Starspot->fitalpha,c->Starspot->fiti,c->Starspot->fitchi,c->Starspot->fitpsi0,c->Starspot->fitmconst,&c->Starspot->chisq[lc2],c->Starspot->correctlc,c->Starspot->omodel,outname);
	  else
	    fitstarspot_amoeba(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],&c->Starspot->period[lc][0],&c->Starspot->a[lc2],&c->Starspot->b[lc2],&c->Starspot->alpha[lc2],&c->Starspot->inclination[lc2],&c->Starspot->chi[lc2],&c->Starspot->psi0[lc2],&c->Starspot->mconst[lc2],c->Starspot->fitP,c->Starspot->fita,c->Starspot->fitb,c->Starspot->fitalpha,c->Starspot->fiti,c->Starspot->fitchi,c->Starspot->fitpsi0,c->Starspot->fitmconst,&c->Starspot->chisq[lc2],c->Starspot->correctlc,c->Starspot->omodel,outname);
	}
      break;

    case CNUM_STATS:
      /* Calculate light curve statistics */
      RunStatsCommand(p, lc, lc2, c->Stats);
      break;

    case CNUM_BLS:
      /* Perform BLS on the light curves */
      RunBLSCommand(p, c->Bls, lc2, lc, thisindex, lc2);
      break;
      /*if(p->NJD[lc2] > 1)
	{
	  if(c->Bls->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname2,"%s/%s%s",c->Bls->modeloutdir,&p->lcnames[lc][i2],c->Bls->modelsuffix);
	    }
	  if(c->Bls->ophcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname3,"%s/%s%s",c->Bls->ophcurveoutdir,&p->lcnames[lc][i2],c->Bls->ophcurvesuffix);
	    }
	  if(c->Bls->ojdcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname4,"%s/%s%s",c->Bls->ojdcurveoutdir,&p->lcnames[lc][i2],c->Bls->ojdcurvesuffix);
	    }*/
	  /* First check to see that the u/v vectors are large enough */
	  /*if(c->Bls->sizeuv[lc2] == 0)
	    {
	      c->Bls->sizeuv[lc2] = p->NJD[lc2];
	      if((c->Bls->u[lc2] = (double *) malloc(c->Bls->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->Bls->v[lc2] = (double *) malloc(c->Bls->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  else if(c->Bls->sizeuv[lc2] < p->NJD[lc2])
	    {
	      c->Bls->sizeuv[lc2] = p->NJD[lc2];
	      free(c->Bls->u[lc2]);
	      free(c->Bls->v[lc2]);
	      if((c->Bls->u[lc2] = (double *) malloc(c->Bls->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->Bls->v[lc2] = (double *) malloc(c->Bls->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }

	  if(c->Bls->operiodogram)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->Bls->outdir,&p->lcnames[lc][i2],c->Bls->suffix);
	    }
	  c->Bls->fmin[lc2] = dmax((1./(p->t[lc2][p->NJD[lc2]-1] - p->t[lc2][0])),1./c->Bls->maxper);
	  if(!c->Bls->freqsteptype) {
	    c->Bls->nf2[lc2] = floor((((1./c->Bls->minper) - c->Bls->fmin[lc2])/c->Bls->df)+1.);
	  } else if(c->Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_PERIOD) {
	    c->Bls->nf2[lc2] = floor((((1./c->Bls->fmin[lc2]) - c->Bls->minper)/c->Bls->df)+1.);
	  } else if(c->Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_LOGPERIOD) {
	    c->Bls->nf2[lc2] = floor(((log(1./c->Bls->fmin[lc2]) - log(c->Bls->minper))/c->Bls->df)+1.);
	  }*/
	  /* Now either run bls using the fixed q range or the fixed stellar radius range */
	  /*if(c->Bls->nf2[lc2] > 0 && c->Bls->nbins > 0 && c->Bls->Npeak > 0) {
	    if(!c->Bls->rflag)
	      {
		eebls(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->Bls->u[lc2],c->Bls->v[lc2],c->Bls->nf2[lc2],c->Bls->fmin[lc2],c->Bls->df,c->Bls->nbins,c->Bls->qmin,c->Bls->qmax,
#ifdef PARALLEL
		      c->Bls->p[lc2]
#else
		      c->Bls->p
#endif
		      ,c->Bls->Npeak,c->Bls->bper[lc2],c->Bls->bt0[lc2],c->Bls->bpow[lc2],c->Bls->sde[lc2],c->Bls->snval[lc2],c->Bls->depth[lc2],c->Bls->qtran[lc2],c->Bls->i1[lc2],c->Bls->i2[lc2],c->Bls->i1_ph[lc2],c->Bls->i2_ph[lc2],c->Bls->chisqrplus[lc2],&c->Bls->chisqrminus[lc2],&c->Bls->bperpos[lc2],&c->Bls->meanmagval[lc2], c->Bls->timezone, c->Bls->fraconenight[lc2], c->Bls->operiodogram, outname, c->Bls->omodel, outname2, c->Bls->correctlc,p->ascii, c->Bls->nt[lc2], c->Bls->Nt[lc2], c->Bls->Nbefore[lc2], c->Bls->Nafter[lc2], c->Bls->rednoise[lc2], c->Bls->whitenoise[lc2], c->Bls->sigtopink[lc2], c->Bls->fittrap, c->Bls->qingress[lc2], c->Bls->OOTmag[lc2], c->Bls->ophcurve, outname3, c->Bls->phmin, c->Bls->phmax, c->Bls->phstep, c->Bls->ojdcurve, outname4, c->Bls->jdstep, c->Bls->nobinnedrms, c->Bls->freqsteptype, c->Bls->adjust_qmin_mindt, c->Bls->reduce_nb, c->Bls->reportharmonics, c->Bls, lc2, lc, c->Bls->usemask, c->Bls->maskvar);
	      }
	    else
	      {
		eebls_rad(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->Bls->u[lc2],c->Bls->v[lc2],c->Bls->nf2[lc2],c->Bls->fmin[lc2],c->Bls->df,c->Bls->nbins,c->Bls->rmin,c->Bls->rmax,
#ifdef PARALLEL
			  c->Bls->p[lc2]
#else
			  c->Bls->p
#endif
			  ,c->Bls->Npeak,c->Bls->bper[lc2],c->Bls->bt0[lc2],c->Bls->bpow[lc2],c->Bls->sde[lc2],c->Bls->snval[lc2],c->Bls->depth[lc2],c->Bls->qtran[lc2],c->Bls->i1[lc2],c->Bls->i2[lc2],c->Bls->i1_ph[lc2],c->Bls->i2_ph[lc2],c->Bls->chisqrplus[lc2],&c->Bls->chisqrminus[lc2],&c->Bls->bperpos[lc2],&c->Bls->meanmagval[lc2], c->Bls->timezone, c->Bls->fraconenight[lc2], c->Bls->operiodogram,outname, c->Bls->omodel, outname2, c->Bls->correctlc,p->ascii, c->Bls->nt[lc2], c->Bls->Nt[lc2], c->Bls->Nbefore[lc2], c->Bls->Nafter[lc2], c->Bls->rednoise[lc2], c->Bls->whitenoise[lc2], c->Bls->sigtopink[lc2], c->Bls->fittrap, c->Bls->qingress[lc2], c->Bls->OOTmag[lc2], c->Bls->ophcurve, outname3, c->Bls->phmin, c->Bls->phmax, c->Bls->phstep, c->Bls->ojdcurve, outname4, c->Bls->jdstep, c->Bls->nobinnedrms,c->Bls->freqsteptype, c->Bls->adjust_qmin_mindt, c->Bls->reduce_nb, c->Bls->reportharmonics, c->Bls, lc2, lc, c->Bls->usemask, c->Bls->maskvar);
	      }
	  } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve is either too short, or an invalid set of parameter options were supplied to BLS.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	} else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc, p->lcnames[lc]);
	    }
      }
      break;*/

    case CNUM_FIXPERBLS:
      /* Perform BLS on the light curves */
      if(c->BlsFixPer->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname2,"%s/%s%s",c->BlsFixPer->modeloutdir,&p->lcnames[lc][i2],c->BlsFixPer->modelsuffix);
	}
      if(p->NJD[lc2] > 1) {
	/* First check to see that the u/v vectors are large enough */
	if(c->BlsFixPer->sizeuv[lc2] == 0)
	  {
	    c->BlsFixPer->sizeuv[lc2] = p->NJD[lc2];
	    if((c->BlsFixPer->u[lc2] = (double *) malloc(c->BlsFixPer->sizeuv[lc2] * sizeof(double))) == NULL ||
	       (c->BlsFixPer->v[lc2] = (double *) malloc(c->BlsFixPer->sizeuv[lc2] * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	else if(c->BlsFixPer->sizeuv[lc2] < p->NJD[lc2])
	  {
	    c->BlsFixPer->sizeuv[lc2] = p->NJD[lc2];
	    free(c->BlsFixPer->u[lc2]);
	    free(c->BlsFixPer->v[lc2]);
	    if((c->BlsFixPer->u[lc2] = (double *) malloc(c->BlsFixPer->sizeuv[lc2] * sizeof(double))) == NULL ||
	       (c->BlsFixPer->v[lc2] = (double *) malloc(c->BlsFixPer->sizeuv[lc2] * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	
	/* Find the period if we're getting it from a previous command */
	if(c->BlsFixPer->pertype == PERTYPE_AOV)
	  {
	    i1 = c->BlsFixPer->lastaovindex;
	    if(c[i1-thisindex].cnum == CNUM_AOV)
	      c->BlsFixPer->period[lc2][0] = c[i1-thisindex].Aov->peakperiods[lc2][0];
	    else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	      c->BlsFixPer->period[lc2][0] = c[i1-thisindex].AovHarm->peakperiods[lc2][0];
	  }
	else if(c->BlsFixPer->pertype == PERTYPE_LS)
	  {
	    i1 = c->BlsFixPer->lastlsindex;
	    c->BlsFixPer->period[lc2][0] = c[i1-thisindex].Ls->peakperiods[lc2][0];
	  }
	else if(c->BlsFixPer->pertype == PERTYPE_FIX)
	  {
	    c->BlsFixPer->period[lc2][0] = c->BlsFixPer->perfix;
	  }
	else if(c->BlsFixPer->pertype == PERTYPE_FIXCOLUMN)
	  {
	    getoutcolumnvalue(c->BlsFixPer->linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &(c->BlsFixPer->period[lc2][0]));
	  }
	else if(c->BlsFixPer->pertype == PERTYPE_EXPR)
	  {
	    c->BlsFixPer->period[lc2][0] = EvaluateExpression(lc, lc2, 0, c->BlsFixPer->perexpr);
	  }


	if(p->NJD[lc2] > 1)
	  {
	    /* Now either run bls using the fixed q range or the fixed stellar radius range */
	    if(!c->BlsFixPer->rflag)
	      {
		if(c->BlsFixPer->pertype != PERTYPE_SPECIFIED)
		  eeblsfixper(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixPer->u[lc2],c->BlsFixPer->v[lc2],c->BlsFixPer->nbins,c->BlsFixPer->qmin,c->BlsFixPer->qmax,&c->BlsFixPer->period[lc2][0],&c->BlsFixPer->bt0[lc2],&c->BlsFixPer->bpow[lc2],&c->BlsFixPer->depth[lc2],&c->BlsFixPer->qtran[lc2],&c->BlsFixPer->i1[lc2],&c->BlsFixPer->i2[lc2],&c->BlsFixPer->i1_ph[lc2],&c->BlsFixPer->i2_ph[lc2],&c->BlsFixPer->chisqrplus[lc2],&c->BlsFixPer->chisqrminus[lc2],&c->BlsFixPer->meanmagval[lc2], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc2], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc2], &c->BlsFixPer->Nt[lc2], &c->BlsFixPer->Nbefore[lc2], &c->BlsFixPer->Nafter[lc2], &c->BlsFixPer->rednoise[lc2], &c->BlsFixPer->whitenoise[lc2], &c->BlsFixPer->sigtopink[lc2], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc2], &c->BlsFixPer->OOTmag[lc2], NULL, lc2, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
		else
		  eeblsfixper(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixPer->u[lc2],c->BlsFixPer->v[lc2],c->BlsFixPer->nbins,c->BlsFixPer->qmin,c->BlsFixPer->qmax,&c->BlsFixPer->period[lc][0],&c->BlsFixPer->bt0[lc2],&c->BlsFixPer->bpow[lc2],&c->BlsFixPer->depth[lc2],&c->BlsFixPer->qtran[lc2],&c->BlsFixPer->i1[lc2],&c->BlsFixPer->i2[lc2],&c->BlsFixPer->i1_ph[lc2],&c->BlsFixPer->i2_ph[lc2],&c->BlsFixPer->chisqrplus[lc2],&c->BlsFixPer->chisqrminus[lc2],&c->BlsFixPer->meanmagval[lc2], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc2], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc2], &c->BlsFixPer->Nt[lc2], &c->BlsFixPer->Nbefore[lc2], &c->BlsFixPer->Nafter[lc2], &c->BlsFixPer->rednoise[lc2], &c->BlsFixPer->whitenoise[lc2], &c->BlsFixPer->sigtopink[lc2], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc2], &c->BlsFixPer->OOTmag[lc2], NULL, lc2, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
	      }
	    else
	      {
		if(c->BlsFixPer->pertype != PERTYPE_SPECIFIED)
		  eeblsfixper_rad(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixPer->u[lc2],c->BlsFixPer->v[lc2],c->BlsFixPer->nbins,c->BlsFixPer->rmin,c->BlsFixPer->rmax,&c->BlsFixPer->period[lc2][0],&c->BlsFixPer->bt0[lc2],&c->BlsFixPer->bpow[lc2],&c->BlsFixPer->depth[lc2],&c->BlsFixPer->qtran[lc2],&c->BlsFixPer->i1[lc2],&c->BlsFixPer->i2[lc2],&c->BlsFixPer->i1_ph[lc2],&c->BlsFixPer->i2_ph[lc2],&c->BlsFixPer->chisqrplus[lc2],&c->BlsFixPer->chisqrminus[lc2],&c->BlsFixPer->meanmagval[lc2], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc2], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc2], &c->BlsFixPer->Nt[lc2], &c->BlsFixPer->Nbefore[lc2], &c->BlsFixPer->Nafter[lc2], &c->BlsFixPer->rednoise[lc2], &c->BlsFixPer->whitenoise[lc2], &c->BlsFixPer->sigtopink[lc2], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc2], &c->BlsFixPer->OOTmag[lc2], NULL, lc2, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
		else
		  eeblsfixper_rad(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixPer->u[lc2],c->BlsFixPer->v[lc2],c->BlsFixPer->nbins,c->BlsFixPer->rmin,c->BlsFixPer->rmax,&c->BlsFixPer->period[lc][0],&c->BlsFixPer->bt0[lc2],&c->BlsFixPer->bpow[lc2],&c->BlsFixPer->depth[lc2],&c->BlsFixPer->qtran[lc2],&c->BlsFixPer->i1[lc2],&c->BlsFixPer->i2[lc2],&c->BlsFixPer->i1_ph[lc2],&c->BlsFixPer->i2_ph[lc2],&c->BlsFixPer->chisqrplus[lc2],&c->BlsFixPer->chisqrminus[lc2],&c->BlsFixPer->meanmagval[lc2], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc2], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc2], &c->BlsFixPer->Nt[lc2], &c->BlsFixPer->Nbefore[lc2], &c->BlsFixPer->Nafter[lc2], &c->BlsFixPer->rednoise[lc2], &c->BlsFixPer->whitenoise[lc2], &c->BlsFixPer->sigtopink[lc2], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc2], &c->BlsFixPer->OOTmag[lc2], NULL, lc2, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
	      }
	  }
      } else {
	if(!p->quiet_mode) {
	  fprintf(stderr,"Warning: skipping -BLSFixPer command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc, p->lcnames[lc]);
	}
      }
      break;

    case CNUM_BLSFIXDURTC:
      /* Perform BLS with fixed transit duration and epoch on the light curves */
      if(p->NJD[lc2] > 1)
	{
	  if(c->BlsFixDurTc->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname2,"%s/%s%s",c->BlsFixDurTc->modeloutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->modelsuffix);
	    }
	  if(c->BlsFixDurTc->ophcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname3,"%s/%s%s",c->BlsFixDurTc->ophcurveoutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->ophcurvesuffix);
	    }
	  if(c->BlsFixDurTc->ojdcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname4,"%s/%s%s",c->BlsFixDurTc->ojdcurveoutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->ojdcurvesuffix);
	    }
	  /* First check to see that the u/v vectors are large enough */
	  if(c->BlsFixDurTc->sizeuv[lc2] == 0)
	    {
	      c->BlsFixDurTc->sizeuv[lc2] = p->NJD[lc2];
	      if((c->BlsFixDurTc->u[lc2] = (double *) malloc(c->BlsFixDurTc->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->BlsFixDurTc->v[lc2] = (double *) malloc(c->BlsFixDurTc->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  else if(c->BlsFixDurTc->sizeuv[lc2] < p->NJD[lc2])
	    {
	      c->BlsFixDurTc->sizeuv[lc2] = p->NJD[lc2];
	      free(c->BlsFixDurTc->u[lc2]);
	      free(c->BlsFixDurTc->v[lc2]);
	      if((c->BlsFixDurTc->u[lc2] = (double *) malloc(c->BlsFixDurTc->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->BlsFixDurTc->v[lc2] = (double *) malloc(c->BlsFixDurTc->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  
	  if(c->BlsFixDurTc->operiodogram)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->BlsFixDurTc->outdir,&p->lcnames[lc][i2],c->BlsFixDurTc->suffix);
	    }
	  if(c->BlsFixDurTc->durtype == PERTYPE_FIX)
	    {
	      c->BlsFixDurTc->inputdur[lc2] = c->BlsFixDurTc->fixdur;
	      d1 = c->BlsFixDurTc->inputdur[lc2];
	    }
	  else if(c->BlsFixDurTc->durtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->BlsFixDurTc->fixdur_linkedcolumn, lc2, lc, 
			      VARTOOLS_TYPE_DOUBLE, 
			      &(c->BlsFixDurTc->inputdur[lc2]));
	    d1 = c->BlsFixDurTc->inputdur[lc2];
	  } else {
	    d1 = c->BlsFixDurTc->inputdur[lc];
	  }
	  if(c->BlsFixDurTc->TCtype == PERTYPE_FIX)
	    {
	      c->BlsFixDurTc->inputTC[lc2] = c->BlsFixDurTc->fixTC;
	      d2 = c->BlsFixDurTc->fixTC;
	    }
	  else if(c->BlsFixDurTc->TCtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->BlsFixDurTc->fixTC_linkedcolumn, lc2, lc, 
			      VARTOOLS_TYPE_DOUBLE, 
			      &(c->BlsFixDurTc->inputTC[lc2]));
	    d2 = c->BlsFixDurTc->inputTC[lc2];
	  }
	  else {
	    d2 = c->BlsFixDurTc->inputTC[lc];
	  }
	  if(c->BlsFixDurTc->fixdepth) {
	    if(c->BlsFixDurTc->depthtype == PERTYPE_FIX) {
	      c->BlsFixDurTc->inputdepth[lc2] = c->BlsFixDurTc->fixdepthval;
	      d3 = c->BlsFixDurTc->fixdepthval;
	    }
	    else if(c->BlsFixDurTc->depthtype == PERTYPE_FIXCOLUMN) {
	      getoutcolumnvalue(c->BlsFixDurTc->fixdepth_linkedcolumn, lc2, lc, 
				VARTOOLS_TYPE_DOUBLE, 
				&(c->BlsFixDurTc->inputdepth[lc2]));
	      d3 = c->BlsFixDurTc->inputdepth[lc2];
	    }
	    else {
	      d3 = c->BlsFixDurTc->inputdepth[lc];
	    }
	    if(c->BlsFixDurTc->qgresstype == PERTYPE_FIX) {
	      c->BlsFixDurTc->inputqgress[lc2] = c->BlsFixDurTc->qgressval;
	      d4 = c->BlsFixDurTc->qgressval;
	    }
	    else if(c->BlsFixDurTc->qgresstype == PERTYPE_FIXCOLUMN) {
	      getoutcolumnvalue(c->BlsFixDurTc->fixqgress_linkedcolumn, lc2, lc, 
				VARTOOLS_TYPE_DOUBLE, 
				&(c->BlsFixDurTc->inputqgress[lc2]));
	      d4 = c->BlsFixDurTc->inputqgress[lc2];
	    }
	    else {
	      d4 = c->BlsFixDurTc->inputqgress[lc];
	    }
	  }
	  c->BlsFixDurTc->fmin[lc2] = dmax((2./(p->t[lc2][p->NJD[lc2]-1] - p->t[lc2][0])),1./c->BlsFixDurTc->maxper);
	  c->BlsFixDurTc->nf2[lc2] = floor((((1./c->BlsFixDurTc->minper) - c->BlsFixDurTc->fmin[lc2])/c->BlsFixDurTc->df)+1.);
	  if(c->BlsFixDurTc->nf2[lc2] > 0 && c->BlsFixDurTc->Npeak > 0) {
	    /* Now either run bls  */
	    eeblsfixdurtc(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixDurTc->u[lc2],c->BlsFixDurTc->v[lc2],d2,d1,c->BlsFixDurTc->fixdepth,d3,d4,c->BlsFixDurTc->nf2[lc2],c->BlsFixDurTc->fmin[lc2],c->BlsFixDurTc->df,
#ifdef PARALLEL
			  c->BlsFixDurTc->p[lc2]
#else
			  c->BlsFixDurTc->p
#endif
			  ,c->BlsFixDurTc->Npeak,c->BlsFixDurTc->bper[lc2],c->BlsFixDurTc->bt0[lc2],c->BlsFixDurTc->bpow[lc2],c->BlsFixDurTc->sde[lc2],c->BlsFixDurTc->snval[lc2],c->BlsFixDurTc->depth[lc2],c->BlsFixDurTc->qtran[lc2],c->BlsFixDurTc->chisqrplus[lc2],&c->BlsFixDurTc->chisqrminus[lc2],&c->BlsFixDurTc->bperpos[lc2],&c->BlsFixDurTc->meanmagval[lc2], c->BlsFixDurTc->timezone, c->BlsFixDurTc->fraconenight[lc2], c->BlsFixDurTc->operiodogram, outname, c->BlsFixDurTc->omodel, outname2, c->BlsFixDurTc->correctlc,p->ascii, c->BlsFixDurTc->nt[lc2], c->BlsFixDurTc->Nt[lc2], c->BlsFixDurTc->Nbefore[lc2], c->BlsFixDurTc->Nafter[lc2], c->BlsFixDurTc->rednoise[lc2], c->BlsFixDurTc->whitenoise[lc2], c->BlsFixDurTc->sigtopink[lc2], c->BlsFixDurTc->fittrap, c->BlsFixDurTc->qingress[lc2], c->BlsFixDurTc->OOTmag[lc2], c->BlsFixDurTc->ophcurve, outname3, c->BlsFixDurTc->phmin, c->BlsFixDurTc->phmax, c->BlsFixDurTc->phstep, c->BlsFixDurTc->ojdcurve, outname4, c->BlsFixDurTc->jdstep, lc2, lc, c->BlsFixDurTc->usemask, c->BlsFixDurTc->maskvar);
	  } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLSFixDurTc command index %d for light curve number: %d, filename: %s. The light curve is either too short, or an invalid set of parameter options were supplied to BLS.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	} else {
	if(!p->quiet_mode) {
	  fprintf(stderr,"Warning: skipping -BLSFixDurTc command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc, p->lcnames[lc]);
	}
      }
      break;

    case CNUM_BLSFIXPERDURTC:
      /* Perform BLS with fixed period, transit duration and epoch on the light curves */
      if(p->NJD[lc2] > 1)
	{
	  if(c->BlsFixPerDurTc->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname2,"%s/%s%s",c->BlsFixPerDurTc->modeloutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->modelsuffix);
	    }
	  if(c->BlsFixPerDurTc->ophcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname3,"%s/%s%s",c->BlsFixPerDurTc->ophcurveoutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->ophcurvesuffix);
	    }
	  if(c->BlsFixPerDurTc->ojdcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname4,"%s/%s%s",c->BlsFixPerDurTc->ojdcurveoutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->ojdcurvesuffix);
	    }
	  /* First check to see that the u/v vectors are large enough */
	  if(c->BlsFixPerDurTc->sizeuv[lc2] == 0)
	    {
	      c->BlsFixPerDurTc->sizeuv[lc2] = p->NJD[lc2];
	      if((c->BlsFixPerDurTc->u[lc2] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->BlsFixPerDurTc->v[lc2] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  else if(c->BlsFixPerDurTc->sizeuv[lc2] < p->NJD[lc2])
	    {
	      c->BlsFixPerDurTc->sizeuv[lc2] = p->NJD[lc2];
	      free(c->BlsFixPerDurTc->u[lc2]);
	      free(c->BlsFixPerDurTc->v[lc2]);
	      if((c->BlsFixPerDurTc->u[lc2] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[lc2] * sizeof(double))) == NULL ||
		 (c->BlsFixPerDurTc->v[lc2] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[lc2] * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  
	  if(c->BlsFixPerDurTc->pertype == PERTYPE_FIX)
	    {
	      c->BlsFixPerDurTc->inputper[lc2] = c->BlsFixPerDurTc->fixper;
	      d0 = c->BlsFixPerDurTc->inputper[lc2];
	    }
	  else if(c->BlsFixPerDurTc->pertype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->BlsFixPerDurTc->fixper_linkedcolumn, lc2, lc, 
			      VARTOOLS_TYPE_DOUBLE, 
			      &(c->BlsFixPerDurTc->inputper[lc2]));
	    d0 = c->BlsFixPerDurTc->inputper[lc2];
	  } else {
	    d0 = c->BlsFixPerDurTc->inputper[lc];
	  }
	  if(c->BlsFixPerDurTc->durtype == PERTYPE_FIX)
	    {
	      c->BlsFixPerDurTc->inputdur[lc2] = c->BlsFixPerDurTc->fixdur;
	      d1 = c->BlsFixPerDurTc->inputdur[lc2];
	    }
	  else if(c->BlsFixPerDurTc->durtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->BlsFixPerDurTc->fixdur_linkedcolumn, lc2, lc, 
			      VARTOOLS_TYPE_DOUBLE, 
			      &(c->BlsFixPerDurTc->inputdur[lc2]));
	    d1 = c->BlsFixPerDurTc->inputdur[lc2];
	  } else {
	    d1 = c->BlsFixPerDurTc->inputdur[lc];
	  }
	  if(c->BlsFixPerDurTc->TCtype == PERTYPE_FIX)
	    {
	      c->BlsFixPerDurTc->inputTC[lc2] = c->BlsFixPerDurTc->fixTC;
	      d2 = c->BlsFixPerDurTc->fixTC;
	    }
	  else if(c->BlsFixPerDurTc->TCtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->BlsFixPerDurTc->fixTC_linkedcolumn, lc2, lc, 
			      VARTOOLS_TYPE_DOUBLE, 
			      &(c->BlsFixPerDurTc->inputTC[lc2]));
	    d2 = c->BlsFixPerDurTc->inputTC[lc2];
	  }
	  else {
	    d2 = c->BlsFixPerDurTc->inputTC[lc];
	  }
	  if(c->BlsFixPerDurTc->fixdepth) {
	    if(c->BlsFixPerDurTc->depthtype == PERTYPE_FIX) {
	      c->BlsFixPerDurTc->inputdepth[lc2] = c->BlsFixPerDurTc->fixdepthval;
	      d3 = c->BlsFixPerDurTc->fixdepthval;
	    }
	    else if(c->BlsFixPerDurTc->depthtype == PERTYPE_FIXCOLUMN) {
	      getoutcolumnvalue(c->BlsFixPerDurTc->fixdepth_linkedcolumn, lc2, lc, 
				VARTOOLS_TYPE_DOUBLE, 
				&(c->BlsFixPerDurTc->inputdepth[lc2]));
	      d3 = c->BlsFixPerDurTc->inputdepth[lc2];
	    }
	    else {
	      d3 = c->BlsFixPerDurTc->inputdepth[lc];
	    }
	    if(c->BlsFixPerDurTc->qgresstype == PERTYPE_FIX) {
	      c->BlsFixPerDurTc->inputqgress[lc2] = c->BlsFixPerDurTc->qgressval;
	      d4 = c->BlsFixPerDurTc->qgressval;
	    }
	    else if(c->BlsFixPerDurTc->qgresstype == PERTYPE_FIXCOLUMN) {
	      getoutcolumnvalue(c->BlsFixPerDurTc->fixqgress_linkedcolumn, lc2, lc, 
				VARTOOLS_TYPE_DOUBLE, 
				&(c->BlsFixPerDurTc->inputqgress[lc2]));
	      d4 = c->BlsFixPerDurTc->inputqgress[lc2];
	    }
	    else {
	      d4 = c->BlsFixPerDurTc->inputqgress[lc];
	    }
	  }
	  eeblsfixperdurtc(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->BlsFixPerDurTc->u[lc2],c->BlsFixPerDurTc->v[lc2],d0,d2,d1,c->BlsFixPerDurTc->fixdepth,d3,d4,&(c->BlsFixPerDurTc->depth[lc2]),&(c->BlsFixPerDurTc->qtran[lc2]),&(c->BlsFixPerDurTc->chisqrplus[lc2]),&(c->BlsFixPerDurTc->meanmagval[lc2]), c->BlsFixPerDurTc->timezone, &(c->BlsFixPerDurTc->fraconenight[lc2]), c->BlsFixPerDurTc->omodel, outname2, c->BlsFixPerDurTc->correctlc, &(c->BlsFixPerDurTc->nt[lc2]), &(c->BlsFixPerDurTc->Nt[lc2]), &(c->BlsFixPerDurTc->Nbefore[lc2]), &(c->BlsFixPerDurTc->Nafter[lc2]), &(c->BlsFixPerDurTc->rednoise[lc2]), &(c->BlsFixPerDurTc->whitenoise[lc2]), &(c->BlsFixPerDurTc->sigtopink[lc2]), c->BlsFixPerDurTc->fittrap, &(c->BlsFixPerDurTc->qingress[lc2]), &(c->BlsFixPerDurTc->OOTmag[lc2]), c->BlsFixPerDurTc->ophcurve, outname3, c->BlsFixPerDurTc->phmin, c->BlsFixPerDurTc->phmax, c->BlsFixPerDurTc->phstep, c->BlsFixPerDurTc->ojdcurve, outname4, c->BlsFixPerDurTc->jdstep, lc2, lc, c->BlsFixPerDurTc->usemask, c->BlsFixPerDurTc->maskvar);
	} else {
	if(!p->quiet_mode) {
	  fprintf(stderr,"Warning: skipping -BLSFixPerDurTc command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc, p->lcnames[lc]);
	}
      }
      break;

    case CNUM_SOFTENEDTRANSIT:
      /* Fit a Softened Transit model to the light curve and remove it if we're doing that */
      if(c->SoftenedTransit->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->SoftenedTransit->modeloutdir,&p->lcnames[lc][i2],c->SoftenedTransit->modelsuffix);
	}
      if(c->SoftenedTransit->frombls)
	{
	  /* Get the starting parameters from BLS if that's what we're doing */
	  i1 = c->SoftenedTransit->lastblsindex;
	  c->SoftenedTransit->period0 = c[i1 - thisindex].Bls->bper[lc2][0];
	  c->SoftenedTransit->T00 = c[i1 - thisindex].Bls->bt0[lc2][0];
	  /*
	  if(c[i1 - thisindex].Bls->i2[lc2][0] > c[i1 - thisindex].Bls->i1[lc2][0])
	    c->SoftenedTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].Bls->i1[lc2][0] + c[i1 - thisindex].Bls->i2[lc2][0])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].Bls->nbins);
	  else
	    c->SoftenedTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].Bls->i1[lc2][0] + 1 - c[i1 - thisindex].Bls->i2[lc2][0])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].Bls->nbins);
	  */
	  c->SoftenedTransit->eta0 = c[i1 - thisindex].Bls->qtran[lc2][0];
	  c->SoftenedTransit->delta0 = c[i1 - thisindex].Bls->depth[lc2][0];
	  c->SoftenedTransit->mconst0 = -1;
	}
      else if(c->SoftenedTransit->fromblsfixper)
	{
	  /* Get the starting parameters from BLSFixPer if that's what we're doing */
	  i1 = c->SoftenedTransit->lastblsfixperindex;
	  c->SoftenedTransit->period0 = c[i1 - thisindex].BlsFixPer->period[lc2][0];
	  c->SoftenedTransit->T00 = c[i1 - thisindex].BlsFixPer->bt0[lc2];

	  /*
	  if(c[i1 - thisindex].BlsFixPer->i2[lc2] > c[i1 - thisindex].BlsFixPer->i1[lc2])
	    c->SoftenedTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc2] + c[i1 - thisindex].BlsFixPer->i2[lc2])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	  else
	    c->SoftenedTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc2] + 1 - c[i1 - thisindex].BlsFixPer->i2[lc2])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	  */
	  c->SoftenedTransit->eta0 = c[i1 - thisindex].BlsFixPer->qtran[lc2];
	  c->SoftenedTransit->delta0 = c[i1 - thisindex].BlsFixPer->depth[lc2];
	  c->SoftenedTransit->mconst0 = -1;
	}
      if(c->SoftenedTransit->dokillharm)
	{
	  /* Get the parameters for removing the harmonic function */
	  if(c->SoftenedTransit->pertype == PERTYPE_AOV)
	    {
	      i1 = c->SoftenedTransit->lastaovindex;
	      if(c[i1-thisindex].cnum == CNUM_AOV)
		c->SoftenedTransit->per_harm_out[lc2] = c[i1-thisindex].Aov->peakperiods[lc2][0];
	      else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		c->SoftenedTransit->per_harm_out[lc2] = c[i1 - thisindex].AovHarm->peakperiods[lc2][0];
	    }
	  else if(c->SoftenedTransit->pertype == PERTYPE_LS)
	    {
	      i1 = c->SoftenedTransit->lastlsindex;
	      c->SoftenedTransit->per_harm_out[lc2] = c[i1 - thisindex].Ls->peakperiods[lc2][0];
	    }
	  else if(c->SoftenedTransit->pertype == PERTYPE_BLS)
	    {
	      i1 = c->SoftenedTransit->lastblsindex;
	      c->SoftenedTransit->per_harm_out[lc2] = c[i1 - thisindex].Bls->bper[lc2][0];
	    }
	  else if(c->SoftenedTransit->pertype == PERTYPE_FIX)
	    {
	      c->SoftenedTransit->per_harm_out[lc2] = c->SoftenedTransit->per_harm;
	    }
	  else if(c->SoftenedTransit->pertype == PERTYPE_SPECIFIED)
	    {
	      c->SoftenedTransit->per_harm_out[lc2] = c->SoftenedTransit->per_harm_spec[lc];
	    }
	}
      c->SoftenedTransit->period[lc2] = c->SoftenedTransit->period0;
      c->SoftenedTransit->T0[lc2] = c->SoftenedTransit->T00;
      c->SoftenedTransit->eta[lc2] = c->SoftenedTransit->eta0;
      c->SoftenedTransit->cval[lc2] = c->SoftenedTransit->cval0;
      c->SoftenedTransit->delta[lc2] = c->SoftenedTransit->delta0;
      c->SoftenedTransit->mconst[lc2] = c->SoftenedTransit->mconst0;
      if(p->NJD[lc2] > 1)
	{
	  if(c->SoftenedTransit->dokillharm)
	    fitsoftened_transit(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],&c->SoftenedTransit->period[lc2],&c->SoftenedTransit->T0[lc2],&c->SoftenedTransit->eta[lc2],&c->SoftenedTransit->cval[lc2],&c->SoftenedTransit->delta[lc2],&c->SoftenedTransit->mconst[lc2],c->SoftenedTransit->fitephem,c->SoftenedTransit->fiteta,c->SoftenedTransit->fitcval,c->SoftenedTransit->fitdelta,c->SoftenedTransit->fitmconst,&c->SoftenedTransit->chisq[lc2],c->SoftenedTransit->correctlc,c->SoftenedTransit->omodel,outname,c->SoftenedTransit->dokillharm,c->SoftenedTransit->nharm,c->SoftenedTransit->nsubharm,c->SoftenedTransit->per_harm_out[lc2],c->SoftenedTransit->subharmA[lc2],c->SoftenedTransit->subharmB[lc2],c->SoftenedTransit->harmA[lc2],c->SoftenedTransit->harmB[lc2],&c->SoftenedTransit->fundA[lc2],&c->SoftenedTransit->fundB[lc2]);
	  else
	    fitsoftened_transit(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],&c->SoftenedTransit->period[lc2],&c->SoftenedTransit->T0[lc2],&c->SoftenedTransit->eta[lc2],&c->SoftenedTransit->cval[lc2],&c->SoftenedTransit->delta[lc2],&c->SoftenedTransit->mconst[lc2],c->SoftenedTransit->fitephem,c->SoftenedTransit->fiteta,c->SoftenedTransit->fitcval,c->SoftenedTransit->fitdelta,c->SoftenedTransit->fitmconst,&c->SoftenedTransit->chisq[lc2],c->SoftenedTransit->correctlc,c->SoftenedTransit->omodel,outname,c->SoftenedTransit->dokillharm,0,0,0.,NULL,NULL,NULL,NULL,NULL,NULL);
	}
      break;

    case CNUM_MANDELAGOLTRANSIT:
      /* Fit a Mandel and Agol (2002) Transit model to the light curve and remove it if we're doing that */
      if(c->MandelAgolTransit->omodel)
	{
	  if(!strncmp(c->MandelAgolTransit->modeloutdir,"-",1) && strlen(c->MandelAgolTransit->modeloutdir) == 1)
	    {
	      outname[0] = '-';
	      outname[1] = '\0';
	    }
	  else
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->MandelAgolTransit->modeloutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->modelsuffix);
	    }
	}
      if(c->MandelAgolTransit->ophcurve)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname3,"%s/%s%s",c->MandelAgolTransit->ophcurveoutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->ophcurvesuffix);
	}
      if(c->MandelAgolTransit->ojdcurve)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname4,"%s/%s%s",c->MandelAgolTransit->ojdcurveoutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->ojdcurvesuffix);
	}
      if(c->MandelAgolTransit->frombls)
	{
	  /* Get the starting parameters from BLS if that's what we're doing */
	  i1 = c->MandelAgolTransit->lastblsindex;
	  c->MandelAgolTransit->P0 = c[i1 - thisindex].Bls->bper[lc2][0];
	  c->MandelAgolTransit->T00 = c[i1 - thisindex].Bls->bt0[lc2][0];
	  /*
	    if(c[i1 - thisindex].Bls->i2[lc2][0] > c[i1 - thisindex].Bls->i1[lc2][0])
	    c->MandelAgolTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].Bls->i1[lc2][0] + c[i1 - thisindex].Bls->i2[lc2][0])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].Bls->nbins);
	  else
	    c->MandelAgolTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].Bls->i1[lc2][0] + 1 - c[i1 - thisindex].Bls->i2[lc2][0])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].Bls->nbins);
	  */

	  c->MandelAgolTransit->a0 = 1./c[i1 - thisindex].Bls->qtran[lc2][0]/M_PI;
	  c->MandelAgolTransit->r0 = sqrt(c[i1 - thisindex].Bls->depth[lc2][0]);
	  c->MandelAgolTransit->mconst0 = -1;
	  c->MandelAgolTransit->inputinclterm = 1;
	  c->MandelAgolTransit->bimpact0 = 0.2;
	  //c->MandelAgolTransit->sin_i0 = 0.5*(1.0 + getminsini(c->MandelAgolTransit->a0, 0., 0., c->MandelAgolTransit->r0));
	  c->MandelAgolTransit->e0 = 0.;
	  c->MandelAgolTransit->omega0 = 0.;
	}
      else if(c->MandelAgolTransit->fromblsfixper)
	{
	  /* Get the starting parameters from BLSFixPer if that's what we're doing */
	  /* Get the starting parameters from BLS if that's what we're doing */
	  i1 = c->MandelAgolTransit->lastblsfixperindex;
	  c->MandelAgolTransit->P0 = c[i1 - thisindex].BlsFixPer->period[lc2][0];
	  c->MandelAgolTransit->T00 = c[i1 - thisindex].BlsFixPer->bt0[lc2];
	  /*
	  if(c[i1 - thisindex].BlsFixPer->i2[lc2] > c[i1 - thisindex].BlsFixPer->i1[lc2])
	    c->MandelAgolTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc2] + c[i1 - thisindex].BlsFixPer->i2[lc2])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	  else
	    c->MandelAgolTransit->T00 = p->t[lc2][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc2] + 1 - c[i1 - thisindex].BlsFixPer->i2[lc2])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	  */

	  c->MandelAgolTransit->a0 = 1./c[i1 - thisindex].BlsFixPer->qtran[lc2]/M_PI;
	  c->MandelAgolTransit->r0 = sqrt(c[i1 - thisindex].BlsFixPer->depth[lc2]);
	  c->MandelAgolTransit->mconst0 = -1;
	  c->MandelAgolTransit->inputinclterm = 1;
	  c->MandelAgolTransit->bimpact0 = 0.2;
	  //c->MandelAgolTransit->sin_i0 = 0.5*(1.0 + getminsini(c->MandelAgolTransit->a0, 0., 0., c->MandelAgolTransit->r0));
	  c->MandelAgolTransit->e0 = 0.;
	  c->MandelAgolTransit->omega0 = 0.;
	}
      c->MandelAgolTransit->period[lc2] = c->MandelAgolTransit->P0;
      c->MandelAgolTransit->T0[lc2] = c->MandelAgolTransit->T00;
      c->MandelAgolTransit->r[lc2] = c->MandelAgolTransit->r0;
      c->MandelAgolTransit->a[lc2] = c->MandelAgolTransit->a0;
      c->MandelAgolTransit->e[lc2] = c->MandelAgolTransit->e0;
      c->MandelAgolTransit->omega[lc2] = c->MandelAgolTransit->omega0;
      if(c->MandelAgolTransit->inputinclterm) {
	c->MandelAgolTransit->bimpact[lc2] = c->MandelAgolTransit->bimpact0;
	c->MandelAgolTransit->inc[lc2] = c->MandelAgolTransit->bimpact0*(1. + c->MandelAgolTransit->e0*cos(c->MandelAgolTransit->omega0))/(1. - c->MandelAgolTransit->e0*c->MandelAgolTransit->e0)/(c->MandelAgolTransit->a0);
	c->MandelAgolTransit->inc[lc2] = 180.*acos(c->MandelAgolTransit->inc[lc2])/M_PI;
      }
      else {
	c->MandelAgolTransit->inc[lc2] = c->MandelAgolTransit->inc0;
	c->MandelAgolTransit->bimpact[lc2] = cos(c->MandelAgolTransit->inc0*M_PI/180.)*(1. - c->MandelAgolTransit->e0*c->MandelAgolTransit->e0)*(c->MandelAgolTransit->a0)/(1. + c->MandelAgolTransit->e0*cos(c->MandelAgolTransit->omega0));
      }
      c->MandelAgolTransit->sin_i[lc2] = c->MandelAgolTransit->sin_i0;
      c->MandelAgolTransit->ldcoeffs[lc2][0] = c->MandelAgolTransit->ldcoeffs0[lc2];
      c->MandelAgolTransit->ldcoeffs[lc2][1] = c->MandelAgolTransit->ldcoeffs0[1];
      c->MandelAgolTransit->ldcoeffs[lc2][2] = c->MandelAgolTransit->ldcoeffs0[2];
      c->MandelAgolTransit->ldcoeffs[lc2][3] = c->MandelAgolTransit->ldcoeffs0[3];
      c->MandelAgolTransit->mconst[lc2] = c->MandelAgolTransit->mconst0;
      c->MandelAgolTransit->gamma[lc2] = c->MandelAgolTransit->gamma0;
      c->MandelAgolTransit->K[lc2] = c->MandelAgolTransit->K0;
      if(p->NJD[lc2] > 1)
	{
	  fitmandelagoltransit_amoeba(p, p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],&c->MandelAgolTransit->period[lc2],&c->MandelAgolTransit->T0[lc2],&c->MandelAgolTransit->r[lc2],&c->MandelAgolTransit->a[lc2],&c->MandelAgolTransit->inc[lc2],&c->MandelAgolTransit->bimpact[lc2],&c->MandelAgolTransit->e[lc2],&c->MandelAgolTransit->omega[lc2],&c->MandelAgolTransit->mconst[lc2],c->MandelAgolTransit->type,c->MandelAgolTransit->ldcoeffs[lc2],c->MandelAgolTransit->fitephem,c->MandelAgolTransit->fitr,c->MandelAgolTransit->fita,c->MandelAgolTransit->fitinclterm,c->MandelAgolTransit->fite,c->MandelAgolTransit->fitomega,c->MandelAgolTransit->fitmconst,c->MandelAgolTransit->fitldcoeffs, &c->MandelAgolTransit->chisq[lc2],c->MandelAgolTransit->correctlc,c->MandelAgolTransit->omodel,outname, c->MandelAgolTransit->fitRV, c->MandelAgolTransit->RVinputfile, c->MandelAgolTransit->RVmodeloutfile, &c->MandelAgolTransit->K[lc2], &c->MandelAgolTransit->gamma[lc2], c->MandelAgolTransit->fitK, c->MandelAgolTransit->fitgamma, c->MandelAgolTransit->refititer, c->MandelAgolTransit->ophcurve, outname3, c->MandelAgolTransit->phmin, c->MandelAgolTransit->phmax, c->MandelAgolTransit->phstep, c->MandelAgolTransit->ojdcurve, outname4, c->MandelAgolTransit->jdstep, c->MandelAgolTransit->modelvarname, c->MandelAgolTransit->modelvar, lc2);
	}
      break;

    case CNUM_MICROLENS:
      /* Fit a Microlens model to the light curve */
      if(c->MicroLens->omodel)
	{
	  if(!strncmp(c->MicroLens->modeloutdir,"-",1) && strlen(c->MicroLens->modeloutdir) == 1)
	    {
	      outname[0] = '-';
	      outname[1] = '\0';
	    }
	  else
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->MicroLens->modeloutdir,&p->lcnames[lc][i2],c->MicroLens->modelsuffix);
	    }
	}
      microlens(p->t[lc2],p->mag[lc2],p->sig[lc2],p->NJD[lc2],lc,c->MicroLens,outname,&c->MicroLens->f0[lc2],&c->MicroLens->f1[lc2],&c->MicroLens->u0[lc2],&c->MicroLens->t0[lc2],&c->MicroLens->tmax[lc2],&c->MicroLens->chi2_[lc2]);
      break;

    case CNUM_SYSREM:
      /* There must be a bug in the program if we're calling sysrem in single light curve process mode! */
      error(ERR_CODEERROR);
      break;

    case CNUM_TFA:
      /* Run the trend filtering algorithm */
      if(p->NJD[lc2] > 1)
	{
	  if(c->TFA->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->TFA->model_outdir,&p->lcnames[lc][i2],c->TFA->model_suffix);
	    }
	  if(c->TFA->ocoeff)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname2,"%s/%s%s",c->TFA->coeff_outdir,&p->lcnames[lc][i2],c->TFA->coeff_suffix);
	    }
	  if(p->matchstringid)
	    detrend_tfa(p, c->TFA, p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->TFA->lcx[lc][0], c->TFA->lcy[lc][0], p->lcnames[lc], outname2, c->TFA->ocoeff, c->TFA->correctlc, c->TFA->omodel, outname, &c->TFA->ave_out[lc2], &c->TFA->rms_out[lc2], p->matchstringid, p->stringid[lc2], p->stringid_idx[lc2], lc2);
	  else
	    detrend_tfa(p, c->TFA, p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->TFA->lcx[lc][0], c->TFA->lcy[lc][0], p->lcnames[lc], outname2, c->TFA->ocoeff, c->TFA->correctlc, c->TFA->omodel, outname, &c->TFA->ave_out[lc2], &c->TFA->rms_out[lc2], 0, NULL, NULL, lc2);
	}
      break;

    case CNUM_TFA_SR:
      /* Run the trend filtering algorithm in SR mode */
      if(p->NJD[lc2] > 1)
	{
	  if(c->TFA_SR->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->TFA_SR->model_outdir,&p->lcnames[lc][i2],c->TFA_SR->model_suffix);
	    }
	  if(c->TFA_SR->ocoeff)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname2,"%s/%s%s",c->TFA_SR->coeff_outdir,&p->lcnames[lc][i2],c->TFA_SR->coeff_suffix);
	    }
	  if((c->TFA_SR->use_bin && c->TFA_SR->use_period) || (c->TFA_SR->use_harm && c->TFA_SR->use_period))
	    {
	      if(c->TFA_SR->pertype == PERTYPE_AOV)
		{
		  i1=c->TFA_SR->lastindex;
		  if(c[i1-thisindex].cnum == CNUM_AOV)
		    c->TFA_SR->periods[lc2][0] = c[i1-thisindex].Aov->peakperiods[lc2][0];
		  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		    c->TFA_SR->periods[lc2][0] = c[i1-thisindex].AovHarm->peakperiods[lc2][0];
		  d1 = c->TFA_SR->periods[lc2][0];
		}
	      else if(c->TFA_SR->pertype == PERTYPE_LS)
		{
		  i1 = c->TFA_SR->lastindex;
		  c->TFA_SR->periods[lc2][0] = c[i1-thisindex].Ls->peakperiods[lc2][0];
		  d1 = c->TFA_SR->periods[lc2][0];
		}
	      else if(c->TFA_SR->pertype == PERTYPE_BLS)
		{
		  i1 = c->TFA_SR->lastindex;
		  c->TFA_SR->periods[lc2][0] = c[i1-thisindex].Bls->bper[lc2][0];
		  d1 = c->TFA_SR->periods[lc2][0];
		}
	      else if(c->TFA_SR->pertype == PERTYPE_SPECIFIED)
		{
		  d1 = c->TFA_SR->periods[lc][0];
		}
	      else if(c->TFA_SR->pertype == PERTYPE_FIX)
		{
		  d1 = c->TFA_SR->fixperiod;
		}
	    }
	  else if(c->TFA_SR->use_harm)
	    d1 = p->t[lc2][p->NJD[lc2]-1] - p->t[lc2][0];
	  else
	    d1 = 1.0;
	  if(!c->TFA_SR->use_bin && !c->TFA_SR->use_harm)
	    s1 = c->TFA_SR->signalfilenames[lc];
	  else
	    s1 = NULL;
	  if(p->matchstringid)
	    detrend_tfa_sr(p, c->TFA_SR, p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->TFA_SR->lcx[lc][0], c->TFA_SR->lcy[lc][0], p->lcnames[lc], outname2, c->TFA_SR->ocoeff, c->TFA_SR->correctlc, c->TFA_SR->omodel, outname, &c->TFA_SR->ave_out[lc2], &c->TFA_SR->rms_out[lc2], d1, s1, p->matchstringid, p->stringid[lc2], p->stringid_idx[lc2], lc2, lc2);
	  else
	    detrend_tfa_sr(p, c->TFA_SR, p->NJD[lc2], p->t[lc2], p->mag[lc2], p->sig[lc2], c->TFA_SR->lcx[lc][0], c->TFA_SR->lcy[lc][0], p->lcnames[lc], outname2, c->TFA_SR->ocoeff, c->TFA_SR->correctlc, c->TFA_SR->omodel, outname, &c->TFA_SR->ave_out[lc2], &c->TFA_SR->rms_out[lc2], d1, s1, 0, NULL, NULL, lc2, lc2);
	}
      break;

    case CNUM_DFTCLEAN:
      /* Run DFT/CLEAN on the light curve */
      if(p->NJD[lc2] > 1)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  dodftclean(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],lc2,c->Dftclean,&p->lcnames[lc][i2],p->ascii, lc2, lc, c->Dftclean->usemask, c->Dftclean->maskvar);
	}
      break;

    case CNUM_BINLC:
      /* Bin the light curve */
      binlc(p,c->Binlc,lc2,lc);
      break;

    case CNUM_MATCHCOMMAND:
      /* Match the light curve to an external datafile */
      RunMatchCommand(p, c->MatchCommand, lc, lc2);
      break;

    case CNUM_MEDIANFILTER:
      /* median filter the light curve */
      medianfilter(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],c->MedianFilter->time,c->MedianFilter->usemean,c->MedianFilter->replace);
      break;

    case CNUM_PHASE:
      /* Replace the time coordinate of a light curve with its phase */
      if(p->NJD[lc2] > 1)
	{
	  if(c->Phase->pertype == PERTYPE_AOV)
	    {
	      i1 = c->Phase->lastaovindex;
	      if(c[i1-thisindex].cnum == CNUM_AOV) {
		d1 = c[i1-thisindex].Aov->peakperiods[lc2][0];
		c->Phase->period[lc2][0] = d1;
	      }
	      else if(c[i1-thisindex].cnum == CNUM_HARMAOV) {
		d1 = c[i1-thisindex].AovHarm->peakperiods[lc2][0];
		c->Phase->period[lc2][0] = d1;
	      }
	    }
	  else if(c->Phase->pertype == PERTYPE_LS)
	    {
	      i1 = c->Phase->lastlsindex;
	      d1 = c[i1-thisindex].Ls->peakperiods[lc2][0];
	      c->Phase->period[lc2][0] = d1;
	    }
	  else if(c->Phase->pertype == PERTYPE_BLS)
	    {
	      i1 = c->Phase->lastblsindex;
	      d1 = c[i1-thisindex].Bls->bper[lc2][0];
	      c->Phase->period[lc2][0] = d1;
	    }
	  else if(c->Phase->pertype == PERTYPE_SPECIFIED)
	    {
	      d1 = c->Phase->period[lc][0];
	    }
	  else if(c->Phase->pertype == PERTYPE_FIX)
	    {
	      d1 = c->Phase->fixperiod;
	      c->Phase->period[lc2][0] = d1;
	    }
	  else if(c->Phase->pertype == PERTYPE_FIXCOLUMN)
	    {
	      getoutcolumnvalue(c->Phase->period_linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &d1);
	      c->Phase->period[lc2][0] = d1;
	    }
	  if(c->Phase->t0type == PERTYPE_AUTOFIND) {
	    d2 = 0.; i2 = 0;
	  }
	  else {
	    i2 = 1;
	    if(c->Phase->t0type == PERTYPE_BLS) {
	      i1 = c->Phase->lastblsindex;
	      /* Get Tc for the transit */
	      if(c[i1-thisindex].Bls->i1[lc2][0] > c[i1-thisindex].Bls->i2[lc2][0]) {
		d2 = p->t[lc2][0] + c[i1-thisindex].Bls->bper[lc2][0]*0.5*(c[i1-thisindex].Bls->i1[lc2][0]+1.+c[i1-thisindex].Bls->i2[lc2][0]+(1./((double)c[i1-thisindex].Bls->nbins_val[lc2])));
	      }
	      else {
		d1 = p->t[lc2][0] + c[i1-thisindex].Bls->bper[lc2][0]*0.5*(c[i1-thisindex].Bls->i1[lc2][0]+c[i1-thisindex].Bls->i2[lc2][0]+(1./((double)c[i1-thisindex].Bls->nbins_val[lc2])));
	      }
	      /* adjust so that Tc has phase phaseTc */
	      d2 = d2 - c->Phase->phaseTc*d1;
	    }
	    else if(c->Phase->t0type == PERTYPE_SPECIFIED) {
	      d2 = c->Phase->T0[lc][0];
	    }
	    else if(c->Phase->t0type == PERTYPE_FIX) {
	      d2 = c->Phase->fixT0;
	    }
	    else if(c->Phase->t0type == PERTYPE_FIXCOLUMN) {
	      getoutcolumnvalue(c->Phase->T0_linkedcolumn, lc2, lc, VARTOOLS_TYPE_DOUBLE, &d2);
	    }
	  }

	  phaselc(p->NJD[lc2],p->t[lc2],p->mag[lc2],p->sig[lc2],d1,i2,d2,c->Phase->phasevarname,c->Phase->phasevar,lc2,c->Phase->startphase);
	  sortlcbytime(p->NJD[lc2], p->t[lc2], lc2, p);
	  /*
	  if(c->Phase->pertype != PERTYPE_SPECIFIED)
	    phaselc(p->NJD[0],p->t[0],p->mag[0],p->sig[0],c->Phase->period[0][0]);
	  else
	    phaselc(p->NJD[0],p->t[0],p->mag[0],p->sig[0],c->Phase->period[lc][0]);
	  */
	}
      break;

    case CNUM_OUTPUTLCS:
      /* Write out the light curves in their present form */
      DoOutputLightCurve(p, c->Outputlcs, lc, lc2);
      /*
      if(p->listflag)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  if(!c->Outputlcs->useformat)
	    sprintf(outname,"%s/%s",c->Outputlcs->outdir,&p->lcnames[lc][i2]);
	  else
	    {
	      sprintf(outname,"%s/",c->Outputlcs->outdir);
	      i1=strlen(outname);
	      i3=0;
	      while(c->Outputlcs->format[i3] != '\0')
		{
		  if(c->Outputlcs->format[i3] != '%')
		    {
		      outname[i1] = c->Outputlcs->format[i3];
		      i1++;
		      outname[i1] = '\0';
		      i3++;
		    }
		  else
		    {
		      i3++;
		      if(c->Outputlcs->format[i3] == 's')
			{
			  i3++;
			  sprintf(&outname[i1],"%s",&p->lcnames[lc][i2]);
			  i1 = strlen(outname);
			}
		      else if(c->Outputlcs->format[i3] == 'd')
			{
			  i3++;
			  sprintf(&outname[i1],"%d",lc+1);
			  i1 = strlen(outname);
			}
		      else if(c->Outputlcs->format[i3] == '0')
			{
			  i3++;
			  tmpstring[0] = '%';
			  tmpstring[1] = '0';
			  i4 = 2;
			  while(c->Outputlcs->format[i3] >= '1' && c->Outputlcs->format[i3] <= '9')
			    {
			      tmpstring[i4] = c->Outputlcs->format[i3];
			      i4++;
			      i3++;
			    }
			  if(c->Outputlcs->format[i3] != 'd')
			    error(ERR_INVALIDOUTPUTFORMAT);
			  i3++;
			  tmpstring[i4] = 'd';
			  i4++;
			  tmpstring[i4] = '\0';
			  sprintf(&outname[i1],tmpstring,lc+1);
			  i1 = strlen(outname);
			}
		      else if(c->Outputlcs->format[i3] == '%')
			{
			  i3++;
			  outname[i1] = '%';
			  i1++;
			  outname[i1] = '\0';
			}
		      else
			error(ERR_INVALIDOUTPUTFORMAT);
		    }
		}
	    }
	  writelightcurves(p, lc2, lc, outname, c->Outputlcs->usecolumnformat, c->Outputlcs->Nvar, c->Outputlcs->variables, c->Outputlcs->printfformats);
	}
      else if(p->fileflag)
	{
	  writelightcurves(p, lc2, lc, c->Outputlcs->outdir, c->Outputlcs->usecolumnformat, c->Outputlcs->Nvar, c->Outputlcs->variables, c->Outputlcs->printfformats);
	  }*/
      break;

#ifdef DYNAMICLIB
    case CNUM_USERCOMMAND:
      RunUserCommand(p,c,lc,lc2);
      break;

#ifdef _HAVE_PYTHON
    case CNUM_PYTHON:
      RunPythonCommand(p, lc, lc2, lc2, c->PythonCommand);
      break;
#endif

#ifdef _HAVE_R
    case CNUM_R:
      RunRCommand(p, lc, lc2, lc2, c->RCommand);
      break;
#endif

#endif

    case CNUM_IF:
      break;

    default:
      error(ERR_CODEERROR);
      break;
    }
}


void ProcessCommandAll(ProgramData *p, Command *c, int thisindex)
{
  double d0, d1, d2, d3, d4;
  double *ers1, *ers2, *ers3;
  int i1, i2, i3, i, lc, i4;
  char *s1;
  double *d1ptr, *d2ptr, *d3ptr;
  char outname[MAXLEN], outname2[MAXLEN], outname3[MAXLEN], outname4[MAXLEN],
    tmpstring[256];
  /* This routine runs a command on all of the light curves */

  /* Sort the light curves in time, and merge unequal values if needed */
  if(c->require_sort || c->require_distinct) {
    for(lc=0; lc<p->Nlcs;lc++) {
      i1 = sortlcbytime(p->NJD[lc], p->t[lc], lc, p);
      if(c->require_distinct && i1) {
	mergeequallctimes(p, lc);
      }
    }
  }

  switch(c->cnum)
    {

    case CNUM_SORTLC:
      /* Sort the light curve by a vector */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	
	if(!c->SortLC->issortvar && !c->SortLC->isreverse) {
	  i1 = sortlcbytime(p->NJD[lc], p->t[lc], lc, p);
	}
	else if(!c->SortLC->issortvar && c->SortLC->isreverse) {
	  i1 = sortlcbytime_rev(p->NJD[lc], p->t[lc], lc, p);
	}
	else if(!c->SortLC->isreverse) {
	  switch(c->SortLC->sortdtype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    i1 = sortlcbyvardbl(p->NJD[lc], (*((double ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_INT:
	    i1 = sortlcbyvarint(p->NJD[lc], (*((int ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    i1 = sortlcbyvarstring(p->NJD[lc], (*((char ****) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    i1 = sortlcbyvarfloat(p->NJD[lc], (*((float ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    i1 = sortlcbyvarchar(p->NJD[lc], (*((char ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    i1 = sortlcbyvarlong(p->NJD[lc], (*((long ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  default:
	    error(ERR_CODEERROR);
	    break;
	  }
	} else {
	  switch(c->SortLC->sortdtype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    i1 = sortlcbyvardbl_rev(p->NJD[lc], (*((double ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_INT:
	    i1 = sortlcbyvarint_rev(p->NJD[lc], (*((int ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    i1 = sortlcbyvarstring_rev(p->NJD[lc], (*((char ****) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    i1 = sortlcbyvarfloat_rev(p->NJD[lc], (*((float ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    i1 = sortlcbyvarchar_rev(p->NJD[lc], (*((char ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    i1 = sortlcbyvarlong_rev(p->NJD[lc], (*((long ***) (c->SortLC->sortvar->dataptr)))[lc], lc, p);
	    break;
	  default:
	    error(ERR_CODEERROR);
	    break;
	  }
	}
      }
      break;

    case CNUM_DIFFFLUXTOMAG:
      /* Convert from isis differential flux to magnitudes */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	difffluxtomag(p->t[lc],p->mag[lc],p->sig[lc],p->NJD[lc],c->DiffFluxtomag->magstar[lc][0],c->DiffFluxtomag->mag_constant1, c->DiffFluxtomag->offset);
      }
      break;

    case CNUM_FLUXTOMAG:
      /* Convert from isis differential flux to magnitudes */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	fluxtomag(p->t[lc],p->mag[lc],p->sig[lc],p->NJD[lc],c->Fluxtomag->mag_constant1, c->Fluxtomag->offset);
      }
      break;

    case CNUM_EXPRESSION:
      /* Evaluate an analytic expression */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	RunExpressionCommand(lc, lc, p, c->ExpressionCommand);
      }
      break;

    case CNUM_LINFIT:
      /* Fit a model that is linear in its free parameters to the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	DoLinfit(p, c->Linfit, lc, lc);
      }
      break;

    case CNUM_NONLINFIT:
      /* Fit a model that is nonlinear in its free parameters to the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	DoNonlinfit(p, c->Nonlinfit, lc, lc);
      }
      break;

    case CNUM_WWZ:
      /* Run the WWZ Transform */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	DoWWZ(p, c->WWZ, lc, lc);
      }
      break;

    case CNUM_FINDBLENDS:
      /* Find variability blends */
      /* Set the period and potential variable x, y coordinates*/
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(c->FindBlends->pertype == PERTYPE_FIXCOLUMN)
	    {
	      getoutcolumnvalue(c->FindBlends->linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &(c->FindBlends->periods[lc][0]));
	    }
	  else if(c->FindBlends->pertype == PERTYPE_FIX)
	    {
	      c->FindBlends->periods[lc][0] = c->FindBlends->fixperiod;
	    }
	  getoutcolumnvalue(c->FindBlends->linkedcolumn_varname, lc, lc, VARTOOLS_TYPE_STRING, &(c->FindBlends->varnames[lc][0]), MAXLEN);
	  c->FindBlends->varx[lc] = c->FindBlends->varxyin[lc][0];
	  c->FindBlends->vary[lc] = c->FindBlends->varxyin[lc][1];
	}
      findblends(p->Nlcs,p->NJD,p->t,p->mag,p->sig,c->FindBlends);
      break;

    case CNUM_CLIP:
      /* Clip the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	c->Clip->Nclip[lc] = sigclip(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &d1, &d2, &d3, &i1, c->Clip->sigclip, c->Clip->iter, lc, p, c->Clip->niter, c->Clip->usemedian, c->Clip->markclip, c->Clip->clipvar, c->Clip->noinitmark);
      }
      break;

    case CNUM_CONVERTTIME:
      /* Perform a time conversion */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	converttime(p->NJD[lc], p->t[lc], lc, lc, c->ConvertTime, p);
      }
      break;

    case CNUM_ENSEMBLERESCALESIG:
      /* Perform Ensemble sigma rescaling */
      if((ers1 = (double *) malloc(p->Nlcs * sizeof(double))) == NULL ||
	 (ers2 = (double *) malloc(p->Nlcs * sizeof(double))) == NULL ||
	 (ers3 = (double *) malloc(p->Nlcs * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      i2 = 0;
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  d1 = rms(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&d1,&d2,&i1,c->Ensemblerescalesig->usemask,c->Ensemblerescalesig->maskvar,lc,lc);
	  c->Ensemblerescalesig->chi2_old[lc] = chi2(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&d3,&i1,c->Ensemblerescalesig->usemask,c->Ensemblerescalesig->maskvar,lc,lc);
	  if(i1 > 1 && c->Ensemblerescalesig->chi2_old[lc] > 0. && d2 > 0. && d2 < RMSTHYCUT && d1 > 0.)
	    {
	      c->Ensemblerescalesig->chi2_old[lc] /= (double) (i1 - 1);
	      ers1[i2] = d2*d2;
	      ers2[i2] = d2*d2*(c->Ensemblerescalesig->chi2_old[lc]);
	      ers3[i2] = 4.0 * d1 * d1 / sqrt(i1);
	      if(ers3[i2] > 0. && !isnan(ers3[i2]))
		i2++;
	    }
	}
      if(i2 < 5)
	error(ERR_NOTENOUGHSTARS_ERS);
      c->Ensemblerescalesig->a = ((mean2(i2,ers1,ers2,ers3)*mean1(i2,ers3))-(mean(i2,ers2,ers3)*mean(i2,ers1,ers3)))/((mean1(i2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean(i2,ers1,ers3)));
      c->Ensemblerescalesig->b = ((mean(i2,ers2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean2(i2,ers1,ers2,ers3)))/((mean1(i2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean(i2,ers1,ers3)));
      while((i3 = purge_bad(&i2,ers1,ers2,ers3,c->Ensemblerescalesig->a,c->Ensemblerescalesig->b,c->Ensemblerescalesig->erssigclip,1)) > 0)
	{
	  c->Ensemblerescalesig->a = ((mean2(i2,ers1,ers2,ers3)*mean1(i2,ers3))-(mean(i2,ers2,ers3)*mean(i2,ers1,ers3)))/((mean1(i2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean(i2,ers1,ers3)));
	  c->Ensemblerescalesig->b = ((mean(i2,ers2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean2(i2,ers1,ers2,ers3)))/((mean1(i2,ers3)*mean2(i2,ers1,ers1,ers3))-(mean(i2,ers1,ers3)*mean(i2,ers1,ers3)));
	}
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  c->Ensemblerescalesig->chi2_old[lc] = chi2(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&d3,&i1,c->Ensemblerescalesig->usemask,c->Ensemblerescalesig->maskvar,lc,lc);
	  if(i1 > 1)
	    c->Ensemblerescalesig->chi2_old[lc] /= (double) (i1 - 1);
	  else
	    c->Ensemblerescalesig->chi2_old[lc] = -1.;

	  rescalesigma_linear(p->NJD[lc], p->sig[lc], c->Ensemblerescalesig->a, c->Ensemblerescalesig->b);
	  c->Ensemblerescalesig->chi2_new[lc] = chi2(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&d3,&i1,c->Ensemblerescalesig->usemask,c->Ensemblerescalesig->maskvar,lc,lc);
	  if(i1 > 1)
	    c->Ensemblerescalesig->chi2_new[lc] /= (double) (i1 - 1);
	  else
	    c->Ensemblerescalesig->chi2_new[lc] = -1.;

	  c->Ensemblerescalesig->rescalefactor[lc] = sqrt(c->Ensemblerescalesig->chi2_new[lc] / c->Ensemblerescalesig->chi2_old[lc]);
	}
      free(ers1);
      free(ers2);
      free(ers3);
      break;

    case CNUM_HARMONICFILTER:
      /* Apply a fourier filter to the light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	doHarmonicFilter(p, c->HarmonicFilter, lc, lc);
      break;

    case CNUM_RESAMPLE:
      /* Resample the times of observation in the light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	DoResample(p, c->Resample, lc, lc);
      break;

    case CNUM_RESCALESIG:
      /* Rescale sigma for a light curve so that it has chi2 = 1 */
      /* First get the old chi2 value */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  c->Rescalesig->chi2_old[lc] = chi2(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &d1, &i1, c->Rescalesig->usemask, c->Rescalesig->maskvar, lc, lc);
	  if(i1 > 1)
	    {
	      c->Rescalesig->chi2_old[lc] /= (double) (i1 - 1);
	    }
	  else
	    c->Rescalesig->chi2_old[lc] = -1.;
	  /* Rescale sigma */
	  rescalesigma_chi2(p->NJD[lc], p->sig[lc], c->Rescalesig->chi2_old[lc]);
	  /* Get the new chi2 value */
	  c->Rescalesig->chi2_new[lc] = chi2(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &d1, &i2, c->Rescalesig->usemask, c->Rescalesig->maskvar, lc, lc);
	  if(i2 > 1)
	    {
	      c->Rescalesig->chi2_new[lc] /= (double) (i2 - 1);
	    }
	  else
	    c->Rescalesig->chi2_new[lc] = -1.;
	  /* Calculate the rescale factor */
	  c->Rescalesig->rescalefactor[lc] = sqrt(c->Rescalesig->chi2_new[lc] / c->Rescalesig->chi2_old[lc]);
	}
      break;

    case CNUM_CHI2_NOBIN:
      /* calculate chi2 without binning */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  c->Chi2_NoBin->chi2val[lc] = chi2(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &c->Chi2_NoBin->wtave[lc], &i1, c->Chi2_NoBin->usemask, c->Chi2_NoBin->maskvar, lc, lc);
	  if(i1 > 1)
	    c->Chi2_NoBin->chi2val[lc] /= (double) (i1 - 1);
	  else
	    c->Chi2_NoBin->chi2val[lc] = -1.;
	}
      break;

    case CNUM_SAVELC:
      /* Save the light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	dosavelc(p, c->Savelc, lc, lc);
      break;

    case CNUM_COPYLC:
      /* This command cannot be used with readallformat */
      error(ERR_READALL_ANDCOPYLC);
      break;

    case CNUM_RESTORELC:
      /* Restore the saved light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	dorestorelc(p, c[c->Restorelc->saveindex - thisindex].Savelc, c->Restorelc, lc, lc, lc);
      break;

    case CNUM_RESTRICTTIMES:
      /* Filter points from the light curves based on the times of
	 observation */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDRANGE) {
	  if(c->RestrictTimes->minJDtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->RestrictTimes->minJD_linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &(c->RestrictTimes->minJD[lc]));
	  }
	  else if(c->RestrictTimes->minJDtype == PERTYPE_FIX) {
	    c->RestrictTimes->minJD[lc] = c->RestrictTimes->minJDfixval;
	  }
	  else if(c->RestrictTimes->minJDtype == PERTYPE_EXPR) {
	    d1 = EvaluateExpression(lc, lc, 0, c->RestrictTimes->minJDexpr);
	    c->RestrictTimes->minJD[lc] = d1;
	  }
	  if(c->RestrictTimes->maxJDtype == PERTYPE_FIXCOLUMN) {
	    getoutcolumnvalue(c->RestrictTimes->maxJD_linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &(c->RestrictTimes->maxJD[lc]));
	  }
	  else if(c->RestrictTimes->maxJDtype == PERTYPE_FIX) {
	    c->RestrictTimes->maxJD[lc] = c->RestrictTimes->maxJDfixval;
	  }
	  else if(c->RestrictTimes->maxJDtype == PERTYPE_EXPR) {
	    d2 = EvaluateExpression(lc, lc, 0, c->RestrictTimes->maxJDexpr);
	    c->RestrictTimes->maxJD[lc] = d2;
	  }
	  RestrictTimes_JDrange_apply(p->NJD[lc], p->t[lc], lc, p,
				      c->RestrictTimes,
				      c->RestrictTimes->minJD[lc],
				      c->RestrictTimes->maxJD[lc],
				      c->RestrictTimes->exclude,
				      c->RestrictTimes->markrestrict,
				      c->RestrictTimes->markvar, 
				      c->RestrictTimes->noinitmark);
	}
	else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDLIST) {
	  RestrictTimes_JDlist_apply(p->NJD[lc], p->t[lc], lc, p,
				     c->RestrictTimes,
				     c->RestrictTimes->JD_restrictlist,
				     c->RestrictTimes->N_restrictlist,
				     c->RestrictTimes->exclude,
				     c->RestrictTimes->markrestrict,
				     c->RestrictTimes->markvar, 
				     c->RestrictTimes->noinitmark);
	}
	else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_IMAGELIST) {
	  RestrictTimes_imagelist_apply(p->NJD[lc], p->stringid[lc], 
					p->stringid_idx[lc], lc, p,
					c->RestrictTimes,
					c->RestrictTimes->image_restrictlist,
					c->RestrictTimes->image_restrictlist_indx,
					c->RestrictTimes->N_restrictlist,
					c->RestrictTimes->exclude,
					c->RestrictTimes->markrestrict,
					c->RestrictTimes->markvar, 
					c->RestrictTimes->noinitmark);
	}
	else if(c->RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_EXPR) {
	  RestrictTimes_expr_apply(p, c->RestrictTimes, lc, lc,
				   c->RestrictTimes->markrestrict,
				   c->RestrictTimes->markvar, 
				   c->RestrictTimes->noinitmark);
	}
      }
      break;

    case CNUM_RESTORETIMES:
      /* Restore the saved light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	RestoreTimes(p, c->RestoreTimes, lc, lc);
      break;

    case CNUM_CHI2_BIN:
      /* calculate chi2 with binning */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  for(i=0; i < c->Chi2_Bin->Nbin ; i++)
	    {
	      c->Chi2_Bin->chi2binvals[lc][i] = binnedchi2(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->Chi2_Bin->bintimes[i], &c->Chi2_Bin->wtavebin[lc][i], &i1, c->Chi2_Bin->usemask, c->Chi2_Bin->maskvar, lc, lc);
	      if(i1 > 1)
		c->Chi2_Bin->chi2binvals[lc][i] /= (double) (i1 - 1);
	      else
		c->Chi2_Bin->chi2binvals[lc][i] = -1.;
	    }
	}
      break;

    case CNUM_CHANGEERROR:
      /* Replace the formal errors in a light curve with the RMS */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  c->Changeerror->rmsval[lc] = changeerror(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &c->Changeerror->ave[lc], &c->Changeerror->ngood[lc], c->Changeerror->usemask, c->Changeerror->maskvar, lc, lc);
	}
      break;

    case CNUM_CHANGEVARIABLE:
      /* Switch the time, mag, sig, or string-id variable */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  DoChangeVariable(p, c->Changevariable, lc);
	}
      break;

    case CNUM_RMS_NOBIN:
      /* calculate RMS without binning */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  c->RMS_NoBin->rmsval[lc] = rms(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &c->RMS_NoBin->ave[lc], &c->RMS_NoBin->rmsthy[lc], &c->RMS_NoBin->ngood[lc], c->RMS_NoBin->usemask, c->RMS_NoBin->maskvar, lc, lc);
	}
      break;

    case CNUM_RMS_BIN:
      /* Calculate RMS with binning */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  for(i=0; i < c->RMS_Bin->Nbin ; i++)
	    {
	      c->RMS_Bin->rmsbinvals[lc][i] = binnedrms(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->RMS_Bin->bintimes[i], &d1, &c->RMS_Bin->rmsthybin[lc][i], &i1, c->RMS_Bin->usemask, c->RMS_Bin->maskvar, lc, lc);
	    }
	}
      break;

    case CNUM_JSTET:
      /* Calculate JSTET */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  getJstet(p->NJD[lc], c->Jstet->Jstet_time, c->Jstet->wkmax, p->t[lc], p->mag[lc], p->sig[lc], &d1, &c->Jstet->jst[lc], &c->Jstet->kur[lc], &c->Jstet->lst[lc], lc, lc, c->Jstet->usemask, c->Jstet->maskvar);
	}
      break;

    case CNUM_ALARM:
      /* Calculate the Alarm */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  c->Alarm->alarmvals[lc] = doalarm(p->NJD[lc], p->mag[lc], p->sig[lc],
					 lc, lc, c->Alarm->usemask,
					 c->Alarm->maskvar);
	}
      break;

    case CNUM_ADDFITSKEYWORD:
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	Run_AddFitsKeyword_Command(p, c->AddFitsKeyword, lc, lc);
	break;
      }

    case CNUM_ADDNOISE:
      /* Add time-correlated noise to the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	addnoise(p, c->AddNoise, lc, lc);
	break;
      }

    case CNUM_AUTOCORR:
      /* Calculate the auto-correlation */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  i2 = 0;
	  i1 = 0;
	  while(p->lcnames[lc][i1] != '\0')
	    {
	      if(p->lcnames[lc][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",c->Autocorr->outdir,&p->lcnames[lc][i2],c->Autocorr->suffix);
	  autocorrelation(p->t[lc], p->mag[lc], p->sig[lc], p->NJD[lc], c->Autocorr->start, c->Autocorr->stop, c->Autocorr->step, outname, lc, lc, c->Autocorr->usemask, c->Autocorr->maskvar);
	}
      break;

    case CNUM_AOV:
      /* Calculate the AoV with phase binning */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  RunAOVCommand(p, c, c->Aov, lc, lc, thisindex);
	}
      break;

    case CNUM_HARMAOV:
      /* Calculate the AoV with harmonics */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  RunAOVHarmCommand(p, c, c->AovHarm, lc, lc, thisindex);
	}
      break;

    case CNUM_LS:
      /* Calculate the Lomb-Scargle Periodogram */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  RunLombScargleCommand(p, c->Ls, c, lc, lc, thisindex);
	}
      break;

#ifdef _HAVE_GSL
    case CNUM_FFT:
      /* Calculate light curve statistics */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  RunFFTCommand(p, lc, lc, c->FFT);
	}
      break;
#endif

    case CNUM_GETLSAMPTHRESH:
      /* Get the amplitude scale-factor for which the signal just passes */
      /* The LS threshhold */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->GetLSAmpThresh->pertype == PERTYPE_LS)
		{
		  i1 = c->GetLSAmpThresh->lastlsindex;
		  c->GetLSAmpThresh->period[lc][0] = c[i1-thisindex].Ls->peakperiods[lc][0];
		}
	      if(c->GetLSAmpThresh->harm_specsigflag)
		{
		  if(gnu_getline(&(c->GetLSAmpThresh->line),&(c->GetLSAmpThresh->line_size),c->GetLSAmpThresh->listfile) < 0)
		    {
		      error2(ERR_GETLSAMPTHRESH_FILETOSHORT,c->GetLSAmpThresh->listfilename);
		    }
		  sscanf(c->GetLSAmpThresh->line,"%s %lf",c->GetLSAmpThresh->filename,c->GetLSAmpThresh->amp);
		  if((c->GetLSAmpThresh->infile = fopen(c->GetLSAmpThresh->filename,"r")) == NULL)
		    error2(ERR_FILENOTFOUND,c->GetLSAmpThresh->filename);
		  c->GetLSAmpThresh->sizesigfile = 0;
		  while(gnu_getline(&(c->GetLSAmpThresh->line),&(c->GetLSAmpThresh->line_size),c->GetLSAmpThresh->infile) >= 0)
		    c->GetLSAmpThresh->sizesigfile = c->GetLSAmpThresh->sizesigfile + 1;
		  if(c->GetLSAmpThresh->sizesigfile != p->NJD[lc])
		    error2(ERR_SIGFILEWRONGLENGTH,c->GetLSAmpThresh->filename);
		}
	      getlsampthresh(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->GetLSAmpThresh->period[lc][0],c->GetLSAmpThresh->harm_specsigflag,c->GetLSAmpThresh->infile,c->GetLSAmpThresh->Nsubharm,c->GetLSAmpThresh->Nharm,c->GetLSAmpThresh->minPer,c->GetLSAmpThresh->thresh,&c->GetLSAmpThresh->ampthresh_scale[lc],&c->GetLSAmpThresh->amp[lc],c->GetLSAmpThresh->use_orig_ls);
	    }
	}
      break;

    case CNUM_DECORR:
      /* Decorrelate the light curves */
      /* First get the lc average */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  d1 = rms(p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], &d2, &d3, &i1, 0, NULL, lc, lc);
	  if(c->Decorr->omodel)
	    {
	      i3 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i3] != '\0')
		{
		  if(p->lcnames[lc][i3] == '/')
		    i2 = i3 + 1;
		  i3++;
		}
	      sprintf(outname,"%s/%s%s",c->Decorr->modeloutdir,&p->lcnames[lc][i2],c->Decorr->modelsuffix);
	    }
	  /* Do the decorrelation only if there is at least 1 degree of freedom left over */
	  if(i1 >= c->Decorr->N_decorrterms_total + 1)
	    {
	      docorr(p->mag[lc], p->sig[lc], p->NJD[lc], c->Decorr->N_decorrterms, c->Decorr->decorr_terms[lc], c->Decorr->order, c->Decorr->b[lc], c->Decorr->b_err[lc], d2, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc);

	      if(c->Decorr->correctlc)
		{
		  magcorr(((void *) p->t[lc]),VARTOOLS_TYPE_DOUBLE,p->mag[lc], p->sig[lc], p->NJD[lc], c->Decorr->N_decorrterms, c->Decorr->decorr_terms[lc], c->Decorr->order, c->Decorr->b[lc], &c->Decorr->chi2val[lc], &d3, d2, c->Decorr->omodel, outname, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc);
		}
	      else
		magcorr_chi2only(p->t[lc],p->mag[lc],p->sig[lc], p->NJD[lc], c->Decorr->N_decorrterms,c->Decorr->decorr_terms[lc],c->Decorr->order,c->Decorr->b[lc],&c->Decorr->chi2val[lc], &d3, d2, c->Decorr->omodel, outname, c->Decorr->zeropointterm, c->Decorr->usemask, c->Decorr->maskvar, lc, lc);
	      c->Decorr->chi2val[lc] /= (i1 - c->Decorr->N_decorrterms_total);
	    }
	  else
	    {
	      for(i2=0;i2<c->Decorr->N_decorrterms_total;i2++)
		{
		  c->Decorr->b[lc][i2] = -1.;
		  c->Decorr->b_err[lc][i2] =  -1.;
		}
	      c->Decorr->chi2val[lc] = -1.;
	    }
	}
      break;

    case CNUM_KILLHARM:
      /* Remove harmonics */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->Killharm->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->Killharm->modeloutdir,&p->lcnames[lc][i2],c->Killharm->modelsuffix);
	    }

	  if(c->Killharm->pertype == PERTYPE_AOV)
	    {
	      c->Killharm->Nper = 1;
	      i1=c->Killharm->lastaovindex;
	      if(c[i1-thisindex].cnum == CNUM_AOV)
		c->Killharm->periods[lc][0] = c[i1-thisindex].Aov->peakperiods[lc][0];
	      else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		c->Killharm->periods[lc][0] = c[i1-thisindex].AovHarm->peakperiods[lc][0];
	    }
	  else if(c->Killharm->pertype == PERTYPE_LS)
	    {
	      c->Killharm->Nper = 1;
	      i1 = c->Killharm->lastlsindex;
	      c->Killharm->periods[lc][0] = c[i1-thisindex].Ls->peakperiods[lc][0];
	    }
	  else if(c->Killharm->pertype == PERTYPE_BOTH)
	    {
	      c->Killharm->Nper = 2;
	      i1=c->Killharm->lastaovindex;
	      i2=c->Killharm->lastlsindex;
	      if(c[i1-thisindex].cnum == CNUM_AOV)
		c->Killharm->periods[lc][0] = c[i1-thisindex].Aov->peakperiods[lc][0];
	      else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		c->Killharm->periods[lc][0] = c[i1-thisindex].AovHarm->peakperiods[lc][0];
	      c->Killharm->periods[lc][1] = c[i2-thisindex].Ls->peakperiods[lc][0];
	    }
	  else if(c->Killharm->pertype == PERTYPE_INJECTHARM)
	    {
	      c->Killharm->Nper = 1;
	      i1 = c->Killharm->lastaovindex;
	      c->Killharm->periods[lc][0] = c[i1-thisindex].Injectharm->periodinject[lc];
	    }
	  else if(c->Killharm->pertype == PERTYPE_FIX)
	    {
	      for(i1 = 0; i1 < c->Killharm->Nper; i1++)
		{
		  c->Killharm->periods[lc][i1] = c->Killharm->fixedperiods[i1];
		}
	    }
	  if(p->NJD[lc] > 1)
	    dokillharms(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->Killharm->Nper,c->Killharm->periods[lc],c->Killharm->Nsubharm,c->Killharm->Nharm,c->Killharm->subharmA[lc],c->Killharm->subharmB[lc],c->Killharm->harmA[lc],c->Killharm->harmB[lc],c->Killharm->fundA[lc],c->Killharm->fundB[lc],&c->Killharm->mean[lc],c->Killharm->omodel,outname,c->Killharm->amp[lc], c->Killharm->fitonly, c->Killharm->outtype, c->Killharm->clip);
	}
      break;

    case CNUM_INJECTHARM:
      /* Inject a harmonic series */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->Injectharm->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->Injectharm->modeloutdir,&p->lcnames[lc][i2],c->Injectharm->modelsuffix);
	    }
	  doinjectharm(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],lc,lc,c->Injectharm,outname);
	}
      break;

    case CNUM_INJECTTRANSIT:
      /* Inject a transit model */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->Injecttransit->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->Injecttransit->modeloutdir,&p->lcnames[lc][i2],c->Injecttransit->modelsuffix);
	    }
	  doinjecttransit(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],lc,lc,c->Injecttransit,outname);
	}
      break;

    case CNUM_STARSPOT:
      /* Fit a starspot model to the light curve and remove it if we're doing that */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->Starspot->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->Starspot->modeloutdir,&p->lcnames[lc][i2],c->Starspot->modelsuffix);
	    }

	  if(c->Starspot->pertype == PERTYPE_AOV)
	    {
	      i1 = c->Starspot->lastaovindex;
	      if(c[i1-thisindex].cnum == CNUM_AOV)
		c->Starspot->period[lc][0] = c[i1-thisindex].Aov->peakperiods[lc][0];
	      else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		c->Starspot->period[lc][0] = c[i1-thisindex].AovHarm->peakperiods[lc][0];
	    }
	  else if(c->Starspot->pertype == PERTYPE_LS)
	    {
	      i1 = c->Starspot->lastlsindex;
	      c->Starspot->period[lc][0] = c[i1-thisindex].Ls->peakperiods[lc][0];
	    }
	  else if(c->Starspot->pertype == PERTYPE_FIXCOLUMN)
	    {
	      getoutcolumnvalue(c->Starspot->linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &(c->Starspot->period[lc][0]));
	    }
	  else if(c->Starspot->pertype == PERTYPE_FIX)
	    {
	      c->Starspot->period[lc][0] = c->Starspot->fixedperiod;
	    }
	  c->Starspot->a[lc] = c->Starspot->a0;
	  c->Starspot->b[lc] = c->Starspot->b0;
	  c->Starspot->alpha[lc] = c->Starspot->alpha0;
	  c->Starspot->inclination[lc] = c->Starspot->inclination0;
	  c->Starspot->chi[lc] = c->Starspot->chi0;
	  c->Starspot->psi0[lc] = c->Starspot->psi00;
	  c->Starspot->mconst[lc] = c->Starspot->mconst0;
	  if(p->NJD[lc] > 1)
	    fitstarspot_amoeba(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&c->Starspot->period[lc][0],&c->Starspot->a[lc],&c->Starspot->b[lc],&c->Starspot->alpha[lc],&c->Starspot->inclination[lc],&c->Starspot->chi[lc],&c->Starspot->psi0[lc],&c->Starspot->mconst[lc],c->Starspot->fitP,c->Starspot->fita,c->Starspot->fitb,c->Starspot->fitalpha,c->Starspot->fiti,c->Starspot->fitchi,c->Starspot->fitpsi0,c->Starspot->fitmconst,&c->Starspot->chisq[lc],c->Starspot->correctlc,c->Starspot->omodel,outname);
	}
      break;

    case CNUM_STATS:
      /* Calculate light curve statistics */
      for(lc=0; lc < p->Nlcs; lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	RunStatsCommand(p, lc, lc, c->Stats);
      }
      break;

    case CNUM_BLS:
      /* Perform BLS on the light curves */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  RunBLSCommand(p, c->Bls, lc, lc, thisindex, 0);
	}
      break;
	  /* First check to see that the u/v vectors are large enough */
	  /*if(p->NJD[lc] > 1)
	    {
	      if(c->Bls->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->Bls->modeloutdir,&p->lcnames[lc][i2],c->Bls->modelsuffix);
		}
	      if(c->Bls->ophcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname3,"%s/%s%s",c->Bls->ophcurveoutdir,&p->lcnames[lc][i2],c->Bls->ophcurvesuffix);
		}
	      if(c->Bls->ojdcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname4,"%s/%s%s",c->Bls->ojdcurveoutdir,&p->lcnames[lc][i2],c->Bls->ojdcurvesuffix);
		}
	      
	      if(c->Bls->sizeuv[0] == 0)
		{
		  c->Bls->sizeuv[0] = p->NJD[lc];
		  if((c->Bls->u[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->Bls->v[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      else if(c->Bls->sizeuv[0] < p->NJD[lc])
		{
		  c->Bls->sizeuv[0] = p->NJD[lc];
		  free(c->Bls->u[0]);
		  free(c->Bls->v[0]);
		  if((c->Bls->u[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->Bls->v[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      
	      if(c->Bls->operiodogram)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->Bls->outdir,&p->lcnames[lc][i2],c->Bls->suffix);
		}

	      c->Bls->fmin[lc] = dmax((2./(p->t[lc][p->NJD[lc]-1] - p->t[lc][0])),1./c->Bls->maxper);
	      if(!c->Bls->freqsteptype) {
		c->Bls->nf2[lc] = floor((((1./c->Bls->minper) - c->Bls->fmin[lc])/c->Bls->df)+1.);
	      } else if(c->Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_PERIOD) {
		c->Bls->nf2[lc] = floor((((1./c->Bls->fmin[lc]) - c->Bls->minper)/c->Bls->df)+1.);
	      } else if(c->Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_LOGPERIOD) {
		c->Bls->nf2[lc] = floor(((log(1./c->Bls->fmin[lc]) - log(c->Bls->minper))/c->Bls->df)+1.);
	      }*/
	      /* Now either run bls using the fixed q range or the fixed stellar radius range */
	      /*if(c->Bls->nf2[lc] > 0 && c->Bls->nbins > 0 && c->Bls->Npeak > 0) {
		if(!c->Bls->rflag)
		  {
		    eebls(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->Bls->u[0],c->Bls->v[0],c->Bls->nf2[lc],c->Bls->fmin[lc],c->Bls->df,c->Bls->nbins,c->Bls->qmin,c->Bls->qmax,
#ifdef PARALLEL
			  c->Bls->p[0]
#else
			  c->Bls->p
#endif
			  ,c->Bls->Npeak,c->Bls->bper[lc],c->Bls->bt0[lc],c->Bls->bpow[lc],c->Bls->sde[lc],c->Bls->snval[lc],c->Bls->depth[lc],c->Bls->qtran[lc],c->Bls->i1[lc],c->Bls->i2[lc],c->Bls->i1_ph[lc],c->Bls->i2_ph[lc],c->Bls->chisqrplus[lc],&c->Bls->chisqrminus[lc],&c->Bls->bperpos[lc],&c->Bls->meanmagval[lc],c->Bls->timezone,c->Bls->fraconenight[lc],c->Bls->operiodogram,outname,c->Bls->omodel,outname2,c->Bls->correctlc,p->ascii, c->Bls->nt[lc], c->Bls->Nt[lc], c->Bls->Nbefore[lc], c->Bls->Nafter[lc], c->Bls->rednoise[lc], c->Bls->whitenoise[lc], c->Bls->sigtopink[lc], c->Bls->fittrap, c->Bls->qingress[lc], c->Bls->OOTmag[lc], c->Bls->ophcurve, outname3, c->Bls->phmin, c->Bls->phmax, c->Bls->phstep, c->Bls->ojdcurve, outname4, c->Bls->jdstep, c->Bls->nobinnedrms,c->Bls->freqsteptype, c->Bls->adjust_qmin_mindt, c->Bls->reduce_nb, c->Bls->reportharmonics, c->Bls, lc, lc, c->Bls->usemask, c->Bls->maskvar);
		  }
		else
		  {
		    eebls_rad(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->Bls->u[0],c->Bls->v[0],c->Bls->nf2[lc],c->Bls->fmin[lc],c->Bls->df,c->Bls->nbins,c->Bls->rmin,c->Bls->rmax,
#ifdef PARALLEL
			      c->Bls->p[0]
#else
			      c->Bls->p
#endif
			      ,c->Bls->Npeak,c->Bls->bper[lc],c->Bls->bt0[lc],c->Bls->bpow[lc],c->Bls->sde[lc],c->Bls->snval[lc],c->Bls->depth[lc],c->Bls->qtran[lc],c->Bls->i1[lc],c->Bls->i2[lc],c->Bls->i1_ph[lc],c->Bls->i2_ph[lc],c->Bls->chisqrplus[lc],&c->Bls->chisqrminus[lc],&c->Bls->bperpos[lc],&c->Bls->meanmagval[lc],c->Bls->timezone,c->Bls->fraconenight[lc],c->Bls->operiodogram,outname,c->Bls->omodel,outname2,c->Bls->correctlc,p->ascii, c->Bls->nt[lc], c->Bls->Nt[lc], c->Bls->Nbefore[lc], c->Bls->Nafter[lc], c->Bls->rednoise[lc], c->Bls->whitenoise[lc], c->Bls->sigtopink[lc], c->Bls->fittrap, c->Bls->qingress[lc], c->Bls->OOTmag[lc], c->Bls->ophcurve, outname3, c->Bls->phmin, c->Bls->phmax, c->Bls->phstep, c->Bls->ojdcurve, outname4, c->Bls->jdstep, c->Bls->nobinnedrms, c->Bls->freqsteptype, c->Bls->adjust_qmin_mindt, c->Bls->reduce_nb, c->Bls->reportharmonics, c->Bls, lc, lc, c->Bls->usemask, c->Bls->maskvar);
		  }
	      } else {
		if(!p->quiet_mode) {
		  fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve is either too short, or an invalid set of parameter options were supplied to BLS.\n", thisindex, lc, p->lcnames[lc]);
		}
	      }
	    } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	}
      break;*/

    case CNUM_FIXPERBLS:
      /* Perform BLS on the light curves */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->BlsFixPer->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->BlsFixPer->modeloutdir,&p->lcnames[lc][i2],c->BlsFixPer->modelsuffix);
		}
	      /* First check to see that the u/v vectors are large enough */
	      if(c->BlsFixPer->sizeuv[0] == 0)
		{
		  c->BlsFixPer->sizeuv[0] = p->NJD[lc];
		  if((c->BlsFixPer->u[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixPer->v[0] = (double *) malloc(c->Bls->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      else if(c->BlsFixPer->sizeuv[0] < p->NJD[lc])
		{
		  c->BlsFixPer->sizeuv[0] = p->NJD[lc];
		  free(c->BlsFixPer->u[0]);
		  free(c->BlsFixPer->v[0]);
		  if((c->BlsFixPer->u[0] = (double *) malloc(c->BlsFixPer->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixPer->v[0] = (double *) malloc(c->BlsFixPer->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      
	      /* Find the period if we're getting it from a previous command */
	      if(c->BlsFixPer->pertype == PERTYPE_AOV)
		{
		  i1 = c->BlsFixPer->lastaovindex;
		  if(c[i1-thisindex].cnum == CNUM_AOV)
		    c->BlsFixPer->period[lc][0] = c[i1-thisindex].Aov->peakperiods[lc][0];
		  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		    c->BlsFixPer->period[lc][0] = c[i1-thisindex].AovHarm->peakperiods[lc][0];
		}
	      else if(c->BlsFixPer->pertype == PERTYPE_LS)
		{
		  i1 = c->BlsFixPer->lastlsindex;
		  c->BlsFixPer->period[lc][0] = c[i1-thisindex].Ls->peakperiods[lc][0];
		}
	      else if(c->BlsFixPer->pertype == PERTYPE_FIX)
		{
		  c->BlsFixPer->period[lc][0] = c->BlsFixPer->perfix;
		}
	      else if(c->BlsFixPer->pertype == PERTYPE_FIXCOLUMN)
		{
		  getoutcolumnvalue(c->BlsFixPer->linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &(c->BlsFixPer->period[lc][0]));
		}
	      else if(c->BlsFixPer->pertype == PERTYPE_EXPR)
		{
		  c->BlsFixPer->period[lc][0] = EvaluateExpression(lc, lc, 0, c->BlsFixPer->perexpr);
		}
	      

	      /* Now either run bls using the fixed q range or the fixed stellar radius range */
	      if(!c->BlsFixPer->rflag)
		{
		  eeblsfixper(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->BlsFixPer->u[0],c->BlsFixPer->v[0],c->BlsFixPer->nbins,c->BlsFixPer->qmin,c->BlsFixPer->qmax,&c->BlsFixPer->period[lc][0],&c->BlsFixPer->bt0[lc],&c->BlsFixPer->bpow[lc],&c->BlsFixPer->depth[lc],&c->BlsFixPer->qtran[lc],&c->BlsFixPer->i1[lc],&c->BlsFixPer->i2[lc],&c->BlsFixPer->i1_ph[lc],&c->BlsFixPer->i2_ph[lc],&c->BlsFixPer->chisqrplus[lc],&c->BlsFixPer->chisqrminus[lc],&c->BlsFixPer->meanmagval[lc], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc], &c->BlsFixPer->Nt[lc], &c->BlsFixPer->Nbefore[lc], &c->BlsFixPer->Nafter[lc], &c->BlsFixPer->rednoise[lc], &c->BlsFixPer->whitenoise[lc], &c->BlsFixPer->sigtopink[lc], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc], &c->BlsFixPer->OOTmag[lc], NULL, lc, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
		}
	      else
		{
		  eeblsfixper_rad(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->BlsFixPer->u[0],c->BlsFixPer->v[0],c->BlsFixPer->nbins,c->BlsFixPer->rmin,c->BlsFixPer->rmax,&c->BlsFixPer->period[lc][0],&c->BlsFixPer->bt0[lc],&c->BlsFixPer->bpow[lc],&c->BlsFixPer->depth[lc],&c->BlsFixPer->qtran[lc],&c->BlsFixPer->i1[lc],&c->BlsFixPer->i2[lc],&c->BlsFixPer->i1_ph[lc],&c->BlsFixPer->i2_ph[lc],&c->BlsFixPer->chisqrplus[lc],&c->BlsFixPer->chisqrminus[lc],&c->BlsFixPer->meanmagval[lc], c->BlsFixPer->timezone, &c->BlsFixPer->fraconenight[lc], c->BlsFixPer->omodel, outname2, c->BlsFixPer->correctlc,p->ascii, &c->BlsFixPer->nt[lc], &c->BlsFixPer->Nt[lc], &c->BlsFixPer->Nbefore[lc], &c->BlsFixPer->Nafter[lc], &c->BlsFixPer->rednoise[lc], &c->BlsFixPer->whitenoise[lc], &c->BlsFixPer->sigtopink[lc], c->BlsFixPer->fittrap, &c->BlsFixPer->qingress[lc], &c->BlsFixPer->OOTmag[lc], NULL, lc, lc, c->BlsFixPer->usemask, c->BlsFixPer->maskvar);
		}
	    } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLSFixPer command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLSFixPer.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	}
      break;

    case CNUM_BLSFIXDURTC:
      /* Perform BLS with fixed transit duration and epoch on the light curves */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->BlsFixDurTc->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->BlsFixDurTc->modeloutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->modelsuffix);
		}
	      if(c->BlsFixDurTc->ophcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname3,"%s/%s%s",c->BlsFixDurTc->ophcurveoutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->ophcurvesuffix);
		}
	      if(c->BlsFixDurTc->ojdcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname4,"%s/%s%s",c->BlsFixDurTc->ojdcurveoutdir,&p->lcnames[lc][i2],c->BlsFixDurTc->ojdcurvesuffix);
		}
	      /* First check to see that the u/v vectors are large enough */
	      if(c->BlsFixDurTc->sizeuv[0] == 0)
		{
		  c->BlsFixDurTc->sizeuv[0] = p->NJD[lc];
		  if((c->BlsFixDurTc->u[0] = (double *) malloc(c->BlsFixDurTc->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixDurTc->v[0] = (double *) malloc(c->BlsFixDurTc->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      else if(c->BlsFixDurTc->sizeuv[0] < p->NJD[lc])
		{
		  c->BlsFixDurTc->sizeuv[0] = p->NJD[lc];
		  free(c->BlsFixDurTc->u[0]);
		  free(c->BlsFixDurTc->v[0]);
		  if((c->BlsFixDurTc->u[0] = (double *) malloc(c->BlsFixDurTc->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixDurTc->v[0] = (double *) malloc(c->BlsFixDurTc->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      
	      if(c->BlsFixDurTc->operiodogram)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->BlsFixDurTc->outdir,&p->lcnames[lc][i2],c->BlsFixDurTc->suffix);
		}
	      if(c->BlsFixDurTc->durtype == PERTYPE_FIX)
		{
		  c->BlsFixDurTc->inputdur[lc] = c->BlsFixDurTc->fixdur;
		  d1 = c->BlsFixDurTc->inputdur[lc];
		}
	      else if(c->BlsFixDurTc->durtype == PERTYPE_FIXCOLUMN) {
		getoutcolumnvalue(c->BlsFixDurTc->fixdur_linkedcolumn, lc, lc, 
				  VARTOOLS_TYPE_DOUBLE, 
				  &(c->BlsFixDurTc->inputdur[lc]));
		d1 = c->BlsFixDurTc->inputdur[lc];
	      } else {
		d1 = c->BlsFixDurTc->inputdur[lc];
	      }
	      if(c->BlsFixDurTc->TCtype == PERTYPE_FIX)
		{
		  c->BlsFixDurTc->inputTC[lc] = c->BlsFixDurTc->fixTC;
		  d2 = c->BlsFixDurTc->fixTC;
		}
	      else if(c->BlsFixDurTc->TCtype == PERTYPE_FIXCOLUMN) {
		getoutcolumnvalue(c->BlsFixDurTc->fixTC_linkedcolumn, lc, lc, 
				  VARTOOLS_TYPE_DOUBLE, 
				  &(c->BlsFixDurTc->inputTC[lc]));
		d2 = c->BlsFixDurTc->inputTC[lc];
	      }
	      else {
		d2 = c->BlsFixDurTc->inputTC[lc];
	      }
	      if(c->BlsFixDurTc->fixdepth) {
		if(c->BlsFixDurTc->depthtype == PERTYPE_FIX) {
		  c->BlsFixDurTc->inputdepth[lc] = c->BlsFixDurTc->fixdepthval;
		  d3 = c->BlsFixDurTc->fixdepthval;
		}
		else if(c->BlsFixDurTc->depthtype == PERTYPE_FIXCOLUMN) {
		  getoutcolumnvalue(c->BlsFixDurTc->fixdepth_linkedcolumn, lc, lc, 
				    VARTOOLS_TYPE_DOUBLE, 
				    &(c->BlsFixDurTc->inputdepth[lc]));
		  d3 = c->BlsFixDurTc->inputdepth[lc];
		}
		else {
		  d3 = c->BlsFixDurTc->inputdepth[lc];
		}
		if(c->BlsFixDurTc->qgresstype == PERTYPE_FIX) {
		  c->BlsFixDurTc->inputqgress[lc] = c->BlsFixDurTc->qgressval;
		  d4 = c->BlsFixDurTc->qgressval;
		}
		else if(c->BlsFixDurTc->qgresstype == PERTYPE_FIXCOLUMN) {
		  getoutcolumnvalue(c->BlsFixDurTc->fixqgress_linkedcolumn, lc, lc, 
				    VARTOOLS_TYPE_DOUBLE, 
				    &(c->BlsFixDurTc->inputqgress[lc]));
		  d4 = c->BlsFixDurTc->inputqgress[lc];
		}
		else {
		  d4 = c->BlsFixDurTc->inputqgress[lc];
		}
	      }
	      c->BlsFixDurTc->fmin[lc] = dmax((2./(p->t[lc][p->NJD[lc]-1] - p->t[lc][0])),1./c->BlsFixDurTc->maxper);
	      c->BlsFixDurTc->nf2[lc] = floor((((1./c->BlsFixDurTc->minper) - c->BlsFixDurTc->fmin[lc])/c->BlsFixDurTc->df)+1.);
	      if(c->BlsFixDurTc->nf2[lc] > 0 && c->BlsFixDurTc->Npeak > 0) {
		/* Now either run bls  */
		eeblsfixdurtc(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->BlsFixDurTc->u[lc],c->BlsFixDurTc->v[lc],d2,d1,c->BlsFixDurTc->fixdepth,d3,d4,c->BlsFixDurTc->nf2[lc],c->BlsFixDurTc->fmin[lc],c->BlsFixDurTc->df,
#ifdef PARALLEL
			      c->BlsFixDurTc->p[lc]
#else
			      c->BlsFixDurTc->p
#endif
			      ,c->BlsFixDurTc->Npeak,c->BlsFixDurTc->bper[lc],c->BlsFixDurTc->bt0[lc],c->BlsFixDurTc->bpow[lc],c->BlsFixDurTc->sde[lc],c->BlsFixDurTc->snval[lc],c->BlsFixDurTc->depth[lc],c->BlsFixDurTc->qtran[lc],c->BlsFixDurTc->chisqrplus[lc],&c->BlsFixDurTc->chisqrminus[lc],&c->BlsFixDurTc->bperpos[lc],&c->BlsFixDurTc->meanmagval[lc], c->BlsFixDurTc->timezone, c->BlsFixDurTc->fraconenight[lc], c->BlsFixDurTc->operiodogram, outname, c->BlsFixDurTc->omodel, outname2, c->BlsFixDurTc->correctlc,p->ascii, c->BlsFixDurTc->nt[lc], c->BlsFixDurTc->Nt[lc], c->BlsFixDurTc->Nbefore[lc], c->BlsFixDurTc->Nafter[lc], c->BlsFixDurTc->rednoise[lc], c->BlsFixDurTc->whitenoise[lc], c->BlsFixDurTc->sigtopink[lc], c->BlsFixDurTc->fittrap, c->BlsFixDurTc->qingress[lc], c->BlsFixDurTc->OOTmag[lc], c->BlsFixDurTc->ophcurve, outname3, c->BlsFixDurTc->phmin, c->BlsFixDurTc->phmax, c->BlsFixDurTc->phstep, c->BlsFixDurTc->ojdcurve, outname4, c->BlsFixDurTc->jdstep, lc, lc, c->BlsFixDurTc->usemask, c->BlsFixDurTc->maskvar);
	      } else {
		if(!p->quiet_mode) {
		  fprintf(stderr,"Warning: skipping -BLSFixDurTc command index %d for light curve number: %d, filename: %s. The light curve is either too short, or an invalid set of parameter options were supplied to BLSFixDurTc.\n", thisindex, lc, p->lcnames[lc]);
		}
	      }
	    } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLSFixDurTc command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLSFixDurTc.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	}
      break;

    case CNUM_BLSFIXPERDURTC:
      /* Perform BLS with fixed transit duration and epoch on the light curves */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->BlsFixPerDurTc->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->BlsFixPerDurTc->modeloutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->modelsuffix);
		}
	      if(c->BlsFixPerDurTc->ophcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname3,"%s/%s%s",c->BlsFixPerDurTc->ophcurveoutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->ophcurvesuffix);
		}
	      if(c->BlsFixPerDurTc->ojdcurve)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname4,"%s/%s%s",c->BlsFixPerDurTc->ojdcurveoutdir,&p->lcnames[lc][i2],c->BlsFixPerDurTc->ojdcurvesuffix);
		}
	      /* First check to see that the u/v vectors are large enough */
	      if(c->BlsFixPerDurTc->sizeuv[0] == 0)
		{
		  c->BlsFixPerDurTc->sizeuv[0] = p->NJD[lc];
		  if((c->BlsFixPerDurTc->u[0] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixPerDurTc->v[0] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      else if(c->BlsFixPerDurTc->sizeuv[0] < p->NJD[lc])
		{
		  c->BlsFixPerDurTc->sizeuv[0] = p->NJD[lc];
		  free(c->BlsFixPerDurTc->u[0]);
		  free(c->BlsFixPerDurTc->v[0]);
		  if((c->BlsFixPerDurTc->u[0] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[0] * sizeof(double))) == NULL ||
		     (c->BlsFixPerDurTc->v[0] = (double *) malloc(c->BlsFixPerDurTc->sizeuv[0] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      
	      if(c->BlsFixPerDurTc->pertype == PERTYPE_FIX)
		{
		  c->BlsFixPerDurTc->inputper[lc] = c->BlsFixPerDurTc->fixper;
		  d0 = c->BlsFixPerDurTc->inputper[lc];
		}
	      else if(c->BlsFixPerDurTc->pertype == PERTYPE_FIXCOLUMN) {
		getoutcolumnvalue(c->BlsFixPerDurTc->fixper_linkedcolumn, lc, lc, 
				  VARTOOLS_TYPE_DOUBLE, 
				  &(c->BlsFixPerDurTc->inputper[lc]));
		d0 = c->BlsFixPerDurTc->inputper[lc];
	      } else {
		d0 = c->BlsFixPerDurTc->inputper[lc];
	      }
	      if(c->BlsFixPerDurTc->durtype == PERTYPE_FIX)
		{
		  c->BlsFixPerDurTc->inputdur[lc] = c->BlsFixPerDurTc->fixdur;
		  d1 = c->BlsFixPerDurTc->inputdur[lc];
		}
	      else if(c->BlsFixPerDurTc->durtype == PERTYPE_FIXCOLUMN) {
		getoutcolumnvalue(c->BlsFixPerDurTc->fixdur_linkedcolumn, lc, lc, 
				  VARTOOLS_TYPE_DOUBLE, 
				  &(c->BlsFixPerDurTc->inputdur[lc]));
		d1 = c->BlsFixPerDurTc->inputdur[lc];
	      } else {
		d1 = c->BlsFixPerDurTc->inputdur[lc];
	      }
	      if(c->BlsFixPerDurTc->TCtype == PERTYPE_FIX)
		{
		  c->BlsFixPerDurTc->inputTC[lc] = c->BlsFixPerDurTc->fixTC;
		  d2 = c->BlsFixPerDurTc->fixTC;
		}
	      else if(c->BlsFixPerDurTc->TCtype == PERTYPE_FIXCOLUMN) {
		getoutcolumnvalue(c->BlsFixPerDurTc->fixTC_linkedcolumn, lc, lc, 
				  VARTOOLS_TYPE_DOUBLE, 
				  &(c->BlsFixPerDurTc->inputTC[lc]));
		d2 = c->BlsFixPerDurTc->inputTC[lc];
	      }
	      else {
		d2 = c->BlsFixPerDurTc->inputTC[lc];
	      }
	      if(c->BlsFixPerDurTc->fixdepth) {
		if(c->BlsFixPerDurTc->depthtype == PERTYPE_FIX) {
		  c->BlsFixPerDurTc->inputdepth[lc] = c->BlsFixPerDurTc->fixdepthval;
		  d3 = c->BlsFixPerDurTc->fixdepthval;
		}
		else if(c->BlsFixPerDurTc->depthtype == PERTYPE_FIXCOLUMN) {
		  getoutcolumnvalue(c->BlsFixPerDurTc->fixdepth_linkedcolumn, lc, lc, 
				    VARTOOLS_TYPE_DOUBLE, 
				    &(c->BlsFixPerDurTc->inputdepth[lc]));
		  d3 = c->BlsFixPerDurTc->inputdepth[lc];
		}
		else {
		  d3 = c->BlsFixPerDurTc->inputdepth[lc];
		}
		if(c->BlsFixPerDurTc->qgresstype == PERTYPE_FIX) {
		  c->BlsFixPerDurTc->inputqgress[lc] = c->BlsFixPerDurTc->qgressval;
		  d4 = c->BlsFixPerDurTc->qgressval;
		}
		else if(c->BlsFixPerDurTc->qgresstype == PERTYPE_FIXCOLUMN) {
		  getoutcolumnvalue(c->BlsFixPerDurTc->fixqgress_linkedcolumn, lc, lc, 
				    VARTOOLS_TYPE_DOUBLE, 
				    &(c->BlsFixPerDurTc->inputqgress[lc]));
		  d4 = c->BlsFixPerDurTc->inputqgress[lc];
		}
		else {
		  d4 = c->BlsFixPerDurTc->inputqgress[lc];
		}
	      }
	      eeblsfixperdurtc(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->BlsFixPerDurTc->u[lc],c->BlsFixPerDurTc->v[lc],d0,d2,d1,c->BlsFixPerDurTc->fixdepth,d3,d4,&(c->BlsFixPerDurTc->depth[lc]),&(c->BlsFixPerDurTc->qtran[lc]),&(c->BlsFixPerDurTc->chisqrplus[lc]),&c->BlsFixPerDurTc->meanmagval[lc], c->BlsFixPerDurTc->timezone, &(c->BlsFixPerDurTc->fraconenight[lc]), c->BlsFixPerDurTc->omodel, outname2, c->BlsFixPerDurTc->correctlc, &(c->BlsFixPerDurTc->nt[lc]), &(c->BlsFixPerDurTc->Nt[lc]), &(c->BlsFixPerDurTc->Nbefore[lc]), &(c->BlsFixPerDurTc->Nafter[lc]), &(c->BlsFixPerDurTc->rednoise[lc]), &(c->BlsFixPerDurTc->whitenoise[lc]), &(c->BlsFixPerDurTc->sigtopink[lc]), c->BlsFixPerDurTc->fittrap, &(c->BlsFixPerDurTc->qingress[lc]), &(c->BlsFixPerDurTc->OOTmag[lc]), c->BlsFixPerDurTc->ophcurve, outname3, c->BlsFixPerDurTc->phmin, c->BlsFixPerDurTc->phmax, c->BlsFixPerDurTc->phstep, c->BlsFixPerDurTc->ojdcurve, outname4, c->BlsFixPerDurTc->jdstep, lc, lc, c->BlsFixPerDurTc->usemask, c->BlsFixPerDurTc->maskvar);
	    } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLSFixPerDurTc command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLSFixPerDurTc.\n", thisindex, lc, p->lcnames[lc]);
	    }
	  }
	}
      break;

    case CNUM_SOFTENEDTRANSIT:
      /* Fit a Softened Transit model to the light curve and remove it if we're doing that */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->SoftenedTransit->omodel)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname,"%s/%s%s",c->SoftenedTransit->modeloutdir,&p->lcnames[lc][i2],c->SoftenedTransit->modelsuffix);
	    }
	  if(c->SoftenedTransit->frombls)
	    {
	      /* Get the starting parameters from BLS if that's what we're doing */
	      i1 = c->SoftenedTransit->lastblsindex;
	      c->SoftenedTransit->period0 = c[i1 - thisindex].Bls->bper[lc][0];
	      c->SoftenedTransit->T00 = c[i1 - thisindex].Bls->bt0[lc][0];
	      /*
	      if(c[i1 - thisindex].Bls->i2[lc][0] > c[i1 - thisindex].Bls->i1[lc][0])
		c->SoftenedTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].Bls->i1[lc][0] + c[i1 - thisindex].Bls->i2[lc][0])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].Bls->nbins);
	      else
		c->SoftenedTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].Bls->i1[lc][0] + 1 - c[i1 - thisindex].Bls->i2[lc][0])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].Bls->nbins);
	      */
	      c->SoftenedTransit->eta0 = c[i1 - thisindex].Bls->qtran[lc][0];
	      c->SoftenedTransit->delta0 = c[i1 - thisindex].Bls->depth[lc][0];
	      c->SoftenedTransit->mconst0 = -1;
	    }
	  else if(c->SoftenedTransit->fromblsfixper)
	    {
	      /* Get the starting parameters from BLSFixPer if that's what we're doing */
	      i1 = c->SoftenedTransit->lastblsfixperindex;
	      c->SoftenedTransit->period0 = c[i1 - thisindex].BlsFixPer->period[lc][0];
	      c->SoftenedTransit->T00 = c[i1 - thisindex].BlsFixPer->bt0[lc];
	      /*
	      if(c[i1 - thisindex].BlsFixPer->i2[lc] > c[i1 - thisindex].BlsFixPer->i1[lc])
		c->SoftenedTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc] + c[i1 - thisindex].BlsFixPer->i2[lc])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	      else
		c->SoftenedTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc] + 1 - c[i1 - thisindex].BlsFixPer->i2[lc])*0.5*c->SoftenedTransit->period0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	      */
	      c->SoftenedTransit->eta0 = c[i1 - thisindex].BlsFixPer->qtran[lc];
	      c->SoftenedTransit->delta0 = c[i1 - thisindex].BlsFixPer->depth[lc];
	      c->SoftenedTransit->mconst0 = -1;
	    }
	  if(c->SoftenedTransit->dokillharm)
	    {
	      /* Get the parameters for removing the harmonic function */
	      if(c->SoftenedTransit->pertype == PERTYPE_AOV)
		{
		  i1 = c->SoftenedTransit->lastaovindex;
		  if(c[i1-thisindex].cnum == CNUM_AOV)
		    c->SoftenedTransit->per_harm_out[lc] = c[i1 - thisindex].Aov->peakperiods[lc][0];
		  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
		    c->SoftenedTransit->per_harm_out[lc] = c[i1 - thisindex].AovHarm->peakperiods[lc][0];
		}
	      else if(c->SoftenedTransit->pertype == PERTYPE_LS)
		{
		  i1 = c->SoftenedTransit->lastlsindex;
		  c->SoftenedTransit->per_harm_out[lc] = c[i1 - thisindex].Ls->peakperiods[lc][0];
		}
	      else if(c->SoftenedTransit->pertype == PERTYPE_BLS)
		{
		  i1 = c->SoftenedTransit->lastblsindex;
		  c->SoftenedTransit->per_harm_out[lc] = c[i1 - thisindex].Bls->bper[lc][0];
		}
	      else if(c->SoftenedTransit->pertype == PERTYPE_FIX)
		{
		  c->SoftenedTransit->per_harm_out[lc] = c->SoftenedTransit->per_harm;
		}
	      else if(c->SoftenedTransit->pertype == PERTYPE_SPECIFIED)
		{
		  c->SoftenedTransit->per_harm_out[lc] = c->SoftenedTransit->per_harm_spec[lc];
		}
	    }
	  c->SoftenedTransit->period[lc] = c->SoftenedTransit->period0;
	  c->SoftenedTransit->T0[lc] = c->SoftenedTransit->T00;
	  c->SoftenedTransit->eta[lc] = c->SoftenedTransit->eta0;
	  c->SoftenedTransit->cval[lc] = c->SoftenedTransit->cval0;
	  c->SoftenedTransit->delta[lc] = c->SoftenedTransit->delta0;
	  c->SoftenedTransit->mconst[lc] = c->SoftenedTransit->mconst0;
	  if(p->NJD[lc] > 1)
	    {
	      if(c->SoftenedTransit->dokillharm)
		fitsoftened_transit(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&c->SoftenedTransit->period[lc],&c->SoftenedTransit->T0[lc],&c->SoftenedTransit->eta[lc],&c->SoftenedTransit->cval[lc],&c->SoftenedTransit->delta[lc],&c->SoftenedTransit->mconst[lc],c->SoftenedTransit->fitephem,c->SoftenedTransit->fiteta,c->SoftenedTransit->fitcval,c->SoftenedTransit->fitdelta,c->SoftenedTransit->fitmconst,&c->SoftenedTransit->chisq[lc],c->SoftenedTransit->correctlc,c->SoftenedTransit->omodel,outname,c->SoftenedTransit->dokillharm,c->SoftenedTransit->nharm,c->SoftenedTransit->nsubharm,c->SoftenedTransit->per_harm_out[lc],c->SoftenedTransit->subharmA[lc],c->SoftenedTransit->subharmB[lc],c->SoftenedTransit->harmA[lc],c->SoftenedTransit->harmB[lc],&c->SoftenedTransit->fundA[lc],&c->SoftenedTransit->fundB[lc]);
	      else
		fitsoftened_transit(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&c->SoftenedTransit->period[lc],&c->SoftenedTransit->T0[lc],&c->SoftenedTransit->eta[lc],&c->SoftenedTransit->cval[lc],&c->SoftenedTransit->delta[lc],&c->SoftenedTransit->mconst[lc],c->SoftenedTransit->fitephem,c->SoftenedTransit->fiteta,c->SoftenedTransit->fitcval,c->SoftenedTransit->fitdelta,c->SoftenedTransit->fitmconst,&c->SoftenedTransit->chisq[lc],c->SoftenedTransit->correctlc,c->SoftenedTransit->omodel,outname,c->SoftenedTransit->dokillharm,0,0,0.,NULL,NULL,NULL,NULL,NULL,NULL);
	    }
	}
      break;

    case CNUM_MANDELAGOLTRANSIT:
      /* Fit a Mandel and Agol (2002) Transit model to the light curve and remove it if we're doing that */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->MandelAgolTransit->omodel)
	    {
	      if(!strncmp(c->MandelAgolTransit->modeloutdir,"-",1) && strlen(c->MandelAgolTransit->modeloutdir) == 1)
		{
		  outname[0] = '-';
		  outname[1] = '\0';
		}
	      else
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->MandelAgolTransit->modeloutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->modelsuffix);
		}
	    }
	  if(c->MandelAgolTransit->ophcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname3,"%s/%s%s",c->MandelAgolTransit->ophcurveoutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->ophcurvesuffix);
	    }
	  if(c->MandelAgolTransit->ojdcurve)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      sprintf(outname4,"%s/%s%s",c->MandelAgolTransit->ojdcurveoutdir,&p->lcnames[lc][i2],c->MandelAgolTransit->ojdcurvesuffix);
	    }
	  if(c->MandelAgolTransit->frombls)
	    {
	      /* Get the starting parameters from BLS if that's what we're doing */
	      i1 = c->MandelAgolTransit->lastblsindex;
	      c->MandelAgolTransit->P0 = c[i1 - thisindex].Bls->bper[lc][0];
	      c->MandelAgolTransit->T00 = c[i1 - thisindex].Bls->bt0[lc][0];
	      /*
	      if(c[i1 - thisindex].Bls->i2[lc][0] > c[i1 - thisindex].Bls->i1[lc][0])
		c->MandelAgolTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].Bls->i1[lc][0] + c[i1 - thisindex].Bls->i2[lc][0])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].Bls->nbins);
	      else
		c->MandelAgolTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].Bls->i1[lc][0] + 1 - c[i1 - thisindex].Bls->i2[lc][0])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].Bls->nbins);
	      */
	      c->MandelAgolTransit->a0 = 1./c[i1 - thisindex].Bls->qtran[lc][0]/M_PI;
	      c->MandelAgolTransit->r0 = sqrt(c[i1 - thisindex].Bls->depth[lc][0]);
	      c->MandelAgolTransit->mconst0 = -1;
	      c->MandelAgolTransit->inputinclterm = 1;
	      c->MandelAgolTransit->bimpact0 = 0.2;
	      //c->MandelAgolTransit->sin_i0 = 0.5*(1.0 + getminsini(c->MandelAgolTransit->a0, 0., 0., c->MandelAgolTransit->r0));
	      c->MandelAgolTransit->e0 = 0.;
	      c->MandelAgolTransit->omega0 = 0.;
	    }
	  else if(c->MandelAgolTransit->fromblsfixper)
	    {
	      /* Get the starting parameters from BLSFixPer if that's what we're doing */
	      /* Get the starting parameters from BLS if that's what we're doing */
	      i1 = c->MandelAgolTransit->lastblsfixperindex;
	      c->MandelAgolTransit->P0 = c[i1 - thisindex].BlsFixPer->period[lc][0];
	      c->MandelAgolTransit->T00 = c[i1 - thisindex].BlsFixPer->bt0[lc];
	      /*
	      if(c[i1 - thisindex].BlsFixPer->i2[lc] > c[i1 - thisindex].BlsFixPer->i1[lc])
		c->MandelAgolTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc] + c[i1 - thisindex].BlsFixPer->i2[lc])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	      else
		c->MandelAgolTransit->T00 = p->t[lc][0] + (double)(c[i1 - thisindex].BlsFixPer->i1[lc] + 1 - c[i1 - thisindex].BlsFixPer->i2[lc])*0.5*c->MandelAgolTransit->P0 / (double)(c[i1 - thisindex].BlsFixPer->nbins);
	      */
	      c->MandelAgolTransit->a0 = 1./c[i1 - thisindex].BlsFixPer->qtran[lc]/M_PI;
	      c->MandelAgolTransit->r0 = sqrt(c[i1 - thisindex].BlsFixPer->depth[lc]);
	      c->MandelAgolTransit->mconst0 = -1;
	      c->MandelAgolTransit->bimpact0 = 0.2;
	      c->MandelAgolTransit->inputinclterm = 1;
	      //c->MandelAgolTransit->sin_i0 = 0.5*(1.0 + getminsini(c->MandelAgolTransit->a0, 0., 0., c->MandelAgolTransit->r0));
	      c->MandelAgolTransit->e0 = 0.;
	      c->MandelAgolTransit->omega0 = 0.;
	    }
	  c->MandelAgolTransit->period[lc] = c->MandelAgolTransit->P0;
	  c->MandelAgolTransit->T0[lc] = c->MandelAgolTransit->T00;
	  c->MandelAgolTransit->r[lc] = c->MandelAgolTransit->r0;
	  c->MandelAgolTransit->a[lc] = c->MandelAgolTransit->a0;
	  if(c->MandelAgolTransit->inputinclterm) {
	    c->MandelAgolTransit->bimpact[lc] = c->MandelAgolTransit->bimpact0;
	    c->MandelAgolTransit->inc[lc] = c->MandelAgolTransit->bimpact0*(1. + c->MandelAgolTransit->e0*cos(c->MandelAgolTransit->omega0))/(1. - c->MandelAgolTransit->e0*c->MandelAgolTransit->e0)/(c->MandelAgolTransit->a0);
	    c->MandelAgolTransit->inc[lc] = 180.*acos(c->MandelAgolTransit->inc[lc])/M_PI;
	  }
	  else {
	    c->MandelAgolTransit->inc[lc] = c->MandelAgolTransit->inc0;
	    c->MandelAgolTransit->bimpact[lc] = cos(c->MandelAgolTransit->inc0*M_PI/180.)*(1. - c->MandelAgolTransit->e0*c->MandelAgolTransit->e0)*(c->MandelAgolTransit->a0)/(1. + c->MandelAgolTransit->e0*cos(c->MandelAgolTransit->omega0));
	  }
	  //c->MandelAgolTransit->sin_i[lc] = c->MandelAgolTransit->sin_i0;
	  c->MandelAgolTransit->e[lc] = c->MandelAgolTransit->e0;
	  c->MandelAgolTransit->omega[lc] = c->MandelAgolTransit->omega0;
	  c->MandelAgolTransit->ldcoeffs[lc][0] = c->MandelAgolTransit->ldcoeffs0[0];
	  c->MandelAgolTransit->ldcoeffs[lc][1] = c->MandelAgolTransit->ldcoeffs0[1];
	  c->MandelAgolTransit->ldcoeffs[lc][2] = c->MandelAgolTransit->ldcoeffs0[2];
	  c->MandelAgolTransit->ldcoeffs[lc][3] = c->MandelAgolTransit->ldcoeffs0[3];
	  c->MandelAgolTransit->gamma[lc] = c->MandelAgolTransit->gamma0;
	  c->MandelAgolTransit->K[lc] = c->MandelAgolTransit->K0;
	  c->MandelAgolTransit->mconst[lc] = c->MandelAgolTransit->mconst0;
	  if(p->NJD[lc] > 1)
	    {
	      fitmandelagoltransit_amoeba(p, p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],&c->MandelAgolTransit->period[lc],&c->MandelAgolTransit->T0[lc],&c->MandelAgolTransit->r[lc],&c->MandelAgolTransit->a[lc],&c->MandelAgolTransit->inc[lc],&c->MandelAgolTransit->bimpact[lc],&c->MandelAgolTransit->e[lc],&c->MandelAgolTransit->omega[lc],&c->MandelAgolTransit->mconst[lc],c->MandelAgolTransit->type,c->MandelAgolTransit->ldcoeffs[lc],c->MandelAgolTransit->fitephem,c->MandelAgolTransit->fitr,c->MandelAgolTransit->fita,c->MandelAgolTransit->fitinclterm,c->MandelAgolTransit->fite,c->MandelAgolTransit->fitomega,c->MandelAgolTransit->fitmconst,c->MandelAgolTransit->fitldcoeffs, &c->MandelAgolTransit->chisq[lc],c->MandelAgolTransit->correctlc,c->MandelAgolTransit->omodel,outname, c->MandelAgolTransit->fitRV, c->MandelAgolTransit->RVinputfile, c->MandelAgolTransit->RVmodeloutfile, &c->MandelAgolTransit->K[lc], &c->MandelAgolTransit->gamma[lc], c->MandelAgolTransit->fitK, c->MandelAgolTransit->fitgamma, c->MandelAgolTransit->refititer, c->MandelAgolTransit->ophcurve, outname3, c->MandelAgolTransit->phmin, c->MandelAgolTransit->phmax, c->MandelAgolTransit->phstep, c->MandelAgolTransit->ojdcurve, outname4, c->MandelAgolTransit->jdstep, c->MandelAgolTransit->modelvarname, c->MandelAgolTransit->modelvar, lc);
	    }
	}
      break;

    case CNUM_MICROLENS:
      /* Fit a Microlens model to the light curve */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(c->MicroLens->omodel)
	    {
	      if(!strncmp(c->MicroLens->modeloutdir,"-",1) && strlen(c->MicroLens->modeloutdir) == 1)
		{
		  outname[0] = '-';
		  outname[1] = '\0';
		}
	      else
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->MicroLens->modeloutdir,&p->lcnames[lc][i2],c->MicroLens->modelsuffix);
		}
	    }
	  microlens(p->t[lc],p->mag[lc],p->sig[lc],p->NJD[lc],lc,c->MicroLens,outname,&c->MicroLens->f0[lc],&c->MicroLens->f1[lc],&c->MicroLens->u0[lc],&c->MicroLens->t0[lc],&c->MicroLens->tmax[lc],&c->MicroLens->chi2_[lc]);
	}
      break;


    case CNUM_SYSREM:
      /* Run the SYSREM detrending algorithm */
      do_sysrem(p, c->Sysrem, p->Nlcs, p->NJD, p->t, p->mag, p->sig, p->lcnames,p->matchstringid,p->stringid,p->stringid_idx);
      break;

    case CNUM_TFA:
      /* Run the trend filtering algorithm */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->TFA->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->TFA->model_outdir,&p->lcnames[lc][i2],c->TFA->model_suffix);
		}
	      if(c->TFA->ocoeff)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->TFA->coeff_outdir,&p->lcnames[lc][i2],c->TFA->coeff_suffix);
		}
	      if(p->matchstringid)
		detrend_tfa(p, c->TFA, p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->TFA->lcx[lc][0], c->TFA->lcy[lc][0], p->lcnames[lc], outname2, c->TFA->ocoeff, c->TFA->correctlc, c->TFA->omodel, outname, &c->TFA->ave_out[lc], &c->TFA->rms_out[lc],p->matchstringid,p->stringid[lc],p->stringid_idx[lc], 0);
	      else
		detrend_tfa(p, c->TFA, p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->TFA->lcx[lc][0], c->TFA->lcy[lc][0], p->lcnames[lc], outname2, c->TFA->ocoeff, c->TFA->correctlc, c->TFA->omodel, outname, &c->TFA->ave_out[lc], &c->TFA->rms_out[lc],0,NULL,NULL, 0);
	    }
	}
      break;

    case CNUM_TFA_SR:
      /* Run the trend filtering algorithm in SR mode */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[0] > 1)
	    {
	      if(c->TFA_SR->omodel)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname,"%s/%s%s",c->TFA_SR->model_outdir,&p->lcnames[lc][i2],c->TFA_SR->model_suffix);
		}
	      if(c->TFA_SR->ocoeff)
		{
		  i1 = 0;
		  i2 = 0;
		  while(p->lcnames[lc][i1] != '\0')
		    {
		      if(p->lcnames[lc][i1] == '/')
			i2 = i1 + 1;
		      i1++;
		    }
		  sprintf(outname2,"%s/%s%s",c->TFA_SR->coeff_outdir,&p->lcnames[lc][i2],c->TFA_SR->coeff_suffix);
		}
	      if((c->TFA_SR->use_bin && c->TFA_SR->use_period) || (c->TFA_SR->use_harm && c->TFA_SR->use_period))
		{
		  if(c->TFA_SR->pertype == PERTYPE_AOV)
		    {
		      i1=c->TFA_SR->lastindex;
		      if(c[i1-thisindex].cnum == CNUM_AOV)
			c->TFA_SR->periods[lc][0] = c[i1-thisindex].Aov->peakperiods[lc][0];
		      else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
			c->TFA_SR->periods[lc][0] = c[i1-thisindex].AovHarm->peakperiods[lc][0];
		      d1 = c->TFA_SR->periods[lc][0];
		    }
		  else if(c->TFA_SR->pertype == PERTYPE_LS)
		    {
		      i1 = c->TFA_SR->lastindex;
		      c->TFA_SR->periods[lc][0] = c[i1-thisindex].Ls->peakperiods[lc][0];
		      d1 = c->TFA_SR->periods[lc][0];
		    }
		  else if(c->TFA_SR->pertype == PERTYPE_BLS)
		    {
		      i1 = c->TFA_SR->lastindex;
		      c->TFA_SR->periods[lc][0] = c[i1-thisindex].Bls->bper[lc][0];
		      d1 = c->TFA_SR->periods[lc][0];
		    }
		  else if(c->TFA_SR->pertype == PERTYPE_SPECIFIED)
		    {
		      d1 = c->TFA_SR->periods[lc][0];
		    }
		  else if(c->TFA_SR->pertype == PERTYPE_FIX)
		    {
		      d1 = c->TFA_SR->fixperiod;
		    }
		}
	      else if(c->TFA_SR->use_harm)
		d1 = p->t[lc][p->NJD[lc]-1] - p->t[lc][0];
	      else
		d1 = 1.0;
	      if(!c->TFA_SR->use_bin && !c->TFA_SR->use_harm)
		s1 = c->TFA_SR->signalfilenames[lc];
	      else
		s1 = NULL;
	      if(p->matchstringid)
		detrend_tfa_sr(p, c->TFA_SR, p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->TFA_SR->lcx[lc][0], c->TFA_SR->lcy[lc][0], p->lcnames[lc], outname2, c->TFA_SR->ocoeff, c->TFA_SR->correctlc, c->TFA_SR->omodel, outname, &c->TFA_SR->ave_out[lc], &c->TFA_SR->rms_out[lc], d1, s1,p->matchstringid,p->stringid[lc],p->stringid_idx[lc], lc, 0);
	      else
		detrend_tfa_sr(p, c->TFA_SR, p->NJD[lc], p->t[lc], p->mag[lc], p->sig[lc], c->TFA_SR->lcx[lc][0], c->TFA_SR->lcy[lc][0], p->lcnames[lc], outname2, c->TFA_SR->ocoeff, c->TFA_SR->correctlc, c->TFA_SR->omodel, outname, &c->TFA_SR->ave_out[lc], &c->TFA_SR->rms_out[lc], d1, s1,0,NULL,NULL, lc, 0);
	    }
	}
      break;

    case CNUM_DFTCLEAN:
      /* Run DFT/CLEAN on the light curves */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      dodftclean(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],lc,c->Dftclean,&p->lcnames[lc][i2],p->ascii, lc, lc, c->Dftclean->usemask, c->Dftclean->maskvar);
	    }
	}
      break;

    case CNUM_BINLC:
      /* Bin the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	binlc(p, c->Binlc, lc, lc);
      }
      break;

    case CNUM_MATCHCOMMAND:
      /* Match the light curve to an external datafile */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	RunMatchCommand(p, c->MatchCommand, lc, lc);
      }
      break;

    case CNUM_MEDIANFILTER:
      /* median filter the light curve */
      for(lc=0;lc<p->Nlcs;lc++) {
	if(p->isifcommands) {
	  if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	    SkipCommand(p, c, thisindex, lc, lc);
	    continue;
	  }
	}
	medianfilter(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],c->MedianFilter->time,c->MedianFilter->usemean,c->MedianFilter->replace);
      }
      break;

    case CNUM_PHASE:
      /* Replace the time coordinate of a light curve with its phase */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  if(p->NJD[lc] > 1)
	    {
	      if(c->Phase->pertype == PERTYPE_AOV)
		{
		  i1 = c->Phase->lastaovindex;
		  if(c[i1-thisindex].cnum == CNUM_AOV) {
		    d1 = c[i1-thisindex].Aov->peakperiods[lc][0];
		    c->Phase->period[lc][0] = d1;
		  }
		  else if(c[i1-thisindex].cnum == CNUM_HARMAOV) {
		    d1 = c[i1-thisindex].AovHarm->peakperiods[lc][0];
		    c->Phase->period[lc][0] = d1;
		  }
		}
	      else if(c->Phase->pertype == PERTYPE_LS)
		{
		  i1 = c->Phase->lastlsindex;
		  d1 = c[i1-thisindex].Ls->peakperiods[lc][0];
		  c->Phase->period[lc][0] = d1;
		}
	      else if(c->Phase->pertype == PERTYPE_BLS)
		{
		  i1 = c->Phase->lastblsindex;
		  d1 = c[i1-thisindex].Bls->bper[lc][0];
		  c->Phase->period[lc][0] = d1;
		}
	      else if(c->Phase->pertype == PERTYPE_SPECIFIED)
		{
		  d1 = c->Phase->period[lc][0];
		}
	      else if(c->Phase->pertype == PERTYPE_FIX)
		{
		  d1 = c->Phase->fixperiod;
		  c->Phase->period[lc][0] = d1;
		}
	      else if(c->Phase->pertype == PERTYPE_FIXCOLUMN)
		{
		  getoutcolumnvalue(c->Phase->period_linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &d1);
		  c->Phase->period[lc][0] = d1;
		}
	      if(c->Phase->t0type == PERTYPE_AUTOFIND) {
		d2 = 0.; i2 = 0;
	      }
	      else {
		i2 = 1;
		if(c->Phase->t0type == PERTYPE_BLS) {
		  i1 = c->Phase->lastblsindex;
		  /* Get Tc for the transit */
		  if(c[i1-thisindex].Bls->i1[lc][0] > c[i1-thisindex].Bls->i2[lc][0]) {
		    d2 = p->t[lc][0] + c[i1-thisindex].Bls->bper[lc][0]*0.5*(c[i1-thisindex].Bls->i1[lc][0]+1.+c[i1-thisindex].Bls->i2[lc][0]+(1./((double)c[i1-thisindex].Bls->nbins_val[lc])));
		  }
		  else {
		    d1 = p->t[lc][0] + c[i1-thisindex].Bls->bper[lc][0]*0.5*(c[i1-thisindex].Bls->i1[lc][0]+c[i1-thisindex].Bls->i2[lc][0]+(1./((double)c[i1-thisindex].Bls->nbins_val[lc])));
		  }
		  /* adjust so that Tc has phase phaseTc */
		  d2 = d2 - c->Phase->phaseTc*d1;
		}
		else if(c->Phase->t0type == PERTYPE_SPECIFIED) {
		  d2 = c->Phase->T0[lc][0];
		}
		else if(c->Phase->t0type == PERTYPE_FIX) {
		  d2 = c->Phase->fixT0;
		}
		else if(c->Phase->t0type == PERTYPE_FIXCOLUMN) {
		  getoutcolumnvalue(c->Phase->T0_linkedcolumn, lc, lc, VARTOOLS_TYPE_DOUBLE, &d2);
		}
	      }

	      phaselc(p->NJD[lc],p->t[lc],p->mag[lc],p->sig[lc],d1,i2,d2,c->Phase->phasevarname,c->Phase->phasevar,lc,c->Phase->startphase);
	      sortlcbytime(p->NJD[lc], p->t[lc], lc, p);

	    }
	}
      break;

    case CNUM_OUTPUTLCS:
      /* Write out the light curves in their present form */
      for(lc=0;lc<p->Nlcs;lc++)
	{
	  if(p->isifcommands) {
	    if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
	      SkipCommand(p, c, thisindex, lc, lc);
	      continue;
	    }
	  }
	  DoOutputLightCurve(p, c->Outputlcs, lc, lc);
	  /*
	  if(p->listflag)
	    {
	      i1 = 0;
	      i2 = 0;
	      while(p->lcnames[lc][i1] != '\0')
		{
		  if(p->lcnames[lc][i1] == '/')
		    i2 = i1 + 1;
		  i1++;
		}
	      if(!c->Outputlcs->useformat)
		sprintf(outname,"%s/%s",c->Outputlcs->outdir,&p->lcnames[lc][i2]);
	      else
		{
		  sprintf(outname,"%s/",c->Outputlcs->outdir);
		  i1=strlen(outname);
		  i3=0;
		  while(c->Outputlcs->format[i3] != '\0')
		    {
		      if(c->Outputlcs->format[i3] != '%')
			{
			  outname[i1] = c->Outputlcs->format[i3];
			  i1++;
			  outname[i1] = '\0';
			  i3++;
			}
		      else
			{
			  i3++;
			  if(c->Outputlcs->format[i3] == 's')
			    {
			      i3++;
			      sprintf(&outname[i1],"%s",&p->lcnames[lc][i2]);
			      i1 = strlen(outname);
			    }
			  else if(c->Outputlcs->format[i3] == 'd')
			    {
			      i3++;
			      sprintf(&outname[i1],"%d",lc+1);
			      i1 = strlen(outname);
			    }
			  else if(c->Outputlcs->format[i3] == '0')
			    {
			      i3++;
			      tmpstring[0] = '%';
			      tmpstring[1] = '0';
			      i4 = 2;
			      while(c->Outputlcs->format[i3] >= '1' && c->Outputlcs->format[i3] <= '9')
				{
				  tmpstring[i4] = c->Outputlcs->format[i3];
				  i4++;
				  i3++;
				}
			      if(c->Outputlcs->format[i3] != 'd')
				error(ERR_INVALIDOUTPUTFORMAT);
			      i3++;
			      tmpstring[i4] = 'd';
			      i4++;
			      tmpstring[i4] = '\0';
			      sprintf(&outname[i1],tmpstring,lc+1);
			      i1 = strlen(outname);
			    }
			  else if(c->Outputlcs->format[i3] == '%')
			    {
			      i3++;
			      outname[i1] = '%';
			      i1++;
			      outname[i1] = '\0';
			    }
			  else
			    error(ERR_INVALIDOUTPUTFORMAT);
			}
		    }
		}
	      writelightcurves(p, lc, lc, outname, c->Outputlcs->usecolumnformat, c->Outputlcs->Nvar, c->Outputlcs->variables, c->Outputlcs->printfformats);
	    }
	  else if(p->fileflag)
	    {
	      writelightcurves(p, lc, lc, c->Outputlcs->outdir, c->Outputlcs->usecolumnformat, c->Outputlcs->Nvar, c->Outputlcs->variables, c->Outputlcs->printfformats);
	      }*/
	}
      break;

#ifdef DYNAMICLIB
    case CNUM_USERCOMMAND:
      if(p->UserLib[c->UserCommand->libnum].RequireReadAll) {
	RunUserCommand_all_lcs(p,c);
      }
      else {
	for(lc=0;lc<p->Nlcs;lc++)
	  {
	    if(p->isifcommands) {
	      if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
		SkipCommand(p, c, thisindex, lc, lc);
		continue;
	      }
	    }
	    RunUserCommand(p,c,lc,lc);
	  }
      }
      break;
#endif

#ifdef _HAVE_PYTHON
    case CNUM_PYTHON:
      if(c->PythonCommand->RequireReadAll) {
	RunPythonCommand_all_lcs(p, c->PythonCommand);
	if(!c->PythonCommand->iscontinueprocess &&
	   c->PythonCommand->Nchildren == 0) {
	  StopRunningPythonCommand(p, 0, c->PythonCommand);
	}
	else if(c->PythonCommand->iscontinueprocess) {
	  if(((_PythonCommand **)((_PythonCommand *)c->PythonCommand->continueprocesscommandptr)->childcommandptrs)[(((_PythonCommand *)c->PythonCommand->continueprocesscommandptr)->Nchildren)-1] == c->PythonCommand) {
	    StopRunningPythonCommand(p, 0, c->PythonCommand);
	  }
	}
      } else {
	for(lc = 0; lc < p->Nlcs; lc++)
	  {
	    if(p->isifcommands) {
	      if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
		SkipCommand(p, c, thisindex, lc, lc);
		continue;
	      }
	    }
	    RunPythonCommand(p, lc, lc, 0, c->PythonCommand);
	  }
	if(!c->PythonCommand->iscontinueprocess &&
	   c->PythonCommand->Nchildren == 0) {
	  StopRunningPythonCommand(p, 0, c->PythonCommand);
	}
	else if(c->PythonCommand->iscontinueprocess) {
	  if(((_PythonCommand **)((_PythonCommand *)c->PythonCommand->continueprocesscommandptr)->childcommandptrs)[(((_PythonCommand *)c->PythonCommand->continueprocesscommandptr)->Nchildren)-1] == c->PythonCommand) {
	    StopRunningPythonCommand(p, 0, c->PythonCommand);
	  }
	}
      }
      break;
#endif

#ifdef _HAVE_R
    case CNUM_R:
      if(c->RCommand->RequireReadAll) {
	RunRCommand_all_lcs(p, c->RCommand);
	if(!c->RCommand->iscontinueprocess &&
	   c->RCommand->Nchildren == 0) {
	  StopRunningRCommand(p, 0, c->RCommand);
	}
	else if(c->RCommand->iscontinueprocess) {
	  if(((_RCommand **)((_RCommand *)c->RCommand->continueprocesscommandptr)->childcommandptrs)[(((_RCommand *)c->RCommand->continueprocesscommandptr)->Nchildren)-1] == c->RCommand) {
	    StopRunningRCommand(p, 0, c->RCommand);
	  }
	}
      } else {
	for(lc = 0; lc < p->Nlcs; lc++)
	  {
	    if(p->isifcommands) {
	      if(!TestIf(p->IfStack[lc], p, c, lc, lc) || p->skipfaillc[lc]) {
		SkipCommand(p, c, thisindex, lc, lc);
		continue;
	      }
	    }
	    RunRCommand(p, lc, lc, 0, c->RCommand);
	  }
	if(!c->RCommand->iscontinueprocess &&
	   c->RCommand->Nchildren == 0) {
	  StopRunningRCommand(p, 0, c->RCommand);
	}
	else if(c->RCommand->iscontinueprocess) {
	  if(((_RCommand **)((_RCommand *)c->RCommand->continueprocesscommandptr)->childcommandptrs)[(((_RCommand *)c->RCommand->continueprocesscommandptr)->Nchildren)-1] == c->RCommand) {
	    StopRunningRCommand(p, 0, c->RCommand);
	  }
	}
      }
      break;
#endif

    default:
      error(ERR_CODEERROR);
      break;
    }
}

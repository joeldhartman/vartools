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

#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

#define DEF_QUICKSELECT(type) \
type quickselect_##type (int k, int n, type *arr) \
{ \
  int i, ir, j,m, mid; \
  type a, temp; \
  k--; \
  m = 0; \
  ir = n-1; \
  for(;;) { \
    if(ir<=m+1) { \
      if(ir == m+1 && (arr)[ir] < (arr)[m]) {	\
	SWAP((arr)[m],(arr)[ir]);	  \
      } \
      return (arr)[k];				\
    } else { \
      mid=(m+ir) >> 1; \
      SWAP((arr)[mid],(arr)[m+1]);		\
      if((arr)[m] > (arr)[ir]) {		\
	SWAP((arr)[m],(arr)[ir]);		\
      } \
      if((arr)[m+1] > (arr)[ir]) {		\
	SWAP((arr)[m+1],(arr)[ir]);		\
      } \
      if((arr)[m] > (arr)[m+1]) {		\
	SWAP((arr)[m],(arr)[m+1]);		\
      } \
      i=m+1; \
      j=ir; \
      a=arr[m+1]; \
      for(;;) { \
	do i++; while((arr)[i] < a);		\
	do j--; while((arr)[j] > a);		\
	if(j < i) break; \
	SWAP((arr)[i],(arr)[j]);		\
      } \
      (arr)[m+1]=(arr)[j];			\
      (arr)[j]=a;				\
      if(j >= k) ir=j-1; \
      if(j <=k) m=i; \
    } \
  } \
}

DEF_QUICKSELECT(double)
DEF_QUICKSELECT(float)
DEF_QUICKSELECT(int)
DEF_QUICKSELECT(short)
DEF_QUICKSELECT(long)
DEF_QUICKSELECT(char)

double quickselect(int k, int n, double *arr)
/* Select the k smallest element from the vector 'data' of size n in
   order n time. data is re-arranged by this process. This algorithm
   is taken from Numerical Recipes in C, Press et al. 
   The name has been changed from select to quickselect, the datatypes
   and array indexing have also been changed.
*/
{
  return quickselect_double(k, n, arr);
}

/* The version of the Quickselect routine below includes a correction for weights following Rauh and Arce, 2010, in ICIP(2010) p 105-108 in this case select the element X for which the sum of all weights of elements less than X is equal to W0. */
#define DEF_QUICKSELECT_WEIGHT(type) \
  type quickselect_weight_##type (double W0, int n, type *arr, double *weight, double Wl, double Wr) \
{ \
  int i, ir, j,m, mid; \
  type a, temp; \
  double wSumLeft = 0.0; \
  double wLeft = Wl; \
  double wRight = Wr; \
  m = 0; \
  ir = n-1; \
  for(;;) { \
    if(ir <= m) { \
       return arr[m]; \
    } \
    if(ir == m+1) { \
      if(arr[m] > arr[ir]) { \
        SWAP(arr[m],arr[ir]); \
        SWAP(weight[m],weight[ir]); \
      } \
      if( wLeft + weight[m] >= W0 ) { \
        return arr[m]; \
      } else { \
        return arr[ir]; \
      } \
    } \
    mid=(m+ir) >> 1; \
    if(arr[mid] > arr[ir]) { \
      SWAP(arr[mid],arr[ir]); \
      SWAP(weight[mid],weight[ir]); \
    } \
    if(arr[m] > arr[ir]) { \
      SWAP(arr[m],arr[ir]); \
      SWAP(weight[m],weight[ir]); \
    } \
    if(arr[mid] > arr[m]) { \
      SWAP(arr[mid], arr[m]); \
      SWAP(weight[mid], weight[m]); \
    } \
    SWAP(arr[mid], arr[m+1]); \
    SWAP(weight[mid], weight[m+1]); \
    a = arr[m]; \
    i = m + 1; \
    j = ir; \
    for(;;) { \
      do { \
        wSumLeft += weight[i]; \
        i++; \
      } while(a > arr[i]); \
      do { \
	  j--; \
      } while(arr[j] > a); \
      if(j < i) break; \
      SWAP(arr[i],arr[j]); \
      SWAP(weight[i],weight[j]); \
    } \
    SWAP(arr[m],arr[j]); \
    SWAP(weight[m],weight[j]); \
    while(i > j+1) { \
      i--; \
      wSumLeft -= weight[i]; \
    } \
    if( (wLeft + wSumLeft) <= W0) { \
      if( (wLeft + wSumLeft + weight[j]) > W0) { \
        return arr[j]; \
      } \
      m = j+1; \
      wLeft += weight[j] + wSumLeft; \
    } else { \
      ir = j-1; \
      wRight = 2*W0 - wSumLeft - wLeft; \
    } \
    wSumLeft = 0.0; \
  } \
}

DEF_QUICKSELECT_WEIGHT(double)
DEF_QUICKSELECT_WEIGHT(float)
DEF_QUICKSELECT_WEIGHT(int)
DEF_QUICKSELECT_WEIGHT(short)
DEF_QUICKSELECT_WEIGHT(long)
DEF_QUICKSELECT_WEIGHT(char)

double quickselect_weight(double W0, int n, double *arr, double *weight, double wl, double wr)
/* Select the smallest element from the vector 'data' of size n for
   which the sum of weights of all elements smaller than the element,
   plus the element's own weight, is greater than W0. This is done in
   order n time. The elements of arr and weight are re-ordered in this
   process. wl is the sum of weights already discarded to the left, wr the
   sum already discarded to the right.
*/
{
  return quickselect_weight_double(W0, n, arr, weight, wl, wr);
}


#define DEF_GETWEIGHTEDMEAN(type) \
type getweightedmean_##type (int n, type *data, double *sig){ \
  int i; \
  double tmp; \
  double num = (double) n; \
  long double mean_val1 = 0, mean_val2 = 0; \
  for (i=0;i<n;i++){ \
    tmp = 1./sig[i]/sig[i]; \
    mean_val1 += (long double) data[i]*tmp; \
    mean_val2 += (long double) tmp; \
  } \
  return (type) (mean_val1/mean_val2); \
}

DEF_GETWEIGHTEDMEAN(double)
DEF_GETWEIGHTEDMEAN(float)
DEF_GETWEIGHTEDMEAN(int)
DEF_GETWEIGHTEDMEAN(long)
DEF_GETWEIGHTEDMEAN(short)    

double getweightedmean(int n, double *data, double *sig){
  return getweightedmean_double(n, data, sig);
}

#define DEF_GETMEAN(type) \
type getmean_##type (int n, type *data){ \
  int i; \
  double num = (double) n; \
  long double mean_val = 0; \
  for (i=0;i<n;i++){ \
    mean_val += (double) (data[i]/ num); \
  } \
  return (type) mean_val; \
}

DEF_GETMEAN(double)
DEF_GETMEAN(float)
DEF_GETMEAN(int)
DEF_GETMEAN(long)
DEF_GETMEAN(short)    

double getmean(int n, double *data){
  return getmean_double(n, data);
}

#define DEF_MEDIAN(type) \
type median_##type (int n, type *data) \
{ \
  type temp; \
  type *data2; \
  if(n == 1) \
    return(data[0]); \
  if((data2 = (type *) malloc(n * sizeof(type))) == NULL) \
    error(ERR_MEMALLOC); \
  memcpy(data2,data,n*sizeof(type)); \
  if(n % 2 == 0) { \
    temp = quickselect_##type ((n/2),n,data2); \
    temp += quickselect_##type ((n/2)+1,n,data2); \
    temp = temp/2.0; \
  } \
  else { \
    temp = quickselect_##type ((n/2)+1,n,data2); \
  } \
  free(data2); \
  return(temp); \
}

DEF_MEDIAN(double)
DEF_MEDIAN(float)
DEF_MEDIAN(int)
DEF_MEDIAN(long)
DEF_MEDIAN(short)    

double median(int n, double *data)
{
  return median_double(n, data);
}

#define DEF_MEDIAN_NOCOPY(type) \
type median_nocopy_##type (int n, type *data) \
{ \
  type temp; \
  if(n == 1) \
    return(data[0]); \
  if(n % 2 == 0) { \
    temp = quickselect_##type ((n/2),n,data); \
    temp += quickselect_##type ((n/2)+1,n,data); \
    temp = temp/2.0; \
  } \
  else { \
    temp = quickselect_##type ((n/2)+1,n,data); \
  } \
  return(temp); \
}
DEF_MEDIAN_NOCOPY(double)
DEF_MEDIAN_NOCOPY(float)
DEF_MEDIAN_NOCOPY(int)
DEF_MEDIAN_NOCOPY(long)
DEF_MEDIAN_NOCOPY(short)    


double median_nocopy(int n, double *data)
{
  return median_nocopy_double(n, data);
}

#define DEF_MEDIAN_WEIGHT(type) \
type median_weight_##type (int n, type *data, double *err)	\
{ \
  type temp; \
  type *data2; \
  double *weight; \
  double wsum, W0; \
  int i; \
  if(n == 1) \
    return(data[0]); \
  if((data2 = (type *) malloc(n * sizeof(type))) == NULL || \
     (weight = (double *) malloc(n * sizeof(double))) == NULL)	\
    error(ERR_MEMALLOC); \
  memcpy(data2,data,n*sizeof(type)); \
  wsum = 0.0; \
  for(i=0;i<n;i++){ \
     if(err[i] > 0) {weight[i]=1./err[i]/err[i];} \
     else weight[i] = 0.0; \
     wsum += weight[i]; \
  } \
  if(wsum == 0.0){ \
    temp = median_nocopy_##type (n,data2);	\
  } else { \
     W0 = wsum*0.5; \
     temp = quickselect_weight_##type (W0,n,data2,weight,0.0,0.0);	\
  } \
  free(data2); \
  free(weight); \
  return(temp); \
}

DEF_MEDIAN_WEIGHT(double)
DEF_MEDIAN_WEIGHT(float)
DEF_MEDIAN_WEIGHT(int)
DEF_MEDIAN_WEIGHT(long)
DEF_MEDIAN_WEIGHT(short)    

double median_weight(int n, double *data, double *err)
{
  return median_weight_double(n, data, err);
}

#define DEF_MEDIAN_WEIGHT_NOCOPY(type) \
type median_weight_nocopy_##type (int n, type *data, double *err)	\
{ \
  type temp; \
  double *weight; \
  double wsum, W0; \
  int i; \
  if(n == 1) \
    return(data[0]); \
  if((weight = (double *) malloc(n * sizeof(double))) == NULL)	\
    error(ERR_MEMALLOC); \
  wsum = 0.0; \
  for(i=0;i<n;i++){ \
     if(err[i] > 0) {weight[i]=1./err[i]/err[i];} \
     else weight[i] = 0.0; \
     wsum += weight[i]; \
  } \
  if(wsum == 0.0){ \
    temp = median_nocopy_##type (n,data);	\
  } else { \
     W0 = wsum*0.5; \
     temp = quickselect_weight_##type (W0,n,data,weight,0.0,0.0);	\
  } \
  free(weight); \
  return(temp); \
}

DEF_MEDIAN_WEIGHT_NOCOPY(double)
DEF_MEDIAN_WEIGHT_NOCOPY(float)
DEF_MEDIAN_WEIGHT_NOCOPY(int)
DEF_MEDIAN_WEIGHT_NOCOPY(long)
DEF_MEDIAN_WEIGHT_NOCOPY(short)    

double median_weight_nocopy(int n, double *data, double *err)
{
  return median_weight_nocopy_double(n, data, err);
}

#define DEF_MEDMEDDEV(type) \
type medmeddev_##type (int n, type *data) \
{ \
  type medval1; \
  int i; \
  type *data2; \
  type temp; \
  if(n == 1) \
    return (type) 0.0;					  \
  if((data2 = (type *) malloc(n * sizeof(type))) == NULL) \
    error(ERR_MEMALLOC); \
  memcpy(data2,data,n*sizeof(type)); \
  medval1 = median_nocopy_##type (n,data2); \
  for(i=0;i<n;i++) { \
    data2[i] = data2[i] - medval1; \
    if(data2[i] < 0) \
      data2[i] = -data2[i]; \
  } \
  temp = median_nocopy_##type (n,data2); \
  free(data2); \
  return temp; \
}

DEF_MEDMEDDEV(double)
DEF_MEDMEDDEV(float)
DEF_MEDMEDDEV(int)
DEF_MEDMEDDEV(long)
DEF_MEDMEDDEV(short)    

double medmeddev(int n, double *data)
{
  return medmeddev_double(n, data);
}

#define DEF_MAD(type) \
type MAD_##type (int n, type *data) \
/* Compute the MAD statistic = 1.483*(median absolute deviation from \
   the median). For a gaussian stddev = 1.483*medmeddev, hence the factor \
   of 1.483 here. \
*/ \
 \
{ \
  return (type) (1.483*((double) medmeddev_##type (n, data)));	\
}

DEF_MAD(double)
DEF_MAD(float)
DEF_MAD(int)
DEF_MAD(long)
DEF_MAD(short)

double MAD (int n, double *data)
/* Compute the MAD statistic = 1.483*(median absolute deviation from
   the median). For a gaussian stddev = 1.483*medmeddev, hence the factor
   of 1.483 here.
*/
{
 return 1.483*medmeddev_double (n, data);
}


#define DEF_MEDDEV(type) \
type meddev_##type(int n, type *data) \
/* Compute the deviation about the median, this is basically the standard \
   deviation, but calculated with respect to the median rather than the \
   mean */ \
{ \
  type medval1, var2, var3, dn; \
  int i; \
  if(n <= 1) \
    return 0.0; \
  medval1 = median_##type(n, data); \
  var2 = 0.; \
  dn = (type) n - 1; \
  for(i=0; i < n; i++) { \
    var3 = data[i] - medval1; \
    var2 += var3*var3 / (type) dn; \
  } \
  return((type) (sqrt((double)var2))); \
}

DEF_MEDDEV(double)
DEF_MEDDEV(float)
DEF_MEDDEV(int)
DEF_MEDDEV(long)
DEF_MEDDEV(short)

double meddev(int n, double *data)
{
  return meddev_double(n, data);
}

#define DEF_STDDEV(type) \
type stddev_##type (int n, type *data) \
{ \
  int i; \
  type var1, var2, var3, dn; \
  var1 = getmean_##type (n, data); \
  var2 = 0.; \
  dn = (type) n - 1; \
  if(n > 1) \
    { \
      for(i=0;i<n;i++){ \
	var3 = data[i] - var1; \
	var2 += var3*var3 / (type) dn; \
      } \
      return((type) (sqrt((double) var2))); \
    } \
  else \
    return(0.); \
}

DEF_STDDEV(double)
DEF_STDDEV(float)
DEF_STDDEV(int)
DEF_STDDEV(long)
DEF_STDDEV(short)

double stddev(int n, double *data)
{
  return stddev_double(n, data);
}

#define DEF_KURTOSIS(type) \
type kurtosis_##type (int n, type *data) \
{ \
  int i; \
  long double var1 = 0, var2 = 0, var3 = 0, var4 = 0, var5; \
  double dat;						    \
  for(i=0;i<n;i++){ \
    dat = (double) data[i]; \
    var1 += (long double) (dat/ (double) n); \
    var2 += (long double) (dat*dat) / (double) n; \
    var3 += (long double) (dat*dat*dat) / (double) n; \
    var4 += (long double) (dat*dat*dat*dat) / (double) n; \
  } \
  /* var5 = stddev_##type(n,data); */					\
  return((type) ((-3*(var1*var1*var1*var1)+6*(var1*var1*var2)-4*(var1*var3)+var4)/(var2*var2 - 2*var1*var1*var2 + var1*var1*var1*var1))); \
}

DEF_KURTOSIS(double)
DEF_KURTOSIS(float)
DEF_KURTOSIS(int)
DEF_KURTOSIS(long)
DEF_KURTOSIS(short)

double kurtosis(int n, double *data) {
  return kurtosis_double(n, data);
}

#define DEF_SKEWNESS(type) \
type skewness_##type (int n, type *data) \
{ \
  int i; \
  long double var1 = 0, var2 = 0, var3 = 0, var5; \
  double dat; \
  for(i=0;i<n;i++){ \
    dat = (double) data[i]; \
    var1 += (long double) (dat/(double) n); \
    var2 += (long double) (dat*dat)/(double)n; \
    var3 += (long double) (dat*dat*dat)/(double)n; \
  } \
  var5 = pow((double)(var2 - var1*var1),1.5); \
  return (type)((var3 - 3*var1*var2 + 2*var1*var1*var1)/var5); \
}

DEF_SKEWNESS(double)
DEF_SKEWNESS(float)
DEF_SKEWNESS(int)
DEF_SKEWNESS(long)
DEF_SKEWNESS(short)

double skewness(int n, double *data) {
  return skewness_double(n, data);
}


#define DEF_PERCENTILE(type) \
type percentile_##type (int n, type *data, double pct) \
{ \
  int i, n1, n2; \
  double N_; \
  double dN; \
  type temp; \
  type *data2; \
  \
  if(n == 1) \
    return(data[0]);  \
  \
  if((data2 = (type *) malloc(n * sizeof(type))) == NULL) \
    error(ERR_MEMALLOC); \
  memcpy(data2,data,n*sizeof(type)); \
  \
  \
  dN = 1.0 / ((double) n); \
  \
  if(pct <= (dN/2.0)) { \
    temp = data2[0]; \
    for(i = 1; i < n; i++) { \
      if(data[i] < temp) \
	temp = data2[i]; \
    } \
    free(data2); \
    return temp; \
  } \
  if(pct >= 100. - (dN/2.0)) { \
    temp = data2[0]; \
    for(i = 1; i < n; i++) { \
      if(data2[i] > temp) \
	temp = data2[i]; \
    } \
    free(data2); \
    return temp; \
  } \
  \
  N_ = pct*((double) (n-1))/100. - dN/2.0; \
  \
  n1 = floor((double) N_); \
  n2 = ceil((double) N_); \
  if(n1 > n-1) n1 = n-1; \
  if(n2 > n-1) n2 = n-1; \
  if(n1 < 0) n1 = 0; \
  if(n2 < 0) n2 = 0; \
  \
  if(n1 != n2) { \
    temp = quickselect_##type (n1+1,n,data2); \
    temp += (type) ((N_ - (double) n1)*((double) (quickselect_##type (n2+1,n,data2) - temp))); \
  } else { \
    temp = quickselect_##type (n1+1,n,data2); \
  } \
  free(data2); \
  return temp; \
}

DEF_PERCENTILE(double)
DEF_PERCENTILE(float)
DEF_PERCENTILE(int)
DEF_PERCENTILE(long)
DEF_PERCENTILE(short)

double percentile (int n, double *data, double pct) 
{
  return percentile_double(n, data, pct);
}

#define DEF_PERCENTILE_NOCOPY(type) \
type percentile_nocopy_##type (int n, type *data, double pct) \
{ \
  int i, n1, n2; \
  double N_; \
  double dN; \
  type temp; \
  \
  if(n == 1) \
    return(data[0]); \
  \
  dN = 1.0 / ((double) n); \
  \
  if(pct <= (dN/2.0)) { \
    temp = data[0]; \
    for(i = 1; i < n; i++) { \
      if(data[i] < temp) \
	temp = data[i]; \
    } \
    return temp; \
  } \
  if(pct >= 100. - (dN/2.0)) { \
    temp = data[0]; \
    for(i = 1; i < n; i++) { \
      if(data[i] > temp) \
	temp = data[i]; \
    } \
    return temp; \
  } \
  \
  N_ = pct*((double) (n-1))/100. - dN/2.0; \
  \
  n1 = floor((double) N_); \
  n2 = ceil((double) N_); \
  if(n1 > n-1) n1 = n-1; \
  if(n2 > n-1) n2 = n-1; \
  if(n1 < 0) n1 = 0; \
  if(n2 < 0) n2 = 0; \
  \
  if(n1 != n2) { \
    temp = quickselect_##type (n1+1,n,data); \
    temp += (type) ((N_ - (double) n1)*(double) (quickselect_##type(n2+1,n,data) - temp)); \
  } else { \
    temp = quickselect_##type (n1+1,n,data); \
  } \
  return temp; \
}

DEF_PERCENTILE_NOCOPY(double)
DEF_PERCENTILE_NOCOPY(float)
DEF_PERCENTILE_NOCOPY(int)
DEF_PERCENTILE_NOCOPY(long)
DEF_PERCENTILE_NOCOPY(short)

double percentile_nocopy (int n, double *data, double pct) 
{
  return percentile_nocopy_double(n, data, pct);
}

#define DEF_PERCENTILE_WEIGHT(type) \
type percentile_weight_##type (int n, type *data, double *err, double pct) \
{ \
  int i, n1, n2; \
  double N_; \
  double dN; \
  type temp; \
  type *data2; \
  double *weight; \
  double wsum, wmin, wmax, W0; \
  \
  if(n == 1) \
    return(data[0]); \
  \
  if((data2 = (type *) malloc(n * sizeof(type))) == NULL || \
     (weight = (double *) malloc(n * sizeof(double))) == NULL) \
    error(ERR_MEMALLOC); \
  memcpy(data2,data,n*sizeof(type)); \
  \
  wsum = 0.0; \
  for(i=0; i < n; i++) { \
    if(err[i] > 0) { \
      weight[i] = 1./err[i]/err[i]; \
    } \
    else weight[i] = 0.0; \
    wsum += weight[i]; \
  } \
  \
  if(wsum == 0) { \
    temp = percentile_nocopy_##type (n, data2, pct); \
    free(weight); \
    return temp; \
  } \
  \
  W0 = pct*wsum/100.0; \
  \
  temp = quickselect_weight_##type (W0, n, data2, weight, 0.0, 0.0); \
  \
  free(data2); \
  free(weight); \
  return temp; \
  \
}

DEF_PERCENTILE_WEIGHT(double)
DEF_PERCENTILE_WEIGHT(float)
DEF_PERCENTILE_WEIGHT(int)
DEF_PERCENTILE_WEIGHT(long)
DEF_PERCENTILE_WEIGHT(short)

double percentile_weight (int n, double *data, double *err, double pct) 
{
  return percentile_weight_double(n, data, err, pct);
}

#define DEF_PERCENTILE_WEIGHT_NOCOPY(type) \
type percentile_weight_nocopy_##type (int n, type *data, double *err, double pct) \
{ \
  int i, n1, n2; \
  double N_; \
  double dN; \
  type temp; \
  double *weight; \
  double wsum, wmin, wmax, W0; \
  \
  if(n == 1) \
    return(data[0]);  \
  \
  if((weight = (double *) malloc(n * sizeof(double))) == NULL) \
    error(ERR_MEMALLOC); \
  \
  wsum = 0.0; \
  for(i=0; i < n; i++) { \
    if(err[i] > 0) { \
      weight[i] = 1./err[i]/err[i]; \
    } \
    else weight[i] = 0.0; \
    wsum += weight[i]; \
  } \
  \
  if(wsum == 0) { \
    temp = percentile_nocopy_##type (n, data, pct); \
    free(weight); \
    return temp; \
  } \
  \
  W0 = pct*wsum/100.0; \
  \
  temp = quickselect_weight_##type (W0, n, data, weight, 0.0, 0.0); \
  \
  free(weight); \
  return temp; \
  \
}

DEF_PERCENTILE_WEIGHT_NOCOPY(double)
DEF_PERCENTILE_WEIGHT_NOCOPY(float)
DEF_PERCENTILE_WEIGHT_NOCOPY(int)
DEF_PERCENTILE_WEIGHT_NOCOPY(long)
DEF_PERCENTILE_WEIGHT_NOCOPY(short)

double percentile_weight_nocopy (int n, double *data, double *err, double pct) 
{
  return percentile_weight_nocopy_double(n, data, err, pct);
}

#define DEF_GETMAXIMUM(type) \
type getmaximum_##type (int n, type *data) \
{ \
  int j; \
  type temp; \
  temp = data[0]; \
  for(j=1; j < n; j++) { \
    if(data[j] > temp) \
      temp = data[j]; \
  } \
  return temp; \
}

DEF_GETMAXIMUM(double)
DEF_GETMAXIMUM(float)
DEF_GETMAXIMUM(int)
DEF_GETMAXIMUM(long)
DEF_GETMAXIMUM(short)

double getmaximum(int n, double *data) {
  return getmaximum_double(n, data);
}

#define DEF_GETMINIMUM(type) \
type getminimum_##type (int n, type *data) \
{ \
  int j; \
  type temp; \
  temp = data[0]; \
  for(j=1; j < n; j++) { \
    if(data[j] < temp) \
      temp = data[j]; \
  } \
  return temp; \
}

DEF_GETMINIMUM(double)
DEF_GETMINIMUM(float)
DEF_GETMINIMUM(int)
DEF_GETMINIMUM(long)
DEF_GETMINIMUM(short)

double getminimum(int n, double *data) {
  return getminimum_double(n, data);
}

#define DEF_GETSUM(type) \
type getsum_##type (int n, type *data) \
{ \
  int j; \
  type temp = 0; \
  for(j=0; j < n; j++) { \
    temp += data[j]; \
  } \
  return temp; \
}

DEF_GETSUM(double)
DEF_GETSUM(float)
DEF_GETSUM(int)
DEF_GETSUM(long)
DEF_GETSUM(short)

double getsum(int n, double *data) {
  return getsum_double(n, data);
}
   

void RunStatsCommand(ProgramData *p, int lcindex, int threadindex, _Stats *s)
{
  int i, j, k, Npct;
  double *tmpdata = NULL, *tmpweight = NULL;
  if(p->NJD[threadindex] <= 0) {
    for(i=0, k=0; i < s->Nvar; i++) {
      for(j=0; j < s->Nstats; j++, k++) {
	s->statsout[threadindex][k] = 0.0;
      }
    }
    return;
  }
  if((tmpdata = (double *) malloc(p->NJD[threadindex]*sizeof(double))) == NULL) {
    error(ERR_MEMALLOC);
  }
  for(i = 0, k=0; i < s->Nvar; i++) {
    if(s->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC) {
      error(ERR_BADVARIABLETYPE_STATSCOMMAND);
    }
    for(j=0; j < p->NJD[threadindex]; j++) {
      tmpdata[j] = EvaluateVariable_Double(lcindex, threadindex, j, s->vars[i]);
    }
    Npct = 0;
    for(j = 0; j < s->Nstats; j++, k++) {
      switch(s->statstocalc[j]) {
      case VARTOOLS_STATSTYPE_MEAN:
	s->statsout[threadindex][k] = getmean(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_WEIGHTEDMEAN:
	s->statsout[threadindex][k] = getweightedmean(p->NJD[threadindex], tmpdata, p->sig[threadindex]);
	break;
      case VARTOOLS_STATSTYPE_MEDIAN:
	s->statsout[threadindex][k] = median(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT:
	s->statsout[threadindex][k] = median_weight(p->NJD[threadindex], tmpdata, p->sig[threadindex]);
	break;
      case VARTOOLS_STATSTYPE_STDDEV:
	s->statsout[threadindex][k] = stddev(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_MEDDEV:
	s->statsout[threadindex][k] = meddev(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_MEDMEDDEV:
	s->statsout[threadindex][k] = medmeddev(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_MAD:
	s->statsout[threadindex][k] = MAD(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_KURTOSIS:
	s->statsout[threadindex][k] = kurtosis(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_SKEWNESS:
	s->statsout[threadindex][k] = skewness(p->NJD[threadindex], tmpdata);
	break;
      case VARTOOLS_STATSTYPE_PERCENTILE:
	s->statsout[threadindex][k] = percentile(p->NJD[threadindex], 
							tmpdata,
							s->pctval[Npct]);
	Npct++;
	break;
      case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT:
	s->statsout[threadindex][k] = percentile_weight(p->NJD[threadindex], 
							tmpdata,
							       p->sig[threadindex],
							s->pctval[Npct]);
	Npct++;
	break;
      case VARTOOLS_STATSTYPE_MAXIMUM:
	s->statsout[threadindex][k] = getmaximum(p->NJD[threadindex],tmpdata);
	break;
      case VARTOOLS_STATSTYPE_MINIMUM:
	s->statsout[threadindex][k] = getminimum(p->NJD[threadindex],tmpdata);
	break;
      case VARTOOLS_STATSTYPE_SUM:
	s->statsout[threadindex][k] = getsum(p->NJD[threadindex],tmpdata);
	break;
      default:
	error(ERR_CODEERROR);
      }
    }
  }
  if(tmpdata != NULL)
    free(tmpdata);
}

int ParseStatsCommand(int *iret, int argc, char **argv, ProgramData *p, _Stats *s)
{
  int i, j, k, i1, i2, Nvar, Nstat, lentmp, stattype, Npct = 0;
  double pctval;
  char *tmpstring;
  i = *iret;
  if(i >= argc)
    return 1;
  Nvar = 1;
  j = 0;
  while(argv[i][j] != '\0') {
    if(argv[i][j] == ',')
      Nvar++;
    j++;
  }
  if((s->varnames = (char **) malloc(Nvar * sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);
  i1 = 0;
  i2 = 0;
  for(k = 0; k < Nvar; k++) {
    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
      i2++;
    }
    if((s->varnames[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    for(j=i1; j < i2; j++) {
      s->varnames[k][j-i1] = argv[i][j];
    }
    s->varnames[k][j-i1] = '\0';
    i1 = i2+1;
    i2 = i2+1;
  }

  s->Nvar = Nvar;
  if((s->vars = (_Variable **) malloc(Nvar * sizeof(_Variable *))) == NULL)
    error(ERR_MEMALLOC);

  i++;
  if(i >= argc) {
    *iret = i; return 1;
  }

  if((tmpstring = (char *) malloc(256 * sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  lentmp = 256;

  Nstat = 1;
  j = 0;
  while(argv[i][j] != '\0') {
    if(argv[i][j] == ',')
      Nstat++;
    j++;
  }
  if((s->statstocalc = (int *) malloc(Nstat * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  i1 = 0;
  i2 = 0;
  for(k = 0; k < Nstat; k++) {
    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
      i2++;
    }

    if((i2 - i1 + 1) > lentmp) {
      lentmp *= 2;
      if((tmpstring = (char *) realloc(tmpstring, lentmp*sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
    }
    for(j=i1; j < i2; j++) {
      tmpstring[j-i1] = argv[i][j];
    }
    tmpstring[j-i1] = '\0';

    if(!strcmp(tmpstring,"mean")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MEAN;
    }
    else if(!strcmp(tmpstring,"weightedmean")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_WEIGHTEDMEAN;
    }
    else if(!strcmp(tmpstring,"median")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MEDIAN;
    }
    else if(!strcmp(tmpstring,"wmedian")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MEDIAN_WEIGHT;
    }
    else if(!strcmp(tmpstring,"stddev")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_STDDEV;
    }
    else if(!strcmp(tmpstring,"meddev")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MEDDEV;
    }
    else if(!strcmp(tmpstring,"medmeddev")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MEDMEDDEV;
    }
    else if(!strcmp(tmpstring,"MAD")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MAD;
    }
    else if(!strcmp(tmpstring,"kurtosis")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_KURTOSIS;
    }
    else if(!strcmp(tmpstring,"skewness")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_SKEWNESS;
    }
    else if(sscanf(tmpstring,"wpct%lf",&pctval) == 1) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT;
      if(!Npct) {
	s->pctval = (double *) malloc(sizeof(double));
      } else {
	s->pctval = (double *) realloc(s->pctval, (Npct + 1)*sizeof(double));
      }
      s->pctval[Npct] = pctval;
      Npct++;
    }
    else if(sscanf(tmpstring,"pct%lf",&pctval) == 1) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_PERCENTILE;
      if(!Npct) {
	s->pctval = (double *) malloc(sizeof(double));
      } else {
	s->pctval = (double *) realloc(s->pctval, (Npct + 1)*sizeof(double));
      }
      s->pctval[Npct] = pctval;
      Npct++;
    }
    else if(!strcmp(tmpstring,"max")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MAXIMUM;
    }
    else if(!strcmp(tmpstring,"min")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_MINIMUM;
    }
    else if(!strcmp(tmpstring,"sum")) {
      s->statstocalc[k] = VARTOOLS_STATSTYPE_SUM;
    }
    else {
      error2(ERR_INVALIDSTATISTIC,tmpstring);
    }

    i1 = i2+1;
    i2 = i2+1;
  }

  s->Nstats = Nstat;

  s->Nstatstot = s->Nvar * s->Nstats;

  free(tmpstring);	

  *iret = i;
  return 0;
}

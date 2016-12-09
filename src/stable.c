/*
 *  stable : A Library of Functions for Stable Distributions
 *  Copyright (C) 1998, 1999, 2000, 2001 P. Lambert and J.K. Lindsey
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 *  void stable(int *n, double *y, double *beta, double *alpha, int *npt,
 *	    double *up, double *eps, int *type, int *err, double *ffy)
 *  void pstable(int *n, double *y, double *beta, double *alpha,
 *	    double *eps, int *err, double *ffy)
 *
 *  DESCRIPTION
 *
 *    Functions to compute the density and cdf of a stable distribution
 *
 */

#include <math.h>
#include <stddef.h>
#include "R.h"

static int nn;
static double alphai, etai, setai, cetai, yi, yyi;

static double fcn1(double s){
  double sa;
  sa=pow(s,alphai);
  return(cos(-yi*s+sa*setai)*exp(-sa*cetai));}

static double fcn2(double s){
  double sa;
  sa=pow(s,-alphai);
  return(cos(-yi/s+sa*setai)*exp(-sa*cetai)/(s*s));}

static double fcn3(double s){
  double sa;
  sa=pow(s,alphai);
  return((sin(yyi*s-sa*setai)/s)*exp(-sa*cetai));}

static double fcn4(double s){
  double sa;
  sa=pow(s,-alphai);
  return((sin(yyi/s-sa*setai)*s)*exp(-sa*cetai)/(s*s));}

/* integration routines: stripped down version of Romberg integration
   from rmutil library */

static void interp(double x[], double fx[], double *f, double *df)
{
  int i, j, ni=0;
  double diff1, diff2, tmp1, tmp2, lim1, lim2, tab1[5], tab2[5];
 
  tmp1=fabs(x[0]);
  for(i=0;i<5;i++){
    tmp2=fabs(x[i]);
    if(tmp2<tmp1){
      ni=i;
      tmp1=tmp2;}
    tab1[i]=fx[i];
    tab2[i]=fx[i];}
  *f=fx[ni--];
  for(j=0;j<4;j++){
    for(i=0;i<=4-j;i++){
      lim1=x[i];
      lim2=x[i+j+1];
      diff1=tab1[i+1]-tab2[i];
      diff2=lim1-lim2;
      if(diff2==0.0)return;
      diff2=diff1/diff2;
      tab2[i]=lim2*diff2;
      tab1[i]=lim1*diff2;}
    *df=2*ni<(2-j)?tab1[ni+1]:tab2[ni--];
    *f+=*df;}}

#define FCN(x) ((*fcn)(x))

static double evalfn(double (*fcn)(double), int n)
{
  double x, nn, tmpsum, pnt1, pnt2;
  static double sum;
  int i,j;

  if (n==1){
    sum=FCN(0.5);
    return(sum);}
  else {
    for(i=1,j=1;j<n-1;j++) i*=3;
    nn=i;
    pnt1=1/(3.0*nn);
    pnt2=2.0*pnt1;
    x=0.5*pnt1;
    tmpsum=0.0;
    for(j=1;j<=i;j++){
      tmpsum+=FCN(x);
      x+=pnt2;
      tmpsum+=FCN(x);
      x+=pnt1;}
    sum=(sum+tmpsum/nn)/3.0;
    return(sum);}}

static double romberg(double (*fcn)(double), double eps)
{
  int j,j1;
  double sum,errsum=0,x[16],fx[16];

  x[0]=1.0;
  for(j=0;j<16;j++){
    j1=j+1;
    fx[j]=evalfn(fcn,j1);
    if(j1>=5){
      interp(&x[j1-5],&fx[j1-5],&sum,&errsum);
      if(fabs(errsum)<eps*fabs(sum))return(sum);}
    x[j1]=x[j]/9.0;
    fx[j1]=fx[j];}
  sum=R_NaN;
  return(sum);}

/* density of a stable distribution */
void stable(int *n, double *y, double *beta, double *alpha, int *npt,
	    double *up, double *eps, int *type, int *err, double *ffy)
{
  int i, j;
  double h, s, *eta, *seta, *ceta, *sa;
  *err=0;
  eta=(double*)R_alloc((size_t)(*n),sizeof(double));
  seta=(double*)R_alloc((size_t)(*n),sizeof(double));
  ceta=(double*)R_alloc((size_t)(*n),sizeof(double));
  sa=(double*)R_alloc((size_t)(*n),sizeof(double));
  nn=*n;
  if(!eta||!seta||!ceta||!sa){
    *err=1;
    return;}
  for(i=0;i<*n;i++){
    ffy[i]=0.0;
    eta[i]=beta[i]*(1.0-fabs(1.0-alpha[i]))*M_PI/2.0;
    seta[i]=sin(eta[i]);
    ceta[i]=cos(eta[i]);}
  if(*type==1){
    *npt=*npt-*npt%2;
    h=*up/ *npt;
	   for(j=0;j<=*npt;j++){
	     s=(*npt-j)*h;
	     for(i=0;i<*n;i++){
	       sa[i]=pow(s,alpha[i]);
	       ffy[i]=ffy[i]+(4-2*(j%2==0)-(j==1||j==*npt))*cos(-y[i]*s+sa[i]*seta[i])*exp(-sa[i]*ceta[i]);}}
	   for(i=0;i<*n;i++)ffy[i]=ffy[i]*h/3.0/M_PI;}
  else {
    for(i=0;i<*n;i++){
      alphai=alpha[i];
      yi=y[i];
      setai=seta[i];
      cetai=ceta[i];
      ffy[i]=(romberg(fcn1, *eps)+romberg(fcn2, *eps))/M_PI;}}}

/* cdf of a stable distribution */
void pstable(int *n, double *y, double *beta, double *alpha,
	     double *eps, int *err, double *ffy)
{
  int i;
  *err=0;
  nn=*n;
  for(i=0;i<*n;i++){
    ffy[i]=0.0;
    etai=beta[i]*(1.0-fabs(1.0-alpha[i]))*M_PI/2.0;
    setai=sin(etai);
    cetai=cos(etai);
    alphai=alpha[i];
    yyi=y[i];
    if(etai==0.&&yyi==0)
      ffy[i]=0.5;
    else ffy[i]=0.5+(romberg(fcn3, *eps)+romberg(fcn4, *eps))/M_PI;}}


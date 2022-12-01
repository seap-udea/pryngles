#include "cpixx_dev.c"
//#include "cpixx.c"

#define VERBOSE 0
//////////////////////////////////////////////////////////////
// MAIN PROGRAM
//////////////////////////////////////////////////////////////
int main(void)
{
  int i,j,k,m;

  int shape[3];
  double xmu[MAX_MUS];
  double rfou[MAX_MAT*MAX_MUS][MAX_MUS][MAX_FOU];
  double rtra[MAX_MAT*MAX_MUS][MAX_MUS][MAX_FOU];
  char filename[MAX_STRING];
  int nmat,nmugs,nfou;
  int npix;
  double phi[MAX_PIX],beta[MAX_PIX],theta0[MAX_PIX],theta[MAX_PIX],apix[MAX_PIX],Sarr[MAX_PIX][MAX_MAT+1];
  char tfile[MAX_STRING];
  FILE *g,*f;
  double tmp;
  double Spixx[4];
  double dif,difmax=-1e100;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST READ FOURIER
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Planet
  sprintf(filename,"fou_gasplanet.dat");read_fourier(filename,shape,xmu,rfou,rtra);
  //Ring
  //sprintf(filename,"fou_ring_0_4_0_8.dat");read_fourier(filename,shape,xmu,rfou,rtra);
  //*/  

  //Get size of information
  nmat=shape[0];
  nmugs=shape[1];
  nfou=shape[2];
  printf("nmat = %d, nmugs = %d, nfou = %d\n",shape[0],shape[1],shape[2]);

  //*
  //Print rfou matrix
  sprintf(tfile,"test-fou-%s",filename);
  printf("Saving %s\n",tfile);
  g=fopen(tfile,"w");
  for(i=0;i<nmugs*nmat;i++){
    for(j=0;j<nmugs;j++){
      for(m=0;m<nfou;m++){
	fprintf(g,"%d %d %d %.16e\n",i,j,m,rfou[i][j][m]);
      }
    }
  }
  fclose(g);
  sprintf(tfile,"test-tra-%s",filename);
  printf("Saving %s\n",tfile);
  g=fopen(tfile,"w");
  for(i=0;i<nmugs*nmat;i++){
    for(j=0;j<nmugs;j++){
      for(m=0;m<nfou;m++){
	fprintf(g,"%d %d %d %.16e\n",i,j,m,rtra[i][j][m]);
      }
    }
  }
  fclose(g);
  //*/
  exit(0);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST REFLECTION SINGLE VALUE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*
  //Test values: planet
  npix=1;
  phi[0]=1.0448451569439827;
  beta[0]=3.069394277348945;
  theta0[0]=0.04990329026929557;
  theta[0]=0.02509670973070432;
  apix[0]=9.432328787795567e-05;
  reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);

  for(i=0;i<npix;i++){
    for(j=0;j<nmat+1;j++){
      printf("%.16e ",Sarr[i][j]);
    }
    printf("\n");
  }
  //*/

  npix=1;

  //Test values: forward
  /*
  phi[0]=1.490116119384765625e-08;
  beta[0]=0.000000000000000000e+00;
  theta0[0]=4.999999999999998890e-01;
  theta[0]=5.000000000000000000e-01;
  apix[0]=1.163314390931110409e-04;
  reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);
  printf("Expected:\n6.193646058775441171e-06 -3.046132218287162406e-07 -9.759550180589223642e-15 4.918156751904256135e-02\n");
  //*/
  
  //Test values: backward
  //*
  phi[0]=1.601029385538801364e+00;
  beta[0]=1.601029385538801364e+00;
  theta0[0]=1.744974835125044643e-02;
  theta[0]=5.000000000000000000e-01;
  apix[0]=1.163314390931110409e-04;
  //reflection(npix,phi,beta,theta0,theta,shape,xmu,rtra,apix,Sarr);
  reflection(npix,phi,beta,theta0,theta,shape,xmu,rtra,apix,Sarr);
  printf("Expected:\n1.688771016436214060e-07 -1.485273114913191718e-08 2.674513729889046925e-09 8.936444506313186154e-02\n");
  //*/

  for(i=0;i<npix;i++){
    for(j=0;j<nmat+1;j++){
      printf("%.16e ",Sarr[i][j]);
    }
    printf("\n");
  }
  //*/

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST REFLECTION (ALL VALUES)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*
  //Planet
  char line[MAX_STRING];
  
  sprintf(filename,"fou_gasplanet.dat");
  read_fourier(filename,shape,xmu,rfou,rtra);
  nmat=shape[0];
  nmugs=shape[1];
  nfou=shape[2];
  printf("nmat = %d, nmugs = %d, nfou = %d\n",shape[0],shape[1],shape[2]);
  f=fopen("planet-interpolation.mat","r");
  double tmp;
  double Spixx[4];
  double dif,difmax=-1e100;
	       
  npix=1;
  i=0;
  while(fgets(line,sizeof(line),f)){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	   &phi[0],&beta[0],&theta0[0],&theta[0],
	   &tmp,&tmp,&tmp,
	   &apix[0],
	   &Spixx[0],&Spixx[1],&Spixx[2],&Spixx[3]);
    //Compute with new routine
    reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);

    if(VERBOSE){
      printf("phi = %e, beta = %e, theta0 = %e, theta = %e, apix = %e\n",
	     phi[0],beta[0],theta0[0],theta[0],apix[0]);
      printf("S (pixx) = %e, %e, %e, %e\n",
	     Spixx[0],Spixx[1],Spixx[2],Spixx[3]);
      printf("S (cpixx) = %e, %e, %e, %e\n",
	     Sarr[0][0],Sarr[0][1],Sarr[0][2],Sarr[0][3]);
    }  
    //Compute difference
    dif=0;
    for(j=0;j<nmat+1;j++){
      dif+=fabs(Spixx[j]-Sarr[0][j]);
    }
    if(VERBOSE) printf("Difference: %e\n",dif);
    difmax=dif>difmax?dif:difmax;
    //if(i>100) break;
    if((i%50000)==0)
      printf("We are in %d\n",i);
    i++;
  }  
  printf("Maximum difference: %e\n",difmax);
  //*/

  /*
  //Ring forward
  char line[MAX_STRING];
  
  sprintf(filename,"fou_ring_0_4_0_8.dat");
  read_fourier(filename,shape,xmu,rfou,rtra);
  nmat=shape[0];
  nmugs=shape[1];
  nfou=shape[2];
  printf("nmat = %d, nmugs = %d, nfou = %d\n",shape[0],shape[1],shape[2]);
  FILE *f=fopen("ring_forward-interpolation.mat","r");
  double tmp;
  double Spixx[4];
  double dif,difmax=-1e100;
	       
  npix=1;
  i=0;
  while(fgets(line,sizeof(line),f)){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	   &phi[0],&beta[0],&theta0[0],&theta[0],
	   &tmp,&tmp,&tmp,
	   &apix[0],
	   &Spixx[0],&Spixx[1],&Spixx[2],&Spixx[3]);
    //Compute with new routine
    reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);

    if(VERBOSE){
      printf("phi = %e, beta = %e, theta0 = %e, theta = %e, apix = %e\n",
	     phi[0],beta[0],theta0[0],theta[0],apix[0]);
      printf("S (pixx) = %e, %e, %e, %e\n",
	     Spixx[0],Spixx[1],Spixx[2],Spixx[3]);
      printf("S (cpixx) = %e, %e, %e, %e\n",
	     Sarr[0][0],Sarr[0][1],Sarr[0][2],Sarr[0][3]);
    }  
    //Compute difference
    dif=0;
    for(j=0;j<nmat+1;j++){
      dif+=fabs(Spixx[j]-Sarr[0][j]);
    }
    if(VERBOSE) printf("Difference: %e\n",dif);
    difmax=dif>difmax?dif:difmax;
    //if(i>100) break;
    if((i%50000)==0)
      printf("We are in %d\n",i);
    i++;
    //if(i>100) break;
  }  
  printf("Maximum difference: %e\n",difmax);
  //*/

  //*
  //Ring backward
  char line[MAX_STRING];
  
  sprintf(filename,"fou_ring_0_4_0_8.dat");
  read_fourier(filename,shape,xmu,rfou,rtra);
  nmat=shape[0];
  nmugs=shape[1];
  nfou=shape[2];
  printf("nmat = %d, nmugs = %d, nfou = %d\n",shape[0],shape[1],shape[2]);
  f=fopen("ring_back-interpolation.mat","r");

  npix=1;
  i=0;
  while(fgets(line,sizeof(line),f)){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	   &phi[0],&beta[0],&theta0[0],&theta[0],
	   &tmp,&tmp,&tmp,
	   &apix[0],
	   &Spixx[0],&Spixx[1],&Spixx[2],&Spixx[3]);
    //Compute with new routine
    reflection(npix,phi,beta,theta0,theta,shape,xmu,rtra,apix,Sarr);

    if(VERBOSE){
      printf("phi = %e, beta = %e, theta0 = %e, theta = %e, apix = %e\n",
	     phi[0],beta[0],theta0[0],theta[0],apix[0]);
      printf("S (pixx) = %e, %e, %e, %e\n",
	     Spixx[0],Spixx[1],Spixx[2],Spixx[3]);
      printf("S (cpixx) = %e, %e, %e, %e\n",
	     Sarr[0][0],Sarr[0][1],Sarr[0][2],Sarr[0][3]);
    }  
    //Compute difference
    dif=0;
    for(j=0;j<nmat+1;j++){
      dif+=fabs(Spixx[j]-Sarr[0][j]);
    }
    if(VERBOSE) printf("Difference: %e\n",dif);
    difmax=dif>difmax?dif:difmax;
    //if(i>100) break;
    if((i%50000)==0)
      printf("We are in %d\n",i);
    i++;
    //if(i>100) break;
  }  
  printf("Maximum difference: %e\n",difmax);
  //*/

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST BRACK
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*
  FILE *f=fopen("xmu.mat","r");
  double xmu[21];

  for(i=0;i<20;i++){
    fscanf(f,"%lf",&xmu[i]);
    printf("%.17lf\n",xmu[i]);
  }
  */
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST SPLINE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*
  #define n 10
  double x[n],y[n],y2[n];
  for(i=0;i<n;i++){
    x[i]=i+1;
    y[i]=x[i]*x[i];
  }
  spline(x,y,n,y2);
  for(i=0;i<n;i++){
    printf("%d %.17lf %.17lf %.17lf\n",i,x[i],y[i],y2[i]);
  }
  
  double xv,yv;
  for(i=1;i<n;i++){
    xv=i+0.3;
    yv=splint(x,y,y2,n,xv);
    printf("%.17lf %.17lf\n",xv,yv);
  }
  //*/
}

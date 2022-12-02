#include "cpixx.c"

#define VERBOSE 0
#define MAX_MAT 5
#define MAX_PIX 1000000

enum OBJS {PLANET,RING_REFLECTION,RING_TRANSMISSION};

//////////////////////////////////////////////////////////////
// MAIN PROGRAM
//////////////////////////////////////////////////////////////
int main(void)
{
  int type;
  //type=PLANET;
  //type=RING_REFLECTION;
  type=RING_TRANSMISSION;

  FILE *f;

  int i,j,k,m;
  int nmat,nmugs,nfou;
  int npix;

  char filename[MAX_STRING];
  char line[MAX_STRING];

  struct FourierCoefficients F;
  
  double check_sum;
  double *phi,*beta,*theta0,*theta,*apix;
  double **Sarr;
  double **Spixx;

  double tmp;
  double dif,difmax=-1e100;
  int qreflection;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // READ FOURIER
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("================================================================================\n");
  printf("READ TESTING:\n");
  printf("================================================================================\n");
  if(type==PLANET){
    sprintf(filename,"fou_gasplanet.dat");
    f=fopen("planet-interpolation.mat","r");
    qreflection=1;
  }
  if(type==RING_REFLECTION){
    sprintf(filename,"fou_ring_0_4_0_8.dat");
    f=fopen("ring_forward-interpolation.mat","r");
    qreflection=1;
  }
  if(type==RING_TRANSMISSION){
    sprintf(filename,"fou_ring_0_4_0_8.dat");
    f=fopen("ring_back-interpolation.mat","r");
    qreflection=0;
  }
  check_sum=read_fourier(filename,&F);
  nmat=F.nmat;nmugs=F.nmugs;nfou=F.nfou;
  printf("Main: nmat = %d, numgs = %d, nfou = %d\n",nmat,nmugs,nfou);
  printf("Checksum = %.16e\n",check_sum);
  
  /*
    Checksum of a Fourier coefficient file

        python -c "import numpy as np;data=np.loadtxt('fou_gasplanet.dat',skiprows=43);print(data[:,3:].sum())"

    where 43 is the starting point of data. For the ring Fourier coefficients use:

        python -c "import numpy as np;data=np.loadtxt('fou_ring_0_4_0_8.dat',skiprows=22);print(data[:,3:].sum())"
   */

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // CHECK INTERPOLATION SINGLE VALUE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("================================================================================\n");
  printf("REFLECTION TEST (SINGLE VALUE):\n");
  printf("================================================================================\n");
  double **rmat;

  npix=1;
  Sarr=zeros_matrix(npix,nmat+1);
  phi=zeros_vector(npix);
  beta=zeros_vector(npix);
  theta0=zeros_vector(npix);
  theta=zeros_vector(npix);
  apix=zeros_vector(npix);
  
  if(type==PLANET){
    phi[0]=1.0448451569439827;
    beta[0]=3.069394277348945;
    theta0[0]=0.04990329026929557;
    theta[0]=0.02509670973070432;
    apix[0]=9.432328787795567e-05;
    printf("Expected values:\n4.597857424902560283e-07 2.251229972872198058e-07 1.834400800563127439e-09 4.896421313424954569e-01\n");
  }
  if(type==RING_REFLECTION){
    phi[0]=1.490116119384765625e-08;
    beta[0]=0.000000000000000000e+00;
    theta0[0]=4.999999999999998890e-01;
    theta[0]=5.000000000000000000e-01;
    apix[0]=1.163314390931110409e-04;
    printf("Expected:\n6.193646058775441171e-06 -3.046132218287162406e-07 -9.759550180589223642e-15 4.918156751904256135e-02\n");
  }
  if(type==RING_TRANSMISSION){
    phi[0]=1.601029385538801364e+00;
    beta[0]=1.601029385538801364e+00;
    theta0[0]=1.744974835125044643e-02;
    theta[0]=5.000000000000000000e-01;
    apix[0]=1.163314390931110409e-04;
    printf("Expected:\n1.688771016436214060e-07 -1.485273114913191718e-08 2.674513729889046925e-09 8.936444506313186154e-02\n");
  }
  reflection(F,qreflection,npix,phi,beta,theta0,theta,apix,Sarr);
  printf("Calculated values:\n");
  for(i=0;i<npix;i++){
    for(j=0;j<nmat+1;j++){
      printf("%.16e ",Sarr[i][j]);
    }
    printf("\n");
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // CHECK INTERPOLATION VALUES MASSIVELY
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ç
  printf("================================================================================\n");
  printf("REFLECTION TEST (MULTIPLE VALUE):\n");
  printf("================================================================================\n");
  printf("nmat = %d, nmugs = %d, nfou = %d\n",nmat,nmugs,nfou);

  i=0;
  phi=zeros_vector(MAX_PIX);
  beta=zeros_vector(MAX_PIX);
  theta0=zeros_vector(MAX_PIX);
  theta=zeros_vector(MAX_PIX);
  apix=zeros_vector(MAX_PIX);
  Spixx=zeros_matrix(MAX_PIX,MAX_MAT+1);
  while(fgets(line,sizeof(line),f)){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	   &phi[i],&beta[i],&theta0[i],&theta[i],
	   &tmp,&tmp,&tmp,
	   &apix[i],
	   &Spixx[i][0],&Spixx[i][1],&Spixx[i][2],&Spixx[i][3]);
    i++;
  }
  npix=i;
  Sarr=zeros_matrix(npix,nmat+1);
  printf("Comparing %d values...\n",npix);
  reflection(F,qreflection,npix,phi,beta,theta0,theta,apix,Sarr);
  
  for(i=0;i<npix;i++){
    dif=0;
    for(j=0;j<nmat+1;j++){
      dif+=fabs(Spixx[i][j]-Sarr[i][j]);
    }
    difmax=dif>difmax?dif:difmax;
  }
  printf("Maximum difference: %e\n",difmax);

  return 0;
}

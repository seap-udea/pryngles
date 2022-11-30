#include "cpixx_dev.c"

//////////////////////////////////////////////////////////////
// MAIN PROGRAM
//////////////////////////////////////////////////////////////
int main(void)
{
  int i;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // TEST BRACK
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE *f=fopen("xmu.mat","r");
  double xmu[21];

  for(i=0;i<20;i++){
    fscanf(f,"%lf",&xmu[i]);
    printf("%.17lf\n",xmu[i]);
  }
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
  */
}

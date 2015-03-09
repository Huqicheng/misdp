#include <stdio.h>

extern "C"{
  void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
}

#define D 4

int main()
{
  FILE *out;

  double  A[D*D];
 
  int M, N, LDA, IPIV[D], INFO;

 int i, j;

 N = D;
 M = D;
 LDA = M;

A[0+0*D]=1, A[0+D*1]=2, A[0+2*D]=3, A[0+D*3]=4;
A[1+D*0]=2, A[1+D*1]=3, A[1+D*2]=4, A[1+D*3]=1;
A[2+D*0]=3, A[2+D*1]=4, A[2+D*2]=1, A[2+D*3]=2;
A[3+D*0]=4, A[3+D*1]=1, A[3+D*2]=2, A[3+D*3]=3;




out = stderr;

fprintf(out, "A before call:\n");
for(i=0;i<D; i++)
  {
    for(j=0;j<D; j++)
      fprintf(out, "%f\t", A[i+D*j]);
    fprintf(out, "\n");
  }


 dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);


fprintf(out, "INFO: %d\n", INFO);

fprintf(out, "A after call:\n");
for(i=0;i<D; i++)
  {
    for(j=0;j<D; j++)
      fprintf(out, "%f\t", A[i+D*j]);
    fprintf(out, "\n");
  }


fprintf(out, "VECR:\n");
 double piii=1;
 for(i=0;i<D; i++)
   {
     j=i;
     fprintf(out, "%f\t", A[i+D*j]);
     piii=piii*A[i+D*j];
  }
 fprintf(out, "%f\n",piii);
 fprintf(out, "P:\n");
 for(j=0;j<D; j++)
   fprintf(out, "%d\n", IPIV[j]);




return 0;
}

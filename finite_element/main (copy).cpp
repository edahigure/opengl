#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//compile as  gcc main.cpp -o finite
 
double* new_double(int N)
{
   return (double*)malloc( (size_t)(N*sizeof(double)));
}

int* new_int(int N)
{
   return (int*)malloc( (size_t)(N*sizeof(int)));
}


int **new_matrix_int( long nrh,  long nch)
/* allocate a double matrix with subscript range m[0..nrh][0..nch] */
{
	long i, nrow=nrh+1,ncol=nch+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow)*sizeof(int*)));
	if (!m) fprintf(stderr,"allocation failure");

	/* allocate rows and set pointers to them */
	m[0]=(int *) malloc((size_t)((nrow*ncol)*sizeof(int)));
	if (!m[0]) fprintf(stderr,"allocation failure");

	for(i=1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


double **new_matrix_double( long nrh,  long nch)
/* allocate a double matrix with subscript range m[0..nrh][0..nch] */
{
	long i, nrow=nrh+1,ncol=nch+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow)*sizeof(double*)));
	if (!m) fprintf(stderr,"allocation failure");

	/* allocate rows and set pointers to them */
	m[0]=(double *) malloc((size_t)((nrow*ncol)*sizeof(double)));
	if (!m[0]) fprintf(stderr,"allocation failure");

	for(i=1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}



struct element 
{

  int Id;   
  int i,  j,  k;
  int ei, ej, ek;
  int si, sj, sk;
  double xV, yV;
} ;

struct node 
{
  int Id;   
  double x, y;
  int mark;
} ;




typedef struct element ELEMENT;
typedef struct node NODE;



ELEMENT* new_element(int N)
{
   return (ELEMENT*)malloc( (size_t)(N*sizeof(ELEMENT)));
}
 
NODE* new_node(int N)
{
   return (NODE*)malloc( (size_t)(N*sizeof(NODE)));
}



int Ne,Nn;

ELEMENT *ele;
NODE *nod;

int **nod_ele;
int *N_nod_len;
double **K;
double *F;
double *FF;




void getK(void)
{

   double ai[3],bi[3],ci[3],Ae;
   double xi,yi;
   double xj,yj;
   double xk,yk;

   int ni,nj,nk;
   double fac_1=0.0;
   double temp;

   double KK[3][3];


   fprintf(stderr,"Here Ok\n");


   for(int k=0;k<3;k++){
   for(int i=0;i<3;i++){
   
      KK[k][i]=0.0;
     
   }}


   for (int J=0;J<Ne;J++){
      
      xi= nod[ ele[J].i ].x ;
      yi= nod[ ele[J].i ].y ;

      xj= nod[ ele[J].j ].x ;
      yj= nod[ ele[J].j ].y ;

      xk= nod[ ele[J].k ].x ;
      yk= nod[ ele[J].k ].y ;


      ai[0]= xj*yk-xk*yj;
      bi[0]= yj-yk;
      ci[0]= xk-xj;


      Ae=0.5*(  (xj*yk-xk*yj) - xi* (yk-yj) + yi*(xk-xj) ) ;


      xi= nod[ ele[J].j ].x ;
      yi= nod[ ele[J].j ].y ;

      xj= nod[ ele[J].k ].x ;
      yj= nod[ ele[J].k ].y ;

      xk= nod[ ele[J].i ].x ;
      yk= nod[ ele[J].i ].y ;


      ai[1]= xj*yk-xk*yj;
      bi[1]= yj-yk;
      ci[1]= xk-xj;


      xi= nod[ ele[J].k ].x ;
      yi= nod[ ele[J].k ].y ;

      xj= nod[ ele[J].i ].x ;
      yj= nod[ ele[J].i ].y ;

      xk= nod[ ele[J].j ].x ;
      yk= nod[ ele[J].j ].y ;


      ai[2]= xj*yk-xk*yj;
      bi[2]= yj-yk;
      ci[2]= xk-xj;


      ni=ele[J].i;
      nj=ele[J].j;
      nk=ele[J].k;

 
 
      KK[0][0]=( bi[0] * bi[0] + ci[0] * ci[0] )*1.0/(4.0*Ae);
      KK[0][1]=( bi[0] * bi[1] + ci[0] * ci[1] )*1.0/(4.0*Ae);
      KK[0][2]=( bi[0] * bi[2] + ci[0] * ci[2] )*1.0/(4.0*Ae);

      KK[1][0]=( bi[1] * bi[0] + ci[1] * ci[0] )*1.0/(4.0*Ae);
      KK[1][1]=( bi[1] * bi[1] + ci[1] * ci[1] )*1.0/(4.0*Ae);
      KK[1][2]=( bi[1] * bi[2] + ci[1] * ci[2] )*1.0/(4.0*Ae);

      KK[2][0]=( bi[2] * bi[0] + ci[2] * ci[0] )*1.0/(4.0*Ae);
      KK[2][1]=( bi[2] * bi[1] + ci[2] * ci[1] )*1.0/(4.0*Ae);
      KK[2][2]=( bi[2] * bi[2] + ci[2] * ci[2] )*1.0/(4.0*Ae);



      K[ni][ni]  = K[ni][ni]  + KK[0][0];
      K[ni][nj]  = K[ni][nj]  + KK[0][1];
      K[ni][nk]  = K[ni][nk]  + KK[0][2];


      K[nj][ni]  = K[nj][ni]  + KK[1][0];
      K[nj][nj]  = K[nj][nj]  + KK[1][1];
      K[nj][nk]  = K[nj][nk]  + KK[1][2];

      K[nk][ni]  = K[nk][ni]  + KK[2][0];
      K[nk][nj]  = K[nk][nj]  + KK[2][1];
      K[nk][nk]  = K[nk][nk]  + KK[2][2];

   }

}






// local assembly of the load vector


void getF(void)
{

   double ai[3],bi[3],ci[3],Ae;
   double xi,yi;
   double xj,yj;
   double xk,yk;
   double xc,yc;
   double xd,yd;


   for (int J=0;J<Ne;J++){

      xi= nod[ ele[J].i ].x ;
      yi= nod[ ele[J].i ].y ;

      xj= nod[ ele[J].j ].x ;
      yj= nod[ ele[J].j ].y ;

      xk= nod[ ele[J].k ].x ;
      yk= nod[ ele[J].k ].y ;


      Ae=0.5*(  (xj*yk-xk*yj) - xi* (yk-yj) + yi*(xk-xj) ) ;


      FF[ele[J].i]= FF[ele[J].i]+ Ae/12.0*(2*F[ele[J].i] +   F[ele[J].j]+  F[ele[J].k]) ;
      FF[ele[J].j]= FF[ele[J].j]+ Ae/12.0*(  F[ele[J].i] + 2*F[ele[J].j]+  F[ele[J].k]) ;
      FF[ele[J].k]= FF[ele[J].k]+ Ae/12.0*(  F[ele[J].i] +   F[ele[J].j]+2*F[ele[J].k]) ;

   }

// the load vector as in your article
      FF[1]=0.0104;
      FF[2]=0.0417;
      FF[3]=0.0208;
      FF[4]=0.0417;
      FF[5]=0.0104;
      FF[6]=0.0417;
      FF[7]=0.0417;
      FF[8]=0.0833;
      FF[9]=0.0417;
      FF[10]=0.0417;
      FF[11]=0.0208;
      FF[12]=0.0833;
      FF[13]=0.0417;
      FF[14]=0.0833;
      FF[15]=0.0208;
      FF[16]=0.0417;
      FF[17]=0.0417;
      FF[18]=0.0833;
      FF[19]=0.0417;
      FF[20]=0.0417;
      FF[21]=0.0104;
      FF[22]=0.0417;
      FF[23]=0.0208;
      FF[24]=0.0417;
      FF[25]=0.0104;
 
   for(int i=1;i<=Nn;i++)
     FF[i]=-FF[i];


}




//for solving the linear system

#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=new_double(n+1);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) fprintf(stderr,"Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free(vv);
}
#undef TINY



void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


double * F_jordan(double *x)
{
     double *y;
     y=new_double(Nn+1);
     for(int k=1;k<=Nn;k++){
        y[k]=0.0;
        for(int j=1;j<=Nn;j++){

         if(k!=j) 
         y[k] = y[k] - K[k][j]*x[j];

        }
         y[k]= y[k]+ FF[k] ;
     }
     return y;
}


double * F_gauss_seidel(double *X)
{
     double *x;
     x=new_double(Nn+1);
     double sum_1;
     double sum_2;

     for(int i=1;i<=Nn;i++){
        
        sum_1=0.0;
        sum_2=0.0;

        for( int j=1;j<= i-1;j++)
        sum_1 = sum_1 + K[i][j]*x[j];

        for( int j=i+1;j<= Nn;j++)
        sum_2 = sum_2 + K[i][j]*X[j];
   
 
        x[i]= 1.0/K[i][i]*(-sum_1-sum_2+FF[i]);

     }

     return x;
}


//the main code


void GetM_sub(double **M,int n,int k,double **Mp)
{
   int ip,jp;

   ip=0;
   for(int i=1;i<n;i++){
      jp=0;
      for(int j=0;j<n;j++){
            if( j!= k){
                Mp[ip][jp]=M[i][j];
                jp++;
            }
            
      }
      ip++;
   }

}


double det(double **M,int n)
{

   double sum;
   double prod;
   double res;
   double **Mp;




   Mp=new_matrix_double(n-1,n-1);
   for(int i=0;i<n-1;i++){
   for(int j=0;j<n-1;j++){
       Mp[i][j]=0.0;
   }}


   
   prod=1.0;
 
   if(n==2){
       res=M[0][0]*M[1][1]-M[1][0]*M[0][1];
       free(Mp);
       return res;
   }
   else{
       for(int k=0;k<n;k++){
           if(M[0][k]!=0){
              GetM_sub(M,n,k,Mp);
              sum=sum + prod*M[0][k]*det(Mp,n-1);
              prod=prod*(-1);      
           }
       
       }
       free(Mp);
       return sum;       
   }  
}


double *Solve(double **M,double *B, int n)
{
   double *X;
   double delta;
   double **Mp;
   X=new_double(n);

   delta=det(M,n);

   fprintf(stderr,"det: %f\n",delta);

   Mp=new_matrix_double(n,n);

   
   for(int k=0;k<n;k++){

      for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
          Mp[i][j]=M[i][j];
      }}

      for(int j=0;j<n;j++){
        Mp[j][k]=B[j];
      }   
      X[k]=det(Mp,n)/delta;
       
      fprintf(stderr,"X[%d]: %f\n",k,X[k]);
 
      
   }

     
   return X;
    

}



int main(int argc, char *argv[])
{

   int dum,cont;
   FILE *fd;



   double d;
   int n,*indx;


   double **A;
   double *B;
   double *sol;

   double delta;






   //opening the elements file

   fd=fopen("example.e","r");
   fscanf(fd,"%d\n",&Ne);

   fprintf(stderr,"%d\n",Ne); 
   
   ele=new_element(Ne);

   fprintf(stderr,"elements\n");

   for(int  k=0;k<Ne;k++)
   {
      fscanf(fd,"%d %d %d %d \n",
      &ele[k].Id,&ele[k].i,&ele[k].j,&ele[k].k);

   }
   fclose(fd);

   //opening the nodes file   

   fd=fopen("example.n","r");
   fscanf(fd,"%d\n",&Nn);
   fprintf(stderr,"%d\n",Nn); 


   A= new_matrix_double(Nn+1,Nn+1);
   B=new_double(Nn+1);
   indx=new_int(Nn+1);


   nod=new_node(Nn+1);





   fprintf(stderr,"nodes\n");
   for(int  k=1;k<=Nn;k++)
   {
      dum=fscanf(fd     ,"%d %lf %lf %d\n", &nod[k].Id,&nod[k].x ,&nod[k].y,&nod[k].mark ) ;
      fprintf(stderr,"%d %lf %lf %d\n", nod[k].Id,nod[k].x ,nod[k].y,nod[k].mark ) ;
   }   

   fclose(fd);



   K=new_matrix_double(Nn+1,Nn+1);
   F=new_double(Nn+1);
   FF=new_double(Nn+1);


   
   for(int k=1;k<=Nn;k++){
      FF[k]=0;
      F[k]=-1.0; 
   }


   for(int k=1;k<=Nn;k++){
   for(int i=1;i<=Nn;i++){
   
      K[k][i]=0.0;
     
   }}
  
   getK();   //calculates the stiffness  Matrix


   fprintf(stderr,"stiffness  Matrix\n");
   for(int k=1;k<=Nn;k++){  
      for(int i=1;i<=Nn;i++){
         fprintf(stderr,"K(%d,%d)=%f;\n",k,i,K[k][i] );;
      }
   }

   getF();       //calculates the load vector

   fprintf(stderr,"load vector\n");
   for(int i=1;i<=Nn;i++)
      fprintf(stderr,"FF(%d)=%e;\n",i,FF[i]);        





   


   for(int i=1;i<=Nn;i++){   
   for(int j=1;j<=Nn;j++){
      A[i][j]=K[i][j];
    
   }}


   ludcmp(A,Nn,indx,&d);

   delta=1;
   for(int i=1;i<=Nn;i++){   
     delta=delta*A[i][i];    
     fprintf(stderr,"A : %f \n",A[i][i]);
   }
   

   for (int jj=1;jj<=Nn;jj++){

      for(int i=1;i<=Nn;i++){   
      for(int j=1;j<=Nn;j++){
         A[i][j]=K[i][j];    
      }}

   
      for(int i=1;i<=Nn;i++){
         A[i][jj]=K[i][jj];    
      }

      ludcmp(A,Nn,indx,&d);

      delta=1;
      for(int i=1;i<=Nn;i++){   
        delta=delta*A[i][i];    
      }
      fprintf(stderr,"delta :%d %f \n",jj,delta);

   }



    
//   sol=Solve(A,B,Nn);
//   fprintf(stderr,"det=%f",det(A,3) );

   

   


}






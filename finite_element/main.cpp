#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//compile as  gcc main.cpp -o finite  -lm
 
double* new_double(int N)
{
   return (double*)malloc( (size_t)(N*sizeof(double)));
}

int* new_int(int N)
{
   return (int*)malloc( (size_t)(N*sizeof(int)));
}


int **new_matrix_int(long nrh, long nch)
/* allocate an int matrix with subscript range m[0..nrh][0..nch] */
{
    long i, nrow = nrh + 1, ncol = nch + 1;
    int **m;

    /* allocate pointers to rows */
    m = (int **)malloc((size_t)(nrow * sizeof(int*)));
    if (!m) {
        fprintf(stderr, "Error: Cannot allocate matrix pointers (%ld rows)\n", nrow);
        return NULL;
    }

    /* allocate rows and set pointers to them */
    m[0] = (int *)calloc((size_t)(nrow * ncol), sizeof(int)); // calloc INICIALIZA A CERO
    if (!m[0]) {
        fprintf(stderr, "Error: Cannot allocate matrix data (%ld x %ld elements)\n", nrow, ncol);
        free(m); // ¡LIBERAR memoria ya asignada!
        return NULL;
    }

    for(i = 1; i <= nrh; i++) {
        m[i] = m[i-1] + ncol;
    }

    return m;
}

double **new_matrix_double(long nrh, long nch)
{
    long i, nrow = nrh + 1, ncol = nch + 1;
    double **m;

    // Asignar array de punteros
    m = (double **)malloc(nrow * sizeof(double*));
    if (!m) {
        fprintf(stderr, "Error: Cannot allocate matrix pointers (%ld rows)\n", nrow);
        return NULL;
    }

    // Asignar memoria contigua para todos los datos
    m[0] = (double *)calloc(nrow * ncol, sizeof(double)); // calloc INICIALIZA A CERO
    if (!m[0]) {
        fprintf(stderr, "Error: Cannot allocate matrix data (%ld x %ld = %ld elements)\n", 
                nrow, ncol, nrow * ncol);
        free(m); // ¡IMPORTANTE! Liberar lo ya asignado
        return NULL;
    }

    // Configurar punteros a filas
    for(i = 1; i <= nrh; i++) {
        m[i] = m[i-1] + ncol;
    }

    return m;
}

void free_matrix_int(int **m)
{
    if (m != NULL) {
        if (m[0] != NULL) {
            free(m[0]);
        }
        free(m);
    }
}

void free_matrix_double(double **m) {
    if (m) {
        if (m[0]) free(m[0]);
        free(m);
    }
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
    ELEMENT *elem = (ELEMENT*)calloc(N, sizeof(ELEMENT)); // calloc INICIALIZA A CERO
    if (!elem) {
        fprintf(stderr, "Error: Cannot allocate %d elements\n", N);
        return NULL;
    }
    return elem;
}

NODE* new_node(int N) 
{
    NODE *nodes = (NODE*)calloc(N, sizeof(NODE)); // calloc INICIALIZA A CERO
    if (!nodes) {
        fprintf(stderr, "Error: Cannot allocate %d nodes\n", N);
        return NULL;
    }
    return nodes;
}


void free_element(ELEMENT *elem) {
    if (elem) free(elem);
}

void free_node(NODE *node_ptr) {
    if (node_ptr) free(node_ptr);
}

int Ne = 0, Nn = 0;
ELEMENT *ele = NULL;
NODE *nod = NULL;
int **nod_ele = NULL;
int *N_nod_len = NULL;
double **K = NULL;
double *F = NULL;
double *FF = NULL;

void applyBoundaryConditions(void)
{
    for(int i = 1; i <= Nn; i++) {
        if(nod[i].mark != 0) {  // Nodo de frontera
            // Condición Dirichlet: u = 0
            for(int j = 1; j <= Nn; j++) {
                K[i][j] = 0.0;
                K[j][i] = 0.0;
            }
            K[i][i] = 1.0;
            FF[i] = 0.0;
        }
    }
}

void getK(void)
{
    double bi[3], ci[3], Ae;
    double x[3], y[3];
    int nodes[3];
    double KK[3][3];

    for (int J = 0; J < Ne; J++) {
        // Get node coordinates
        nodes[0] = ele[J].i;
        nodes[1] = ele[J].j;
        nodes[2] = ele[J].k;
        
        x[0] = nod[nodes[0]].x; y[0] = nod[nodes[0]].y;
        x[1] = nod[nodes[1]].x; y[1] = nod[nodes[1]].y;
        x[2] = nod[nodes[2]].x; y[2] = nod[nodes[2]].y;

        // Area calculation (with absolute value)
        Ae = 0.5 * fabs((x[1]*y[2] - x[2]*y[1]) - x[0]*(y[2] - y[1]) + y[0]*(x[2] - x[1]));
        
        // Correct shape function derivatives
        bi[0] = y[1] - y[2]; ci[0] = x[2] - x[1];
        bi[1] = y[2] - y[0]; ci[1] = x[0] - x[2];
        bi[2] = y[0] - y[1]; ci[2] = x[1] - x[0];

        // Element stiffness matrix
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                KK[i][j] = (bi[i] * bi[j] + ci[i] * ci[j]) / (4.0 * Ae);
            }
        }

        // Assembly to global matrix
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                K[nodes[i]][nodes[j]] += KK[i][j];
            }
        }
    }
}





// local assembly of the load vector


void getF(void)
{
    // -------------------------------------------------------------------
    // PARTE 1: Ensamblar vector de carga por el método de elementos finitos
    // -------------------------------------------------------------------
    // F[] debe contener los valores de la función fuente f(x,y) evaluada en los nodos
    // Usamos integración nodal con la fórmula: ∫fφᵢ ≈ (Ae/12)*(2fᵢ + fⱼ + fₖ)
    // -----------------

    // Inicializar a cero
    for(int i = 1; i <= Nn; i++) {
        FF[i] = 0.0;
    }
    
    for (int J = 0; J < Ne; J++) {
        int i = ele[J].i, j = ele[J].j, k = ele[J].k;
        
        double xi = nod[i].x, yi = nod[i].y;
        double xj = nod[j].x, yj = nod[j].y;
        double xk = nod[k].x, yk = nod[k].y;
        
        double Ae = 0.5 * fabs((xj*yk - xk*yj) - xi*(yk - yj) + yi*(xk - xj));
        
        // Integración nodal (cuadratura)
        FF[i] += Ae * (2.0 * F[i] + F[j] + F[k]) / 12.0;
        FF[j] += Ae * (F[i] + 2.0 * F[j] + F[k]) / 12.0;
        FF[k] += Ae * (F[i] + F[j] + 2.0 * F[k]) / 12.0;
    }

    for (int i = 1; i <= Nn; i++) {
        FF[i] = -FF[i];  // ← APLICAR EL SIGNO NEGATIVO
    }
    

}



//for solving the linear system

#define TINY 1.0e-20  // ← SIN punto y coma

void ludcmp(double **a, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = new_double(n + 1);
    *d = 1.0;
    
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++) {
            if ((temp = fabs(a[i][j])) > big) big = temp;
        }
        if (big == 0.0) {
            fprintf(stderr, "Singular matrix in routine ludcmp\n");
            free(vv);
            exit(EXIT_FAILURE);  // ← Terminar programa
        }
        vv[i] = 1.0 / big;
    }
    
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }
        
        big = 0.0;
        imax = j;  // ← Inicializar imax
        
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
            
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        
        indx[j] = imax;
        
        if (a[j][j] == 0.0) {
            a[j][j] = TINY;
        }
        
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++) {
                a[i][j] *= dum;
            }
        }
    }
    
    free(vv);
}



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
    FILE *fd;
    
    // 1. Leer malla
    fd = fopen("example.e", "r");
    fscanf(fd, "%d\n", &Ne);
    ele = new_element(Ne);
    
    for(int k = 0; k < Ne; k++) {
        fscanf(fd, "%d %d %d %d\n", &ele[k].Id, &ele[k].i, &ele[k].j, &ele[k].k);
		  fprintf(stderr, "%d %d %d %d\n", ele[k].Id, ele[k].i, ele[k].j, ele[k].k);
    }
    fclose(fd);
    
    fd = fopen("example.n", "r");
    fscanf(fd, "%d\n", &Nn);
    nod = new_node(Nn+1);
    
    for(int k = 1; k <= Nn; k++) {
        fscanf(fd, "%d %lf %lf %d\n", &nod[k].Id, &nod[k].x, &nod[k].y, &nod[k].mark);
    	  fprintf(stderr, "%d %lf %lf %d\n", nod[k].Id, nod[k].x, nod[k].y, nod[k].mark);
    }
    fclose(fd);
    
    // 2. Inicializar matrices
    K = new_matrix_double(Nn+1, Nn+1);
    F = new_double(Nn+1);
    FF = new_double(Nn+1);
    
    // Inicializar a cero
    for(int k = 1; k <= Nn; k++) {
        F[k] = -1.0;  // f(x,y) = -1 para Poisson: -∇²u = -1
        FF[k] = 0.0;
        for(int i = 1; i <= Nn; i++) {
            K[k][i] = 0.0;
        }
    }
    
    // 3. Ensamblar sistema
    getK();  // Ensamblar matriz de rigidez
    getF();  // Ensamblar vector de carga
    
    // 4. Aplicar condiciones de frontera (¡IMPORTANTE!)
    applyBoundaryConditions();
    
    // 5. Resolver sistema
    int *indx = new_int(Nn+1);
    double det;
    
    // Hacer copia de K para preservar original
    double **K_copy = new_matrix_double(Nn+1, Nn+1);
    for(int i = 1; i <= Nn; i++) {
        for(int j = 1; j <= Nn; j++) {
            K_copy[i][j] = K[i][j];
        }
    }
    
    // Descomposición LU
    ludcmp(K_copy, Nn, indx, &det);
    
    // Hacer copia de FF para preservar original
    double *solution = new_double(Nn+1);
    for(int i = 1; i <= Nn; i++) {
        solution[i] = FF[i];
    }
    
    // Resolver
    lubksb(K_copy, Nn, indx, solution);
    
    // 6. Imprimir resultados
    printf("Solución del sistema Poisson -∇²u = -1:\n");
    FILE *fd_out;
    fd_out = fopen("node_f.dat","w");
        fprintf(fd_out,"%d\n", Nn);
    for(int i = 1; i <= Nn; i++) {
        fprintf(fd_out,"%f\n", solution[i]);
    }
    fclose(fd_out);
    
    // 7. Liberar memoria
    free_matrix_double(K_copy);
    free(indx);
    free(solution);
    free_matrix_double(K);
    free(F);
    free(FF);
    free_element(ele);
    free_node(nod);
    
    return 0;
}




#include "util.h"


void convert2Engf(float z,float *z1,int *n)
{
	
   double zt=z;
   int nexp=0;
       		 
   if(fabs(z)>10)
   {
   	  zt=z;
      do
	  {
		 zt=zt/10;
		 nexp++;
	   }while(fabs(zt)>10);
	   
  	}
    if(fabs(z)<1)
    {
   	  zt=z;
      do
	  {
		 zt=zt*10;
		 nexp--;		 
	   }while(fabs(zt)<1);
	   
  	}
   *z1=zt;
   *n=nexp;  	  	     
}


double my_sin(double x) {
	return 1/x+sin(x);
}


double ejemplo_func(double x) {
    return x * x + sin(x);
}


int load_i(FILE *in, int *numb)
{
 char dum, dummy[128];
 int ok=-1;

 for(;;)
  {ok=fscanf(in,"%s", dummy);
   if(dummy[0]=='#' && strlen(dummy)>1 && dummy[strlen(dummy)-1]=='#') {}
   else if(dummy[0]=='#') {do{ok=fscanf(in,"%c", &dum);} while(dum!='#');}
   else                   {*numb=atoi(dummy); break;} }
 return  ok;
}

int load_d(FILE *in, double *numb)
{
 char dum, dummy[128];
 bool ok=false;

 for(;;)
  {ok=fscanf(in,"%s", dummy);
   if(dummy[0]=='#' && strlen(dummy)>1 && dummy[strlen(dummy)-1]=='#') {}
   else if(dummy[0]=='#') {do{ok=fscanf(in,"%c", &dum);} while(dum!='#');}
   else                   {*numb=atof(dummy); break;} }
 return  ok;   
}

double *newvec(long nh) {
    double *v = (double *)malloc(nh * sizeof(double));
    if (v == NULL) {
        // Manejo de error (ej: imprimir mensaje y salir)
        fprintf(stderr, "Error: No se pudo asignar memoria.\n");
        exit(1);
    }
    return v;  // Retorna v[0..nh-1]
}

char * newstr(int N)
{
  return (char*) malloc((N)*sizeof(char) );
}

// Función que evalúa f(x) en [xi, xf] y guarda en archivo
void save_func(
    double (*func)(double), // Puntero a función
    double xi,             // Inicio del rango
    double xf,             // Fin del rango
    int puntos,           // Número de puntos
    char *filename  // Nombre del archivo
) {
    double paso = (xf - xi) / (puntos - 1);
    double x, y;
    double y_min = FLT_MAX;  // Inicializar con el máximo valor posible
    double y_max = -FLT_MAX; // Inicializar con el mínimo valor posible

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error al abrir el archivo.\n");
        return;
    }


	 fprintf(file, "%d\n", puntos);

    for (int i = 0; i < puntos; i++) {
        x = xi + i * paso;
        y = func(x);

        // Guardar en el archivo
        fprintf(file, "%lf\t%lf\n", x, y);


    }

    fclose(file);

}




 

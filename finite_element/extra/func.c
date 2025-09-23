#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

}

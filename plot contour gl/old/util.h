#ifndef __UTIL_H__
#define __UTIL_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h> 

double my_sin(double x);
double ejemplo_func(double x);
void save_func(
    double (*func)(double), // Puntero a función
    double xi,             // Inicio del rango
    double xf,             // Fin del rango
    int puntos,           // Número de puntos
    char *filename  // Nombre del archivo
);

char * newstr(int N);

int load_i(FILE *in, int *numb);
int load_d(FILE *in, double *numb);
double *newvec(long nh) ;

void convert2Engf(float z,float *z1,int *n);


#endif

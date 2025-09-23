#ifndef MESH_READER_H
#define MESH_READER_H

#define MAX_NODES 10000
#define MAX_ELEMENTS 10000
#define MAX_SIDES 10000
#define GREAT 1e10
#define SMALL 1e-10
#define OFF -1

#include <GL/gl.h>

typedef struct {
    int index;
    double x, y;
    int mark;
} Node;

typedef struct {
    int index;
    int i, j, k;
    int ei, ej, ek;
    int si, sj, sk;
    double xv, yv;
    int material;
} Element;

typedef struct {
    int index;
    int a, b, c, d;
    int ea, eb;
    int mark;
    double s;
} Side;

// VARIABLES GLOBALES (añadir 'extern')
extern Node* nodes;
extern Element* elements; 
extern Side* sides;
extern int Nn, Ne, Ns;
extern double f_min, f_max;

// Funciones
int load_mesh_data(const char* base_name);
void free_mesh_data();
void initDisplayLists();
double min(double a, double b);
double max(double a, double b);

#endif

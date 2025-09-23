#include "mesh_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/gl.h>

// VARIABLES GLOBALES
Node* nod = NULL;
Element* ele = NULL;
Side* sides = NULL;
int Nn = 0, Ne = 0, Ns = 0;
double f_min = 0.0, f_max = 0.0;

// Funciones auxiliares
double min(double a, double b) { return a < b ? a : b; }
double max(double a, double b) { return a > b ? a : b; }

void remove_comments(char* line) {
    char* comment = strchr(line, '#');
    if (comment) *comment = '\0';
}

int load_mesh_data(const char* base_name) {
    char filename[256];
    
    // Construir nombres de archivos
    snprintf(filename, sizeof(filename), "%s.n", base_name);
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error abriendo archivo de nodos: %s\n", filename);
        return 0;
    }
    
    char line[256];
    
    // Leer número de nodos
    fgets(line, sizeof(line), file);
    remove_comments(line);
    Nn = atoi(line);
    nod = (Node*)malloc(Nn * sizeof(Node));
    
    // Leer nodos
    for (int i = 0; i < Nn; i++) {
        fgets(line, sizeof(line), file);
        remove_comments(line);
        sscanf(line, "%d: %lf %lf %d", 
               &nod[i].Id, &nod[i].x, &nod[i].y, &nod[i].mark);
    }
    fclose(file);
    
    // Leer elementos
    snprintf(filename, sizeof(filename), "%s.e", base_name);
    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error abriendo archivo de elementos: %s\n", filename);
        return 0;
    }
    
    fgets(line, sizeof(line), file);
    remove_comments(line);
    Ne = atoi(line);
    ele = (Element*)malloc(Ne * sizeof(Element));
    
    for (int i = 0; i < Ne; i++) {
        fgets(line, sizeof(line), file);
        remove_comments(line);
        sscanf(line, "%d: %d %d %d %d %d %d %d %d %d %lf %lf %d",
               &ele[i].Id, &ele[i].i, &ele[i].j, &ele[i].k,
               &ele[i].ei, &ele[i].ej, &ele[i].ek,
               &ele[i].si, &ele[i].sj, &ele[i].sk,
               &ele[i].xv, &ele[i].yv, &ele[i].material);
    }
    fclose(file);
    
    // Leer lados
    snprintf(filename, sizeof(filename), "%s.s", base_name);
    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error abriendo archivo de lados: %s\n", filename);
        return 0;
    }
    
    fgets(line, sizeof(line), file);
    remove_comments(line);
    Ns = atoi(line);
    sides = (Side*)malloc(Ns * sizeof(Side));
    
    for (int i = 0; i < Ns; i++) {
        fgets(line, sizeof(line), file);
        remove_comments(line);
        sscanf(line, "%d: %d %d %d %d %d",
               &sides[i].Id, &sides[i].c, &sides[i].d,
               &sides[i].ea, &sides[i].eb, &sides[i].mark);
    }
    fclose(file);
    
    return 1;
}

void free_mesh_data() {
    if (nod) free(nod);
    if (ele) free(ele);
    if (sides) free(sides);
    nod = NULL;
    ele = NULL;
    sides = NULL;
    Nn = Ne = Ns = 0;
}

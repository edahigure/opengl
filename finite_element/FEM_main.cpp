#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TINY 1.0e-20

// Funciones de asignación de memoria
double* new_double(int N) {
    double *ptr = (double*)malloc(N * sizeof(double));
    if (!ptr) {
        fprintf(stderr, "Error: Cannot allocate %d doubles\n", N);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

int* new_int(int N) {
    int *ptr = (int*)malloc(N * sizeof(int));
    if (!ptr) {
        fprintf(stderr, "Error: Cannot allocate %d ints\n", N);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

double **new_matrix_double(int nrh, int nch) {
    int i, nrow = nrh, ncol = nch;
    double **m = (double **)malloc(nrow * sizeof(double*));
    if (!m) {
        fprintf(stderr, "Error: Cannot allocate matrix pointers (%d rows)\n", nrow);
        return NULL;
    }
    m[0] = (double *)calloc(nrow * ncol, sizeof(double));
    if (!m[0]) {
        fprintf(stderr, "Error: Cannot allocate matrix data (%d x %d)\n", nrow, ncol);
        free(m);
        return NULL;
    }
    for (i = 1; i < nrow; i++) {
        m[i] = m[i-1] + ncol;
    }
    return m;
}

void free_matrix_double(double **m) {
    if (m) {
        if (m[0]) free(m[0]);
        free(m);
    }
}

typedef struct {
    int Id;
    int i, j, k;
    int ei, ej, ek;
    int si, sj, sk;
    double xV, yV;
} ELEMENT;

typedef struct {
    int Id;
    double x, y;
    int mark;
} NODE;

ELEMENT* new_element(int N) {
    ELEMENT *elem = (ELEMENT*)calloc(N, sizeof(ELEMENT));
    if (!elem) {
        fprintf(stderr, "Error: Cannot allocate %d elements\n", N);
        return NULL;
    }
    return elem;
}

NODE* new_node(int N) {
    NODE *nodes = (NODE*)calloc(N, sizeof(NODE));
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
double **K = NULL;
double *F = NULL;
double *FF = NULL;

void applyBoundaryConditions(void) {
    for (int i = 0; i < Nn; i++) {
        if (nod[i].mark != 0) {
            for (int j = 0; j < Nn; j++) {
                K[i][j] = 0.0;
                K[j][i] = 0.0;
            }
            K[i][i] = 1.0;
            FF[i] = 0.0;
        }
    }
}

void getK(void) {
    double bi[3], ci[3], Ae;
    double x[3], y[3];
    int nodes[3];
    double KK[3][3];

    for (int J = 0; J < Ne; J++) {
        nodes[0] = ele[J].i - 1;
        nodes[1] = ele[J].j - 1;
        nodes[2] = ele[J].k - 1;
        
        x[0] = nod[nodes[0]].x; y[0] = nod[nodes[0]].y;
        x[1] = nod[nodes[1]].x; y[1] = nod[nodes[1]].y;
        x[2] = nod[nodes[2]].x; y[2] = nod[nodes[2]].y;

        Ae = 0.5 * fabs((x[1]*y[2] - x[2]*y[1]) - x[0]*(y[2] - y[1]) + y[0]*(x[2] - x[1]));
        
        bi[0] = y[1] - y[2]; ci[0] = x[2] - x[1];
        bi[1] = y[2] - y[0]; ci[1] = x[0] - x[2];
        bi[2] = y[0] - y[1]; ci[2] = x[1] - x[0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                KK[i][j] = (bi[i] * bi[j] + ci[i] * ci[j]) / (4.0 * Ae);
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                K[nodes[i]][nodes[j]] += KK[i][j];
            }
        }
    }
}

void getF(void) {
    for (int i = 0; i < Nn; i++) {
        FF[i] = 0.0;
    }
    
    for (int J = 0; J < Ne; J++) {
        int i = ele[J].i - 1, j = ele[J].j - 1, k = ele[J].k - 1;
        
        double xi = nod[i].x, yi = nod[i].y;
        double xj = nod[j].x, yj = nod[j].y;
        double xk = nod[k].x, yk = nod[k].y;
        
        double Ae = 0.5 * fabs((xj*yk - xk*yj) - xi*(yk - yj) + yi*(xk - xj));
        
        FF[i] += Ae * (2.0 * F[i] + F[j] + F[k]) / 12.0;
        FF[j] += Ae * (F[i] + 2.0 * F[j] + F[k]) / 12.0;
        FF[k] += Ae * (F[i] + F[j] + 2.0 * F[k]) / 12.0;
    }
}

void ludcmp(double **a, int n, int *indx, double *d) {
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv = new_double(n);
    *d = 1.0;
    
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            if ((temp = fabs(a[i][j])) > big) big = temp;
        }
        if (big == 0.0) {
            fprintf(stderr, "Singular matrix in routine ludcmp\n");
            free(vv);
            exit(EXIT_FAILURE);
        }
        vv[i] = 1.0 / big;
    }
    
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            sum = a[i][j];
            for (k = 0; k < i; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }
        
        big = 0.0;
        imax = j;
        
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
            
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        
        if (j != imax) {
            for (k = 0; k < n; k++) {
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
        
        if (j != n - 1) {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i < n; i++) {
                a[i][j] *= dum;
            }
        }
    }
    
    free(vv);
}

void lubksb(double **a, int n, int *indx, double b[]) {
    int i, ii = 0, ip, j;
    double sum;

    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii - 1; j < i; j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i + 1;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

int main(int argc, char *argv[]) {
    FILE *fd;
    
    // 1. Leer malla
    fd = fopen("example.e", "r");
    if (!fd) {
        fprintf(stderr, "Error: Cannot open example.e\n");
        return 1;
    }
    fscanf(fd, "%d\n", &Ne);
    ele = new_element(Ne);
    for (int k = 0; k < Ne; k++) {
        fscanf(fd, "%d %d %d %d\n", &ele[k].Id, &ele[k].i, &ele[k].j, &ele[k].k);
    }
    fclose(fd);
    
    fd = fopen("example.n", "r");
    if (!fd) {
        fprintf(stderr, "Error: Cannot open example.n\n");
        return 1;
    }
    fscanf(fd, "%d\n", &Nn);
    nod = new_node(Nn);
    for (int k = 0; k < Nn; k++) {
        fscanf(fd, "%d %lf %lf %d\n", &nod[k].Id, &nod[k].x, &nod[k].y, &nod[k].mark);
    }
    fclose(fd);
    
    // 2. Inicializar matrices
    K = new_matrix_double(Nn, Nn);
    F = new_double(Nn);
    FF = new_double(Nn);
    
    for (int k = 0; k < Nn; k++) {
        F[k] = -1.0; 
        FF[k] = 0.0;
        for (int i = 0; i < Nn; i++) {
            K[k][i] = 0.0;
        }
    }
    
    // 3. Ensamblar sistema
    getK();
    getF();
    
    // 4. Aplicar condiciones de frontera
    applyBoundaryConditions();
    
    // 5. Resolver sistema
    int *indx = new_int(Nn);
    double det;
    double **K_copy = new_matrix_double(Nn, Nn);
    for (int i = 0; i < Nn; i++) {
        for (int j = 0; j < Nn; j++) {
            K_copy[i][j] = K[i][j];
        }
    }
    
    double *solution = new_double(Nn);
    for (int i = 0; i < Nn; i++) {
        solution[i] = FF[i];
    }
    
    ludcmp(K_copy, Nn, indx, &det);
    lubksb(K_copy, Nn, indx, solution);
    
    // 6. Imprimir resultados
    FILE *fd_out = fopen("node_f.dat", "w");
    if (!fd_out) {
        fprintf(stderr, "Error: Cannot open node_f.dat\n");
        return 1;
    }
    fprintf(fd_out, "%d\n", Nn);
    for (int i = 0; i < Nn; i++) {
        fprintf(fd_out, "%d %.6f %.6f %d %.6f\n", nod[i].Id, nod[i].x, nod[i].y, nod[i].mark, solution[i]);
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

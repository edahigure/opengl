#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Parámetros
double L = 1.0;          // Lado del dominio cuadrado
int N = 50;              // Número de términos de Fourier
int NX = 5;              // Puntos en x
int NY = 5;              // Puntos en y
double PI = 3.14159265358979323846;

// Coeficiente A_mn
double compute_Amn(int m, int n ) {
    double term_m = (m % 2 == 0) ? 0.0 : 2.0;
    double term_n = (n % 2 == 0) ? 0.0 : 2.0;
    return -4.0 * term_m * term_n  / (m * n * PI * PI * PI * PI * (m * m + n * n));
}

// Solución u(x, y)
double compute_u(double x, double y, double L, int N) {
    double u = 0.0;
    for (int m = 1; m <= N; m++) {
        for (int n = 1; n <= N; n++) {
            double Amn = compute_Amn(m, n);
            u += Amn * sin(m * PI * x / L) * sin(n * PI * y / L);
        }
    }
    return u;
}

int main() {
    // Espaciado de la malla
    double dx = L / (NX - 1);
    double dy = L / (NY - 1);
    
    // Archivo de salida
    FILE *fp = fopen("poisson_analytic.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: No se pudo abrir poisson_analytic.dat\n");
        return 1;
    }
    
    // Escribir encabezado
    fprintf(fp, "%d \n", NX*NY);
    
    // Calcular solución en la malla
    for (int i = 0; i < NX; i++) {
        double x = i * dx;
        for (int j = 0; j < NY; j++) {
            double y = j * dy;
            double u = compute_u(x, y, L, N);
            fprintf(fp, "%f %f %f\n", x, y, u);
        }
    }
    
    fclose(fp);
    printf("Solución escrita en poisson_analytic.dat\n");
    printf("Formato: NX NY\nx y u\n");
    
    return 0;
}

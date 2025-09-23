#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

int main() {
    FILE *file = fopen("circulos.d", "w");
    if (file == NULL) {
        printf("Error al crear archivo\n");
        return 1;
    }
    
    int segments = 17; // Más segmentos para mejor calidad
    double center_x = 0.0, center_y = 0.0;
    double r_outer = 2.0, r_inner = 1.0;
    double spacing_outer = 0.3;
    double spacing_inner = 0.3;
    
    // PUNTOS - los índices EMPIEZAN en 0
    fprintf(file, "%d\n", 2 * segments);
    
    // Círculo exterior - SENTIDO ANTIHORARIO (positivo)
    for (int i = 0; i < segments; i++) {
        double angle = 2 * PI * i / segments;
        double x = center_x + r_outer * cos(angle);
        double y = center_y + r_outer * sin(angle);
        fprintf(file, "%d %.6f %.6f %.2f %d\n", 
                i, x, y, spacing_outer, 1);
    }
    
    // Círculo interior - SENTIDO HORARIO (negativo) PARA HUECO
    for (int i = 0; i < segments; i++) {
        double angle = -2 * PI * i / segments; // NEGATIVO para sentido horario
        double x = center_x + r_inner * cos(angle);
        double y = center_y + r_inner * sin(angle);
        fprintf(file, "%d %.6f %.6f %.2f %d\n", 
                segments + i, x, y, spacing_inner, 2);
    }
    
    // EDGES
    fprintf(file, "%d\n", 2 * segments);
    
    // Edges círculo exterior - ANTIHORARIO
    for (int i = 0; i < segments; i++) {
        int next = (i + 1) % segments;
        fprintf(file, "%d %d %d %d\n", 
                i, i, next, 1);
    }
    
    // Edges círculo interior - HORARIO
    for (int i = 0; i < segments; i++) {
        int current = segments + i;
        int next = segments + ((i + 1) % segments);
        fprintf(file, "%d %d %d %d\n", 
                segments + i, current, next, 2);
    }
    
    // ELEMENTOS, HUECOS y REGIONES
    fprintf(file, "0\n"); // elementos
    fprintf(file, "1\n"); // 1 hueco
    fprintf(file, "1 %.1f %.1f\n", center_x, center_y); // hueco en centro
    fprintf(file, "0\n"); // regiones
    
    fclose(file);
    printf("Archivo circulos_hueco.mesh generado\n");
    printf("Círculo exterior: antihorario\n");
    printf("Círculo interior: horario (para hueco)\n");
    
    return 0;
}

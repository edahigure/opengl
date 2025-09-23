#include "mesh_reader.h"
#include <stdio.h>
#include <string.h> 

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s <archivo_base>\n", argv[0]);
        printf("Ejemplo: %s example\n", argv[0]);
        return 1;
    }
    
    MeshData mesh;
    memset(&mesh, 0, sizeof(MeshData));
    
    if (load_mesh_data(argv[1], &mesh)) {
        print_mesh_summary(&mesh);
        printf("\n✅ Archivos leídos correctamente\n");
    } else {
        printf("\n❌ Error leyendo archivos\n");
        return 1;
    }
    
    return 0;
}

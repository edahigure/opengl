#include <QApplication>
#include <QGLWidget>
#include <QGLFunctions>
#include <QPushButton>
#include <QVBoxLayout>
#include <QDebug>
#include <QMouseEvent>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/gl.h>

#define MAX_NODES 10000
#define MAX_ELEMENTS 10000
#define MAX_SIDES 10000
#define GREAT 1e10
#define SMALL 1e-10
#define OFF -1

// ==================== ESTRUCTURAS ====================
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

typedef struct {
    float r, g, b;
} RGB;

typedef struct {
    int v1, v2, v3;
} Triangle;

// ==================== VARIABLES GLOBALES ====================
Node* nodes = NULL;
Element* elements = NULL;
Side* sides = NULL;
int Nn = 0, Ne = 0, Ns = 0;
double f_min = 0.0, f_max = 0.0;

// ==================== FUNCIONES GLOBALES ====================
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
    nodes = (Node*)malloc(Nn * sizeof(Node));
    
    // Leer nodos
    for (int i = 0; i < Nn; i++) {
        fgets(line, sizeof(line), file);
        remove_comments(line);
        sscanf(line, "%d: %lf %lf %d", 
               &nodes[i].index, &nodes[i].x, &nodes[i].y, &nodes[i].mark);
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
    elements = (Element*)malloc(Ne * sizeof(Element));
    
    for (int i = 0; i < Ne; i++) {
        fgets(line, sizeof(line), file);
        remove_comments(line);
        sscanf(line, "%d: %d %d %d %d %d %d %d %d %d %lf %lf %d",
               &elements[i].index, &elements[i].i, &elements[i].j, &elements[i].k,
               &elements[i].ei, &elements[i].ej, &elements[i].ek,
               &elements[i].si, &elements[i].sj, &elements[i].sk,
               &elements[i].xv, &elements[i].yv, &elements[i].material);
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
               &sides[i].index, &sides[i].c, &sides[i].d,
               &sides[i].ea, &sides[i].eb, &sides[i].mark);
    }
    fclose(file);
    
    return 1;
}

void free_mesh_data() {
    if (nodes) free(nodes);
    if (elements) free(elements);
    if (sides) free(sides);
    nodes = NULL;
    elements = NULL;
    sides = NULL;
    Nn = Ne = Ns = 0;
}

char* newstr(int N) {
    return (char*)malloc(N * sizeof(char));
}

double* new_double(int N) {
    return (double*)malloc(N * sizeof(double));
}

RGB map_to_rgb(double value, double a, double b) {
    RGB color;
    color.r = 1.0f;
    color.g = 1.0f;
    color.b = 1.0f;
    
    if (a == b) return color;
    
    double normalized = (value - a) / (b - a);
    if (normalized < 0.0) normalized = 0.0;
    if (normalized > 1.0) normalized = 1.0;
    
    double mid = 0.5;
    if (normalized < mid) {
        double factor = normalized / mid;
        color.r = 0.0f;
        color.g = (float)factor;
        color.b = 1.0f;
    } else {
        double factor = (normalized - mid) / (1.0 - mid);
        color.r = (float)factor;
        color.g = 1.0f;
        color.b = 1.0f - (float)factor;
    }
    
    return color;
}

void convert2Engf(float z, float *z1, int *n);

// ==================== CLASE GLWidget ====================
class GLWidget : public QGLWidget, protected QGLFunctions {
public:
    GLWidget(QWidget *parent = NULL, char *file_1 = NULL) : QGLWidget(parent) {
        xRot = 0; yRot = 0; zRot = 0;
        xTrans = 0; yTrans = 0; zTrans = -10.0;
        scale = 1.0;
        x1 = -0.2; x2 = 1.1; y1 = -0.2; y2 = 1.1;

        draw_mesh = true;
        draw_voronoi = true;
        draw_marks = true; 
        draw_fill = true;

        file_name_1 = newstr(200);
        snprintf(file_name_1, 200, "%s.n", file_1);
        file_name_2 = newstr(200);  
        snprintf(file_name_2, 200, "%s.e", file_1); 

        // Cargar datos usando variables globales
        load_mesh_data(file_1);
    }

    ~GLWidget() {
        // SOLO liberar lo que se allocó en este objeto
        if (file_name_1) free(file_name_1);
        if (file_name_2) free(file_name_2);
        // NO liberar nodes, elements, sides (son globales)
    }

protected:
    bool draw_mesh, draw_voronoi, draw_marks, draw_fill;
    int drawList;	
    QPoint oldPos;
    float x1, x2, y1, y2;
    float xTrans, yTrans, zTrans;
    float scale;
    float xRot, yRot, zRot;
    char *file_name_1;
    char *file_name_2;

    void initializeGL() {
        initializeGLFunctions();
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    }

    void resizeGL(int w, int h) {
        glViewport(0, 0, (GLint)w, (GLint)h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(x1, x2, y1, y2, 5.0, 15.0);
        glMatrixMode(GL_MODELVIEW);
    }

    void transform() {
        glTranslatef(xTrans, yTrans, zTrans);
        glScalef(scale, scale, scale);
        glRotatef(xRot, 1.0, 0.0, 0.0);
        glRotatef(yRot, 0.0, 1.0, 0.0);
        glRotatef(zRot, 0.0, 0.0, 1.0);
    }

    void writeNumber(float x, float y, float num) {
        QFont FontSmall("Helvetica", 10);
        FontSmall.setPointSize(7);
        QString snum;
        snum.setNum(num, 'e', 1);
        renderText(x, y, 0, snum, FontSmall);
    }

    void axis() {
        float xt, yt, Lx, Ly;
        
        glColor3f(0, 0, 1);
        glBegin(GL_LINES);
        glVertex2f(x1 - xTrans / scale, 0);
        glVertex2f(x2 - xTrans / scale, 0);
        glEnd();
        glBegin(GL_LINES);
        glVertex2f(0, y1 - yTrans / scale);
        glVertex2f(0, y2 - yTrans / scale);
        glEnd();
        
        int nr;
        float x1r, y1r;
        
        convert2Engf(x1, &x1r, &nr);
        xt = floor(x1r * 10) / 10 * pow(10, nr);
        convert2Engf(y1, &y1r, &nr);
        yt = floor(y1r * 10) / 10 * pow(10, nr);
        
        Lx = fabs(x2 - x1);
        Ly = fabs(y2 - y1);

        do {
            glBegin(GL_LINES);
            glVertex2f(xt, -Ly / 100 / scale);
            glVertex2f(xt, Ly / 100 / scale);
            glEnd();
            if (fabs(xt) > 0.000001)
                writeNumber(xt - 4 * Lx / 100 / scale, -Ly / 25 / scale, xt);
            xt = xt + Lx / 10;
        } while (xt < x2 - xTrans / scale);

        do {
            glBegin(GL_LINES);
            glVertex2f(-Lx / 100 / scale, yt);
            glVertex2f(Lx / 100 / scale, yt);
            glEnd();
            if (fabs(yt) > 0.000001)
                writeNumber(-Lx / 10 / scale, yt - Ly / 400 / scale, yt);
            yt = yt + Ly / 10;
        } while (yt < y2 - yTrans / scale);
    }

    void circle(double xc, double yc, double R) {
        double t;
        int N = 100;
        double PI = 4 * atan(1);
        glBegin(GL_LINE_STRIP);
        t = 0;
        do {
            glVertex2f(xc + R * cos(t) / scale, yc + R * sin(t) / scale);
            t += (2 * PI) / N;
        } while (t <= 2.2 * PI);
        glEnd();
    }

    void paintGL() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        transform();
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1, 0.5);

        paintMesh();

        axis();
    }

    void mousePressEvent(QMouseEvent *e) {
        e->accept();
        oldPos = e->pos();
    }

    void mouseReleaseEvent(QMouseEvent *e) {
        e->accept();
        oldPos = e->pos();
    }

    void mouseMoveEvent(QMouseEvent *e) {
        e->accept();
        double dx = e->x() - oldPos.x();
        double dy = e->y() - oldPos.y();
        oldPos = e->pos();

        double rx = dx / width();
        double ry = dy / height();

        if (e->buttons() == Qt::LeftButton) {
            xTrans += fabs(x1 - x2) * rx / scale;
            yTrans -= fabs(y1 - y2) * ry / scale;
            updateGL();
        }
    }

    void wheelEvent(QWheelEvent *e) {
        e->accept();
        float zoomFactor = 1.15f;
        float oldScale = scale;

        if (e->delta() > 0) {
            scale *= zoomFactor;
        } else {
            scale /= zoomFactor;
        }

        if (scale < 0.1f) scale = 0.1f;
        if (scale > 10.0f) scale = 10.0f;

        QPointF pos = e->pos();
        double mouseX = x1 + (x2 - x1) * pos.x() / width();
        double mouseY = y2 - (y2 - y1) * pos.y() / height();
        xTrans += (mouseX - xTrans) * (1.0f - oldScale / scale);
        yTrans += (mouseY - yTrans) * (1.0f - oldScale / scale);

        updateGL();
    }

    void paintMesh()
	 {
		 int n, s, ea, eb;
		 double scl, x, y, xc, yc, xd, yd, x1, y1, x2, y2,
		        xmax = -GREAT, xmin = +GREAT, ymax = -GREAT, ymin = +GREAT;
		 
		 float red[3] = {1,0,0}, green[3] = {0,1,0}, blue[3] = {0,0,1},
		       redgreen[3] = {.5,.5,0}, redblue[3] = {.5,0,.5}, 
		       greenblue[3] = {0,.5,.5},
		       white[3] = {1,1,1};
		 float *color[7] = {red, green, blue, redgreen, redblue, greenblue, white};                
		 
		 qDebug() << "paintMesh";
		 
     
		     glClearColor(0.2f, 0.2f, 0.2f, 1.0f); 
		 // Calcular límites
		 for(n = 0; n < Nn; n++) {
		     if(nodes[n].mark != OFF) {
		         xmin = min(xmin, nodes[n].x);
		         ymin = min(ymin, nodes[n].y);
		         xmax = max(xmax, nodes[n].x);
		         ymax = max(ymax, nodes[n].y);
		     }
		 }
		 scl = min(380.0/(ymax-ymin+SMALL), 700.0/(xmax-xmin+SMALL));
		 
		 // Dibujar malla (si está habilitado)
		 if(draw_mesh) {
		     glColor3fv(blue);
		     
		     for(s = 0; s < Ns; s++) {
		         if(sides[s].mark != OFF) {
		             xc = nodes[sides[s].c].x;
		             yc = nodes[sides[s].c].y;
		             xd = nodes[sides[s].d].x;
		             yd = nodes[sides[s].d].y;
		             
		             glBegin(GL_LINES);
		             glVertex2f(xc, yc);
		             glVertex2f(xd, yd);
		             glEnd();
		         }
		     }
		 }
		 
		 // Dibujar lados con marcas (boundary)
		 for(s = 0; s < Ns; s++) {
		     if(sides[s].mark > 0) { // Lado en el boundary
		         if(draw_marks) {
		             glColor3fv(color[7 % sides[s].mark]);
		         } else {
		             glColor3fv(color[5]);
		         }
		         
		         xc = nodes[sides[s].c].x;
		         yc = nodes[sides[s].c].y;
		         xd = nodes[sides[s].d].x;
		         yd = nodes[sides[s].d].y;
		         
		         glBegin(GL_LINES);
		         glVertex2f(xc, yc);
		         glVertex2f(xd, yd);
		         glEnd();
		     }
		 }
		 
		 // Dibujar nodos con marcas (boundary)
		 if(draw_marks) {
		     for(n = 0; n < Nn; n++) {
		         if(nodes[n].mark > 0) { // Nodo en el boundary
		             glColor3fv(color[7 % (7 - nodes[n].mark)]);
		             x = nodes[n].x;
		             y = nodes[n].y;
		             circle(x, y, 0.02);
		         }
		     }
		 }
		 
		 // Dibujar diagrama de Voronoi (si está habilitado)
		 if(draw_voronoi) {
		     glColor3fv(color[4]);
		     
		     for(s = 0; s < Ns; s++) {
		         if(sides[s].mark != OFF) {
		             if((ea = sides[s].ea) != OFF) {
		                 x1 = elements[ea].xv;
		                 y1 = elements[ea].yv;
		             } else {
		                 x1 = 0.5 * (nodes[sides[s].c].x + nodes[sides[s].d].x);
		                 y1 = 0.5 * (nodes[sides[s].c].y + nodes[sides[s].d].y);
		             }
		             
		             if((eb = sides[s].eb) != OFF) {
		                 x2 = elements[eb].xv;
		                 y2 = elements[eb].yv;
		             } else {
		                 x2 = 0.5 * (nodes[sides[s].c].x + nodes[sides[s].d].x);
		                 y2 = 0.5 * (nodes[sides[s].c].y + nodes[sides[s].d].y);
		             }
		             
		             glBegin(GL_LINES);
		             glVertex2f(x1, y1);
		             glVertex2f(x2, y2);
		             glEnd();
		         }
		     }
		 }
		 
	}

};

// ==================== CLASE MainWindow ====================
class MainWindow : public QWidget {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = NULL, char *file_1 = NULL) : QWidget(parent) {
        file_name_1 = newstr(200);
        strcpy(file_name_1, file_1); 

        QVBoxLayout *layout = new QVBoxLayout(this);
        GLWidget *glWidget = new GLWidget(this, file_name_1);
        glWidget->setMinimumSize(500, 500);
        layout->addWidget(glWidget);
    }

    ~MainWindow() {
        if (file_name_1) free(file_name_1);
        free_mesh_data(); // Liberar datos globales al cerrar
    }

private:
    char *file_name_1;
};

// ==================== MAIN ====================
int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    if (argc < 2) {
        qDebug() << "Uso: ./programa <archivo_base>";
        return 1;
    }
    
    MainWindow mainWindow(NULL, argv[1]);
    mainWindow.setWindowTitle("Visualizador de Malla");
    mainWindow.resize(500, 500);
    mainWindow.show();
    
    return app.exec();
}

// ==================== FUNCIONES AUXILIARES ====================
void convert2Engf(float z, float *z1, int *n) {
    double zt = z;
    int nexp = 0;
    
    if(fabs(z) > 10) {
        do {
            zt = zt / 10;
            nexp++;
        } while(fabs(zt) > 10);
    }
    
    if(fabs(z) < 1) {
        do {
            zt = zt * 10;
            nexp--;
        } while(fabs(zt) < 1);
    }
    
    *z1 = zt;
    *n = nexp;
}

#include "main.moc"

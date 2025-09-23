#include <QApplication>
#include <QGLWidget>
#include <QGLFunctions>
#include <QPushButton>
#include <QVBoxLayout>
#include <QDebug>
#include <QMouseEvent>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>

void convert2Engf(float z,float *z1,int *n);

char * newstr(int N)
{
  return (char*) malloc((N)*sizeof(char) );
}


double function_f(double x, double y) {
    return cos(1/x * 1/y);
}

double* new_double(int N) {
    return (double*)malloc((size_t)(N * sizeof(double)));
}



typedef struct {
    float r, g, b;
} RGB;

typedef struct {
    int index;
    float x, y;
    int mark;
    double solution; // 0 para interior, 1 para frontera
} Node;

typedef struct {
    int Id;
    int v1, v2, v3; // Índices de los nodos (1-based, ajustado a 0-based)
    int ei, ej, ek;
    int si, sj, sk;
    double xv, yv;
    int material;

} Triangle;


// Función map_to_rgb (azul → cian claro → amarillo)
RGB map_to_rgb(double value, double a, double b) {
    RGB color;
    color.r = 1.0f;
    color.g = 1.0f;
    color.b = 1.0f;
    
    if (a == b) {
        return color; // Blanco si a == b
    }
    
    double normalized = (value - a) / (b - a);
    if (normalized < 0.0) normalized = 0.0;
    if (normalized > 1.0) normalized = 1.0;
    
    double mid = 0.5;
    if (normalized < mid) {
        // Azul (0, 0, 1) a cian claro (0, 1, 1)
        double factor = normalized / mid;
        color.r = 0.0f;
        color.g = (float)factor;
        color.b = 1.0f;
    } else {
        // Cian claro (0, 1, 1) a amarillo (1, 1, 0)
        double factor = (normalized - mid) / (1.0 - mid);
        color.r = (float)factor;
        color.g = 1.0f;
        color.b = 1.0f - (float)factor;
    }
    
    return color;
}

// Leer el archivo de nodos
int read_nodes(const char* filename, Node** nodes, int* n) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        qDebug() << "Error: No se pudo abrir el archivo" << filename;
        return 0;
    }
    
    fscanf(file, "%d", n);
    *nodes = (Node*)malloc(*n * sizeof(Node));
    if (!*nodes) {
        qDebug() << "Error: No se pudo asignar memoria para los nodos";
        fclose(file);
        return 0;
    }
    
    for (int i = 0; i < *n; i++) {
        int index, mark;
        float x, y;
		  double solution;	
        if (fscanf(file, "%d %f %f %d %lf", &index, &x, &y, &mark, &solution) != 5) {
            qDebug() << "Error: Formato de archivo incorrecto en la línea" << i + 2;
            free(*nodes);
            fclose(file);
            return 0;
        }
        (*nodes)[i].index = index;
        (*nodes)[i].x = x;
        (*nodes)[i].y = y;
        (*nodes)[i].mark = mark;
		  (*nodes)[i].solution = solution;
    }
    
    fclose(file);
    return 1;
}

// Leer el archivo de triangulación
// Leer el archivo de triangulación
int read_triangles(const char* filename, Triangle** triangles, int* num_triangles) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        qDebug() << "Error: No se pudo abrir el archivo" << filename;
        return 0;
    }
    
    // Leer número de triángulos
    if (fscanf(file, "%d", num_triangles) != 1) {
        qDebug() << "Error: No se pudo leer el número de triángulos";
        fclose(file);
        return 0;
    }
    
    qDebug() << "Número de triángulos:" << *num_triangles;
    
    *triangles = (Triangle*)malloc(*num_triangles * sizeof(Triangle));
    if (!*triangles) {
        qDebug() << "Error: No se pudo asignar memoria para los triángulos";
        fclose(file);
        return 0;
    }
    
    for (int i = 0; i < *num_triangles; i++) {
        int id, v1, v2, v3;
        int ei, ej, ek, si, sj, sk, material;
        double xv, yv;
        
        // Leer todos los 13 campos del archivo
        int elementos_leidos = fscanf(file, "%d %d %d %d %d %d %d %d %d %d %lf %lf %d", 
               &id, &v1, &v2, &v3, &ei, &ej, &ek, &si, &sj, &sk, &xv, &yv, &material);
        
        if (elementos_leidos != 13) {
            qDebug() << "Error: Formato de archivo incorrecto en la línea" << i + 2;
            qDebug() << "Se esperaban 13 campos, se leyeron" << elementos_leidos;
            qDebug() << "Leyendo: id, v1, v2, v3, ei, ej, ek, si, sj, sk, xv, yv, material";
            free(*triangles);
            fclose(file);
            return 0;
        }
        
        // Asignar los vértices del triángulo
        (*triangles)[i].v1 = v1;
        (*triangles)[i].v2 = v2;
        (*triangles)[i].v3 = v3;
        
        // Debug: mostrar los primeros 5 triángulos
        if (i < 5) {
            qDebug() << "Triángulo" << i << ": v1 =" << v1 << ", v2 =" << v2 << ", v3 =" << v3;
        }
    }
    
    fclose(file);
    qDebug() << "Triángulos leídos exitosamente:" << *num_triangles;
    return 1;
}

class GLWidget : public QGLWidget, protected QGLFunctions {
public:
    GLWidget(QWidget *parent = NULL, char *file_1 = NULL, char *file_2 = NULL) : QGLWidget(parent) {
        xRot = 0;
        yRot = 0;
        zRot = 0;
        xTrans = 0;
        yTrans = 0;
        zTrans = -10.0;
        scale = 1.0;
        x1 = -0.2;
        x2 = 1.1;
        y1 = -0.2;
        y2 = 1.1;

        file_name_1 = newstr(200);
        strcpy(file_name_1, file_1 ? file_1 : "");


        file_name_2 = newstr(200);
        strcpy(file_name_2, file_2 ? file_2 : "");

        nodes = NULL;
        num_nodes = 0;
        triangles = NULL;
        num_triangles = 0;
        f_val = NULL;

        // Leer datos y calcular valores de la función
        read_data(file_name_1,file_name_2);
        fprintf(stderr,"Here ok %d\n", num_nodes);
        if (num_nodes > 0) {
            f_val = new_double(num_nodes);
            f_min = 10.0e10;
            f_max = -10.0e10;
            for (int i = 0; i < num_nodes; i++) {
                f_val[i] = nodes[i].solution;
                if (f_val[i] < f_min) f_min = f_val[i];
                if (f_val[i] > f_max) f_max = f_val[i];
            }
        }
    }

    ~GLWidget() {
        if (nodes) free(nodes);
        if (triangles) free(triangles);
        if (f_val) free(f_val);
        if (file_name_1) free(file_name_1);
        if (file_name_2) free(file_name_2);
    }

    void read_data(const char* nodes_filename, const char* triangles_filename) {
        if (!read_nodes(nodes_filename, &nodes, &num_nodes)) {
            nodes = NULL;
            num_nodes = 0;
            return;
        }
        if (!read_triangles(triangles_filename, &triangles, &num_triangles)) {
            triangles = NULL;
            num_triangles = 0;
            free(nodes);
            nodes = NULL;
            num_nodes = 0;
            return;
        }
    }

protected:
    QPoint oldPos;
    float x1, x2, y1, y2;
    Node* nodes;
    int num_nodes;
    Triangle* triangles;
    int num_triangles;
    double* f_val;
    double f_min, f_max;
    float xTrans, yTrans, zTrans;
    float scale;
    float xRot, yRot, zRot;
    char *file_name_1;
    char *file_name_2;

    void initializeGL() {
        initializeGLFunctions();
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Fondo blanco
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
        glVertex2f(x1 - xTrans / scale, 0); // Ajustar ejes al zoom
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
            glVertex2f(xc + R * cos(t) / scale, yc + R * sin(t) / scale); // Ajustar radio al zoom
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

        // Dibujar triángulos coloreados
        if (nodes && triangles && f_val) {
            glBegin(GL_TRIANGLES);
            for (int i = 0; i < num_triangles; i++) {
                int v1 = triangles[i].v1 ;
                int v2 = triangles[i].v2 ;
                int v3 = triangles[i].v3 ;

                RGB color = map_to_rgb(f_val[v1], f_min, f_max);
                glColor3f(color.r, color.g, color.b);
                glVertex2f(nodes[v1].x, nodes[v1].y);

                color = map_to_rgb(f_val[v2], f_min, f_max);
                glColor3f(color.r, color.g, color.b);
                glVertex2f(nodes[v2].x, nodes[v2].y);

                color = map_to_rgb(f_val[v3], f_min, f_max);
                glColor3f(color.r, color.g, color.b);
                glVertex2f(nodes[v3].x, nodes[v3].y);
            }
            glEnd();
				if(0){
            // Dibujar la malla con líneas
		         glDisable(GL_POLYGON_OFFSET_FILL);
		         glColor3f(0.0f, 0.0f, 0.0f);
		         glBegin(GL_LINES);
		         for (int i = 0; i < num_triangles; i++) {
		             int v1 = triangles[i].v1 ;
		             int v2 = triangles[i].v2 ;
		             int v3 = triangles[i].v3 ;

		             glVertex2f(nodes[v1].x, nodes[v1].y);
		             glVertex2f(nodes[v2].x, nodes[v2].y);
		             glVertex2f(nodes[v2].x, nodes[v2].y);
		             glVertex2f(nodes[v3].x, nodes[v3].y);
		             glVertex2f(nodes[v3].x, nodes[v3].y);
		             glVertex2f(nodes[v1].x, nodes[v1].y);
		         }
		         glEnd();

		         // Dibujar nodos como círculos
		         glColor3f(0.0f, 0.0f, 0.0f);
		         for (int i = 0; i < num_nodes; i++) {
		             circle(nodes[i].x, nodes[i].y, 0.03);
		         }

            }
        }

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
            xTrans += fabs(x1 - x2) * rx / scale; // Ajustar traslación al zoom
            yTrans -= fabs(y1 - y2) * ry / scale;
            updateGL();
        }
    }

    void wheelEvent(QWheelEvent *e) {
        e->accept();
        float zoomFactor = 1.15f; // Factor de zoom (1.15 = 15% por paso)
        float oldScale = scale;

        // Obtener el ángulo de desplazamiento (positivo para acercar, negativo para alejar)
        if (e->delta() > 0) {
            scale *= zoomFactor; // Acercar
        } else {
            scale /= zoomFactor; // Alejar
        }

        // Limitar el rango de scale
        if (scale < 0.1f) scale = 0.1f; // Zoom mínimo
        if (scale > 10.0f) scale = 10.0f; // Zoom máximo

        // Ajustar traslación para mantener el punto bajo el ratón fijo
        QPointF pos = e->pos();
        double mouseX = x1 + (x2 - x1) * pos.x() / width();
        double mouseY = y2 - (y2 - y1) * pos.y() / height(); // y invertido en OpenGL
        xTrans += (mouseX - xTrans) * (1.0f - oldScale / scale);
        yTrans += (mouseY - yTrans) * (1.0f - oldScale / scale);

        updateGL();
    }
};

class MainWindow : public QWidget {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = NULL,char * file_1 = NULL,char * file_2 = NULL) : QWidget(parent) {
        // Layout y botón
		  file_name_1 = newstr(200);
        strcpy(file_name_1,file_1); 
		  file_name_2 = newstr(200);
        strcpy(file_name_2,file_2); 

        QVBoxLayout *layout = new QVBoxLayout(this);
        GLWidget *glWidget = new GLWidget(this, file_name_1, file_name_2);
        glWidget->setMinimumSize(500, 500);

        layout->addWidget(glWidget);

    }

    ~MainWindow() {
        if (file_name_1) free(file_name_1);
        if (file_name_2) free(file_name_2);
    }

public slots:


private:
    char *file_name_1;
    char *file_name_2;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    MainWindow mainWindow(NULL,argv[1],argv[2]);
    mainWindow.setWindowTitle("OpenGL + Botón en Qt4");
    mainWindow.resize(500, 500);
    mainWindow.show();
    return app.exec();
}


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


#include "main.moc"

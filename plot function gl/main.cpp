#include <QApplication>
#include <QGLWidget>
#include <QGLFunctions>
#include <QPushButton>
#include <QVBoxLayout>
#include <QDebug>
#include <QMouseEvent>
#include "util.h"



class GLWidget : public QGLWidget, protected QGLFunctions
{
public:
    GLWidget(QWidget *parent = NULL, char *file = NULL) : QGLWidget(parent) {
	 xRot=0;
    yRot=0;
    zRot=0;
    xTrans=0;
    yTrans=0;
    zTrans=-10.0;
    scale=1.0;	
    
	 x1=-5;
    x2=5;
    y1=-5;
    y2=5;

	file_name = newstr(200);
   strcpy(file_name,file);
	}

	void mousePressEvent( QMouseEvent *e )
	{
 	   e->accept();
	   oldPos = e->pos();
	}

	void mouseReleaseEvent( QMouseEvent *e )
	{
		 e->accept();
		 oldPos = e->pos();
	}

	void mouseMoveEvent( QMouseEvent *e )
	{
		 e->accept();
		 double dx = e->x() - oldPos.x();
		 double dy = e->y() - oldPos.y();

		 oldPos = e->pos();

		 double rx = dx / width();
		 double ry = dy / height();

		 if ( e->buttons() == Qt::LeftButton )
		 {
		  	   
		  	   xTrans=xTrans+fabs(x1-x2)*rx;
		  	   yTrans=yTrans-fabs(y1-y2)*ry; 
		        transform();
		        updateGL();
			}

	}
    


protected:

   QPoint oldPos;

   float x1,x2,y1,y2;
	void initializeGL()
	{
	  initializeGLFunctions();
	  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Color de fondo (azul grisáceo)
	}
    

	void resizeGL( int w, int h )
	{
		 glViewport( 0, 0, (GLint)w, (GLint)h );
		 glMatrixMode( GL_PROJECTION );
		 glLoadIdentity();
		 glOrtho( x1, x2, y1, y2, 5.0, 15.0 );
		 glMatrixMode( GL_MODELVIEW );

	}

	void transform()
	{
		 glTranslatef( xTrans, yTrans, zTrans );
		 glScalef( scale, scale, scale );
		 glRotatef( xRot, 1.0, 0.0, 0.0 );
		 glRotatef( yRot, 0.0, 1.0, 0.0 );
		 glRotatef( zRot, 0.0, 0.0, 1.0 );
	}

	void writeNumber( float x,float y,float num)
	{
		QFont FontSmall("Helvetica", 10);
		FontSmall.setPointSize(7);
		QString snum;		    
		snum.setNum(num,'e',1);   
		renderText(x,y,0,snum,FontSmall,2000); 
	}

	void axis(void)
	{
		float xt,yt,Lx,Ly;
		
		glColor3f(0,0,1);
		glBegin(GL_LINES);
		   glVertex2f(x1-xTrans,0);
		   glVertex2f(x2-xTrans,0);
		glEnd();
		glBegin(GL_LINES);
		   glVertex2f(0,y1-yTrans);
		   glVertex2f(0,y2-yTrans);
		glEnd();
	  
		int nr;
		float x1r,y1r;
		
		convert2Engf(x1,&x1r,&nr);   
		xt=floor(x1r*10)/10*pow(10,nr);
		convert2Engf(y1,&y1r,&nr);
		yt=floor(y1r*10)/10*pow(10,nr);
		
		Lx=fabs(x2-x1);
		Ly=fabs(y2-y1);

		do
		{
		   glBegin(GL_LINES);
		      glVertex2f(xt,-Ly/100);
		      glVertex2f(xt,Ly/100);
		   glEnd();
		   if(fabs(xt)>.000001 )
		      writeNumber(xt-4*Lx/100,-Ly/25,xt);
		   xt=xt+Lx/10;
		}while(xt<x2-xTrans);
		do
		{
		   glBegin(GL_LINES);
		      glVertex2f(-Lx/100,yt);
		      glVertex2f(Lx/100,yt);
		   glEnd();
		   if(fabs(yt)>.000001 )
		      writeNumber(-Lx/10,yt-Ly/400,yt);
		   yt=yt+Ly/10;
		}
		while(yt<y2-yTrans);
		
			
	}

	void circle(double xc,double yc,double R)
	{
		 double t;
		 int N=100;
		 double PI;
       PI= 4*atan(1);
		 glBegin(GL_LINE_STRIP);
		 t=0;
		 do
		 {

		      glVertex2f(xc+R*cos(t),yc+R*sin(t));
		      t=t+(2*PI)/N;
		    }while(t<=2.2*PI);
		 glEnd();

	}



	void paintGL()
	{

		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glLoadIdentity();
		transform();
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1, .5);


      glColor3f(0.0f, 0.0f, 1.0f);   
		axis();
    

      double *x,*y;
	   int Npoints;

		FILE *file;
		
		
		if((file=fopen(file_name, "r"))==NULL)
		{fprintf(stderr, "Cannot load file %s !\n", file_name);
		return ;}

		load_i(file, &Npoints);
   	fprintf(stderr,"N: %d\n",Npoints);
   	x = newvec(Npoints);
   	y = newvec(Npoints);

   	for(int i = 0; i < Npoints ; i++) {
      	load_d(file, &x[i]);
      	load_d(file, &y[i]);
		}

		for(int i=0;i<Npoints-1;i++)
		{

         fprintf(stderr,"%d   %f  %f\n",i,x[i],y[i]);

	 	 	glColor3f(0.0f, 0.3f, 0.3f);
			circle(x[i],y[i],(x2-x1)/((double)Npoints*2));
		   glColor3f(1.0f, 0.0f, 0.0f);
         
			glBegin(GL_LINES);
				glVertex2f(x[i],y[i]);
				glVertex2f(x[i+1],y[i+1]);
			glEnd();
		}


	 
    }

private: 

   float xTrans, yTrans, zTrans;
   float scale;
   float xRot, yRot, zRot;
	char *file_name;



};

class MainWindow : public QWidget {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = NULL,char * file = NULL) : QWidget(parent) {
        // Layout y botón
		  file_name = newstr(200);
        strcpy(file_name,file); 
        QVBoxLayout *layout = new QVBoxLayout(this);
        QPushButton *button = new QPushButton("Calculate", this);
        GLWidget *glWidget = new GLWidget(this,file_name);
        glWidget->setMinimumSize(400, 400);

        layout->addWidget(button);
        layout->addWidget(glWidget);

        // Conexión Qt4 (SIGNAL/SLOT)
        connect(button, SIGNAL(clicked()), this, SLOT(calcular()));
    }



public slots:
	void calcular() {

		double xi = 0.0;   // Inicio del rango
		double xf = 10.0;  // Fin del rango
		int puntos = 100; // Número de puntos
	
		save_func(my_sin, xi, xf, puntos, file_name);

		qDebug() << "Botón presionado!";
		// Aquí va tu lógica de cálculo.

	}
private:
	char *file_name;
};

int main(int argc, char *argv[])
{

    QApplication app(argc, argv);
    fprintf(stderr,"file_name %s\n",argv[1]);
    MainWindow mainWindow(NULL,argv[1]);
    mainWindow.setWindowTitle("OpenGL + Botón en Qt4");
    mainWindow.resize(400, 500);
    mainWindow.show();    
    return app.exec();
}

#include "main.moc"  // Necesario para que Qt4 procese los slots


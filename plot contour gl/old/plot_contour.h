
#ifndef __PLOT_CONTOUR_GL__
#define __PLOT_CONTOUR_GL__


#include "glplot.h"
#include "complex.h"
#include <QtGui>
#include "vectord.h"
#include "matrixd.h"




class plot_contour:public glplot
{
   public:
	plot_contour();

	~plot_contour(){};
	
	void    setx1(float );
	void    setx2(float );
	void    sety1(float );
	void    sety2(float );
  
          
   private:
   
   double xc1;
   double xc2;
   double yc1;
   double yc2;
   int   NX;
   int   NY;   
      
	void  initializeGL(void);
	void  initDisplayLists();
	void  paintGL(void);
	GLint list;


};

typedef class plot_contour plot_contour;

#endif

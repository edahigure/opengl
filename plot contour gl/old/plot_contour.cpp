#include "plot_contour.h"

plot_contour::plot_contour()
:glplot(0,0,0,0)
{

	xc1=0.0;
	xc2=1.0;
	yc1=0.0;
	yc2=1.0;
	NX=50;
	NY=50;   

}

void plot_contour::initializeGL()
{  
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glDepthFunc( GL_LESS );
   initDisplayLists();   
}


void plot_contour::paintGL()
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glMatrixMode( GL_PROJECTION );
	 glLoadIdentity();
	 glOrtho( x1, x2, y1, y2, 5.0, 15.0 );
	 glMatrixMode( GL_MODELVIEW );

    glLoadIdentity();
    transform();    
    glCallList(list);  
}


double an(double L,double H,int n)
{
	return -2/(sinh(-n*PI*L/H))*1/(n*PI)*(cos(n*PI)-1);
}

double harmonic_f(double L,double H,double x,double y)
{
   double sum=0;
   for(int n=1;n<=50;n++)
   {
    	sum=sum+an(L,H,n)*sin(n*PI*y/H)*sinh(n*PI/H*(x-L));
    	
   } 
//	printf("%f\n",sum);      
   return sum;    
}


void plot_contour::initDisplayLists()
{
 
	vectord x,y;
   matrixd z;
   float x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
   float n1[3],n2[3],n3[3],n4[3],v1[3],v2[3],v3[3],v4[3];
   x.table(xc1,xc2,fabs(xc2-xc1)/(float)NX);
   y.table(yc1,yc2,fabs(yc2-yc1)/(float)NX);
   z.zeros(x.length(),y.length());

   list = glGenLists(1);
   glNewList(list, GL_COMPILE);      
//   light();   

   glColor3f (1.0, 1.0, 1.0);   
   
   glBegin(GL_LINES);			     
   glVertex3f(0,0,0);
   glVertex3f(1,0,0);         
   glEnd();			         
   
   glBegin(GL_LINES);			     
   glVertex3f(0,0,0);
   glVertex3f(0,1,0);         
   glEnd();			         

   glBegin(GL_LINES);			     
   glVertex3f(0,0,0);
   glVertex3f(0,0,1);         
   glEnd();			         
   

   glBegin(GL_LINES);			     
   glVertex3f(0,0,0);
   glVertex3f(1,0,0);         
   glEnd();			         

	

   glColor3f (0.0, 1.0, 1.0);   
   
   
	for(int i=1;i<=x.length();i++)
	{
		 for(int j=1;j<=y.length();j++)
		 {
           z[i][j]=harmonic_f(1,1,x[i],y[j]);
		 }
	}	    
   
	for(int i=1;i<=x.length()-1;i++)
	{
  	    for(int j=1;j<=y.length()-1;j++)
	    {
	      x1=x[i];
	      y1=y[j];
	      z1=z[i][j];
	      
	      x2=x[i+1];
	      y2=y[j];
	      z2=z[i+1][j];
	      
	      x3=x[i+1];
	      y3=y[j+1];
	      z3=z[i+1][j+1];
	      
	      x4=x[i];
	      y4=y[j+1];
	      z4=z[i][j+1];

			v1[0]=x2-x1;
			v1[1]=y2-y1;
			v1[2]=z2-z1;

			v2[0]=x3-x2;
			v2[1]=y3-y2;
			v2[2]=z3-z2;
				
			v3[0]=x4-x3;
			v3[1]=y4-y3;
			v3[2]=z4-z3;

			v4[0]=x1-x4;
			v4[1]=y1-y4;
			v4[2]=z1-z4;				      

         normcrossprod( v1, v2, n1); 				   
         normcrossprod( v2, v3, n2); 				   
         normcrossprod( v3, v4, n3); 				   
         normcrossprod( v4, v1, n4); 				   
	      	      	       
	    	glBegin(GL_POLYGON);			         
				glColor3f (z1, 0.0, 1-z1);		
			   glNormal3fv(n1);								
				glVertex3f(x1,y1,z1);
		
				glColor3f (z2, 0.0, 1-z2);	
				glNormal3fv(n2);				
				glVertex3f(x2,y2,z2);
		
				glColor3f (z3, 0.0, 1-z3);							
 	 		   glNormal3fv(n3);								
				glVertex3f(x3,y3,z3);															

				glColor3f (z4, 0.0, 1-z4);							
				glNormal3fv(n4);				
				glVertex3f(x4,y4,z4);																			
			glEnd();  			
	    }
   }	      
	glEndList();        
}



void plot_contour::setx1(float xt1)
{
    x1=xt1;    
    initDisplayLists();       
	 paintGL();
}



void plot_contour::setx2(float xt2)
{
	x2=xt2;
   initDisplayLists();   	
	paintGL();
}


void plot_contour::sety1(float yt1)
{  
	y1=yt1;
   initDisplayLists();   	
	paintGL();
}

void plot_contour::sety2(float yt2)
{
   y2=yt2;
   initDisplayLists();      
	paintGL();
}



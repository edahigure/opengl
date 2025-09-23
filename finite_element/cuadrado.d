
#*******************************************************

  Everything enclosed inside the cashes is a comment.

  This input file is used to generate the triangular
  mesh over the following domain:

     3--------------------2
     |                    |
     |    5----6          |
     |    |    |          |
     |    |    |          |
     |    4----7          |
     |                    |
     0--------------------1


  Run EasyMesh with the command:

  > EasyMesh example +fig
 
  if you want to see the results with xfig,
  or with

  > EasyMesh example +dxf

  if you want to see the results with some tool that
  suports Autodesk DXF format

*******************************************************#

#*********
  POINTS
*********#
4 # number of points defining the boundary #

# rectangular domain #
#-------+-----+-----+---------+--------#
# point |  x  |  y  | spacing | marker #
#-------+-----+-----+---------+--------#
   0:     0.0   0.0    0.1       1
   1:     1.0   0.0    0.1       1
   2:     1.0   1.0    0.1       1
   3:     0.0   1.0    0.1       1



#***********
  SEGMENTS
***********#
4 # number of segments #

# domain #
#---------+-----+-----+--------#
# segment | 1st | 2nd | marker #
#---------+-----+-----+--------#
     0:      0     1      1
     1:      1     2      1
     2:      2     3      1
     3:      3     0      1


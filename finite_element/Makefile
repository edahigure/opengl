#
# TransFig makefile
#

all: mesh.epic 

# translation into epic

mesh.tex: mesh.fig Makefile
	fig2dev -L epic mesh.fig mesh.tex
clean::
	rm -f mesh.epic


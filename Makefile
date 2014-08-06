LD=`root-config --libs` 
#For the latest Ubuntu system, we need -Wl,--no-as-needed option to link the libs of ROOT
AS_NEEDED = $(shell echo "void main() {return 0;}" | g++ -Wl,--no-as-needed -xc -o /tmp/a.out - 2>/dev/null && echo as_need)
ifeq ($(AS_NEEDED),as_need)
	OPT=-Wl,--no-as-needed
else
	OPT=
endif
INC=-I`root-config --incdir`

application:applicationMain.o lex.yy.o y.tab.o kMvaApplication.o
	g++ $(OPT) -fopenmp  -o $@ $^ $(LD)
applicationMain.o:applicationMain.cpp
	g++ -fopenmp $(INC) -c $^ -o $@
kMvaApplication.o:kMvaApplication.cpp
	g++ $(INC) -fopenmp -c $^ -o $@
classification:classificationMain.o lex.yy.o y.tab.o kMvaClassification.o iBalltreeDensity.o iBalltree.o
	g++ -fopenmp $(OPT) -o $@ $^ $(LD)
classificationMain.o:classificationMain.cpp
	g++ -fopenmp $(INC) -c $^ -o $@ 
kMvaClassification.o:kMvaClassification.cpp
	g++ $(INC) -fopenmp -c $^ -o $@
iBalltreeDensity.o:iBalltreeDensity.cpp
	g++ $(INC) -fopenmp -c $^ -o $@
iBalltree.o:iBalltree.cpp
	g++ $(INC) -fopenmp -c $^ -o $@	
lex.yy.o:lex.yy.c y.tab.c
	gcc -c lex.yy.c -o lex.yy.o
y.tab.o:y.tab.c
	gcc -c y.tab.c -o y.tab.o
lex.yy.c:mva.l
	flex mva.l
y.tab.c:mva.y
	yacc -d mva.y
clean:
	rm -fr *.o classification application lex.yy* y.tab*


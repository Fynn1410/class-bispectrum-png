INCLUDES =  -I../Include/ -I/opt/local/include -I../Library/Cuba-4.2.1 -I../Class/include  
LIBS =  -L/opt/local/lib  -L../Class/lib -L../Library/Cuba-4.2.1
FLAGS =  -lgsl -lgslcblas -lcuba  -lclass -lm -o
OFILESA = main.o utilities.o cosmology.o 
OFILESB = wnw_split.o IR_res.o ps_halo_1loop.o

limHaloPT: ${OFILESA} ${OFILESB}
	mpicc ${OFILESA}${OFILESB} -fopenmp -Wall $(INCLUDES) $(LIBS) $(FLAGS) limHaloPT     
main.o: main.c 
	mpicc -g -c -fopenmp main.c -o main.o $(INCLUDES) 
utilities.o: utilities.c ../Include/utilities.h
	mpicc -g -c -fopenmp utilities.c -o utilities.o $(INCLUDES)
cosmology.o: cosmology.c ../Include/cosmology.h
	mpicc -g -c -fopenmp cosmology.c -o cosmology.o $(INCLUDES)
wnw_split.o: wnw_split.c ../Include/wnw_split.h
	mpicc -g -c  -fopenmp wnw_split.c -o wnw_split.o $(INCLUDES)
IR_res.o: IR_res.c ../Include/IR_res.h
	mpicc -g -c  -fopenmp IR_res.c -o IR_res.o $(INCLUDES)
ps_halo_1loop.o: ps_halo_1loop.c ../Include/ps_halo_1loop.h
	mpicc -g -c  -fopenmp ps_halo_1loop.c -o ps_halo_1loop.o $(INCLUDES)

clean:
	rm *.o limHaloPT

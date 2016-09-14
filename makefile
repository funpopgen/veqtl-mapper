#Replace with location of static GSL libbrary if you want to make statically linked version
GSL = /usr/local/software/gsl-install/
LDC = ldc
GDC = gdc
DMD = dmd

DSOURCES = src/main.d src/arg_parse.d src/read_data.d src/calculation.d src/run_analysis.d

ldc : ${DSOURCES} src/beta.o
	${LDC} -release -enable-inlining -O -w -oq ${DSOURCES} src/beta.o -L-lgsl -L-lgslcblas -of="bin/VEQM"
	rm -f src/*.o *.o

test : ${DSOURCES} src/beta.o
	${LDC} -d-debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/beta.o -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

static : ${DSOURCES} src/static_beta.o
	${LDC} -release -enable-inlining -O -w -oq -d-version=STATICLINKED ${DSOURCES} src/static_beta.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="bin/VEQM"
	rm -f src/*.o *.o

static_test : ${DSOURCES} src/static_beta.o
	${LDC} -d-debug -g -unittest -w -d-version=STATICLINKED ${DSOURCES} src/static_beta.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

gdc : ${DSOURCES} src/beta.o
	${GDC} -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_VEQM ${DSOURCES} src/beta.o -lgsl -lgslcblas -O3 -o bin/VEQM
	strip bin/VEQM
	rm -f src/*.o

gdc_test : ${DSOURCES} src/beta.o
	${GDC} -fdebug -g -funittest -Werror -Wall ${DSOURCES} -lgsl -lgslcblas src/beta.o -o unittest
	./unittest
	rm -f unittest src/*.o *.o

gdc_static : ${DSOURCES} src/static_beta.o
	${GDC} -frelease -finline-functions -O3 -Werror -Wall -fversion=STATICLINKED ${DSOURCES} src/static_beta.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -O3 -o bin/VEQM
	strip bin/VEQM
	rm -f src/*.o

gdc_static_test : ${DSOURCES} src/static_beta.o
	${GDC} -fdebug -g -funittest -Werror -Wall -fversion=STATICLINKED ${DSOURCES} ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a src/static_beta.o -o unittest
	./unittest
	rm -f unittest src/*.o *.o

dmd : ${DSOURCES} src/beta.o
	${DMD} -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas ${DSOURCES} src/beta.o -ofbin/VEQM
	rm src/*.o

dmd_test : ${DSOURCES} src/beta.o
	${DMD} -debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/beta.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

dmd_static : ${DSOURCES} src/static_beta.o
	${DMD} -O -release -noboundscheck -version=STATICLINKED -inline ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a ${DSOURCES} src/static_beta.o -ofbin/VEQM
	rm src/*.o

dmd_static_test : ${DSOURCES} src/static_beta.o
	${DMD} -debug -g -unittest -w -version=STATICLINKED ${DSOURCES} src/static_beta.o ${GSL}/lib/libgsl.a ${GSL}/lib/libgslcblas.a -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

src/static_beta.o : src/beta.c
	cc -c src/beta.c -I{GSL}/include/ -o src/static_beta.o

.PHONY : test ldc static static_test dmd dmd_test dmd_static dmd_static_test gdc gdc_test gdc_static gdc_static_test clean install

clean :
	rm -f src/*.o *.o bin/VEQM

install : ${DSOURCES} VEQM.1
	cp -v $(shell pwd)/bin/VEQM /usr/local/bin/
	cp -v $(shell pwd)/VEQM.1 /usr/local/man/man1/

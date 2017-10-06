LDC := ldc
DMD := dmd

D_SOURCES := ${wildcard src/*.d}

CHECK_LDC := ${shell command ${LDC} | grep -i 'ldc' 2> /dev/null}
CHECK_DMD := ${shell command ${DMD} | grep -i 'dmd' 2> /dev/null}

ifneq (${CHECK_LDC},)
	COMPILER := ${LDC}
	RELEASE_FLAGS := -Jviews -release -enable-inlining -O -w -oq
	DEBUG_FLAGS := -Jviews -d-debug -g -unittest -w
else
	COMPILER := ${DMD}
	RELEASE_FLAGS := -Jviews -release -inline -O -noboundscheck 
	DEBUG_FLAGS := -Jviews -debug -g -unittest -w
endif

ifeq (${CHECK_LDC},)
ifeq (${CHECK_DMD},)
${error No D compiler found at ${LDC} or ${DMD}}
endif
endif

CLEAN_OBJECTS := rm -f src/*.o bin/*.o *.o

bin/veqtl-mapper :	${D_SOURCES} src/beta.o views/commit
	${COMPILER} ${RELEASE_FLAGS} ${D_SOURCES} src/beta.o -of="bin/veqtl-mapper"
	${CLEAN_OBJECTS}

test :	${D_SOURCES} src/beta.o views/commit
	${COMPILER} ${DEBUG_FLAGS} ${D_SOURCES} src/beta.o -of="unittest"
	./unittest
	${CLEAN_OBJECTS} unittest

views/commit :	${D_SOURCES} src/beta.c
	mkdir -p views
	git rev-parse --short HEAD > views/commit

.PHONY : test clean

clean :
	${CLEAN_OBJECTS} bin/veqtl-mapper

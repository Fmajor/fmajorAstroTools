inpath:=./
outpath:=../lib/
outpath:=$(abspath ${outpath})/
#outpath:=./
CC:=gcc

all: test debug

test:
	${CC} makeheaders.c -o makeheaders ${CFLAGS}
	./makeheaders astroMath.c
	cat astroMath.struct astroMath.h > _astroMath.h
	mv _astroMath.h astroMath.h
	${CC} astroMath.c -fPIC -shared -O3 -o astroMath.o ${CFLAGS}
	${CC} testMath.c astroMath.c -o testMath ${CFLAGS}

run:
	./testMath

debug:
	-mkdir -p ${outpath}
	${CC} ${inpath}makeheaders.c -o ${outpath}makeheaders ${CFLAGS}
	${outpath}makeheaders ${inpath}nrlib.c
	${CC} ${inpath}astroMath.c -fPIC -shared -ggdb -DDEBUG -o ${outpath}astroMath.so ${CFLAGS}
	cp ${outpath}astroMath.so ${outpath}libastroMath.so

.PHONY: clean
clean:
	-rm ${outpath}*.so
	-rm -r ${outpath}*.so.dSYM

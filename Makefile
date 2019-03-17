PROG= SingleCellTools

all: $(PROG)

HTSDIR = htslib-1.9
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
LIBA = src/liba.a

CC       = gcc
CFLAGS   = -Wall -O0 -g
DFLAGS   =
INCLUDES = -Isrc -I$(HTSDIR)/ 
LIBS = -lz -lbz2 -llzma -pthread -lm

all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) src/$< -o src/$@

LIB_OBJ = \
	barcode_list.o \
	bed_lite.o \
	number.o

liba.a: $(LIB_OBJ)
	@-rm -f src/$@
	$(AR) -rc $@ $(LIB_OBJ)
	-$(RANLIB) $@

test: $(HTSLIB) version.h

SingleCellTools: $(HTSLIB)
	$(CC) $(CFLAGS) $(INCLUDES) src/main.c src/bam_anno.c src/bam_count.c src/barcode_list.c src/bed_lite.c src/number.c $(HTSLIB) $(LIBS) -o $@

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS src/*.[ch] $(HTSDIR)/*.c

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

.PHONY:all clean clean-all distclean install lib tags test testclean 

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

LIB_OBJ = src/barcode_list.o src/bed_lite.o src/number.o src/fastq.o src/thread_pool.o src/kson.o src/json_config.o

AOBJ = src/bam_anno.o src/bam_count.o src/bam_pick.o src/fastq_parse_barcode.o src/dyncut.o src/sam2bam.o src/bam_rmdup.o

liba.a: $(LIB_OBJ)
	@-rm -f src/$@
	$(AR) -rcs src/$@ $(LIB_OBJ)

test: $(HTSLIB) version.h

SingleCellTools: $(HTSLIB) liba.a $(AOBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ src/main.c $(AOBJ) src/liba.a $(HTSLIB) $(LIBS)

src/bam_anno.o: src/bam_anno.c
src/bam_count.o: src/bam_count.c
src/bam_pick.o: src/bam_pick.c
src/fastq_parse_barcode.o: src/fastq_parse_barcode.c
src/dyncut.o: src/dyncut.c
src/sam2bam.o: src/sam2bam.c
src/barcode_list.o: src/barcode_list.c
src/bed_lite.o: src/bed_lite.c
src/number.o: src/number.c
src/fastq.o: src/fastq.c
src/thread_pool.o: src/thread_pool.c
src/json_config.o: src/json_config.c
src/kson.o: src/kson.c
src/bam_rmdup.o: src/bam_rmdup.c

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm src/*.o src/liba.a

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS src/*.[ch] $(HTSDIR)/*.c

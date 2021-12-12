PROG= PISA

all: $(PROG)

HTSDIR = third_party/htslib-1.10.2
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSVERSION = $(HTSDIR)/version.h

#FMLDIR = third_party/fermi-lite
#include $(FMLDIR)/fermi.mk
#FMLLIB = $(FMLDIR)/libfml.a

LIBA = src/liba.a

ZLIBDIR= third_party/zlib-1.2.11
include $(ZLIBDIR)/zlib.mk
LIBZ = $(ZLIBDIR)/libz.a

CC       = gcc
CFLAGS   = -Wall -O0 -g -D_FILE_OFFSET_BITS=64
DEBUGFLAGS = -fsanitize=address -fno-omit-frame-pointer -O0 -g
DFLAGS   =
INCLUDES = -Isrc -I$(HTSDIR)/ -I. -I$(ZLIBDIR)
LIBS = -lbz2 -llzma -pthread -lm -lcurl 

#all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

pisa_version.h:
	-rm -f pisa_version.h
	echo '#define PISA_VERSION "$(PACKAGE_VERSION)"' > $@

.SUFFIXES:.c .o

.PHONY:all clean clean-all distclean install lib tags test testclean 

force:

.c.o:
	@printf "Compiling $<...                               \r"
	@$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@ || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

LIB_OBJ = src/barcode_list.o \
	src/bam_anno_vcf.o \
	src/bed.o \
	src/number.o \
	src/fastq.o \
	src/kson.o \
	src/json_config.o \
	src/gtf.o \
	src/region_index.o \
	src/dict.o \
	src/ksa.o \
	src/bam_pool.o \
	src/umi_corr.o \
	src/dict.o \
	src/read_thread.o \
	src/read_tags.o \
	src/sim_search.o \
	src/thread.o \
	src/fragment.o \
	src/compactDNA.o \
	src/bam_region.o \
	src/dna_pool.o \
	src/bam_files.o \
	src/biostring.o \
	src/kthread.o

AOBJ = src/bam_anno.o \
	src/bam_count.o \
	src/bam_pick.o \
	src/sam2bam.o \
	src/bam_attr_count.o \
	src/fastq_sort.o \
	src/fastq_parse_barcode.o \
	src/bam_tag_corr.o \
	src/bam2fq.o \
	src/bam_extract_tags.o \
	src/bam_rmdup.o\
	src/gene_fusion.o \
	src/bam_depth.o \
	src/usage.o

liba.a: $(LIB_OBJ)
	@-rm -f src/$@
	$(AR) -rcs src/$@ $(LIB_OBJ)

test: $(HTSLIB) $(HTSVERSION)

PISA: $(HTSLIB) $(LIBZ) liba.a $(AOBJ) pisa_version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ src/main.c $(AOBJ) src/liba.a $(HTSLIB) $(LIBS) $(LIBZ)

debug: $(HTSLIB) $(LIBZ) liba.a $(AOBJ) pisa_version.h 
	$(CC) $(DEBUGFLAGS) $(INCLUDES) -o PISA src/main.c $(AOBJ) src/liba.a $(HTSLIB) $(LIBS) $(LIBZ)

src/sim_search.o: src/sim_search.c
src/bam_depth.o: src/bam_depth.c
src/bam2fq.o: src/bam2fq.c
src/bam_anno.o: src/bam_anno.c
src/bam_count.o: src/bam_count.c pisa_version.h
src/bam_pick.o: src/bam_pick.c
src/bam_anno_vcf.o: src/bam_anno_vcf.c
src/bam_tag_corr.o: src/bam_tag_corr.c
src/umi_corr.o: src/umi_corr.c
src/fastq_parse_barcode.o: src/fastq_parse_barcode.c
src/fastq_sort.o: src/fastq_sort.c
src/dict.o: src/dict.c
src/sam2bam.o: src/sam2bam.c
src/barcode_list.o: src/barcode_list.c
src/bed.o: src/bed.c
src/number.o: src/number.c
src/gtf.o: src/gtf.c
src/region_index.o: src/region_index.c
src/dict.o: src/dict.c
src/fastq.o: src/fastq.c
src/json_config.o: src/json_config.c
src/kson.o: src/kson.c
src/bam_attr_count.o: src/bam_attr_count.c
src/read_thread.o: src/read_thread.c
src/read_tags.o: src/read_tags.c
src/ksa.o: src/ksa.c
src/bam_pool.o: src/bam_pool.c
src/bam_extract_tags.o: src/bam_extract_tags.c
src/usage.o:src/usage.c
src/bam_rmdup.o:src/bam_rmdup.c
src/dna_pool.o:src/dna_pool.c
src/gene_fusion.o:src/gene_fusion.c
src/bam_files.o:src/bam_files.c
src/biostring.o:src/biostring.c
src/kthread.o:src/kthread.c

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) pisa_version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm src/*.o src/liba.a

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS src/*.[ch] $(HTSDIR)/*.c

PROG= PISA

all: $(PROG)

HTSDIR = third_party/htslib-1.10.2
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSVERSION = $(HTSDIR)/version.h
FMLDIR = third_party/fermi-lite
LIBA = src/liba.a
FMLLIB = $(FMLDIR)/libfml.a

CC       = gcc
CFLAGS   = -Wall -O0 -g -D_FILE_OFFSET_BITS=64
DFLAGS   =
INCLUDES = -Isrc -I$(HTSDIR)/ -I. -I $(FMLDIR)
LIBS = -lbz2 -llzma -pthread -lm -lcurl

#all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

pisa_version.h:
	echo '#define PISA_VERSION "$(PACKAGE_VERSION)"' > $@

.SUFFIXES:.c .o

.PHONY:all clean clean-all distclean install lib tags test testclean 

force:

.c.o:
	@printf "Compiling $<...                               \r"
	@$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@ || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

LIB_OBJ = src/barcode_list.o src/bed_lite.o src/number.o src/fastq.o src/thread_pool.o src/kson.o src/json_config.o src/gtf.o src/dict.o src/seq_merge.o src/ksa.o \
	src/bam_pool.o src/umi_corr.o src/dict.o src/read_thread.o src/read_tags.o src/sim_search.o

AOBJ = src/bam_anno.o \
	src/bam_count.o \
	src/bam_pick.o \
	src/sam2bam.o \
	src/bam_rmdup.o \
	src/bam_attr_count.o \
	src/fastq_sort.o \
	src/fastq_parse_barcode.o \
	src/fastq_segment.o \
	src/fastq_assem.o \
	src/LFR_impute.o \
	src/bam_tag_corr.o \
	src/bam2fq.o \
	src/gene_cov.o \
	src/bam_extract_tags.o \
	src/fragment.o

ASSM_LIB_OBJ =	fermi-lite/bfc.o fermi-lite/bseq.o fermi-lite/bubble.o fermi-lite/htab.o fermi-lite/ksw.o fermi-lite/kthread.o fermi-lite/mag.o fermi-lite/misc.o \
	fermi-lite/mrope.o fermi-lite/rld0.o fermi-lite/rle.o fermi-lite/rope.o fermi-lite/unitig.o

liba.a: $(LIB_OBJ)
	@-rm -f src/$@
	$(AR) -rcs src/$@ $(LIB_OBJ)

libfml.a: $(ASSM_LIB_OBJ)
	@-rm -f $(FMLDIR)/$@
	$(AR) -rcs $(FMLDIR)/$@ $(ASSM_LIB_OBJ)

test: $(HTSLIB) $(HTSVERSION)

PISA: $(HTSLIB) liba.a $(AOBJ) pisa_version.h libfml.a 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ src/main.c $(AOBJ) $(FMLLIB) src/liba.a $(HTSLIB) -lz $(LIBS)

src/sim_search.o: src/sim_search.c
src/fragment.o: src/fragment.c
src/gene_cov.o: src/gene_cov.c
src/bam2fq.o: src/bam2fq.c
src/seq_merge.o: src/seq_merge.c
src/bam_anno.o: src/bam_anno.c
src/bam_count.o: src/bam_count.c
src/bam_pick.o: src/bam_pick.c
src/bam_tag_corr.o: src/bam_tag_corr.c
src/umi_corr.o: src/umi_corr.c
src/fastq_parse_barcode.o: src/fastq_parse_barcode.c
src/fastq_sort.o: src/fastq_sort.c
src/dict.o: src/dict.c
src/sam2bam.o: src/sam2bam.c
src/barcode_list.o: src/barcode_list.c
src/bed_lite.o: src/bed_lite.c
src/number.o: src/number.c
src/gtf.o: src/gtf.c
src/dict.o: src/dict.c
src/fastq.o: src/fastq.c
src/thread_pool.o: src/thread_pool.c
src/json_config.o: src/json_config.c
src/kson.o: src/kson.c
src/bam_rmdup.o: src/bam_rmdup.c
src/bam_attr_count.o: src/bam_attr_count.c
src/fastq_segment.o: src/fastq_segment.c
src/fastq_assem.o: src/fastq_assem.c
src/read_thread.o: src/read_thread.c
src/read_tags.o: src/read_tags.c
src/ksa.o: src/ksa.c
src/bam_pool.o: src/bam_pool.c
src/LFR_impute.o: src/LFR_impute.c
src/bam_extract_tags.o: src/bam_extract_tags.c
fermi-lite/bfc.o: fermi-lite/bfc.c
fermi-lite/bseq.o: fermi-lite/bseq.c
fermi-lite/bubble.o: fermi-lite/bubble.c
fermi-lite/htab.o: fermi-lite/htab.c
fermi-lite/ksw.o: fermi-lite/ksw.c
fermi-lite/kthread.o: fermi-lite/kthread.c
fermi-lite/mag.o: fermi-lite/mag.c
fermi-lite/misc.o: fermi-lite/misc.c
fermi-lite/mrope.o: fermi-lite/mrope.c
fermi-lite/rld0.o: fermi-lite/rld0.c
fermi-lite/rle.o: fermi-lite/rle.c
fermi-lite/rope.o: fermi-lite/rope.c
fermi-lite/unitig.o: fermi-lite/unitig.c

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

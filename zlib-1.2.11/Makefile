LIBA = libz.a

all: $(LIBA)

CC       = gcc
CFLAGS   = -Wall -O0 -g -DHAVE_HIDDEN

.SUFFIXES:.c .o

.PHONY:all

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


#libtool -o libz.a 

LIB_OBJ = adler32.o crc32.o deflate.o infback.o inffast.o inflate.o inftrees.o trees.o zutil.o compress.o uncompr.o gzclose.o gzlib.o gzread.o gzwrite.o

libz.a: $(LIB_OBJ)
	@-rm -f $@
	$(AR) -rcs $@ $(LIB_OBJ)

adler32.o: adler32.c
compress.o: compress.c
crc32.o: crc32.c
deflate.o: deflate.c
gzclose.o: gzclose.c
gzlib.o: gzlib.c
gzread.o: gzread.c
gzwrite.o: gzwrite.c
infback.o: infback.c
inffast.o: inffast.c
inflate.o: inflate.c
inftrees.o: inftrees.c
trees.o: trees.c
uncompr.o: uncompr.c
zutil.o: zutil.c


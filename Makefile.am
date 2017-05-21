AM_CPPFLAGS = -Iinclude -Isfmt/include
AM_CFLAGS = -rdc=true -Xptxas -dlcm=ca
bin_PROGRAMS = mpcgs

mpcgs_SOURCES = src/main.c \
                src/phylip.c \
                src/debug.c \
                src/mpcgs.c \
                src/tree.c \
                src/tree-gpu.cu \
                src/mpcgs-gpu.cu \
                sfmt/SFMT.c

.cu.o: 
	$(CC) ${DEFS} ${AM_CFLAGS} ${CFLAGS} -c -o $@ $< ${AM_CPPFLAGS}

   

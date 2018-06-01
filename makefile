CC ?= gcc

CFLAGS = -std=c99 -DNDEBUG #-g
#CFLAGS += -Wall -pedantic

CFLAGS += -DDISABLED_CONFIG_HEADER

CFLAGS += -O3

SRCDIR = src

PHYCDIR = $(SRCDIR)/phyc

CFLAGS += -I$(PHYCDIR)

OBJDIR = obj
BINDIR = bin

ifeq ($(OS),Windows_NT)
    CFLAGS += -D__WIN32__
else
	ARCH := $(shell uname -s)
	ifeq ($(ARCH),Linux)
		CFLAGS += -D__LINUX__ -D_XOPEN_SOURCE=600
	endif
	ifeq ($(ARCH),Darwin)
		CFLAGS += -D__MACH__
	endif
endif

OBJFILES := $(patsubst $(PHYCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(PHYCDIR)/*.c))

LDLIBS += -lm -lgsl -lgslcblas


MKDIR_P = mkdir -p

.PHONY: directories

all: clean directories physher modelavg bootstrap

directories: ${OBJDIR} ${BINDIR}

${OBJDIR}:
	${MKDIR_P} ${OBJDIR}
${BINDIR}:
	${MKDIR_P} ${BINDIR}


$(OBJDIR)/%.o: $(PHYCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS)
	
physher: ${OBJFILES}
	$(CC) $(CFLAGS) $(LDFLAGS) -o ${BINDIR}/physher $(SRCDIR)/physher.c $(OBJFILES) $(LDLIBS)

modelavg: ${OBJFILES}
	$(CC) $(CFLAGS) $(LDFLAGS) -o ${BINDIR}/modelavg $(SRCDIR)/modelAveraging.c $(OBJFILES) $(LDLIBS)

bootstrap: ${OBJFILES}
	$(CC) $(CFLAGS) $(LDFLAGS) -o ${BINDIR}/bootstrap $(SRCDIR)/bootstrap.c $(OBJFILES) $(LDLIBS)

clean:
	rm -f ${BINDIR}/physher
	rm -f ${BINDIR}/modelavg
	rm -f ${BINDIR}/bootstrap
	rm -f $(OBJDIR)/*.o

	


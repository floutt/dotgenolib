LIBNAME = dotgeno

CC      = cc
AR      = ar
ARFLAGS = rcs

PREFIX  ?= /usr/local
LIBDIR  = $(PREFIX)/lib
INCDIR  = $(PREFIX)/include

CFLAGS  ?= -std=c17 -Wall -Wextra -Wpedantic -O2 -lm
PICFLAGS = -fPIC

SRC    = src/dotgeno.c
OBJ    = build/dotgeno.o
PICOBJ = build/dotgeno.pic.o

STATIC_LIB = lib$(LIBNAME).a
SHARED_LIB = lib$(LIBNAME).so

.PHONY: all static shared clean install uninstall

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB)

build:
	mkdir -p build

$(OBJ): $(SRC) src/dotgeno.h src/khash.h | build
	$(CC) $(CFLAGS) -Isrc -c $< -o $@

$(PICOBJ): $(SRC) src/dotgeno.h src/khash.h | build
	$(CC) $(CFLAGS) $(PICFLAGS) -Isrc -c $< -o $@

$(STATIC_LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB): $(PICOBJ)
	$(CC) -shared -o $@ $^

install: all
	mkdir -p $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(INCDIR)

	cp $(STATIC_LIB) $(DESTDIR)$(LIBDIR)/
	cp $(SHARED_LIB) $(DESTDIR)$(LIBDIR)/
	cp src/dotgeno.h $(DESTDIR)$(INCDIR)/

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/$(STATIC_LIB)
	rm -f $(DESTDIR)$(LIBDIR)/$(SHARED_LIB)
	rm -f $(DESTDIR)$(INCDIR)/dotgeno.h

clean:
	rm -rf build
	rm -f $(STATIC_LIB)
	rm -f $(SHARED_LIB)

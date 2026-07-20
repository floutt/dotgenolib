LIBNAME = dotgeno

CC      = cc
AR      = ar
ARFLAGS = rcs

PREFIX  ?= /usr/local
LIBDIR  = $(PREFIX)/lib
INCDIR  = $(PREFIX)/include

CPPFLAGS = -Isrc
CFLAGS   ?= -std=c17 -Wall -Wextra -Wpedantic -O2
PICFLAGS = -fPIC
LDFLAGS  =
LDLIBS   = -lm

SRC    = src/dotgeno.c
OBJ    = build/dotgeno.o
PICOBJ = build/dotgeno.pic.o

STATIC_LIB = lib$(LIBNAME).a
SHARED_LIB = lib$(LIBNAME).so

# Unit tests
UNITY_SRC = Unity/unity.c
UNITY_OBJ = build/unity.o

TEST_SRC = test/test.c
TEST_OBJ = build/test.o
TEST_BIN = build/test.out

.PHONY: all static shared clean install uninstall test check

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB)

build:
	mkdir -p build

$(OBJ): $(SRC) src/dotgeno.h src/khash.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(PICOBJ): $(SRC) src/dotgeno.h src/khash.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PICFLAGS) -c $< -o $@

$(STATIC_LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB): $(PICOBJ)
	$(CC) -shared $(LDFLAGS) -o $@ $^

$(UNITY_OBJ): $(UNITY_SRC) Unity/unity.h Unity/unity_internals.h | build
	$(CC) $(CFLAGS) -IUnity -c $< -o $@

$(TEST_OBJ): $(TEST_SRC) src/dotgeno.h Unity/unity.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -IUnity -c $< -o $@

$(TEST_BIN): $(TEST_OBJ) $(UNITY_OBJ) $(STATIC_LIB)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test: $(TEST_BIN)
	./$(TEST_BIN)

check: test

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

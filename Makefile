CC = gcc

CFLAGS = -g -c -Wall -Wpedantic -Werror
LDFLAGS = 

SRC = $(wildcard src/*.c)
HDR = $(wildcard src/*.h)
OBJ = $(SRC:src/%.c=bin/obj/%.o)
BIN = bin/exec

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDFLAGS)

bin/obj/%.o: src/%.c $(HDR)
	$(CC) -o $@ $< $(CFLAGS)

clean: 
	rm -rf bin/*
	mkdir bin/obj

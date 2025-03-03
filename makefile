# Variables
#
# variables for compiling rules
SHELL=/bin/bash
CC=gcc
CFLAGS=-Wall -pedantic -std=c17 -g3

# paths for testing/submitting
HW3=/c/cs223/hw3
BIN=/c/cs223/bin


# Compiling
#
# compile Heatmap and Unit (test) targets, acts as default target to "make"
all: Heatmap Unit

# compile Heatmap executable, called by default with all
Heatmap: heatmap.o track.o location.o trackpoint.o string_util.o
	${CC} ${CFLAGS} -o $@ $^ -lm

# compile Unit (test) executable, called by default with all
Unit: track_unit.o track.o location.o trackpoint.o
	${CC} ${CFLAGS} -o $@ $^ -lm

# compile object files
heatmap.o: heatmap.h track.h location.h trackpoint.h
track.o: track.h location.h trackpoint.h
string_util.o: string_util.h
trackpoint.o: trackpoint.h location.h
location.o: location.h
track_unit.o: track.h location.h trackpoint.h

# delete the current executable and object files (run to force make to recompile)
clean:
	rm -f Heatmap Unit *.o


# PWD Testing
#
# run the entire test suite within your current directory
test:
	${HW3}/Tests/test.Heatmap

# run a single test using the command, e.g. "make singletest TEST=t001"
singletest:
	BIN="${BIN}" WHERE="${HW3}/Tests" ${HW3}/Tests/$(TEST)


# Submission and Testing
#
# submit your source code/documentation files (change if you created more files)
submit:
	${BIN}/submit 3 makefile heatmap.c heatmap.h string_util.c string_util.h track.c log

# check which files have been submitted and when
check:
	${BIN}/check 3

# run the test suite on your submitted code
testit:
	${BIN}/testit 3 Heatmap

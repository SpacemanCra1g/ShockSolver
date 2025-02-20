CC=clang++
FLAGS= -g -Wall -Wextra -pedantic -Qunused-arguments -march=native
SRC=src
OBJ=build
OUTDIR=OutputData
SRCS1 =$(wildcard $(SRC)/*.cpp)
SRCS2 =$(wildcard ../Source/src/*.cpp)
SRCS = $(SRCS1) $(SRCS2)
OBJS1 =$(patsubst $(SRC)/%.cpp, $(OBJ)/%.o, $(SRCS1))
OBJS2 =$(patsubst ../Source/src/%.cpp, $(OBJ)/%.o, $(SRCS2))
OBJS = $(OBJS1) $(OBJS2)
LINKERS = -std=c++11 -lcblas -fopenmp -llapacke
TARGET=Main.ex
INC_DIR = -I./build/include
OBJSRC=build/src
OBJI=build/include
LFILE=LINKING_FAILED.txt

# OMP_NUM_THREADS=2
# export OMP_NUM_THREADS

all:$(TARGET)

.PHONY: clean run new every print

$(TARGET): $(LFILE) $(OBJS)
	$(CC) $(FLAGS) $(INC_DIR) $(OBJS) -o $@ $(LINKERS) && rm $(LFILE)

$(OBJ)/%.o: $(OBJSRC)/%.cpp
# $(OBJS):  $(OBJS) $(OUTDIR) $(LFILE)
	$(CC) $(FLAGS) $(INC_DIR) -c $< -o $@ $(LINKERS)

$(LFILE): $(OBJ) $(OUTDIR)
	cp ./include/* ../Source/include/* ./build/include/ &&\
	cp ./src/* ../Source/src/* ./build/src/ &&\
	touch $(LFILE)

$(OBJ):
	mkdir $(OBJ) && mkdir $(OBJSRC) $(OBJI)

$(OUTDIR):
	mkdir $(OUTDIR)

clean:
	rm -f $(LFILE) *.ex OutputData/* && rm -r $(OBJ) $(OUTDIR)
new:
	make clean && make
run:
	make new && ./$(TARGET)

print:
	$(info   VAR is $(SRCS))

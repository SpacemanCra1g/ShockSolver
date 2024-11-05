CC=clang++
FLAGS= -O3 -Wall -Wextra -pedantic -Qunused-arguments -march=native
SRC=src
OBJ=build
OUTDIR=OutputData
SRCS=$(wildcard $(SRC)/*.cpp)
OBJS=$(patsubst $(SRC)/%.cpp, $(OBJ)/%.o, $(SRCS))
LINKERS = -std=c++11 -lcblas -fopenmp -llapacke
TARGET=Main.ex
INC_DIR = -I./include


# OMP_NUM_THREADS=2
# export OMP_NUM_THREADS

all:$(TARGET)

.PHONY: clean run new every

$(TARGET): $(OBJS)
	$(CC) $(FLAGS) $(INC_DIR) $(OBJS) -o $@ $(LINKERS)

$(OBJ)/%.o: $(SRC)/%.cpp $(OBJ) $(OUTDIR)
	$(CC) $(FLAGS) $(INC_DIR) -c $< -o $@ $(LINKERS)


$(OBJ):
	mkdir $(OBJ)

$(OUTDIR):
	mkdir $(OUTDIR)

clean:
	rm -f *.ex OutputData/* && rm -r $(OBJ) $(OUTDIR)
new:
	make clean && make
run:
	make new && ./$(TARGET)

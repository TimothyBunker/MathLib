CC = gcc
CFLAGS = -Wall -Wextra -g -Iinclude -O0

SRCS = src/main.c src/vector.c src/matrix.c src/alg.c src/render.c
LDFLAGS = -lm
OBJS = $(SRCS:.c=.o)
TARGET = mathlib 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f $(OBJS) $(TARGET) ./valgrind_output.txt

valgrind: $(TARGET)
	export DEBUGINFOD_URLS="https://debuginfod.archlinux.org/"; valgrind -s --leak-check=full --track-origins=yes ./$(TARGET) > ./valgrind_output.txt 2>&1 

.PHONY: all clean

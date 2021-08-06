CC		= mpicc
COPT	= -lm

LD		= $(CC)
LDFLAGS	= $(COPT)

all: src.x

src.x:
	$(CC) src.c -o src.x -lm

clean:
	rm src.x

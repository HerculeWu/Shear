CC=mpicc
CFLAGE= -lm
source=shear.c

shear: $(source)
	$(CC) $(CFLAGE) -o shear $(source)

clean:
	rm shear
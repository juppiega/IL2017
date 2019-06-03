FC = gfortran
CFLAGS = -O3

TARGET = IL2017_test

all: $(TARGET)

$(TARGET): IL2017_coeff_mod.o IL2017_mod.o $(TARGET).o
	$(FC) -o $@ IL2017_coeff_mod.o IL2017_mod.o $(TARGET).o

%.o: %.f90
	$(FC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.mod *.o

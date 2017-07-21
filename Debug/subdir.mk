################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../IL2017_coeff_mod.f90 \
../IL2017_mod.f90 

OBJS += \
./IL2017_coeff_mod.o \
./IL2017_mod.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

IL2017_coeff_mod.o: ../IL2017_coeff_mod.f90

IL2017_mod.o: ../IL2017_mod.f90 IL2017_coeff_mod.o



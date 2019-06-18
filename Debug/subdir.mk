################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../PRexAnalyzer.cpp 

OBJS += \
./PRexAnalyzer.o 

CPP_DEPS += \
./PRexAnalyzer.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/newdriver/Source/ROOT_CERN/root-6.14.04/root-6.14.04-install/include" -I"/home/newdriver/Source/CODA/2.6.2/Linux-x86_64/include" -I"/home/newdriver/Source/JLab_libs/analyzer/analyzer-1.6.6/analyzer-1.6.6-install/include" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



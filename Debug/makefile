################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: PRexAnalyzer

# Tool invocations
PRexAnalyzer: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -L"/home/newdriver/Source/ROOT_CERN/root-6.14.04/root-6.14.04-install/lib" -L"/home/newdriver/Source/CODA/2.6.2/Linux-x86_64/lib" -L"/home/newdriver/Source/JLab_libs/analyzer/analyzer-1.6.6/analyzer-1.6.6-install/lib64" -Xlinker -L/usr/local/lib `/home/newdriver/Source/ROOT_CERN/root-6.14.04/root-6.14.04-install/bin/root-config --libs --glibs --cflags --ldflags --auxlibs --auxcflags` -pthread -lm -ldl -o "PRexAnalyzer" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) PRexAnalyzer
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets

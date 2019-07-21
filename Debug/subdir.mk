################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../MotionAlgorithm.c \
../MotionTracking.c \
../NE10_random.c \
../fileTesting.c \
../gsl_tests.c \
../main.c \
../polyfit.c \
../seatest.c \
../test_fft.c \
../test_suite_fft_float32.c \
../test_suite_fir.c \
../testneon.c \
../unit_test_common.c 

CPP_SRCS += \
../BgXmlParser.cpp \
../opencv_tests.cpp 

OBJS += \
./BgXmlParser.o \
./MotionAlgorithm.o \
./MotionTracking.o \
./NE10_random.o \
./fileTesting.o \
./gsl_tests.o \
./main.o \
./opencv_tests.o \
./polyfit.o \
./seatest.o \
./test_fft.o \
./test_suite_fft_float32.o \
./test_suite_fir.o \
./testneon.o \
./unit_test_common.o 

C_DEPS += \
./MotionAlgorithm.d \
./MotionTracking.d \
./NE10_random.d \
./fileTesting.d \
./gsl_tests.d \
./main.d \
./polyfit.d \
./seatest.d \
./test_fft.d \
./test_suite_fft_float32.d \
./test_suite_fir.d \
./testneon.d \
./unit_test_common.d 

CPP_DEPS += \
./BgXmlParser.d \
./opencv_tests.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	arm-linux-gnueabihf-g++ -I/home/echocare/work/system/app/opencv/opencv/build_arm/install/include/opencv4 -I/home/echocare/work/system/app/fftw-3.3.8/__install/include -I/home/echocare/work/system/app/libxml2-2.7.2/include -O0 -g3 -Wall -c -fmessage-length=0 -D__ARM -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	arm-linux-gnueabihf-gcc -I/home/echocare/work/system/app/opencv/opencv/build_arm/install/include/opencv4 -I/home/echocare/work/system/app/fftw-3.3.8/__install/include -I/home/echocare/work/system/app/Ne10/inc -I/home/echocare/work/system/app/install/include -I/home/echocare/work/system/app/opencv/opencv/build_arm/install/include/opencv4/opencv2 -I/home/echocare/work/system/app/libxml2-2.7.2/include -O0 -g3 -Wall -c -fmessage-length=0  -D__ARM -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



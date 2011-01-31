CC=g++
CFLAGS=-g -Wall -O3 -mcpu=cortex-a8 -mfpu=neon -ftree-vectorize -mfloat-abi=softfp -ffast-math -fsingle-precision-constant -fomit-frame-pointer -fno-math-errno -fno-signed-zeros
LDFLAGS=-lmathneon -lm 
INCLUDES=-I/usr/include -I/usr/local/include
SOURCES=bmpfile.c Clock.c FeatureRecognition.cpp ProcessImage.cpp ImageProcessingUtilities_Debug.cpp CannyEdgeDetect.cpp ImageProcessingUtilities_DSP.cpp ImageProcessingUtilities_FPGA.cpp MyBitmap.cpp BitmapProcessCppExe.cpp

		 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=BitmapProcessCppExe

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f ./*.o
	rm -f ./BitmapProcessCppExe
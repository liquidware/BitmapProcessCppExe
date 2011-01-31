#include <stdio.h>
#include <fenv.h>
#include <signal.h>
#include <stdlib.h>

#include "Clock.h"
#include "ProcessImage.h"
#include "bmpfile.h"

#define _GNU_SOURCE 1

//void fpehandler(int sig_num)
//{
//        signal(SIGFPE, fpehandler);
//        printf("SIGFPE: floating point exception occured, exiting.\n");
//        abort();
//}

void enable_runfast()
{
	static const unsigned int x = 0x04086060;
	static const unsigned int y = 0x03000000;
	int r;
	asm volatile (
		"fmrx	%0, fpscr			\n\t"	//r0 = FPSCR
		"and	%0, %0, %1			\n\t"	//r0 = r0 & 0x04086060
		"orr	%0, %0, %2			\n\t"	//r0 = r0 | 0x03000000
		"fmxr	fpscr, %0			\n\t"	//FPSCR = r0
		: "=r"(r)
		: "r"(x), "r"(y)
	);
}

void DrawCross(uint8_t* bits, int width, int height,
		int crossHorz, int crossVert, int crossSize, int R, int G, int B)
{
	// Invert verical coordinate to acount for different row order
	int crossV2 = height - crossVert - 1;

	// Horiz line
	for (int i = crossHorz - crossSize; i < crossHorz + crossSize; i++)
	{
		int index = 3 * (i + crossV2 * width);
		bits[index] = B;		// Blue
		bits[index + 1] = G;	// Green
		bits[index + 2] = R;	// Red
	}
	// Vert line
	for (int i = crossV2 - crossSize; i < crossV2 + crossSize; i++)
	{
		int index = 3 * (crossHorz + i * width);
		bits[index] = B;
		bits[index + 1] = G;
		bits[index + 2] = R;
	}
}

void DrawBorder(uint8_t* bits, int width, int height,
		PNT border[], int borderCount, int R, int G, int B)
{
	for (int i=0; i<borderCount; i++)
	{
		PNT* pt = &(border[i]);
		int index = 3 * (pt->X + (height - pt->Y - 1) * width);
		bits[index] = B;
		bits[index + 1] = G;
		bits[index + 2] = R;
	}
}

int CopyBytes(uint8_t** pb, int size)
{
	int b = **pb;
	(*pb)++;
	int g = **pb;
	(*pb)++;
	int r = **pb;
	(*pb)++;

	if (size == 4)
	{
		int a = **pb;
		(*pb)++;
	}
	// Build greyscale image
	return (int)(r * 0.299 + g * 0.587 + b * 0.114);
}

int CreatePixelArray(uint8_t* in, int** out, int width, int height, int bytesPerPixel)
{
	int pixelSize = width * height;
	*out = new int[pixelSize];
	for (int row=0; row<height; row++)
	{
		for (int col=0; col<width; col++)
		{
			int index1 = ((height - 1 - row) * width + col) * bytesPerPixel;
			uint8_t* pIn = in + index1;
			int val = CopyBytes(&pIn, bytesPerPixel);
			int index2 = row * width + col;
			(*out)[index2] = val;
		}
	}
	return pixelSize;
}

int main(void) {
	int width = 320;
	int height = 240;
	int bytesPerPixel = 3;
	int *intArr;

	//int* intArr;

	clock_start("main");
	enable_runfast();
    //feenableexcept(FE_ALL_EXCEPT);
    //signal(SIGFPE, fpehandler);

    bmpfile_t bmp;
    uint8_t buff[320*240*3];
    bmp_create_from_file("KatieGood.bmp", &bmp, buff);

    int i=0;
    bmpfile_t * bmpsave;
    bmpsave = bmp_create(bmp.dib.width, bmp.dib.height, 24);
    for (int y=0; y< bmp.dib.height; y++) {
    	for(int x=0;x < bmp.dib.width; x++) {
    		rgb_pixel_t pixel;
    		pixel.red=buff[i];
    		pixel.green=buff[i+1];
    		pixel.blue=buff[i+2];
    		bmp_set_pixel(bmpsave, x, y, pixel);
    		i+=3;
    	}
    }
    bmp_save(bmpsave,"input.bmp");

    int arrSize = CreatePixelArray(buff, &intArr, width, height, bytesPerPixel);

	// Find pupil and CR
	PUPIL_CR_OBJECT pupil;
	PUPIL_CR_OBJECT CR[MAX_CR_OBJECTS];
	int number_CR_objects;
	FindPupilCR(width, height, intArr, &pupil, CR, &number_CR_objects);

	printf("Results\n pupil.found %d, number_CR_objects %d\n", pupil.found, number_CR_objects);

	// Draw crosses
	uint8_t* inpBytes = (uint8_t*)buff;
	if (pupil.found)
	{
		DrawCross(inpBytes, width, height, (int)pupil.center_horz, (int)pupil.center_vert, 10, 255, 0, 0);
		DrawBorder(inpBytes, width, height, pupil.ellipse_points, pupil.no_ellipse_points, 255, 0, 0);
	}

	if (number_CR_objects > 0)
	{
		for (int i=0; i<number_CR_objects; i++)
		{
			DrawCross(inpBytes, width, height, (int)CR[i].center_horz, (int)CR[i].center_vert, 10, 100, 255, 0);
		}
	}

	//Save the output
    i=0;
    bmpfile_t * bmpout;
    bmpout = bmp_create(bmp.dib.width, bmp.dib.height, 24);
    for (int y=0; y< bmp.dib.height; y++) {
    	for(int x=0;x < bmp.dib.width; x++) {
    		rgb_pixel_t pixel;
    		pixel.red=buff[i];
    		pixel.green=buff[i+1];
    		pixel.blue=buff[i+2];

    		bmp_set_pixel(bmpout, x, y, pixel);
    		i+=3;
    	}
    }
    bmp_save(bmpout,"output.bmp");

	//delete[] intArr;

	clock_end();
	return 0;
}

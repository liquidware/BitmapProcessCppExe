#include <mathneon.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "Clock.h"
#include "ProcessImage.h"
#include "ImageProcessingUtilities_FPGA.h"

ImageProcessingUtilities_FPGA::ImageProcessingUtilities_FPGA(void)
{
}

ImageProcessingUtilities_FPGA::~ImageProcessingUtilities_FPGA(void)
{
}

/*
struct
{
  int64_t val[2];
} int64x2_t;

struct
{
  int32_t val[2];
} int32x2_t;

struct
{
  float val[4];
} float32x4_t;
*/
#if 0
float32x4_t transform (float32x4_t * matrix, float32x4_t vector)
{
  /* in a perfect world this code would compile into just four instructions */
  float32x4_t result;

  result = vml (matrix[0], vector);
  result = vmla (result, matrix[1], vector);
  result = vmla (result, matrix[2], vector);
  result = vmla (result, matrix[3], vector);

  return result;
}
#endif

//Apply a 3x3 convolution to 9 pixels
// Assuming following format:
// kernel (b):
//  [][][][0]
//  [][][][0]
//  [][][][0]
//
// source (c):
//  [][][][dummy]...320 pixel stride
//  [][][][dummy]...320 pixel stride
//  [][][][dummy]
//
// 1.) result = b*c+a
// 2.) result += b*c+a
// 3.) result += b*c+a
void neon_conv(int32_t a[4],
				   int32_t b[4],
				   int32_t c[4])
{
  int32x4_t ain; //= vdupq_n_s32(0); // clear accumulators
  int32x4_t bin = vld1q_s32(b); //load 4 ints into neon
  int32x4_t cin = vld1q_s32(c); //load 4 ints into neon

  ain = vmlaq_s32(ain, bin, cin); //1.) multiply-accumulate 4 int's in one instruction

  bin = vld1q_s32(b+4); //load
  cin = vld1q_s32(c+320); //load

  ain = vmlaq_s32(ain, bin, cin); //2.) multiply-accumulate 4 int's in one instruction

  bin = vld1q_s32(b+8); //load
  cin = vld1q_s32(c+640); //load

  ain = vmlaq_s32(ain, bin, cin); //3.) multiply-accumulate 4 int's in one instruction

  vst1q_s32(a, ain); //store
}

/**
 * Convolution Example Run Notes: 
 *   //pixel
    sum=58.000000,kernel[0][0],source[1605]
    sum=474.000000,kernel[1][0],source[1606]
    sum=1209.000000,kernel[2][0],source[1607]
    
    sum=1393.000000,kernel[0][1],source[1925]
    sum=2929.000000,kernel[1][1],source[1926]
    sum=5971.000000,kernel[2][1],source[1927]
    
    sum=6601.000000,kernel[0][2],source[2245]
    sum=9721.000000,kernel[1][2],source[2246]
    sum=13698.000000,kernel[2][2],source[2247]
    
    //pixel
    sum=104.000000,kernel[0][0],source[1606]
    sum=524.000000,kernel[1][0],source[1607]
    sum=1259.000000,kernel[2][0],source[1608]
    sum=1643.000000,kernel[0][1],source[1926]
    sum=3515.000000,kernel[1][1],source[1927]
    sum=6323.000000,kernel[2][1],source[1928]
    sum=7163.000000,kernel[0][2],source[2246]
    sum=9685.000000,kernel[1][2],source[2247]
    sum=14031.000000,kernel[2][2],source[2248]

 */
void ImageProcessingUtilities_FPGA::GrayscaleConvolution( int *psw, int *psh, int *wL, int *wR, int *wT, int *wB, int source[], int *ksize, int kernel[5][5], float *gain, int result[] )
{
    //Note: we are assuming that any kernel_element*255 is significantly smaller than integer resolution
    //psw and psh are pointers to width and height of image array; wL, wR, wT, and wB are pntrs to boundaries (row or column numbers) of a window within image array
    //row and col numbers start from 0, so max row number is image height-1 and max col number is image width-1.
    //The image array ("source") is a single subscripted array, where subscript == col_number + (image widthe * row_number).
    //The kernel is a square, double subscript array kernel[col][row].  Column and row numbers start from 0. 
    //The "result" is returned in the same single subscript array form as the "source".  
    bool sstop = false; //debug
    int kw, kh, sw, sh, outrow, outcol, krow, kcol, inrow, incol, OutcolOutrow, IncolInrow, kwOfst, khOfst;
    //float sum, ksum, temp;
    int ksum;
    int temp;
    int sum;
    long maxint = 32767.0;
    //float maxint = 2147483647;
    kw = *ksize;
    kh = *ksize;
    sw = *psw;
    sh = *psh;
    kwOfst = (kw-1)/2; 
    khOfst = (kh-1)/2;
    ksum = 0;
    clock_start("ImageProcessingUtilities_FPGA::GrayscaleConvolution");

    /* Construct a kernel for the neon vector processor */
    int kernel_neon[4][4]  = {
            {kernel[0][0], kernel[1][0], kernel[2][0], 0},
            {kernel[0][1], kernel[1][1], kernel[2][1], 0},
            {kernel[0][2], kernel[1][2], kernel[2][2], 0},
    };

    //get sum of kernel elements for normalization
    for (krow = 0; krow < kw; krow++)
        for (kcol = 0; kcol < kh; kcol++)
        {
            ksum += kernel[kcol][krow];
        }
    if (ksum <= 0) ksum = 1;//protect against divide by zero
    //printf("ksum=%d\n", ksum);

    //zero (kh-1)/2 rows and (kw-1)/2 cols on each edge of  the result array (we cannot compute the convolution on these boarders)
    for (outrow = *wT; outrow < *wT + (khOfst); outrow++)
        for (outcol = *wL; outcol <= *wR; outcol++)
        {
            OutcolOutrow = outcol + (sw * outrow);
            result[OutcolOutrow] = 0;
        }
    for (outrow = *wB; outrow > *wB -(khOfst); outrow--)
        for (outcol = *wL; outcol <= *wR; outcol++)
        {
            OutcolOutrow = outcol + (sw * outrow);
            result[OutcolOutrow] = 0;
        }
    for (outcol =  *wL; outcol < *wL + (kwOfst); outcol++)
        for (outrow = *wT; outrow <= *wB; outrow++)
        {
            OutcolOutrow = outcol + (sw * outrow);
            result[OutcolOutrow] = 0;
        }
    for (outcol = *wR; outcol > *wR - (kwOfst); outcol--)
        for (outrow = *wT; outrow <= *wB; outrow++)
        {
            OutcolOutrow = outcol + (sw * outrow);
            result[OutcolOutrow] = 0;
        }

    int32_t my_sum[12];
    int wT_loc = *wT;
    int wB_loc = *wB;
    int wL_loc = *wL;
    int wR_loc = *wR;

    //slide the convolution kernel over the image (within defined image window) to get each pixel of result array
    for (outrow = *wT + ((kh - 1) / 2); outrow <= *wB - (kh - 1) / 2; outrow++)
    {
        for (outcol = *wL + (kh - 1) / 2; outcol <= *wR - (kw - 1) / 2; outcol++)
        {
            //compute current pixel
            //sum = 0;
#ifdef DEBUG
            if( (outrow == 2) && (outcol == 2) )//debug
            {//debug
                sstop = true; //debug
            }//debug
#endif
#if SLOW_MATH
            for (krow = 0; krow < kh; krow++)
                for (kcol = 0; kcol < kw; kcol++)
                {
                    inrow = outrow - ((kh - 1) / 2) + krow;
                    incol = outcol - ((kw - 1) / 2) + kcol;
                    IncolInrow = incol + (sw * inrow);
                    sum += kernel[kcol][krow] * source[IncolInrow];

                }
#else
            OutcolOutrow = outcol + (sw * outrow);
            IncolInrow = OutcolOutrow - 642;

            neon_conv(&my_sum[0], &kernel_neon[0][0],&source[IncolInrow]);
            sum = my_sum[0] + my_sum[1] + my_sum[2];

            //printf("1 my_sum[0]=%d,my_sum[1]=%d,my_sum[2]=%d,my_sum[3]=%d,neon_sum=%ld\n",
            //		my_sum[0], my_sum[1], my_sum[2], my_sum[3], sum);
#endif
            //temp = (*gain * sum / ksum + 1); //because gain is always setup as 1.0 in this example

            temp = (sum / ksum + 1);

            //if (temp > maxint) temp = maxint; //limiter
            //if (temp < -maxint) temp = -maxint;//limiter
            //OutcolOutrow = outcol + (sw * outrow);
            result[OutcolOutrow] = (int)(temp); //normalize & round to nearest integer
            //printf("IncolInrow=%d,OutcolOutrow=%d,sub=%d\n",IncolInrow,OutcolOutrow, (OutcolOutrow-IncolInrow));
        }
    }
   clock_end();
 }

	void ImageProcessingUtilities_FPGA::GrayHist( int *w, int *h, int *wL, int *wR, int*wT, int *wB, int gray[], int grayhist[], int *grayhist_size )
    { //  int[] GrayHist(int[,] gray)
		int ij;
		clock_start("ImageProcessingUtilities_FPGA::GrayHist");

        for (int k = 0; k < *grayhist_size; k++) {
        	grayhist[k] = 0;
        }

        for (int i = *wL; i <= *wR; i++) {
            for (int j = *wT; j <= *wB; j+=16)
            {
				ij = i + (*w * j);
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+1;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+2;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+3;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+4;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+5;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+6;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+7;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+8;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+9;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+10;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+11;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+12;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+13;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+14;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);

				ij = i + (*w * j)+15;
                if( gray[ij] < *grayhist_size )
					++(grayhist[gray[ij]]);
				else ++(grayhist[ *grayhist_size - 1 ]);


            }
        }
        clock_end();
    }

	void ImageProcessingUtilities_FPGA::ThresholdArray(int Array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *Thresh, int tArray[])
    {
        //Output is tArray
		int ij;
		clock_start("ImageProcessingUtilities_FPGA::ThresholdArray");
        for (int i = *wL; i <= *wR; i++)
            for (int j = *wT; j <= *wB; j++)
            {
				ij = i + (*w * j);
                if (Array[ij] < *Thresh) tArray[ij] = 0;
                else tArray[ij] = Array[ij];
            }
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::MakeEdgeBuf( int EdgeImageArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *indx, int *max_edge_buf_size, PNT EdgeBuf[] )
    {
        //Outputs are EdgeBuf and indx (indx is number of values actually in Buf, as opposed to alloted memory space)
		int ij;
        *indx = 0;
        for (int j = *wT; j <= *wB; j++)
		{
			if( *indx >= *max_edge_buf_size ) break; //don't overrun the buffer length
            for (int i = *wL; i <= *wR; i++)
            { //col number (hor. pixel addr) will incr first. Edge pnts in same row will be sequential in buf.
				ij = i + (*w * j);
                if (EdgeImageArray[ij] > 0)
                {
                    EdgeBuf[*indx].X = i;//column number
                    EdgeBuf[*indx].Y = j;//row number
                    (*indx)++;
					if( *indx >= *max_edge_buf_size ) break; //don't overrun the buffer length
                }
            }
		}
    }
	    void ImageProcessingUtilities_FPGA::MakeThreshBuf( int ThrshArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *indx, int *max_thrsh_buf_size, LBLD_AREA_BUF ThrshBuf[] )
    {
        //Outputs are EdgeBuf and indx (indx is number of values actually in Buf, as opposed to alloted memory space)
		int ij;
		clock_start("ImageProcessingUtilities_FPGA::MakeThreshBuf");
        *indx = 0;
        for (int j = *wT; j <= *wB; j++)
		{
			if( *indx >= *max_thrsh_buf_size ) break; //don't overrun the buffer length
            for (int i = *wL; i <= *wR; i++)
            { //col number (hor. pixel addr) will incr first. Edge pnts in same row will be sequential in buf.
				ij = i + (*w * j);
                if (ThrshArray[ij] > 0)
                {
                    ThrshBuf[*indx].X = i;//column number
                    ThrshBuf[*indx].Y = j;//row number
					ThrshBuf[*indx].label = 0; //initialize region label.
                    (*indx)++;
					if( *indx >= *max_thrsh_buf_size ) break; //don't overrun the buffer length
                }
            }
		}
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::MakeGradientArray( int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int grad_array[] )
    {
    	clock_start("ImageProcessingUtilities_FPGA::MakeGradientArray");
#if 0
        // Output is grad_array
		//H_array and V_array MUST have same dimensions
		int ij;
        int max_1byte_int = 32767;
        int t;



        for (int i = *wL; i <= *wR; i++)
        {
            for (int j = *wT; j <= *wB; j++)
            {
				ij = i + (*width * j);

                t = AbsInt(H_array[ij]) + AbsInt(V_array[ij]);
                if (t > max_1byte_int) {
                	grad_array[ij] = max_1byte_int;
                } else if (t < -max_1byte_int) {
                	grad_array[ij] = -max_1byte_int;
                }
            }
        }
#endif
        clock_end();
    }

	void ImageProcessingUtilities_FPGA::Make8LevelDirectionArray(int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir_array_8[], int *thresh)
    {
        // Output is dir_array_8
		//H_array and V_array MUST have same dimensions
		int ij;
        int dir = 0;
        int grad = 0;
        int loop_count=0;
    	float temp;
    	clock_start("ImageProcessingUtilities_FPGA::Make8LevelDirectionArray");

        for (int i = *wL; i <= *wR; i++) {
            for (int j = *wT; j <= *wB; j++)
            {
            	//loop_count++;
				ij = i + (*width * j);
				grad = AbsInt(H_array[ij]) + AbsInt(V_array[ij]);

                if (grad <= *thresh)
                    dir = 0;
                else
                {
                    if ((H_array[ij] != 0) && (V_array[ij] != 0))// if neither V or H are zero
                    {
                        temp = (float)(V_array[ij]) / (float)(H_array[ij]);

                        if (H_array[ij] > 0) //H is positive
                        {
#ifdef SLOW_MATH
                            if ((temp > -0.414) && (temp <= 0.414)) dir = 1;
                            if ((temp > 0.414) && (temp <= 2.41)) dir = 2;
                            if ((temp > 2.41)) dir = 3;
                            if ((temp <= -2.41)) dir = 7;
                            if ((temp > -2.41) && (temp <= -0.414)) dir = 8;
#endif

                            if ((temp - -0.414) <= 0.414) {
                            	dir = 1;
                            } else if ((temp - 0.414) <= 2.41) {
                            	dir = 2;
                            } else if ((temp > 2.41)) {
                            	dir = 3;
                            } else if ((temp <= -2.41)) {
                            	dir = 7;
                            } else if ((temp - -2.41 <= -0.414)) {
                            	dir = 8;
                            }
                        }
                        else //H is negative
                        {
#ifdef SLOW_MATH
                            if ((temp <= -2.41)) dir = 3;
                            if ((temp > -2.41) && (temp <= -0.414)) dir = 4;
                            if ((temp > -0.414) && (temp <= 0.414)) dir = 5;
                            if ((temp > 0.414) && (temp <= 2.41)) dir = 6;
                            if ((temp > 2.41)) dir = 7;
#endif
                            if ((temp <= -2.41)) {
                            	dir = 3;
                            } else if ((temp - -2.41) <= -0.414) {
                            	dir = 4;
                            } else if ((temp - -0.414) <= 0.414) {
                            	dir = 5;
                            } else if ((temp - 0.414) <= 2.41) {
                            	dir = 6;
                            } else if ((temp > 2.41)) {
                            	dir = 7;
                            }
                        }
                    }
                    else //either H or V array is zero
                    {
                        if (H_array[ij] == 0)
                            if (V_array[ij] > 0) dir = 3;
                            else dir = 7;
                        if (V_array[ij] == 0)
                            if (H_array[ij] > 0) dir = 1;
                            else dir = 5;
                    }
                }
                dir_array_8[ij] = dir;
            }
        }
        //printf("loop_count=%d, minval=%f, maxval=%f, minval_az=%f, minval_bz=%f\n", loop_count, minval, maxval, minval_az, minval_bz);
        clock_end();
    }


	void ImageProcessingUtilities_FPGA::LoopStuffing(int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir_array_8[], int *thresh,
			int grad_array[],
			int dir_array_400[],
			int *Thresh,
			int thresh_grad_array[],
			int ThinGradArray[]
			)
    {
        // Output is dir_array_8
		//H_array and V_array MUST have same dimensions
		int ij;
        int dir = 0;
        int grad = 0;
        int loop_count=0;
        int max_1byte_int = 32767;

    	float temp;
    	clock_start("LoopStuffing");
    	clock_start("MakeGradientArray");
    	clock_start("Make8LevelDirectionArray");
    	clock_start("Make_0_to_400_DirectionArray");
    	clock_start("ThresholdArray");

        for (int i = *wL; i <= *wR; i++)
        {
            for (int j = *wT; j <= *wB; j++)
            {
                int H_arrayij = H_array[ij]; //do it once
                int V_arrayij = V_array[ij]; //do it once

				ij = i + (*width * j);

                //Testing loop stuffing


				//MakeGradientArray
                grad = AbsInt(H_arrayij) + AbsInt(V_arrayij);
                if (grad > max_1byte_int) {
                	grad = max_1byte_int;
                } else if (grad < -max_1byte_int) {
                	grad = -max_1byte_int;
                }
                grad_array[ij] = grad;
				//end MakeGradientArray

                //ThresholdArray
                if (grad < *Thresh) thresh_grad_array[ij] = 0;
                else thresh_grad_array[ij] = grad;
                //end  ThresholdArray

                //Make8LevelDirectionArray
                if (grad <= *thresh)
                    dir = 0;
                else
                {
                    if ((H_arrayij != 0) && (V_arrayij != 0))// if neither V or H are zero
                    {
                        temp = (float)(V_arrayij) / (float)(H_arrayij);

                        if (H_arrayij > 0) //H is positive
                        {
#ifdef SLOW_MATH
                            if ((temp > -0.414) && (temp <= 0.414)) dir = 1;
                            if ((temp > 0.414) && (temp <= 2.41)) dir = 2;
                            if ((temp > 2.41)) dir = 3;
                            if ((temp <= -2.41)) dir = 7;
                            if ((temp > -2.41) && (temp <= -0.414)) dir = 8;
#endif

                            if ((temp - -0.414) <= 0.414) {
                            	dir = 1;
                            } else if ((temp - 0.414) <= 2.41) {
                            	dir = 2;
                            } else if ((temp > 2.41)) {
                            	dir = 3;
                            } else if ((temp <= -2.41)) {
                            	dir = 7;
                            } else if ((temp - -2.41 <= -0.414)) {
                            	dir = 8;
                            }
                        }
                        else //H is negative
                        {
#ifdef SLOW_MATH
                            if ((temp <= -2.41)) dir = 3;
                            if ((temp > -2.41) && (temp <= -0.414)) dir = 4;
                            if ((temp > -0.414) && (temp <= 0.414)) dir = 5;
                            if ((temp > 0.414) && (temp <= 2.41)) dir = 6;
                            if ((temp > 2.41)) dir = 7;
#endif
                            if ((temp <= -2.41)) {
                            	dir = 3;
                            } else if ((temp - -2.41) <= -0.414) {
                            	dir = 4;
                            } else if ((temp - -0.414) <= 0.414) {
                            	dir = 5;
                            } else if ((temp - 0.414) <= 2.41) {
                            	dir = 6;
                            } else if ((temp > 2.41)) {
                            	dir = 7;
                            }
                        }
                    }
                    else //either H or V array is zero
                    {
                        if (H_arrayij == 0)
                            if (V_arrayij > 0) dir = 3;
                            else dir = 7;
                        if (V_arrayij == 0)
                            if (H_arrayij > 0) dir = 1;
                            else dir = 5;
                    }
                }
                dir_array_8[ij] = dir;
                //end Make8LevelDirectionArray



                //Make_0_to_400_DirectionArray
                if ((H_arrayij != 0) && (V_arrayij != 0))
                {
                    if ((V_arrayij > 0) && (H_arrayij > 0) && (V_arrayij >= H_arrayij)) //left top (use H/V)
                    {
                        dir_array_400[ij] = (int)(100.0f * (float)(H_arrayij) / (float)(V_arrayij));
                    }
                    else if ((V_arrayij > 0) && (H_arrayij > 0) && (V_arrayij < H_arrayij)) //left above center (use 2 - V/H)
                    {
                        dir_array_400[ij] = (int)(100.0f * (2.0f - (float)(V_arrayij) / (float)(H_arrayij)));
                    }
                    else if ((V_arrayij < 0) && (H_arrayij > 0) && (-V_arrayij <= H_arrayij)) //left below center (use 2 - V/H)
                    {
                        dir_array_400[ij] = (int)(100.0f * (2.0f - (float)(V_arrayij) / (float)(H_arrayij)));
                    }
                    else if ((V_arrayij < 0) && (H_arrayij > 0) && (-V_arrayij > H_arrayij)) //left bottom (use 4 + H/V)
                    {
                        dir_array_400[ij] = (int)(100.0f * (4.0f + (float)(H_arrayij) / (float)(V_arrayij)));
                    }
                    else if ((V_arrayij > 0) && (H_arrayij < 0) && (V_arrayij >= -H_arrayij)) //right top (use  H/V)
                    {
                        dir_array_400[ij] = (int)(100.0f * (float)(H_arrayij) / (float)(V_arrayij));
                    }
                    else if ((V_arrayij > 0) && (H_arrayij < 0) && (V_arrayij < -H_arrayij)) //right above center (use -2 - V/H)
                    {
                        dir_array_400[ij] = (int)(-100.0f * (2.0f + (float)(V_arrayij) / (float)(H_arrayij)));
                    }
                    else if ((V_arrayij < 0) && (H_arrayij < 0) && (V_arrayij >= H_arrayij)) //right below center (use -2 - V/H)
                    {
                        dir_array_400[ij] = (int)(-100.0f * (2.0f + (float)(V_arrayij) / (float)(H_arrayij)));
                    }
                    else if ((V_arrayij < 0) && (H_arrayij < 0) && (V_arrayij < H_arrayij)) //right bottom (use -4 + H/V)
                    {
                        dir_array_400[ij] = (int)(-100.0f * (4.0f - (float)(H_arrayij) / (float)(V_arrayij)));
                    }
                }
                else // H, V, or both are = 0
                {
                    if ((H_arrayij == 0) && (V_arrayij == 0)) dir_array_400[ij] = 0;
                    else
                    {
                        if (H_arrayij == 0)
                            if (V_arrayij > 0) dir_array_400[ij] = 0;
                            else dir_array_400[ij] = 400;
                        else
                            if (H_arrayij > 0) dir_array_400[ij] = 200;
                            else dir_array_400[ij] = -200;
                    }
                }
                //end Make_0_to_400_DirectionArray


            }
        }
        //printf("loop_count=%d, minval=%f, maxval=%f, minval_az=%f, minval_bz=%f\n", loop_count, minval, maxval, minval_az, minval_bz);
        clock_end();
        clock_end();
        clock_end();
        clock_end();
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::Make_0_to_400_DirectionArray( int H[], int V[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir[] )
    {
    	clock_start("ImageProcessingUtilities_FPGA::Make_0_to_400_DirectionArray");
#if 0
        //Output is dir

		//Bright circular object will be 0 dir at top and + or - 400 at bottom
        //To avoid trig functions, use H/V when V larger and V/H when H larger (region where tan doesn't deviate too much from a linear function).
        //As we go around the left edge of the circle starting from the top,
        //from top of bright circle to 45 deg around left edge, use H/V; from 45 to 135 use 2-V/H; from 135 to 180 use 4+H/V
        //Multiply by 100 to convert range to 0 - 400 (instead of 0 - 4)
        //As we go around the right edge starting from the top, values are negative of those described above.
        int grad = 0;
		int ij;


        for (int i = *wL; i <= *wR; i++)
        {
            for (int j = *wT; j <= *wB; j++)
            {
				ij = i + (*width * j);
                if ((H[ij] != 0) && (V[ij] != 0))
                {
                    if ((V[ij] > 0) && (H[ij] > 0) && (V[ij] >= H[ij])) //left top (use H/V)
                    {
                        dir[ij] = (int)(100.0f * (float)(H[ij]) / (float)(V[ij]));
                    }
                    else if ((V[ij] > 0) && (H[ij] > 0) && (V[ij] < H[ij])) //left above center (use 2 - V/H)
                    {
                        dir[ij] = (int)(100.0f * (2.0f - (float)(V[ij]) / (float)(H[ij])));
                    }
                    else if ((V[ij] < 0) && (H[ij] > 0) && (-V[ij] <= H[ij])) //left below center (use 2 - V/H)
                    {
                        dir[ij] = (int)(100.0f * (2.0f - (float)(V[ij]) / (float)(H[ij])));
                    }
                    else if ((V[ij] < 0) && (H[ij] > 0) && (-V[ij] > H[ij])) //left bottom (use 4 + H/V)
                    {
                        dir[ij] = (int)(100.0f * (4.0f + (float)(H[ij]) / (float)(V[ij])));
                    }
                    else if ((V[ij] > 0) && (H[ij] < 0) && (V[ij] >= -H[ij])) //right top (use  H/V)
                    {
                        dir[ij] = (int)(100.0f * (float)(H[ij]) / (float)(V[ij]));
                    }
                    else if ((V[ij] > 0) && (H[ij] < 0) && (V[ij] < -H[ij])) //right above center (use -2 - V/H)
                    {
                        dir[ij] = (int)(-100.0f * (2.0f + (float)(V[ij]) / (float)(H[ij])));
                    }
                    else if ((V[ij] < 0) && (H[ij] < 0) && (V[ij] >= H[ij])) //right below center (use -2 - V/H)
                    {
                        dir[ij] = (int)(-100.0f * (2.0f + (float)(V[ij]) / (float)(H[ij])));
                    }
                    else if ((V[ij] < 0) && (H[ij] < 0) && (V[ij] < H[ij])) //right bottom (use -4 + H/V)
                    {
                        dir[ij] = (int)(-100.0f * (4.0f - (float)(H[ij]) / (float)(V[ij])));
                    }
                }
                else // H, V, or both are = 0
                {
                    if ((H[ij] == 0) && (V[ij] == 0)) dir[ij] = 0;
                    else
                    {
                        if (H[ij] == 0)
                            if (V[ij] > 0) dir[ij] = 0;
                            else dir[ij] = 400;
                        else
                            if (H[ij] > 0) dir[ij] = 200;
                            else dir[ij] = -200;
                    }
                }
            }
        }
#endif
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::NonMaximalSuppression(int grad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int dir_array[], int ThinGradArray[])
    {


        // Output is ThinGradArray
		// dir_array MUST have same dimensions as grad_array
		int i, j, ij, i1j1, i2j2, grad, i1, i2, j1, j2;

		clock_start("ImageProcessingUtilities_FPGA::NonMaximalSuppression");
        for ( i = *wL; i <= *wR; i++)
            for ( j = *wT; j <= *wB; j++)
            {
				ij = i + (*w * j);
                if ((i == 1) || (j == 1) || (i == *w-1) || (j == *h-1)) ThinGradArray[ij] = 0;
				if ((i == 0) || (j == 0) || (i == *w-2) || (j == *h-2)) ThinGradArray[ij] = 0;
                else
                {
                    grad = grad_array[ij];
                    i1 = 0, i2 = 0, j1 = 0, j2 = 0;
                    switch (dir_array[ij])
                    {
                        case 1: // 0 deg (gradient points right)
                            {
                                i1 = i - 1; //left
                                j1 = j;
                                i2 = i + 1; //right
                                j2 = j;
                                break;
                            }
                        case 5: // 180 deg (gradient points lft)
                            {
                                i2 = i - 1; //left
                                j2 = j;
                                i1 = i + 1; //right
                                j1 = j;
                                break;
                            }
                        case 2: // 45 deg (gradient points lower right)
                            {
                                i2 = i + 1; // right
                                j2 = j + 1; // down
                                i1 = i - 1; // left
                                j1 = j - 1; // up
                                break;
                            }
                        case 6: // -135 deg (gradient points upper left)
                            {
                                i1 = i + 1; // right
                                j1 = j + 1; // down
                                i2 = i - 1; // left
                                j2 = j - 1; // up
                                break;
                            }
                        case 3: // 90 deg (gradient points down)
                            {
                                i2 = i;
                                j2 = j + 1; // down
                                i1 = i;
                                j1 = j - 1; // up
                                break;
                            }
                        case 7: // -90 deg (gradient points up)
                            {
                                i1 = i;
                                j1 = j + 1; // down
                                i2 = i;
                                j2 = j - 1; // up
                                break;
                            }
                        case 4: // 135 deg (gradient points lower left)
                            {
                                i2 = i - 1; // left
                                j2 = j + 1; // down
                                i1 = i + 1; // right
                                j1 = j - 1; // up
                                break;
                            }
                        case 8: // -45 deg (gradient points upper right)
                            {
                                i1 = i - 1; // left
                                j1 = j + 1; // down
                                i2 = i + 1; // right
                                j2 = j - 1; // up
                                break;
                            }
							
                    }
					i1j1 = i1 + (*w * j1);
					i2j2 = i2 + (*w * j2);
                    if ((grad_array[i1j1] > grad) || (grad_array[i2j2] >= grad)) ThinGradArray[ij] = 0;
                    else ThinGradArray[ij] = grad_array[ij];
                }
            }
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::GradHist( int grad[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *hdim,  int hist[] )
    {
        //Output hist
		//int localhist[32768];//debug
		int ij;
		bool sstop = false; //debug
		clock_start("ImageProcessingUtilities_FPGA::GradHist");
        for (int k = 0; k < *hdim; k++) //initialize to zero
		{
			hist[k] = 0;
		}
        for (int i = *wL; i <= *wR; i++) {
            for (int j = *wT; j <= *wB; j+=16)
            {
				ij = i + (*w * j);
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+1;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+2;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+3;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+4;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+5;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+6;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+7;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

                //next 8
				ij = i + (*w * j)+8;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+9;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+10;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+11;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+12;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+13;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+14;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

				ij = i + (*w * j)+15;
                if( grad[ij] < *hdim ) {
					++( hist[grad[ij]] );
				} else {
					++( hist[ *hdim - 1 ] );
				}

            }
        }
        clock_end();
    }

	void ImageProcessingUtilities_FPGA::MakeDualThresholdArray(int ThinGrad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *Thigh, int *Tlow, int DualThresh_array[] )
    {
		int ij;
        for (int i = *wL; i <= *wR; i++)
            for (int j = *wT; j < *wB; j++)
            {
				ij = i + (*w * j);
                if (ThinGrad_array[ij] < *Tlow) DualThresh_array[ij] = 0;
                else
                    if (ThinGrad_array[ij] < *Thigh) DualThresh_array[ij] = 1;
                    else DualThresh_array[ij] = 2;
            }
    }

	void ImageProcessingUtilities_FPGA::DualThreshArrayWithThinLineSuppression(int ThinGrad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, 
		                                                                          int DirArray8[], int *Thigh, int *Tlow, int DualThresh_array[], int *hairwidth )
    {
		int ij;
		bool sstop = false; //debug
		clock_start("ImageProcessingUtilities_FPGA::DualThreshArrayWithThinLineSuppression");
        for (int j = *wT; j <= *wB; j++)
            for (int i = *wL; i <= *wR; i++)
            {
				ij = i + (*w * j);
				if( (i == 4) && (j == 3) )
					sstop = true;
				//if ( (i < 80) || (i > 240) || (j < 60) || (j > 160) ) DualThresh_array[ij] = 0;//debug
				//else//debug
                if (ThinGrad_array[ij] < *Tlow) DualThresh_array[ij] = 0;
                else
					if( !SuppressThinDarkLine(ThinGrad_array, w, h, wL, wR, wT, wB, DirArray8, &i, &j, hairwidth) )
					{
						if (ThinGrad_array[ij] < *Thigh) DualThresh_array[ij] = 1;
						else DualThresh_array[ij] = 2;
					}
					else DualThresh_array[ij] = 0;
            }
        clock_end();
    }

	bool ImageProcessingUtilities_FPGA::SuppressThinDarkLine( int Grad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int DirArray8[], int *i, int *j, int *hairwidth )
	{
		//look in bright to dark grad direction for an edge point with opposite grad direction within hiarwidth pixels away.  
        //If such an edge point exists report suppress==true.  Note that for an FPGA to do this it needs acess to (2*hirwidth +1) lines.

        //Alternate strategy is to just look along current line.  If bright to dark grad points upper right, right or lower right, look right. 
        //If there is a point within "hairwidth" in this direction with grad either upper left, left or lower left, set suppress==true.
        //If bright to dark grad points upper left, left or lower left, look left and set suppress==true if point found with brt to drk grad pointing up-rt,rt, or low-rt.
        //(If grad points directly up or down, don't do anything. Leave suppress==false).  An FPGA would need access only to 2*hairwidth elements on current line.  
		int ij, iijj, ii, /*jj,*/ n;
		bool suppress = false;
		ij = *i + (*w * *j);

        switch( DirArray8[ij])
        {
            case 4:
            case 5:
            case 6: //dark area to right 
                ii = *i;
                for (n = 0; n < *hairwidth; n++)
                {
                    ++ii;
					iijj = ii + (*w * *j);
                    if ( ii > *wR ) break;
                    if (Grad_array[iijj] > 0)
                        if ((DirArray8[iijj] == 2) || (DirArray8[iijj] == 8) || (DirArray8[iijj] == 1)) //dark to left
                        {
                            suppress = true;
                            break;
                        }
                }
                break;
            case 8:
            case 1:
            case 2: //dark to left
                ii = *i;
                for (n = 0; n < *hairwidth; n++)
                {
                    --ii;
					iijj = ii + (*w * *j);
                    if ( ii < *wL ) break;
                    if (Grad_array[iijj] > 0)
                        if ((DirArray8[iijj] > 4) || (DirArray8[iijj] <= 6)) //dark to right
                        {
                            suppress = true;
                            break;
                        }
                }
                break;
        }

		//switch( DirArray8[ij] )
		//{
		//	case 1: // look up
		//		jj = *j;
		//		for( n = 0; n < *hairwidth; n++ )
		//		{
		//			--jj;
		//			iijj = *i + (*w * jj);
		//			if ( DirArray8[iijj] == 5 )//bright to dark direction is down
		//			{
		//				Grad_array[iijj] = 0;
		//				DirArray8[iijj] = 0;
		//				Grad_array[ij] = 0;
		//				DirArray8[ij] = 0;
		//				suppress = true;
		//				break;
		//			}
		//		}
		//		break;
		//	case 2: // look up to left
		//		ii = *i;
		//		jj = *j;
		//		for( n = 0; n < *hairwidth; n++ )
		//		{
		//			--ii;
		//			--jj;
		//			iijj = ii + (*w * jj);
		//			if ( DirArray8[iijj] == 6 )//bright to dark direction is diagonal pointing down right
		//			{
		//				Grad_array[iijj] = 0;
		//				DirArray8[iijj] = 0;
		//				Grad_array[ij] = 0;
		//				DirArray8[ij] = 0;
		//				suppress = true;
		//				break;
		//			}
		//		}
		//		break;
		//	case 3: // look to left
		//		ii = *i;
		//		for( n = 0; n < *hairwidth; n++ )
		//		{
		//			--ii;
		//			iijj = ii + (*w * *j);
		//			if ( DirArray8[iijj] == 7 )//bright to dark direction is right
		//			{
		//				Grad_array[iijj] = 0;
		//				DirArray8[iijj] = 0;
		//				Grad_array[ij] = 0;
		//				DirArray8[ij] = 0;
		//				suppress = true;
		//				break;
		//			}
		//		}
		//}
		return (suppress);
	}

	

    void ImageProcessingUtilities_FPGA::MakeDVBuf( int DedgeArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int DirArray[], int *maxBufSize, int *dvbuf_size, DVBUF dvbuf[] )
    {
        //Outputs are dvbuf and dvbuf_size (dbvug size is actural number of elements put in dvbuf, rather than memory allocated)
		int ij;
		int k = 0;
		clock_start("ImageProcessingUtilities_FPGA::MakeDVBuf");
        for (int j = *wT; j <= *wB; j++ )
		{
            for (int i = *wL; i <= *wR; i++)
            {
				ij = i + (*w * j);
                if( DedgeArray[ij] > 0 )
                {
                    dvbuf[k].X = i;
                    dvbuf[k].Y = j;
                    dvbuf[k].value = DedgeArray[ij];
                    dvbuf[k].dir = DirArray[ij];
                    k++;
                }
				if( k >= *maxBufSize ) 
					break;
            }
			if( k >= *maxBufSize ) 
				break;
		}
        *dvbuf_size = k;
        clock_end();
    }

    void ImageProcessingUtilities_FPGA::CopyDvbuf( DVBUF dvbuf[], int *dvbuf_size, DVBUF dvbuf2[] )
    {
    	clock_start("ImageProcessingUtilities_FPGA::CopyDvbuf");
        //Output is dvbuf2
		for (int i = 0; i < *dvbuf_size; i++)
        {
            dvbuf2[i].X = dvbuf[i].X;
            dvbuf2[i].Y = dvbuf[i].Y;
            dvbuf2[i].value = dvbuf[i].value;
            dvbuf2[i].dir = dvbuf[i].dir;
        }
		clock_end();
    }

    void ImageProcessingUtilities_FPGA::MakeDirBuf( PNT EdgeBuf[], int *size, int DirImageArray[], int DirBuf[], int *Isizeh, int *Isizev )
    {
        //Output DirArray
		//DirArray is int[] array with elements corresponding to each element in EdgeBuf.
        //Each DirArray element contains a direction value (from DirImageArray) for the edge point 
        // specified by the corresponding element of EdgeBuf.
		int ij;
        for (int i = 0; i < *size; i++)
        {
			ij = EdgeBuf[i].X + (*Isizeh * EdgeBuf[i].Y);
            DirBuf[i] = DirImageArray[ij];
        }
    }

	int ImageProcessingUtilities_FPGA::AbsInt(int I)
	{
		if(I<0) return(-I);
		else return(I);
	}
















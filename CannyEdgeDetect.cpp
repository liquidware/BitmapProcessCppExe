#include <stdio.h>
#include "Clock.h"
#include "ProcessImage.h"
#include "CannyEdgeDetect.h"
#include "ImageProcessingUtilities_DSP.h"
#include "ImageProcessingUtilities_FPGA.h"

CannyEdgeDetect::CannyEdgeDetect(void)
{
}

CannyEdgeDetect::~CannyEdgeDetect(void)
{
}

void CannyEdgeDetect::PCrEdgeDetect
	(int *Isizeh, int *Isizev, int Iarray[], EDGE_DETECTION_PARAM_FPGA *pfpga, EDGE_DETECTION_PARAM_DSP *pdsp, OBJECT_RECOGNITION_PARAMETERS *orp_CR, int Skern[5][5], int *SkernSize, int GkernH[5][5], int GkernV[5][5],
	int *GkernSize, PNT Pbuf[],  PNT CRbuf[], int PDirBuf[], int CrDirBuf[], int *pbuf_size, int *CRbuf_size, int grad_H_array[], int grad_V_array[], int CRarray[], LBLD_AREA_BUF CRthreshBuf[], int *cr_threshbuf_size,
	CANNY_EDGE_DETECT_DATA* data)
//Inputs:
//	Isizeh;				horiz eye image size (pixels)
//  Isizev;				vert eye image size (pixels)
//  Iarray;				Eye cam image in form of a doubly subscripted array.  The subscripts are col and row addresses and each element value is the gray scale level. 
//	pfpga;				edge detection parameters that will be used by FPGA
//  pdsp;				edge detection parameters that will be used by DSP
//  Skern;				Kernel for smoothing convolution
//  SkernSize;          dimension of Kernel (single value since kernel array is always square)
//  Smooth;				Switch to use (true) or not use (false) smoothing convolution
//  GkernH;				Kernel for horizontal gradient convolution 
//  GkernV;             KErnel for vertical gradiet convolution
//  GkernSize;          dimension of GkernH and GkernV (both are same size and are square)
//  pupil_rec_param;	image processing parameters used by DPS to find pupil
//  CR_rec_param;		image processing parameters used by DPS to find CR
//
//Ouputs:
//  Pbuf;			List of possible pupil edge points.  Each list element has h (col) and v(row) pixel address 
//  CRbuf;			List of possible CR edge points.
//  PDirBuf;		List of elements corresponding to Pbuf.  Each list element has a direction value
//  CrDirBuf;		List of elements corresponding to CRbuf.  Each list element has a direction value
//  Pbuf_size;		No of list elements in Pbuf and PDirBuf.
//  CRbuf_size;		No of list elements in CRbuf and CrDirBuf.
//  grad_H_array;	Although used internally, this is passed out of routine only for debugging help (Object Recognition Routine will use it to create some arrays that are only for debugging)
//  grad_V_array;	Although used internally, this is passed out of routine only for debugging help (Object Recognition Routine will use it to create some arrays that are only for debugging)
//  CRarray;        original gray scale image or smoothed image, with all values below CR threshold zeroed.  
//	LBLD_AREA_BUF CRthreshBuf[]; List of points with non zero values in CRarray (points that were above CR threshold value). 
//	int *cr_threshbuf_size;		 number of elements in CRthreshBuf[] 
//  CANNY_EDGE_DETECT_DATA* data	pointer to memory storage for temp vars

{
	clock_start("CannyEdgeDetect::PCrEdgeDetect");
//int dummy;
//LARGE_INTEGER t1, t2, t3;
//QueryPerformanceFrequency(&t1);
//QueryPerformanceCounter(&t2);

	int CR_featuretype = BRIGHT;
	int IsizehIsizev = *Isizeh * *Isizev;

	//int DirBuf8[800000];
	//int Ismooth_array[800000];
	//int grad_array[800000];
	//int dir_array_8[800000];
	//int dir_array_400[800000];
	//int thresh_grad_array[800000];
	//int thin_array[800000];
	////int thin_CR[800000];
	//int dual_thresh_array[800000];
	//int grayhist[256];
	//int ghist[32767];
	//int *DirBuf8 = new int[IsizehIsizev];
	//int *Ismooth_array = new int[IsizehIsizev];
	//int *grad_array = new int[IsizehIsizev];
	//int *dir_array_8 = new int[IsizehIsizev];
	//int *dir_array_400 = new int[IsizehIsizev];
	//int *thresh_grad_array = new int[IsizehIsizev];
	//int *thin_array = new int[IsizehIsizev];
	//int *thin_CR = new int[IsizehIsizev];
	//int *dual_thresh_array = new int[IsizehIsizev];
	//int *grayhist = new int[pfpga->grayhist_size];
	//int *ghist = new int[pfpga->grad_hist_size];

	int max_grad, gh_cmx, gh_xpeak, gh_ypeak, gh_bndry, Thigh, Tlow;
	int CRthresh, CRthresh_gray2, maxGray, crh_cmx, crh_xpeak, crh_ypeak, brt_pup_bndry,  dvbuf_size, wL, wR, wT, wB, wl, wr, wt, wb;

	//DVBUF dvbuf[100000];
	//DVBUF dvbuf2[100000];
    //DVBUF *dvbuf = new DVBUF[pfpga->max_pbuf_size]; 
    //DVBUF *dvbuf2 = new DVBUF[pfpga->max_pbuf_size];

//QueryPerformanceCounter(&t3); //timer end point
//double dt1 = t1.QuadPart; //timer 
//double dt2 = t2.QuadPart;//timer
//double dt3 = t3.QuadPart;//timer
//double deltat = (dt3 - dt2) / dt1;//timer


	//image boundaries for smoothing convolution
	wl = pfpga->wL;
	wr = pfpga->wR;
	wt = pfpga->wT;
	wb = pfpga->wB;

	//if use smoothing convolution, grad conv and subsequent operation must assume additional (k-1)/2 line boundary, where k is order of smoothing conv. 
	if( pfpga->smooth )
	{
		wL = pfpga->wL + pfpga->extra_bndry;
		wR = pfpga->wR - pfpga->extra_bndry;
		wT = pfpga->wT + pfpga->extra_bndry;
		wB = pfpga->wB - pfpga->extra_bndry;
	}
	else
	{
		wL = pfpga->wL;
		wR = pfpga->wR;
		wT = pfpga->wT;
		wB = pfpga->wB;
	}

	//int dummy;
	//LARGE_INTEGER t1, t2, t3;
	//QueryPerformanceFrequency(&t1);

//QueryPerformanceCounter(&t2);

    //***** EDGE DETECTION -- FPGA tasks  ***********//

	if( pfpga->smooth )//if we are going to do smoothing convolution
	{
		ImageProcessingUtilities_FPGA::GrayscaleConvolution(Isizeh, Isizev, &wl, &wr, &wt, &wb, Iarray, SkernSize, Skern, &pfpga->SmthCnvl_Scale, data->Ismooth_array);
		ImageProcessingUtilities_FPGA::GrayHist(Isizeh, Isizev, &wl, &wr, &wt, &wb, data->Ismooth_array, data->grayhist, &pfpga->grayhist_size);//grayhist is sent by the FPGA to the DSP
	} else {
		ImageProcessingUtilities_FPGA::GrayHist(Isizeh, Isizev, &wL, &wR, &wT, &wB, Iarray, data->grayhist, &pfpga->grayhist_size);
	}

	//***********Computation that will be done by DSP for use in next field *****************************
	//cmx is center of mass; ypeak is number of pixels in hist peak; xpeak is brightness value at hist peak; maxGray is brightest pixel.
	//brt_pup_bndry will be lowest intensity with at least min_pup_pix pixels (bright pupil should have at least min_pup_pix at its brightest intensity and most of bright pupil is probably just below this brightness boundary).
	ImageProcessingUtilities_DSP::ComputeThreshFromHist( data->grayhist, &pfpga->grayhist_size, &pdsp->CRthresh_hist_fraction, &pdsp->min_pup_pix, &maxGray, &crh_cmx, &crh_xpeak, &crh_ypeak, &brt_pup_bndry, &CRthresh); //note: only brt_pup_bndry and maxGray currently used from this routing.
    if (maxGray > pdsp->MaxGrayLimit) maxGray = pdsp->MaxGrayLimit;
    CRthresh_gray2 = brt_pup_bndry + (int)((double)(maxGray - brt_pup_bndry) * pdsp->frac_dist_brtpup_to_brightest_pixel__CRthresh + 0.5);
	//****************************************************************************************************
    
	//****** do a simple intensity threshold to find potential CR pixels ******//
	if( pfpga->smooth ) { //if we did smoothing, used smoothed image
		ImageProcessingUtilities_FPGA::ThresholdArray( data->Ismooth_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, &CRthresh_gray2, CRarray );
	} else { //otherwise, use original image
		ImageProcessingUtilities_FPGA::ThresholdArray( Iarray, Isizeh, Isizev, &wL, &wR, &wT, &wB, &CRthresh_gray2, CRarray );
	}
	//****** turn this into a buffer list to pass to the DSP
    ImageProcessingUtilities_FPGA::MakeThreshBuf( CRarray, Isizeh, Isizev, &wL, &wR, &wT, &wB, cr_threshbuf_size, &pfpga->max_cr_threshbuf_size, CRthreshBuf ); //make list of points (CRthreshBuf), defined by h and v coord, that are above CR intesity threshold

//QueryPerformanceCounter(&t2);
    //-- make gradient map and direction map
    //******* first compute horizontal and vertical gradients
	if( pfpga->smooth )//if we did smoothing, compute gradients using the smoothed image
	{
		ImageProcessingUtilities_FPGA::GrayscaleConvolution( Isizeh, Isizev, &wL, &wR, &wT, &wB, data->Ismooth_array, GkernSize, GkernH, &pfpga->GradCnvl_Scale, grad_H_array );
		ImageProcessingUtilities_FPGA::GrayscaleConvolution( Isizeh, Isizev, &wL, &wR, &wT, &wB, data->Ismooth_array, GkernSize, GkernV, &pfpga->GradCnvl_Scale, grad_V_array );
	}
	else //otherwise use original image
	{
		ImageProcessingUtilities_FPGA::GrayscaleConvolution( Isizeh, Isizev, &wL, &wR, &wT, &wB, Iarray, GkernSize, GkernH, &pfpga->GradCnvl_Scale, grad_H_array );
		ImageProcessingUtilities_FPGA::GrayscaleConvolution( Isizeh, Isizev, &wL, &wR, &wT, &wB, Iarray, GkernSize, GkernV, &pfpga->GradCnvl_Scale, grad_V_array );
	}
//QueryPerformanceCounter(&t2);
    //******* combine horizontal and vertical gradients to make gradient magnitude and gradient direction maps
    ImageProcessingUtilities_FPGA::LoopStuffing( grad_H_array, grad_V_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_8, &pfpga->grad_thresh_for_dir_array,
    		data->grad_array,
    		data->dir_array_400,
    		&pfpga->grad_thresh,
    		data->thresh_grad_array,
    		data->thin_array);
    //ImageProcessingUtilities_FPGA::MakeGradientArray( grad_H_array, grad_V_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->grad_array );
    //ImageProcessingUtilities_FPGA::Make8LevelDirectionArray( grad_H_array, grad_V_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_8, &pfpga->grad_thresh_for_dir_array);

    //ImageProcessingUtilities_FPGA::Make_0_to_400_DirectionArray( grad_H_array, grad_V_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_400 );
//QueryPerformanceCounter(&t2);
    //****** Threshold the gradient magnitude array (pixels with above thresh gradient retain unchanged grad values while those below are zeroed 
    //cladden testing loop stuffing
    //ImageProcessingUtilities_FPGA::ThresholdArray( data->grad_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, &pfpga->grad_thresh, data->thresh_grad_array ); //grad_thresh usually 1500

    //***** do non-maximal suppression to thin edges (zero gradient for pixels with higher gradient neighbor along gradient direction)  
    //int[,] thin8LDir_array;
    ImageProcessingUtilities_FPGA::NonMaximalSuppression( data->thresh_grad_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_8, data->thin_array );
//QueryPerformanceCounter(&t2);
    //***** Apply dual thresholds to the thin edge array -- we can then either do threshold hysteresis, or just use the lower threshold
    //******* First make a histogram of gradient values
    ImageProcessingUtilities_FPGA::GradHist( data->thin_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, &pfpga->grad_hist_size, data->ghist );

	//*** Determination of Thigh and Tlow will be done by DSP for next field ***//
    //use the grad histogram to determine dual threshold values Thigh and Tlow -- this will actually be done by the DSP for use with the subsequent video field
    ImageProcessingUtilities_DSP::ComputeThreshFromHist( data->ghist, &pfpga->grad_hist_size, &pdsp->high_thresh_fraction, &pdsp->min_pup_pix, &max_grad, &gh_cmx, &gh_xpeak, &gh_ypeak, &gh_bndry, &Thigh );
    Tlow = (int)(pdsp->low_thresh_fraction * (float)Thigh);
	//**************************************************************************//

    //******* Assign pixels with gradients above high thresh (Thigh) = 2, pixels above low thresh (Tlow) =1, others = 0 
    //ImageProcessingUtilities_FPGA::MakeDualThresholdArray( data->thin_array, Isizeh, Isizev, &Thigh, &Tlow, dual_thresh_array );
	ImageProcessingUtilities_FPGA::DualThreshArrayWithThinLineSuppression( data->thin_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_8, &Thigh, &Tlow, data->dual_thresh_array, &pfpga->hair_width );
	
	//Convert the dual value array into a buffer.  The buffer is a structure containing pixel h and v address, threshold value (1 or 2), and direction value for every pixel above lower thresh.
	//Note that if we are not going to do threshold hysteresis, this is a wasted step -- we could make the edge buff and direction buff directly using the lower threshold value.
    ImageProcessingUtilities_FPGA::MakeDVBuf( data->dual_thresh_array, Isizeh, Isizev, &wL, &wR, &wT, &wB, data->dir_array_400, &pfpga->max_pbuf_size, &dvbuf_size, data->dvbuf ); //dvbuf sent to DSP by FPGA
    ImageProcessingUtilities_FPGA::CopyDvbuf(data->dvbuf, &dvbuf_size, data->dvbuf2); //make a copy just so we can look at the original during debug (it will be modified by threshold hysteresis)
    
	//******* Threshold Hysteresis:  If a "1" pixel is next to a "2" pixel, make it a 2, NOTE: This will be done by the DSP if at all
    if (pfpga->use_thresh_hyst)
		ImageProcessingUtilities_DSP::DoThresholdHysteresis2( data->dvbuf, &dvbuf_size);//DSP function -- dvbuf will be sent to DSP by FPGA

    //*** Make edge buf -- array of PNT structures with each element contianing h and v (x and y) coordinates of one edge point ***
    ImageProcessingUtilities_DSP::MakeEdgeBuf_from_DVBUF( data->dvbuf, &dvbuf_size, &pfpga->use_thresh_hyst, &pfpga->max_pbuf_size, pbuf_size, Pbuf ); //DSP function -- use dvbuf to make Pbuf

    //*** Make pupil edge direction buffer (int array with each element specifying gradient direction of corresponding point in PupilEdgeBuf)
	// make high resolution direction buffer for use, by DSP, in feature recognition (see "ImageProcessingUtilities_FPGA::Make_0_to_400_DirectionArray" routine for explantion of direction values).
    ImageProcessingUtilities_FPGA::MakeDirBuf(Pbuf, pbuf_size, data->dir_array_400, PDirBuf, Isizeh, Isizev); //the FPGA will actually use the edge array (rather than edge buf) and dir_array to do this
    // make low reslution direction buffer (8 cardinal directions) that will be needed by DSP only if using edge data (as opposed to labeled region alaysis) to find CR
	if(!orp_CR->use_labeled_area_analysis_for_CR)//we don't need to do this if using labeled region analysis to find CR
		ImageProcessingUtilities_FPGA::MakeDirBuf(Pbuf, pbuf_size, data->dir_array_8, data->DirBuf8, Isizeh, Isizev); //the FPGA will actually use the edge array (rather than edge buf) and dir_array to do this

    //*** Make CR edge buf & dir buf -- edge buf is subset of pupil edge buf.  An edge point must have a pixel CRthresh_dist, in grad direction, that is in CRthresBuf (exceeds CRthresh) to be a CR edge point  ******************//
    if(!orp_CR->use_labeled_area_analysis_for_CR)//we don't need to do this if using labeled region analysis to find CR
		ImageProcessingUtilities_DSP::GrayAreaEdgePoints4( Pbuf, pbuf_size, data->DirBuf8, PDirBuf, Isizeh, Isizev, &wL, &wR, &wT, &wB, CRthreshBuf, cr_threshbuf_size, &pdsp->CRthresh_dist, &CR_featuretype, CRbuf_size, CRbuf, CrDirBuf );
	
//QueryPerformanceCounter(&t3); //timer end point
//double dt1 = t1.QuadPart; //timer 
//double dt2 = t2.QuadPart;//timer
//double dt3 = t3.QuadPart;//timer
//double deltat = (dt3 - dt2) / dt1;//timer
	//WCHAR str[80]; //report time in DBCON when running executable
	//swprintf_s(str, 80, L"deltat=%f", deltat); //report time in DBCON when running executable
	//OutputDebugString(str); //report time in DBCON when running executable

	//*** release temporary memory **********************
	//delete[] DirBuf8;
	//delete[] Ismooth_array;
	//delete[] grad_array;
	//delete[] dir_array_8;
	//delete[] dir_array_400;
	//delete[] thresh_grad_array;
	//delete[] thin_array;
	//delete[] thin_CR;
	//delete[] dual_thresh_array;
	//delete[] grayhist;
	//delete[] ghist;
	//delete[] dvbuf;
	//delete[] dvbuf2;
    clock_end();
}

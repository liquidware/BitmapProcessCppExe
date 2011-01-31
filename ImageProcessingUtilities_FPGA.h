#pragma once

class ImageProcessingUtilities_FPGA
{
public:
	ImageProcessingUtilities_FPGA(void);
	~ImageProcessingUtilities_FPGA(void);
	static void GrayscaleConvolution( int *psw, int *psh, int *wL, int *wR, int *wT, int *wB, int source[], int *ksize, int kernel[5][5], float *gain, int result[] );
	static void GrayHist( int *w, int *h, int *wL, int *wR, int *wT, int *wB, int gray[], int grayhist[], int *grayhist_size );
	static void ThresholdArray(int Array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *Thresh, int tArray[]);
	static void MakeEdgeBuf( int EdgeImageArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *indx, int *max_edge_buf_size, PNT EdgeBuf[] );
	static void MakeThreshBuf( int ThrshArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *indx, int *max_thrsh_buf_size, LBLD_AREA_BUF ThrshBuf[] );
	static void MakeGradientArray( int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int grad_array[] );
	static void Make8LevelDirectionArray(int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir_array_8[], int *thresh);
	static void LoopStuffing(int H_array[], int V_array[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir_array_8[], int *thresh,
			int grad_array[],
			int dir_array_400[],
			int *Thresh,
			int thresh_grad_array[],
			int ThinGradArray[]);
	static void Make_0_to_400_DirectionArray( int H[], int V[], int *width, int *height, int *wL, int *wR, int *wT, int *wB, int dir[] );
	static void NonMaximalSuppression( int grad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int dir_array[], int ThinGradArray[] );
	static void GradHist( int grad[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *hdim,  int hist[] );
	static void MakeDualThresholdArray(int ThinGrad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int *Thigh, int *Tlow, int DualThresh_array[] );
	static void MakeDVBuf( int DedgeArray[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int DirArray[], int *maxBufSize, int *dvbuf_size, DVBUF dvbuf[] );
	static void CopyDvbuf( DVBUF dvbuf[], int *dvbuf_size, DVBUF dvbuf2[] );
	static void MakeDirBuf( PNT EdgeBuf[], int *size, int DirImageArray[], int DirBuf[], int *Isizeh, int *Isizev );
	static int AbsInt(int I);
	static void DualThreshArrayWithThinLineSuppression(int ThinGrad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int DirArray8[], int *Thigh, int *Tlow, int DualThresh_array[], int *hairwidth );
	static bool SuppressThinDarkLine( int Grad_array[], int *w, int *h, int *wL, int *wR, int *wT, int *wB, int DirArray8[], int *i, int *j, int *hairwidth );
};


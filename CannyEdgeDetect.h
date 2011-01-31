#pragma once

class CannyEdgeDetect
{
public:
	CannyEdgeDetect(void);
	~CannyEdgeDetect(void);

	static void PCrEdgeDetect
	(int *Isizeh, int *Isizev, int Iarray[], EDGE_DETECTION_PARAM_FPGA *pfpga, EDGE_DETECTION_PARAM_DSP *pdsp, OBJECT_RECOGNITION_PARAMETERS *orp_CR, int Skern[5][5], int *SkernSize, int GkernH[5][5], int GkernV[5][5],
	int *GkernSize, PNT Pbuf[],  PNT CRbuf[], int PDirBuf[], int CrDirBuf[], int *pbuf_size, int *CRbuf_size, int grad_H_array[], int grad_V_array[], int CRarray[], LBLD_AREA_BUF CRthreshBuf[], int *cr_threshbuf_size,
	CANNY_EDGE_DETECT_DATA* data);
};

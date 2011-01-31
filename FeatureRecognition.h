#ifndef FeatureRecognition_h
#define FeatureRecognition_h

#pragma once
#include "ProcessImage.h"

class FeatureRecognition
{
public:
	FeatureRecognition(void);
	static void PCrRecognition
	( int *Isizeh, int *Isizev, EDGE_DETECTION_PARAM_FPGA *pfpga, EDGE_DETECTION_PARAM_DSP *pdsp, OBJECT_RECOGNITION_PARAMETERS *pupil_rec_param, //INPUTS
		OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, PNT Pbuf[], PNT CRbuf[], EDGE_POINT_STATUS edgepoint[], EDGE_POINT_STATUS CR_edgepoint[],
		int PDirBuf[], int CrDirBuf[], int *Pbuf_size, int *CRbuf_size,	FuzzParam *Pfuzz, int grad_H_array[], int grad_V_array[], int CRarray[],
		bool *pupil_found, bool *CR_found, EDGE_STRING es[], int *nstring,  EDGE_STRING CR_es[], int *nCRstring,                                  //OUTPUTS
		LFT_EDGE lft_edge[], int *no_lft_edges, LFT_EDGE lft_CR_edge[], int *no_lft_CR_edges, RT_EDGE rt_edge[], int *no_rt_edges, RT_EDGE rt_CR_edge[], int *no_rt_CR_edges,
		OBJECT objt[], int *no_of_objects, CR_OBJECT CR_objt[], int *no_of_CR_objects, int *pupil_object_no, int *CR_object_no, 
		PNT pupil_ellipse_points[], int *no_p_elps_pnts, PNT CR_ellipse_points[], int *no_CR_elps_pnts,
		LBLD_AREA_BUF CRthreshBuf[], int *cr_threshbuf_size, AREA crArea[], int *no_cr_areas, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int labld_array[], int validAreaCRs[], int *noValidCRAreas
	);

};

#endif

#include <stdio.h>
#include "Clock.h"
#include "ProcessImage.h"
#include "FeatureRecognition.h"
#include "CannyEdgeDetect.h"
#include "ImageProcessingUtilities_DSP.h"
#include "ImageProcessingUtilities_FPGA.h"
#include "ImageProcessingUtilities_Debug.h"


FeatureRecognition::FeatureRecognition(void)
{
}
void FeatureRecognition::PCrRecognition
	( int *Isizeh, int *Isizev, EDGE_DETECTION_PARAM_FPGA *pfpga, EDGE_DETECTION_PARAM_DSP *pdsp, OBJECT_RECOGNITION_PARAMETERS *pupil_rec_param, //INPUTS
		OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, PNT Pbuf[], PNT CRbuf[], EDGE_POINT_STATUS edgepoint[], EDGE_POINT_STATUS CR_edgepoint[], 
		int PDirBuf[], int CrDirBuf[], int *Pbuf_size, int *CRbuf_size,	FuzzParam *Pfuzz, int grad_H_array[], int grad_V_array[], int CRarray[],
		bool *pupil_found, bool *CR_found, EDGE_STRING es[], int *nstring,  EDGE_STRING CR_es[], int *nCRstring,                                  //OUTPUTS
		LFT_EDGE lft_edge[], int *no_lft_edges, LFT_EDGE lft_CR_edge[], int *no_lft_CR_edges, RT_EDGE rt_edge[], int *no_rt_edges, RT_EDGE rt_CR_edge[], int *no_rt_CR_edges,
		OBJECT objt[], int *no_of_objects, CR_OBJECT CR_objt[], int *no_of_CR_objects, int *pupil_object_no, int * CR_object_no, 
		PNT pupil_ellipse_points[], int *no_p_elps_pnts, PNT CR_ellipse_points[], int *no_CR_elps_pnts, 
		LBLD_AREA_BUF CRthreshBuf[], int *cr_threshbuf_size, AREA crArea[], int *no_cr_areas, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int labld_array[], 
		int validCRs[], int *noValidCRs 
	)
//Inputs:
//	Isezeh;				horiz eye image size (pixels)
//  Isezev;				vert eye image size (pixels)
//	pfpga;				edge detection parameters that will be used by FPGA
//  pdsp;				edge detection parameters that will be used by DSP
//  pupil_rec_param;	image processing parameters used by DPS to find pupil
//  CR_rec_param;		image processing parameters used by DPS to find CR
//  Pbuf;				List of possible pupil edge points.  Each list element has h (col) and v(row) pixel address 
//  CRbuf;				List of possible CR edge points.
//  PDirBuf;			List of elements corresponding to Pbuf.  Each list element has a direction value
//  CrDirBuf;			List of elements corresponding to CRbuf.  Each list element has a direction value
//  Pbuf_size;			No of list elements in Pbuf and PDirBuf.
//  CRbuf_size;			No of list elements in CRbuf and CrDirBuf.
//  Pfuzz;				Struct containing fuzzy logic parameters for selection of best pupil
//	***** CRfuzz not yet included				Struct containing fuzzy logic parameters for selection of best CR
//  grad_H_array;		Although used internally by edge detection routine, this is passed to feature recogniton routine only for debugging help (Feature Recognition Routine may use it to create some arrays that are only for debugging)
//  grad_V_array;		Although used internally by edge detection routine, this is passed to feature recogniton routine only for debugging help (Feature Recognition Routine may use it to create some arrays that are only for debugging)
//  CRarray;            original gray scale image or smoothed image, with all values below CR threshold zeroed. Passed to feature rec routing only to enable array based "labeled region" analysis for testing. DSP will do "labeled region" anal with CR buff. 
//
//Outputs:
//  pupil_found;	true if a pupil was identified
//  CR_found;		true if a CR was identified
//  es;				Array of structure elements. Each element contains a list of continuous points that could be part of a pupil edge.
//  nstring;		Number of pupil strings (no of es elements). 
//  lft_edge;		Array of structure elements. Each element contains a list of edge strings (each designated by es array index) that might form the left edge of a pupil object.
//  no_lft_edges;   Number of possible left pupil edges (no of lft_edge elements).
//  rt_edge;		Array of structure elements. Each element contains a list of edge strings (each designated by es array index) that might form the right edge of a pupil object.
//  no_rt_edges;    Number of possible right pupil edges (no of rt_edge elements).
//  objt;			Array of structure elements. Each element contains a list of left edges and rt edges (each designated by lft_edge or rt_edge array index) that might form a single pupil object. 
//  no_of_objects;	Number of pupil objects (no fo objt elements). 
//  pupil_object_no; index of objt element determined to be the pupil.  	
//  CR_es;          Array of structure elements. Each element contains a list of continuous points that could be part of a CR edge.
//  nCRstring;      Number of CR strings (no of CR_es elements).
//  lft_CR_edge;	Array of structure elements. Each element contains a list of edge strings (each designated by CR_es array index) that might form the left edge of a CR object.
//  no_lft_CR_edges;Number of possible left pupil edges (no of lft_edge elements).
//  rt_CR_edge;		Array of structure elements. Each element contains a list of edge strings (each designated by CR_es array index) that might form the right edge of a CR object.
//  no_rt_CR_edges; Number of possible right pupil edges (no of rt_CR edge elements).
//  CR_objt;		Array of structure elements. Each element contains a list of left edges and rt edges (each designated by lft_CR_edge or rt_CR_edge array index) that might form a single CR object.
//  no_of_CR_objects; Number of CR objects (no of CR_objt elements).
//  CR_object_no;	index of CR_objt element determined to be the CR.
//  pupil_ellipse_points;  list of points on the computed pupil ellipse (use these to draw the ellipse).
//  no_p_elps_pnts;  Number of points in pupil_ellipse_points.
//  CR_ellipse_points;  list of points on the computed CR ellipse (use these to draw the ellipse).
//  no_CR_elps_pnts; Number of points in CR_ellipse_points. 
//	LBLD_AREA_BUF CRthreshBuf[]; List of points with non zero values in CRarray (points that were above CR threshold value). 
//	int *cr_threshbuf_size;		 number of elements in CRthreshBuf[] 
//  crArea;          An array of AREA structure elements.  Each element contains information about a possible CR derived from labeled region analysis.
//  no_cr_areas;     Number of AREA elements in crArea array.
//  labeled_area;    Array of LABELED_AREA structure elements. Each element corresponds to a labeled region and contains a list of points in that region. (Used only within labeled regions computation algorithm.)
//  EquivToLabel;    Array of EquivTable structure elements. Each element corresponds to a labeled region and lists other regions that touch it (and are therefore equivalent).  (Used only within labeled regions computation algorithm.)
//  labld_array;     Copy of CRarray, but with each above threshold pixel set to a region label.  (Used only within labeled regions computation algorithm.)
//  validCRs;        Integer array with list of crArea indecies specifying crArea elements that have the characteristics to be valid CRs.
//  noValidCRs;      Number of elements in validCRs array.  

// NOTE: the only outputs really needed by calling routines are: 
//       For pupil information --
//       pupil_found flag; objt, no_of_objects, pupil_object_no, pupil_ellipse_points and no_p_elps_pnts; 
// 
//       For CR information-- 
//       CR_objt, no_of_CR_objects, CR_object_no, CR_ellipse_points, and no_CR_elps_pnts (if edge data used to find CRs).
//       OR
//       crArea, no_cr_areas, validCRs, and noValidCRs (if labeled region analysis used to find CRs).

// All other outputs are intermediate computations, used by the higher level calling routines only to debugging display purposes.
// All inputs and outputs are pointers to variables "owned" by "ProcessImage.cpp".  Memeory for all of them is allocated in "ProcessImage.cpp". 


{
	int pemax;
	clock_start("FeatureRecognition::PCrRecognition");
	//STRING_ARRAY S[500], crS[500]; 
	//STRNG_DIR_ARRAY D[500], crD[500];
	//STRING_ARRAY Ghv[500], crGhv[500];

	//LARGE_INTEGER t1, t2, t3;

	//QueryPerformanceFrequency(&t1);
	//QueryPerformanceCounter(&t2);

    //**** Decide which sets of edge buffer points form pupil strings (unbroken lines that might be part of a pupil edge) 
	ImageProcessingUtilities_DSP::FindStrings(pupil_rec_param, Pbuf, edgepoint, Pbuf_size, PDirBuf, nstring, &es[0] );

	//QueryPerformanceCounter(&t3);

	//double dt1 = t1.QuadPart;
	//double dt2 = t2.QuadPart;
	//double dt3 = t3.QuadPart;
	//double deltat = (dt3 - dt2) / dt1;
	//WCHAR str[80];
	//swprintf_s(str, 80, L"deltat=%f", deltat);
	//OutputDebugString(str);

    //**** Decide which strings or combinations of strings form possible pupil right edge sections and left edge sections 
    ImageProcessingUtilities_DSP::MakeEdges(pupil_rec_param, es, nstring, lft_edge, rt_edge, no_lft_edges, no_rt_edges );

    //**** Decide which left and right edge sections form potential pupil objects
    ImageProcessingUtilities_DSP::MakeObjects(pupil_rec_param, objt, no_of_objects, lft_edge, rt_edge, no_lft_edges, no_rt_edges, es, nstring, Pbuf, Pbuf_size);

    // decide whether to use labeled region analysis or edge information to identify CRs
    if (CR_rec_param->use_labeled_area_analysis_for_CR) //use labeled region analysis to find areas likely to be CRs and the edge points of these areas (return "AREA crArea" structure for each)
		//we can do labeled region analysis either with a binary image array or a buffer (list of points that are above threshold).
		if( CR_rec_param->use_threshold_buf_for_labeled_area_anal ) //use the buffer (list of points).
			ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis2( CR_rec_param, CRthreshBuf, cr_threshbuf_size, labeled_area, EquivToLabel, 
																	     crArea, no_cr_areas, validCRs, noValidCRs  );
		else //use the full array.
			ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis( CR_rec_param, CRarray, Isizeh, Isizev, &pfpga->wL, &pfpga->wR, &pfpga->wT, &pfpga->wB, 
   	                                                                    labeled_area, EquivToLabel, labld_array, crArea, no_cr_areas, validCRs, noValidCRs  );

	else //find strings and edges to make CR objects and identify CRs (return CR_objt structure)
	{
		//*** Make strings, edges,and objects for CR*****
		ImageProcessingUtilities_DSP::FindStrings( CR_rec_param, CRbuf, CR_edgepoint, CRbuf_size, CrDirBuf, nCRstring, &CR_es[0] ); //not dark object (false)
		ImageProcessingUtilities_DSP::MakeEdges( CR_rec_param, CR_es, nCRstring, lft_CR_edge, rt_CR_edge, no_lft_CR_edges, no_rt_CR_edges );
		ImageProcessingUtilities_DSP::MakeCRObjects( CR_rec_param, CR_objt, no_of_CR_objects, lft_CR_edge, rt_CR_edge, no_lft_CR_edges, no_rt_CR_edges );
	}
//QueryPerformanceCounter(&t2); //timer start
    //***** fit elipse to pupil objects and use fuzzy logic to choose the one that is most likely the real pupil 
    *pupil_object_no = 0;//initialize
    *CR_object_no = 0;
    if ( *no_of_objects > pupil_rec_param->feature_max_no_of_elps_computations) pemax = pupil_rec_param->feature_max_no_of_elps_computations;
    else pemax = *no_of_objects;
    if (pemax > 0)
    {
        for (int n = 0; n < pemax; n++)
        {
            if (ImageProcessingUtilities_DSP::FitPupilEllipse(&objt[n], lft_edge, rt_edge, es, Pbuf, Pbuf_size, pupil_rec_param))
                objt[n].scaled_diam = (2.0f * objt[n].elps_maj_rad) / pupil_rec_param->EyeImageScaleFactor;
            else objt[n].scaled_diam = 0.0;
        }
        *pupil_found = ImageProcessingUtilities_DSP::FindBestPupilObject(objt, no_of_objects, pupil_object_no, Pfuzz );//fuzzy logic
    }
    else *pupil_found = false;
	//QueryPerformanceCounter(&t3); //timer stop
	//double dt1 = t1.QuadPart; //timer
	//double dt2 = t2.QuadPart;//timer
	//double dt3 = t3.QuadPart;//timer
	//double deltat = (dt3 - dt2) / dt1; //timer
	//WCHAR str[80]; //output timer value so it can be read in DBCON 
	//swprintf_s(str, 80, L"deltat=%f", deltat); //output timer value so it can be read in DBCON
	//OutputDebugString(str); //output timer value so it can be read in DBCON

    if (CR_rec_param->use_labeled_area_analysis_for_CR)//see whether we have used labled region analysis or edge information to find identify CRs
	{	//***** fit elipse to CR labeled regions 
		ImageProcessingUtilities_DSP::FitElipseToLabeledAreaCRs( CR_rec_param, crArea, validCRs, noValidCRs );
		if (*noValidCRs > 0) *CR_found = true;
		else *CR_found = false;
	}
	else
	{   //***** fit elipse to CR edge detection objects and choose the one that is most likely the real CR
		ImageProcessingUtilities_DSP::FitElipseToEdgeObjectCRs( CR_objt, CR_rec_param, lft_CR_edge, rt_CR_edge,
		                                                         CR_es, CRbuf, CRbuf_size, no_of_CR_objects, validCRs, noValidCRs );
	    if( *noValidCRs > 0) *CR_found = true;
		else *CR_found = false;
	}

    //**** plot pupil elipse outline
    *no_p_elps_pnts = 0; //initialize
    if (*pupil_found)//get list of points on computed pupil ellipse 
    {
        ImageProcessingUtilities_DSP::plot_ellipse( pupil_ellipse_points, no_p_elps_pnts, &objt[*pupil_object_no].elps_hc, &objt[*pupil_object_no].elps_vc,
                     &objt[*pupil_object_no].elps_maj_rad, &objt[*pupil_object_no].elps_min_rad, &objt[*pupil_object_no].elps_angle);
    }

	//**** plot CR elipse outlines
    if (CR_rec_param->use_labeled_area_analysis_for_CR)
	{
		for( int ii = 0; ii < *noValidCRs; ii++ ) //get list of points on computed CR ellipses
		{
			ImageProcessingUtilities_DSP::plot_ellipse( crArea[ validCRs[ii] ].elpspnt, &crArea[ validCRs[ii] ].no_elps_pnts, &crArea[ validCRs[ii] ].hc, &crArea[ validCRs[ii] ].vc,
						 &crArea[ validCRs[ii] ].elps_maj_rad, &crArea[ validCRs[ii ]].elps_min_rad, &crArea[ validCRs[ii] ].elps_angle);
		}
	}
	else
	{
		*no_CR_elps_pnts = 0;
		if (*CR_found)//get list of points on computed CR ellipse 
		{
			ImageProcessingUtilities_DSP::plot_ellipse( CR_ellipse_points, no_CR_elps_pnts, &CR_objt[*CR_object_no].hc, &CR_objt[*CR_object_no].vc,
						 &CR_objt[*CR_object_no].elps_maj_rad, &CR_objt[*CR_object_no].elps_min_rad, &CR_objt[*CR_object_no].elps_angle);
		}
	}
    clock_end();
}

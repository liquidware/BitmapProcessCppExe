#include <stdio.h>
#include <string.h>
#include "ProcessImage.h"
#include "FeatureRecognition.h"
#include "CannyEdgeDetect.h"
#include "ImageProcessingUtilities_DSP.h"
#include "ImageProcessingUtilities_FPGA.h"
#include "Clock.h"

// top level call
void FindPupilCR
(
// inputs
int Isizeh,	// horiz eye image size (pixels)
int Isizev,	// vert eye image size (pixels)
int Iarray[],// Eye cam image in form of a singly subscripted array. The array index is [column# + ( row# * width)]

// outputs
PUPIL_CR_OBJECT* pupil,	// pupil object returned from analysis
PUPIL_CR_OBJECT CR[],	// array of CR objects from analysis. Must be preallocated with size 5
int* number_CR_objects	// number CRs found. Will be <=5
)
{
	//Get image processing parameters defined in C++ DLL (includes max sizes used for some mem allocations)
	int SkernSize, GkernSize;
	int Skern[5][5];
	int GkernH[5][5];
	int GkernV[5][5];
	EDGE_DETECTION_PARAM_FPGA edp_fpga;
	EDGE_DETECTION_PARAM_DSP edp_dsp;
	OBJECT_RECOGNITION_PARAMETERS orp_P;
	OBJECT_RECOGNITION_PARAMETERS orp_CR;
	FuzzParam Pfuzz;

	clock_start("FindPupilCR");
	GetParameters(&SkernSize, &GkernSize, Skern, GkernH, GkernV, &edp_fpga, &edp_dsp, &orp_P, &orp_CR, &Pfuzz);

	//Allocate memory for variables
	bool pupil_found, CR_found;
	int maxIbuf = Isizeh * Isizev;
	int nstring, no_lft_edges, no_rt_edges, no_of_objects, nCRstring, no_lft_CR_edges, no_rt_CR_edges, no_of_CR_objects, pupil_object_no, CR_object_no, no_p_elps_pnts, no_CR_elps_pnts;
	int Pbufsize, CRbufsize, no_cr_areas, noValidCRAreas, cr_threshbuf_size;

	PNT pupil_ellipse_points[MAX_NO_ELLIPSE_PLOT_PNTS];
	PNT CR_ellipse_points[MAX_NO_ELLIPSE_PLOT_PNTS];

	PNT *Pbuf = new PNT[maxIbuf]; //delete[] Pbuf;
	PNT *CRbuf = new PNT[maxIbuf];
	EDGE_POINT_STATUS *edgepoint = new EDGE_POINT_STATUS[maxIbuf];
	EDGE_POINT_STATUS *CR_edgepoint = new EDGE_POINT_STATUS[maxIbuf];
	int *PDirBuf = new int[maxIbuf];
	int *CrDirBuf = new int[maxIbuf];
	int *grad_H_array= new int[maxIbuf];
	int *grad_V_array = new int[maxIbuf];
	int *CRarray = new int[maxIbuf];
	LFT_EDGE *lft_edge = new LFT_EDGE[orp_P.edge_max_no_edges];
	LFT_EDGE *lft_CR_edge = new LFT_EDGE[orp_CR.edge_max_no_edges];
	RT_EDGE *rt_edge = new RT_EDGE[orp_P.edge_max_no_edges];
	RT_EDGE *rt_CR_edge = new RT_EDGE[orp_CR.edge_max_no_edges];
	OBJECT *objt = new OBJECT[orp_P.obj_max_no_objects];
	CR_OBJECT *CR_objt = new CR_OBJECT[orp_CR.obj_max_no_objects];
	EDGE_STRING *es = new EDGE_STRING[orp_P.string_max_no_strings]; 
	EDGE_STRING *CR_es = new EDGE_STRING[orp_CR.string_max_no_strings];
	AREA *crArea = new AREA[orp_CR.max_labels];
	LABELED_AREA *labeled_area = new LABELED_AREA[orp_CR.max_labels];
	EquivTable *EquivToLabel = new EquivTable[orp_CR.max_labels];
	int *labld_array = new int[maxIbuf]; 
	int *validAreaCRs = new int[orp_CR.max_labels];
	LBLD_AREA_BUF *CRthreshBuf = new LBLD_AREA_BUF[edp_fpga.max_cr_threshbuf_size];
	CANNY_EDGE_DETECT_DATA* data = new CANNY_EDGE_DETECT_DATA();

	//Call routine to do image processing
	DoEyeImageProcessing
		( &Isizeh, &Isizev, Iarray, &SkernSize, Skern, &GkernSize, GkernH, GkernV, &edp_fpga, &edp_dsp, &orp_P, &orp_CR, &Pfuzz,
		&nstring, &no_lft_edges, &no_rt_edges, &no_of_objects, &nCRstring, &no_lft_CR_edges, &no_rt_CR_edges, &no_of_CR_objects, &pupil_object_no, &CR_object_no, &pupil_found, &CR_found,
		Pbuf, CRbuf, edgepoint, CR_edgepoint, &Pbufsize, &CRbufsize, lft_edge, lft_CR_edge, rt_edge, rt_CR_edge, objt, CR_objt, es, CR_es, PDirBuf, CrDirBuf, grad_H_array, grad_V_array, 
		CRarray, pupil_ellipse_points, &no_p_elps_pnts, CR_ellipse_points, &no_CR_elps_pnts, CRthreshBuf, &cr_threshbuf_size, crArea, &no_cr_areas, labeled_area, EquivToLabel, labld_array, validAreaCRs, &noValidCRAreas,
		data);

	if (pupil_found)
	{
		pupil->found = true;
		OBJECT* pupilObj = &(objt[pupil_object_no]);
		pupil->center_horz	= pupilObj->elps_hc;
		pupil->center_vert	= pupilObj->elps_vc;
		pupil->major_radius	= pupilObj->elps_maj_rad;
		pupil->minor_radius	= pupilObj->elps_min_rad;
		pupil->elps_angle	= pupilObj->elps_angle;
		pupil->no_ellipse_points = no_p_elps_pnts;
		for (int i=0; i<no_p_elps_pnts; i++)
			pupil->ellipse_points[i] = pupil_ellipse_points[i];
	}
	else
	{
		memset(pupil, 0, sizeof(PUPIL_CR_OBJECT));
	}

	if (CR_found)
	{
		if (orp_CR.use_labeled_area_analysis_for_CR)
		{
			// Create an array of CRs
			*number_CR_objects = (noValidCRAreas <= MAX_CR_OBJECTS)? noValidCRAreas : MAX_CR_OBJECTS;
			for (int i=0; i<*number_CR_objects; i++)
			{
				int index = validAreaCRs[i];
				AREA* crAr = &(crArea[index]);
				PUPIL_CR_OBJECT* pCR = &(CR[i]);
				pCR->center_horz = crAr->hc;
				pCR->center_vert = crAr->vc;
				pCR->major_radius = crAr->elps_maj_rad;
				pCR->minor_radius = crAr->elps_min_rad;
				pCR->elps_angle = crAr->elps_angle;
				pCR->no_ellipse_points = no_CR_elps_pnts;
				for (int j=0; j<no_CR_elps_pnts; j++)
					pCR->ellipse_points[j] = CR_ellipse_points[j];
			}
			
			// Reset unused members of CR array
			for (int i=noValidCRAreas; i<MAX_CR_OBJECTS; i++)
				memset(&CR[i], 0, sizeof(PUPIL_CR_OBJECT));
		}
		else // use ...
		{
			*number_CR_objects = 1;
			CR_OBJECT* crObj = &(CR_objt[CR_object_no]);
			CR->center_horz = crObj->hc;
			CR->center_vert = crObj->vc;
			CR->major_radius = crObj->elps_maj_rad;
			CR->minor_radius = crObj->elps_min_rad;
			CR->elps_angle = crObj->elps_angle;
		}
	}
	else // CR not found
	{
		*number_CR_objects = 0;
		for (int i=0; i<MAX_CR_OBJECTS; i++)
			memset(&CR[i], 0, sizeof(PUPIL_CR_OBJECT));
	}

	delete[] Pbuf;
	delete[] CRbuf;
	delete[] edgepoint;
	delete[] CR_edgepoint;
	delete[] PDirBuf;
	delete[] CrDirBuf;
	delete[] grad_H_array;
	delete[] grad_V_array;
	delete[] CRarray;
	delete[] lft_edge;
	delete[] lft_CR_edge;
	delete[] rt_edge;
	delete[] rt_CR_edge;
	delete[] objt;
	delete[] CR_objt;
	delete[] es; 
	delete[] CR_es;
	delete[] crArea;
	delete[] labeled_area;
	delete[] EquivToLabel;
	delete[] labld_array; 
	delete[] validAreaCRs;
	delete[] CRthreshBuf;
	delete data;
	clock_end();
}

void DoEyeImageProcessing
(

 // Inputs					   
 int *Isizeh,	// horiz eye image size (pixels)
 int *Isizev,	// vert eye image size (pixels)
 int Iarray[],	// Eye cam image in form of a doubly subscripted array. The subscripts are col and row 
				// addresses and each element value is the gray scale level. 
 int *SkernSize,	// dimension of Kernel (single value since kernel array is always square)
 int Skern[5][5],	// Kernel for smoothing convolution
 int *GkernSize,	// dimension of GkernH and GkernV (both are same size and are square)
 int GkernH[5][5],	// Kernel for horizontal gradient convolution 
 int GkernV[5][5],	// Kernel for vertical gradient convolution  
 EDGE_DETECTION_PARAM_FPGA *edp_fpga,	// edge detection parameters that will be used by FPGA
 EDGE_DETECTION_PARAM_DSP *edp_dsp,		// edge detection parameters that will be used by DSP
 OBJECT_RECOGNITION_PARAMETERS *orp_P,	// object recognition parameters for pupil
 OBJECT_RECOGNITION_PARAMETERS *orp_CR,	// object recognition parameters for CR
 FuzzParam *Pfuzz,	// Struct containing fuzzy logic parameters for selection of best pupil

 // Outputs
 int *nstring,			// Number of pupil strings (no of es elements). 
 int *no_lft_edges,		// Number of possible left pupil edges (no of lft_edge elements).
 int *no_rt_edges,		// Number of possible right pupil edges (no of rt_edge elements).
 int *no_of_objects,	// Number of pupil objects (no fo objt elements). 
 int *nCRstring,		// Number of CR strings (no of CR_es elements).
 int *no_lft_CR_edges,	// Number of possible left pupil edges (no of lft_edge elements). 
 int *no_rt_CR_edges,	// Number of possible right pupil edges (no of rt_CR edge elements).
 int *no_of_CR_objects,	// Number of CR objects (no of CR_objt elements).
 int *pupil_object_no,	// most likely pupil object from objt[]
 int *CR_object_no,		// most likely CR object from CR_objt[]
 bool *pupil_found,		// true if a pupil was identified
 bool *CR_found,		// true if a CR was identified
 PNT Pbuf[],			// List of possible pupil edge points.  Each list element has h (col) and v(row) pixel address 
 PNT CRbuf[],			// List of possible CR edge points
 EDGE_POINT_STATUS edgepoint[], // list corresponding to Pbuf[] elements indicating whether point has been assigned to a string
 EDGE_POINT_STATUS CR_edgepoint[], // same as above for CRbuf[]
 int *Pbufsize,			// No of list elements in Pbuf and PDirBuf.
 int *CRbufsize,		// No of list elements in CRbuf and CrDirBuf.
 LFT_EDGE lft_edge[],	// Array of structure elements. Each element contains a list of edge strings (each designated by es array index) that might form the left edge of a pupil object.
 LFT_EDGE lft_CR_edge[],// Array of structure elements. Each element contains a list of edge strings (each designated by CR_es array index) that might form the left edge of a CR object.
 RT_EDGE rt_edge[],		// same as lft_edge[] for possible right pupil edges
 RT_EDGE rt_CR_edge[],	// same as lft_CR_edge[] for possible right CR edges
 OBJECT objt[],			// array of possible pupil objects
 CR_OBJECT CR_objt[],	// array of possible CR objects
 EDGE_STRING es[],		// possible pupil strings
 EDGE_STRING CR_es[],	// possible CR strings (not used if we chose Labeled Region Analysis for CR)
 int PDirBuf[],			// list containing direction value for each member of Pbuf[]
 int CrDirBuf[],		// same as above for CRbuf[]
 int grad_H_array[],	// Although used internally by edge detection routine, this is passed to feature recogniton routine only for debugging help (Feature Recognition Routine may use it to create some arrays that are only for debugging)
 int grad_V_array[],	// same as above
 int CRarray[],			// original gray scale image or smoothed image, with all values below CR threshold zeroed. Passed to feature rec routing only to enable array based "labeled region" analysis for testing. DSP will do "labeled region" anal with CR buff. 
 PNT pupil_ellipse_points[],// list of points on the computed pupil ellipse (use these to draw the ellipse)
 int *no_p_elps_pnts,	// Number of points in pupil_ellipse_points
 PNT CR_ellipse_points[],// list of points on the computed CR ellipse (use these to draw the ellipse)
 int *no_CR_elps_pnts,	// Number of points in CR_ellipse_points
 LBLD_AREA_BUF CRthreshBuf[], //List of points with non zero values in CRarray (points that were above CR threshold value). 
 int *cr_threshbuf_size,// number of elements in CRthreshBuf[] 
 AREA crArea[],			// An array of AREA structure elements.  Each element contains information about a possible CR derived from labeled region analysis
 int *no_cr_areas,		// Number of AREA elements in crArea array
 LABELED_AREA labeled_area[],// Array of LABELED_AREA structure elements. Each element corresponds to a labeled region and contains a list of points in that region (used only within labeled regions computation algorithm)
 EquivTable EquivToLabel[],// Array of EquivTable structure elements. Each element corresponds to a labeled region and lists other regions that touch it (and are therefore equivalent) (used only within labeled regions computation algorithm)
 int labld_array[],		// Copy of CRarray, but with each above threshold pixel set to a region label (used only within labeled regions computation algorithm)
 int validAreaCRs[],	// Integer array with list of crArea indecies specifying crArea elements that have the characteristics to be valid CRs
 int *noValidCRAreas,	// Number of elements in validCRs array
 CANNY_EDGE_DETECT_DATA* data // storage for temp vars in procedure CannyEdgeDetect
)

 //************ C++ image processing ******************
 //Recieve image array, convolution kernels, edge detection and image processing parameters from C++ calling program.
 //Return flags indicating whether pupil and CR were found.
 //Return pupil and CR object struct arrays, with index of object determined to be the real pupil and CR. 
 //Also return structures and arrays with intermediate edge and object detection information for possible display.  
{
	clock_start("DoEyeImageProcessing");
	//int dummy;
	//LARGE_INTEGER t1, t2, t3;
	//QueryPerformanceFrequency(&t1);

	//edp_fpga->wL = 0;
	//edp_fpga->wR = *Isizeh-1;
	//edp_fpga->wT =0;
	//edp_fpga->wB = *Isizev-1;

	edp_fpga->wL = 5;
	edp_fpga->wR = 290; //*Isizeh - 25;
	edp_fpga->wT = 5;
	edp_fpga->wB = 230; //*Isizev - 10;

	//QueryPerformanceCounter(&t2);
	CannyEdgeDetect::PCrEdgeDetect
		( 
		Isizeh, Isizev, Iarray, edp_fpga, edp_dsp, orp_CR, Skern, SkernSize, GkernH, GkernV, //Inputs
		GkernSize, Pbuf,  CRbuf, PDirBuf, CrDirBuf, Pbufsize, CRbufsize, grad_H_array, grad_V_array, CRarray, CRthreshBuf, cr_threshbuf_size, //Outputs
		data
	);
	//QueryPerformanceCounter(&t2);
	FeatureRecognition::PCrRecognition
		( 
		Isizeh, Isizev, edp_fpga, edp_dsp, orp_P, orp_CR, Pbuf, CRbuf, edgepoint, CR_edgepoint, PDirBuf, CrDirBuf, Pbufsize, CRbufsize, Pfuzz, grad_H_array, grad_V_array, CRarray,//Inputs
		pupil_found, CR_found, es, nstring,  CR_es, nCRstring, lft_edge, no_lft_edges, lft_CR_edge, no_lft_CR_edges, rt_edge, no_rt_edges, rt_CR_edge, no_rt_CR_edges, //Outputs
		objt, no_of_objects, CR_objt, no_of_CR_objects, pupil_object_no, CR_object_no, pupil_ellipse_points, no_p_elps_pnts, CR_ellipse_points, no_CR_elps_pnts, 
		CRthreshBuf, cr_threshbuf_size, crArea, no_cr_areas, labeled_area, EquivToLabel, labld_array, validAreaCRs, noValidCRAreas  
		);
	//QueryPerformanceCounter(&t3);
	//double dt1 = t1.QuadPart;
	//double dt2 = t2.QuadPart;
	//double dt3 = t3.QuadPart;
	//double deltat = (dt3 - dt2) / dt1;
	//dummy = 1;
	////WCHAR str[80];
	////swprintf_s(str, 80, L"deltat=%f", deltat);
	////OutputDebugString(str);

	clock_end();

}


void GetParameters(int *SkernSize, int *GkernSize, int Skern[5][5], int GkernH[5][5], int GkernV[5][5], 
				   EDGE_DETECTION_PARAM_FPGA *edp_fpga, EDGE_DETECTION_PARAM_DSP *edp_dsp, OBJECT_RECOGNITION_PARAMETERS *orp_P, OBJECT_RECOGNITION_PARAMETERS *orp_CR, FuzzParam *Pfuzz)
{
	EDGE_DETECTION_PARAM_FPGA edp_fpga_BrightPupil_1;
	EDGE_DETECTION_PARAM_DSP edp_dsp_BrightPupil_1;
	OBJECT_RECOGNITION_PARAMETERS orp_BrightPupil_1;
	OBJECT_RECOGNITION_PARAMETERS orp_CR_BrightPupil_1;
	FuzzParam Pfuzz_BrightPupil_1;
	int K_Gsmooth_5x5[5][5];
	int K_SobelH_3x3[3][3];
	int K_SobelV_3x3[3][3];
	int K_GgradH_5x5[5][5];
	int K_GgradV_5x5[5][5];

	EDGE_DETECTION_PARAM_FPGA edp_fpga_BrightPupil_2;
	EDGE_DETECTION_PARAM_DSP edp_dsp_BrightPupil_2;
	OBJECT_RECOGNITION_PARAMETERS orp_BrightPupil_2;

	//************** BRIGHT PUPIL 2 Parameters (with Sobel grad conv) ************************************************************//
	//EDGE_DETECTION_PARAM_FPGA
	edp_fpga_BrightPupil_2.pupil_type = BRIGHT; //int -- DARK (0) or BRIGHT (1)
	edp_fpga_BrightPupil_2.smooth = false;//bool -- use smoothing convolultion if true
	edp_fpga_BrightPupil_2.extra_bndry = 2;
	edp_fpga_BrightPupil_2.SmthCnvl_Scale = 1.0;
	edp_fpga_BrightPupil_2.grad_thresh = 50; //int - thresh after computing gradient (std = 1000; older std =1500)
	edp_fpga_BrightPupil_2.grayhist_size = 256; //int -- no elements in grayscale level histogram
	edp_fpga_BrightPupil_2.grad_hist_size = 32767; //int -- usually size of largest signed integer (32767)
	edp_fpga_BrightPupil_2.max_cr_threshbuf_size = 1000;//int
	edp_fpga_BrightPupil_2.grad_thresh_for_dir_array = 0; //int
	edp_fpga_BrightPupil_2.max_pbuf_size = 100000;//int
	edp_fpga_BrightPupil_2.max_CRbuf_size = 1000;//int 
	edp_fpga_BrightPupil_2.GradCnvl_Scale = 1.0f; // float -- usually 0.125 for 5x5 convl
	edp_fpga_BrightPupil_2.use_thresh_hyst = false; //bool -- if false, don't do dual thresh hysterisis, just use lower threshold boundary.
	edp_fpga_BrightPupil_2.hair_width = 3; //width of "thin dark line", in pixels, assumed to be hair or similar  aritfact.

	//EDGE_DETECTION_PARAM_DSP
	edp_dsp_BrightPupil_2.pupil_type = BRIGHT; //int -- dark (0) or bright (1)
	edp_dsp_BrightPupil_2.CRthresh_hist_fraction = 0.9995f; //use to compute CR thresh from histogram (1st pixel brighter than this percent of all pixels)(std = 0.9995) [not currently used]
	edp_dsp_BrightPupil_2.min_pup_pix = 100; //int -- The brightest intensity having at least this many pixels is probably just brighter than the pupil (std = 100)
	edp_dsp_BrightPupil_2.MaxGrayLimit = 225; //int -- max value used as "brighteset pixel" in CR thresh computation (std = 225)
	edp_dsp_BrightPupil_2.frac_dist_brtpup_to_brightest_pixel__CRthresh = 0.7f; //float -- Fraction of dist between probable pupil brightness and brightest pixel. (std = 0.7)
	edp_dsp_BrightPupil_2.high_thresh_fraction = 0.43f; //float -- high thresh (for dual thresh) will be value such that at least this fraction of pixels has lower gradient. (std = 0.5)
	edp_dsp_BrightPupil_2.low_thresh_fraction = 0.33f; //float -- low threshold will be this fraction of high threshold. (std = 0.33)
	edp_dsp_BrightPupil_2.CRthresh_dist = 1; //int -- Require a pixel this far in the grad direction to be above CRthresh for CR edge point. CRthresh is determined from hist. (std = 1)

	//Pupil OBJECT_RECOGNITION_PARAMETERS-- General
	orp_BrightPupil_2.pupil_type = BRIGHT; //int -- dark or bright
	orp_BrightPupil_2.EyeImageScaleFactor = 10.0; //float -- pixels per mm

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for identifying strings
	orp_BrightPupil_2.string_minHbreak = 1;//int
	orp_BrightPupil_2.string_minVbreak = 1;//int
	orp_BrightPupil_2.string_minstringlength = 5;//int --  = 5 for std bright pupil
	orp_BrightPupil_2.string_max_delta_dir = 70;//int --  = 50 for std bright pupil
	orp_BrightPupil_2.string_max_no_strings = 10000;//int
	orp_BrightPupil_2.string_max_points_per_string = 1000;//int
	orp_BrightPupil_2.string_HairWidth = 3; //int

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting and combining strings to form left and right edges
	orp_BrightPupil_2.edge_max_no_edges = 1000;//int
	orp_BrightPupil_2.edge_mindir = 50.0f;//float --  = 100.0 //average gradient should have significant rt or left component for right & lft edges(not straight up or down).
	orp_BrightPupil_2.edge_maxdir = 350.0f;//float -- = 300.0 //(0 is straight dwn, -400 is straight up, 200 is rt, -200 is left)
	orp_BrightPupil_2.edge_minedgelength = 10;//int-- = 15 // shortest string length that can be an "edge"
	orp_BrightPupil_2.edge_maxedgelength = 500;//int -- longest string liength that can be an "edge"
	orp_BrightPupil_2.edge_max_delta_curvature = 30.0f;//float -- = 15.0 //we want consistant curvature (rate of change of direction); ave changes in curvature (delta curvature) should be small
	orp_BrightPupil_2.edge_mincurv = 2.0;//float -- = 2;//previous .01 //we only want to consider strings with reasonable ave curvature
	orp_BrightPupil_2.edge_maxcurv = 10.0;//float -- = 8;//previous .05 //we don't want straight lines or curves that are too tight

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting sets of edges to form objects
	orp_BrightPupil_2.obj_max_no_objects = 50;//int
	orp_BrightPupil_2.obj_max_curv_dif = 4.0f;//float -- = 1; //max difference between left edge ave curvature and right edge ave curvature
	orp_BrightPupil_2.obj_minPD_to_width = 0.5f;//float -- = 0.5; //min for ratio: [estimate of PD from ave curvature] / [dist between lftmost and rtmost points]
	orp_BrightPupil_2.obj_maxPD_to_width = 2.5f;//float -- = 2.0; //max for same ratio
	orp_BrightPupil_2.obj_minwdth = 10;//int
	orp_BrightPupil_2.obj_maxwdth = 200;//int
	orp_BrightPupil_2.obj_minht = 10;//int
	orp_BrightPupil_2.obj_maxht = 200;//int
	orp_BrightPupil_2.obj_neutral_delta_dir = 200;//int
	orp_BrightPupil_2.obj_max_edge_end_point_separation_sqrd = 2.0;//float --
	orp_BrightPupil_2.obj_max_no_objects = 100;//int
	orp_BrightPupil_2.obj_max_edges_per_object = 2;//int

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for choosing best object (object most likely to be feature of interest)
	orp_BrightPupil_2.feature_max_no_of_elps_computations = 10;//int
	orp_BrightPupil_2.feature_max_scaled_diam = 10.0f;//float --

	//************** BRIGHT PUPIL 1 Parameters *************************************************************//
	//EDGE_DETECTION_PARAM_FPGA
	edp_fpga_BrightPupil_1.pupil_type = BRIGHT; //int -- DARK (0) or BRIGHT (1)
	edp_fpga_BrightPupil_1.smooth = false;//bool -- use smoothing convolultion if true
	edp_fpga_BrightPupil_1.extra_bndry = 2;
	edp_fpga_BrightPupil_1.SmthCnvl_Scale = 1.0;
	edp_fpga_BrightPupil_1.grad_thresh = 1000; //int - thresh after computing gradient (std = 1000; older std =1500)
	edp_fpga_BrightPupil_1.grayhist_size = 256; //int -- no elements in grayscale level histogram
	edp_fpga_BrightPupil_1.grad_hist_size = 32767; //int -- usually size of largest signed integer (32767)
	edp_fpga_BrightPupil_1.max_cr_threshbuf_size = 1000;//int
	edp_fpga_BrightPupil_1.grad_thresh_for_dir_array = 0; //int
	edp_fpga_BrightPupil_1.max_pbuf_size = 100000;//int
	edp_fpga_BrightPupil_1.max_CRbuf_size = 1000;//int 
	edp_fpga_BrightPupil_1.GradCnvl_Scale = 0.125f; // float -- usually 0.125 for 5x5 convl
	edp_fpga_BrightPupil_1.use_thresh_hyst = false; //bool -- if false, don't do dual thresh hysterisis, just use lower threshold boundary.
	edp_fpga_BrightPupil_1.hair_width = 3; //width of "thin dark line", in pixels, assumed to be hair or similar  aritfact.

	//EDGE_DETECTION_PARAM_DSP
	edp_dsp_BrightPupil_1.pupil_type = BRIGHT; //int -- dark (0) or bright (1)
	edp_dsp_BrightPupil_1.CRthresh_hist_fraction = 0.9995f; //use to compute CR thresh from histogram (1st pixel brighter than this percent of all pixels)(std = 0.9995) [not currently used]
	edp_dsp_BrightPupil_1.min_pup_pix = 100; //int -- The brightest intensity having at least this many pixels is probably just brighter than the pupil (std = 100)
	edp_dsp_BrightPupil_1.MaxGrayLimit = 225; //int -- max value used as "brighteset pixel" in CR thresh computation (std = 225)
	edp_dsp_BrightPupil_1.frac_dist_brtpup_to_brightest_pixel__CRthresh = 0.7f; //float -- Fraction of dist between probable pupil brightness and brightest pixel. (std = 0.7)
	edp_dsp_BrightPupil_1.high_thresh_fraction = 0.7f; //float -- high thresh (for dual thresh) will be value such that at least this fraction of pixels has lower gradient. (std = 0.5)
	edp_dsp_BrightPupil_1.low_thresh_fraction = 0.33f; //float -- low threshold will be this fraction of high threshold. (std = 0.33)
	edp_dsp_BrightPupil_1.CRthresh_dist = 1; //int -- Require a pixel this far in the grad direction to be above CRthresh for CR edge point. CRthresh is determined from hist. (std = 1)

	//Pupil OBJECT_RECOGNITION_PARAMETERS-- General
	orp_BrightPupil_1.pupil_type = BRIGHT; //int -- dark or bright
	orp_BrightPupil_1.EyeImageScaleFactor = 10.0; //float -- pixels per mm

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for identifying strings
	orp_BrightPupil_1.string_minHbreak = 1;//int
	orp_BrightPupil_1.string_minVbreak = 1;//int
	orp_BrightPupil_1.string_minstringlength = 5;//int --  = 5 for std bright pupil
	orp_BrightPupil_1.string_max_delta_dir = 70;//int --  = 50 for std bright pupil
	orp_BrightPupil_1.string_max_no_strings = 10000;//int
	orp_BrightPupil_1.string_max_points_per_string = 1000;//int
	orp_BrightPupil_1.string_HairWidth = 3; //int

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting and combining strings to form left and right edges
	orp_BrightPupil_1.edge_max_no_edges = 1000;//int
	orp_BrightPupil_1.edge_mindir = 50.0f;//float --  = 100.0 //average gradient should have significant rt or left component for right & lft edges(not straight up or down).
	orp_BrightPupil_1.edge_maxdir = 350.0f;//float -- = 300.0 //(0 is straight dwn, -400 is straight up, 200 is rt, -200 is left)
	orp_BrightPupil_1.edge_minedgelength = 10;//int-- = 15 // shortest string length that can be an "edge"
	orp_BrightPupil_1.edge_maxedgelength = 500;//int -- longest string liength that can be an "edge"
	orp_BrightPupil_1.edge_max_delta_curvature = 30.0f;//float -- = 15.0 //we want consistant curvature (rate of change of direction); ave changes in curvature (delta curvature) should be small
	orp_BrightPupil_1.edge_mincurv = 2.0;//float -- = 2;//previous .01 //we only want to consider strings with reasonable ave curvature
	orp_BrightPupil_1.edge_maxcurv = 10.0;//float -- = 8;//previous .05 //we don't want straight lines or curves that are too tight

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting sets of edges to form objects
	orp_BrightPupil_1.obj_max_no_objects = 50;//int
	orp_BrightPupil_1.obj_max_curv_dif = 2.0f;//float -- = 1; //max difference between left edge ave curvature and right edge ave curvature
	orp_BrightPupil_1.obj_minPD_to_width = 0.5f;//float -- = 0.5; //min for ratio: [estimate of PD from ave curvature] / [dist between lftmost and rtmost points]
	orp_BrightPupil_1.obj_maxPD_to_width = 2.5f;//float -- = 2.0; //max for same ratio
	orp_BrightPupil_1.obj_minwdth = 10;//int
	orp_BrightPupil_1.obj_maxwdth = 200;//int
	orp_BrightPupil_1.obj_minht = 10;//int
	orp_BrightPupil_1.obj_maxht = 200;//int
	orp_BrightPupil_1.obj_neutral_delta_dir = 200;//int
	orp_BrightPupil_1.obj_max_edge_end_point_separation_sqrd = 2.0;//float --
	orp_BrightPupil_1.obj_max_no_objects = 100;//int
	orp_BrightPupil_1.obj_max_edges_per_object = 2;//int

	//Pupil OBJECT_RECOGNITION_PARAMETERS -- parameters for choosing best object (object most likely to be feature of interest)
	orp_BrightPupil_1.feature_max_no_of_elps_computations = 10;//int
	orp_BrightPupil_1.feature_max_scaled_diam = 10.0f;//float --

	//CR OBJECT_RECOGNITION_PARAMETERS-- General
	orp_CR_BrightPupil_1.pupil_type = BRIGHT; //int -- dark or bright
	orp_CR_BrightPupil_1.EyeImageScaleFactor = 10.0f; //float -- pixels per mm
	orp_CR_BrightPupil_1.use_labeled_area_analysis_for_CR = true;
	orp_CR_BrightPupil_1.use_threshold_buf_for_labeled_area_anal = true;

	//CR OBJECT_RECOGNITION_PARAMETERS -- parameters for identifying strings
	orp_CR_BrightPupil_1.string_minHbreak = 1;//int
	orp_CR_BrightPupil_1.string_minVbreak = 1;//int
	orp_CR_BrightPupil_1.string_minstringlength = 2;//int --  = 5 for std bright pupil
	orp_CR_BrightPupil_1.string_max_delta_dir = 300;//int --  = 50 for std bright pupil
	orp_CR_BrightPupil_1.string_max_no_strings = 10000;//int
	orp_CR_BrightPupil_1.string_max_points_per_string = 1000;//int

	//CR OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting and combining strings to form left and right edges
	orp_CR_BrightPupil_1.edge_max_no_edges = 1000;//int
	orp_CR_BrightPupil_1.edge_mindir = 0.0;//float --  = 100.0 //average gradient should have significant rt or left component for right & lft edges(not straight up or down).
	orp_CR_BrightPupil_1.edge_maxdir = 400.0f;//float -- = 300.0 //(0 is straight dwn, -400 is straight up, 200 is rt, -200 is left)
	orp_CR_BrightPupil_1.edge_minedgelength = 2;//int-- = 15 // shortest string length that can be an "edge"
	orp_CR_BrightPupil_1.edge_maxedgelength = 20;//int -- longest string liength that can be an "edge"
	orp_CR_BrightPupil_1.edge_max_delta_curvature = 150.0f;//float -- = 15.0 //we want consistant curvature (rate of change of direction); ave changes in curvature (delta curvature) should be small
	orp_CR_BrightPupil_1.edge_mincurv = 2.0f;//float -- = 2;//previous .01 //we only want to consider strings with reasonable ave curvature
	orp_CR_BrightPupil_1.edge_maxcurv = 150.0f;//float -- = 8;//previous .05 //we don't want straight lines or curves that are too tight

	//CR OBJECT_RECOGNITION_PARAMETERS -- parameters for selecting sets of edges to form objects
	orp_CR_BrightPupil_1.obj_max_no_objects = 100;//int
	orp_CR_BrightPupil_1.obj_max_curv_dif = 1.0;//float -- = 1; //max difference between left edge ave curvature and right edge ave curvature
	orp_CR_BrightPupil_1.obj_minPD_to_width = 0.5f;//float -- = 0.5; //min for ratio: [estimate of PD from ave curvature] / [dist between lftmost and rtmost points]
	orp_CR_BrightPupil_1.obj_maxPD_to_width = 2.0f;//float -- = 2.0; //max for same ratio
	orp_CR_BrightPupil_1.obj_minwdth = 2;//int
	orp_CR_BrightPupil_1.obj_maxwdth = 10;//int
	orp_CR_BrightPupil_1.obj_minht = 2;//int
	orp_CR_BrightPupil_1.obj_maxht = 10;//int
	orp_CR_BrightPupil_1.obj_neutral_delta_dir = 200;//int
	orp_CR_BrightPupil_1.obj_max_edge_end_point_separation_sqrd = 16.0f;//float --
	orp_CR_BrightPupil_1.obj_max_no_objects = 100;//int
	orp_CR_BrightPupil_1.obj_max_edges_per_object = 5;//int

	//CR OBJECT_RECOGNITION_PARAMETERS -- parameters for choosing best object (object most likely to be feature of interest)
	orp_CR_BrightPupil_1.feature_max_no_of_elps_computations = 5;//int
	orp_CR_BrightPupil_1.feature_max_scaled_diam = 1.0f;//float --
	orp_CR_BrightPupil_1.elps_aspect_max = 2.0f; //elipse majaxis/minoraxis.

	//CR OBJECT_RECOGNITION_PARAMETERS -- parametes for labeled region analysis
	orp_CR_BrightPupil_1.max_labels = 100; // must not exceed number of labels in "EquivTable" struct definition (in ProcessImage.h)
	orp_CR_BrightPupil_1.max_no_pnts_per_labeled_area = 1000; // must not exceed number of points in "LabeledAReaBuffer" struct definition (in ProcessImage.h)
	orp_CR_BrightPupil_1.max_no_pnts_per_labeled_area_edge = 100;

	//CR OBJECT_RECOGNITION_PARAMETERS -- parameters for choosing best object from labeled region data
	orp_CR_BrightPupil_1.area_max_area_size = 3.0f;// mm^2
	orp_CR_BrightPupil_1.area_min_area_size = 0.03f; // mm^2
	orp_CR_BrightPupil_1.area_rough_aspect_max = 2.0; //aspect ratio based on top most,btm most, rt most and lft most points
	orp_CR_BrightPupil_1.area_maj_axis_max = 2.0f; //mm
	orp_CR_BrightPupil_1.area_maj_axis_min = 0.2f; //mm

	//Fuzzy Logic Parameters
	Pfuzz_BrightPupil_1.Smoothness.importance = 1;
	Pfuzz_BrightPupil_1.Smoothness.binSize = 2000 / (FUZZ_BINS - 1); // largest smoothness value is probaly 2.0; we will multiply by 1000 before converting to int
	Pfuzz_BrightPupil_1.Smoothness.factor[0] = FUZZ_MAX;        //( FUZZ_MAX, FUZZ_MAX / 2, FUZZ_MAX / 15, FUZZ_MAX / 30, FUZZ_MAX / 127 );//best (highest) fuzz is lowest "smoothness" value
	Pfuzz_BrightPupil_1.Smoothness.factor[1] = FUZZ_MAX / 2;
	Pfuzz_BrightPupil_1.Smoothness.factor[2] = FUZZ_MAX / 15;
	Pfuzz_BrightPupil_1.Smoothness.factor[3] = FUZZ_MAX / 30;
	Pfuzz_BrightPupil_1.Smoothness.factor[4] = FUZZ_MAX / 127;
	Pfuzz_BrightPupil_1.Aspect.importance = 1;
	Pfuzz_BrightPupil_1.Aspect.binSize = 100 / (FUZZ_BINS - 1); // "aspect" will be minor_axis/major_axis * 100.  
	Pfuzz_BrightPupil_1.Aspect.factor[0] =  0;                  //{ 0, FUZZ_MAX / 20, FUZZ_MAX / 10, FUZZ_MAX/5, FUZZ_MAX };
	Pfuzz_BrightPupil_1.Aspect.factor[1] = FUZZ_MAX / 20;
	Pfuzz_BrightPupil_1.Aspect.factor[2] = FUZZ_MAX / 10; 
	Pfuzz_BrightPupil_1.Aspect.factor[3] = FUZZ_MAX/5;
	Pfuzz_BrightPupil_1.Aspect.factor[4] = FUZZ_MAX;
	Pfuzz_BrightPupil_1.PD.importance = 1;
	Pfuzz_BrightPupil_1.PD.binSize = 1000 / (FUZZ_BINS - 1); // PD will be pupil diam in mm * 100.  
	Pfuzz_BrightPupil_1.PD.factor[0] = FUZZ_MAX / 20;        //{ FUZZ_MAX / 20, FUZZ_MAX, FUZZ_MAX, FUZZ_MAX, FUZZ_MAX / 20 };
	Pfuzz_BrightPupil_1.PD.factor[1] = FUZZ_MAX;
	Pfuzz_BrightPupil_1.PD.factor[2] = FUZZ_MAX;
	Pfuzz_BrightPupil_1.PD.factor[3] = FUZZ_MAX;
	Pfuzz_BrightPupil_1.PD.factor[4] = FUZZ_MAX / 20;
	//********************End BRIGHT PUPIL 1 parameters **********************************

	//************Convolution Kernels ****************************************************
	//Gausian 5x5 smoothing kernel
	K_Gsmooth_5x5[0][0] = K_Gsmooth_5x5[4][0] = 1;
	K_Gsmooth_5x5[0][1] = K_Gsmooth_5x5[4][1] = 4;
	K_Gsmooth_5x5[0][2] = K_Gsmooth_5x5[4][2] = 7;
	K_Gsmooth_5x5[0][3] = K_Gsmooth_5x5[4][3] = 4;
	K_Gsmooth_5x5[0][4] = K_Gsmooth_5x5[4][4] = 1;
	K_Gsmooth_5x5[1][0] = K_Gsmooth_5x5[3][0] = 4;
	K_Gsmooth_5x5[1][1] = K_Gsmooth_5x5[3][1] = 16;
	K_Gsmooth_5x5[1][2] = K_Gsmooth_5x5[3][2] = 26;
	K_Gsmooth_5x5[1][3] = K_Gsmooth_5x5[3][3] = 16;
	K_Gsmooth_5x5[1][4] = K_Gsmooth_5x5[3][4] = 4;
	K_Gsmooth_5x5[2][0] = 7;
	K_Gsmooth_5x5[2][1] = 26;
	K_Gsmooth_5x5[2][2] = 41;
	K_Gsmooth_5x5[2][3] = 26;
	K_Gsmooth_5x5[2][4] = 7;

	//Sobel 3x3 Horizontal gradient kernel
	K_SobelH_3x3[0][0] = K_SobelH_3x3[0][2] = -1;
	K_SobelH_3x3[1][0] = K_SobelH_3x3[1][2] = 0;
	K_SobelH_3x3[2][0] = K_SobelH_3x3[2][2] = 1;
	K_SobelH_3x3[0][1] = -2;
	K_SobelH_3x3[1][1] = 0;
	K_SobelH_3x3[2][1] = 2;

	//Sobel 3x3 Vertical gradient kernel
	K_SobelV_3x3[0][0] = K_SobelV_3x3[2][0] = -1;
	K_SobelV_3x3[0][1] = K_SobelV_3x3[2][1] = 0;
	K_SobelV_3x3[0][2] = K_SobelV_3x3[2][2] = 1;
	K_SobelV_3x3[1][0] = -2;
	K_SobelV_3x3[1][1] = 0;
	K_SobelV_3x3[1][2] = 2;

	//Gaussian 5x5 Horizontal gradient kernel
	K_GgradH_5x5[0][0] = K_GgradH_5x5[0][4] = -15;
	K_GgradH_5x5[1][0] = K_GgradH_5x5[1][4] = -35;
	K_GgradH_5x5[2][0] = K_GgradH_5x5[2][4] = 0;
	K_GgradH_5x5[3][0] = K_GgradH_5x5[3][4] = 35;
	K_GgradH_5x5[4][0] = K_GgradH_5x5[4][4] = 15;
	K_GgradH_5x5[0][1] = K_GgradH_5x5[0][3] = -69;
	K_GgradH_5x5[1][1] = K_GgradH_5x5[1][3] = -155;
	K_GgradH_5x5[2][1] = K_GgradH_5x5[2][3] = 0;
	K_GgradH_5x5[3][1] = K_GgradH_5x5[3][3] = 155;
	K_GgradH_5x5[4][1] = K_GgradH_5x5[4][3] = 69;
	K_GgradH_5x5[0][2] = -114;
	K_GgradH_5x5[1][2] = -255;
	K_GgradH_5x5[2][2] = 0;
	K_GgradH_5x5[3][2] = 255;
	K_GgradH_5x5[4][2] = 114;

	//Gaussian 5x5 Vertical gradient kernel
	K_GgradV_5x5[0][0] = K_GgradV_5x5[4][0] = -15;
	K_GgradV_5x5[0][1] = K_GgradV_5x5[4][1] = -35;
	K_GgradV_5x5[0][2] = K_GgradV_5x5[4][2] = 0;
	K_GgradV_5x5[0][3] = K_GgradV_5x5[4][3] = 35;
	K_GgradV_5x5[0][4] = K_GgradV_5x5[4][4] = 15;
	K_GgradV_5x5[1][0] = K_GgradV_5x5[3][0] = -69;
	K_GgradV_5x5[1][1] = K_GgradV_5x5[3][1] = -155;
	K_GgradV_5x5[1][2] = K_GgradV_5x5[3][2] = 0;
	K_GgradV_5x5[1][3] = K_GgradV_5x5[3][3] = 155;
	K_GgradV_5x5[1][4] = K_GgradV_5x5[3][4] = 69;
	K_GgradV_5x5[2][0] = -114;
	K_GgradV_5x5[2][1] = -255;
	K_GgradV_5x5[2][2] = 0;
	K_GgradV_5x5[2][3] = 255;
	K_GgradV_5x5[2][4] = 114;

	//************End Convolution Kernels ****************************************************


	//******* set parameter structures ****************
#ifdef CONV_5BY5
	*edp_fpga = edp_fpga_BrightPupil_1;//use _2 for 3x3 Sobel
	*edp_dsp = edp_dsp_BrightPupil_1;//use _2 for 3x3 Sobel
	*orp_P = orp_BrightPupil_1;//use _2 for 3x3 Sobel
#else
	*edp_fpga = edp_fpga_BrightPupil_2;//use _2 for 3x3 Sobel
	*edp_dsp = edp_dsp_BrightPupil_2;//use _2 for 3x3 Sobel
	*orp_P = orp_BrightPupil_2;//use _2 for 3x3 Sobel
#endif
	*orp_CR =  orp_CR_BrightPupil_1;
	*Pfuzz =  Pfuzz_BrightPupil_1;

	//****** set convolution Kernal arrays ************
#ifdef CONV_5BY5
	//******5x5 grad *****************
	*SkernSize = 5; 
	*GkernSize = 5;
	for( int ii = 0; ii < *SkernSize; ii++)
		for( int jj = 0; jj < *SkernSize; jj++)
			Skern[ii][jj] = K_Gsmooth_5x5[ii][jj];
	for( int ii = 0; ii < *GkernSize; ii++)
		for( int jj = 0; jj < *GkernSize; jj++)
		{
			GkernH[ii][jj] = K_GgradH_5x5[ii][jj];
			GkernV[ii][jj] = K_GgradV_5x5[ii][jj];
		}
#else
		//******3x3 grad *****************
		*SkernSize = 5;
		*GkernSize = 3;
		for( int ii = 0; ii < *SkernSize; ii++)
			for( int jj = 0; jj < *SkernSize; jj++)
				Skern[ii][jj] = K_Gsmooth_5x5[ii][jj];
		for( int ii = 0; ii < *GkernSize; ii++)
			for( int jj = 0; jj < *GkernSize; jj++)
			{
				GkernH[ii][jj] = K_SobelH_3x3[ii][jj];
				GkernV[ii][jj] = K_SobelV_3x3[ii][jj];
			}
#endif
}

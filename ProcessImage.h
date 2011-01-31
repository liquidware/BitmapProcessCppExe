#ifndef ProcessImage_h
#define ProcessImage_h

#pragma once

#define DARK 0
#define BRIGHT 1
#define FUZZ_BINS  5
#define FUZZ_MAX  255
#define INT_MAX 2147483647
#define MAX_NO_ELLIPSE_PLOT_PNTS 2000

#define MAX_I_BUF_SIZE 100000
#define MAX_NO_P_EDGES 1000
#define MAX_NO_CR_EDGES 1000
#define MAX_NO_P_OBJECTS 50
#define MAX_NO_CR_OBJECTS 100
#define MAX_NO_P_STRINGS 10000
#define MAX_NO_CR_STRINGS 10000
#define MAX_PTS_PER_STRING 1000
#define MAX_STRINGS_PER_EDGE 10
#define MAX_EDGES_PER_OBJECT 10
#define MAX_LABELS 100;
#define MAX_NO_PNTS_PER_AREA 1000;
#define MAX_NO_PNTS_PER_CR_EDGE 100;


typedef struct
{
	int pupil_type; //dark or bright
	bool smooth;    //use smoothing convolution if true
	int extra_bndry; // if smooth==true, grad convolution must leave extra (k-1)/2 line boundary from image boundary, where k is the order of the smoothing conv (to avoid detecting edges at image borders).
	float SmthCnvl_Scale;
	int grad_thresh; //thresh after computing gradient (std = 1000; older std =1500)
	int grayhist_size; //no elements in grayscale level histogram
	int grad_hist_size; //largest expected gradient value, usually size of largest signed integer (32767)
	int max_cr_threshbuf_size;
	int grad_thresh_for_dir_array; //usually zero (no thresh used)
	int max_pbuf_size;
	int max_CRbuf_size;
	float GradCnvl_Scale; // usually 0.125 for 5x5 convl
	bool use_thresh_hyst; //if false, don't do dual thresh hysterisis, just use lower threshold boundary.
	int hair_width; //width of "thin dark line" in pixels
	int wL;  //Window Left edge [must be >=0](look for features only subset of camera field of view defined by window) 
	int wR;  //window right edge [must be <= image width-1]
	int wT;  //window top [must be >= 0]
	int wB;  //window bottom [ must be <= image height-1]
}EDGE_DETECTION_PARAM_FPGA;

typedef struct
{
	int pupil_type; //dark or bright
	float CRthresh_hist_fraction; //use to compute CR thresh from histogram (1st pixel brighter than this percent of all pixels)(std = 0.9995) [not currently used]
	int min_pup_pix; //The brightest intensity having at least this many pixels is probably just brighter than the pupil (std = 100)[not currently used].
	int MaxGrayLimit; //max value used as "brighteset pixel" in CR thresh computation (std = 225)
	float frac_dist_brtpup_to_brightest_pixel__CRthresh; //Fraction of dist between probable pupil brightness and brightest pixel. (std = 0.7)
	float high_thresh_fraction; //high thresh (for dual thresh) will be value such that at least this fraction of pixels has lower gradient. (std = 0.5)
	float low_thresh_fraction; //low threshold will be this fraction of high threshold. (std = 0.33)
	int CRthresh_dist; //Require a pixel this far in the grad direction to be above CRthresh for CR edge point. CRthresh is determined from hist. (std = 1)
}EDGE_DETECTION_PARAM_DSP;

typedef struct
{
	int pupil_type; //dark or bright
	float EyeImageScaleFactor; //pixels per mm
	bool use_labeled_area_analysis_for_CR;
	bool use_threshold_buf_for_labeled_area_anal;

	//parameters for identifying strings
	int string_minHbreak;
	int string_minVbreak;
	int string_minstringlength; // = 5 for std bright pupil
	int string_max_delta_dir;// = 50 for std bright pupil
	int string_max_no_strings;
	int string_max_points_per_string;
	int string_HairWidth; //=3 for std bright pupil. If an edge point seems likely to be a point on one edge of a thin hair (thin dark line), don't include it in a string. Value is currently in pixels.  Should be a real measurement unit scaled by PD.

	//parameters for selecting and combining strings to form left and right edges
	int edge_max_strings_per_edge; //currently 1
	int edge_max_no_edges;
	float edge_mindir;// = 100.0 //average gradient should have significant rt or left component for right & lft edges(not straight up or down).
	float edge_maxdir;// = 300.0 //(0 is straight dwn, -400 is straight up, 200 is rt, -200 is left)
	int edge_minedgelength;// = 15 // shortest string length that can be an "edge"
	int edge_maxedgelength;// longest string liength that can be an "edge"
	float edge_max_delta_curvature;// = 15.0 //we want consistant curvature (rate of change of direction); ave changes in curvature (delta curvature) should be small
	float edge_mincurv;// = 2;//previous .01 //we only want to consider strings with reasonable ave curvature
	float edge_maxcurv;// = 8;//previous .05 //we don't want straight lines or curves that are too tight

	//parameters for selecting sets of edges to form objects
	int obj_max_no_objects;
	float obj_max_curv_dif;// = 1; //max difference between left edge ave curvature and right edge ave curvature
	float obj_minPD_to_width;// = 0.5; //min for ratio: [estimate of PD from ave curvature] / [dist between lftmost and rtmost points]
	float obj_maxPD_to_width;// = 2.0; //max for same ratio
	int obj_minwdth;
	int obj_maxwdth;
	int obj_minht;
	int obj_maxht;
	int obj_neutral_delta_dir;
	float obj_max_edge_end_point_separation_sqrd;
	int obj_max_edges_per_object;

	//parameters for choosing acceptable objects and best object (object most likely to be feature of interest)
	int feature_max_no_of_elps_computations;
	int feature_max_no_of_points_in_elps;
	float feature_max_scaled_diam;
	float elps_aspect_max; //elipse majaxis/minoraxis.

	//parmeters for labeled region analysis
	int max_labels; 
	int max_no_pnts_per_labeled_area;
	int max_no_pnts_per_labeled_area_edge;

	//parameters for choosing best object from labeled region data
	float area_max_area_size; // mm^2
	float area_min_area_size; // mm^2
	float area_rough_aspect_max; //aspect ratio based on top most,btm most, rt most and lft most points
	float area_maj_axis_max; //mm
	float area_maj_axis_min; //mm

} OBJECT_RECOGNITION_PARAMETERS;

typedef struct
{
	int X;
	int Y;
}PNT;

typedef struct
{
	int X;
	int Y;
	int value;
	int dir;
}DVBUF;

typedef struct
{
	int X;
	int Y;
	int label;
}LBLD_AREA_BUF;

typedef struct 
{
	//int *edge_buf_index;
	int edge_buf_index[MAX_PTS_PER_STRING];
	bool opened;
	int pixel_lngth; //number of pixels in edge string.
	int hmin;
	int hmax;
	int top;
	int btm;
	int startX;
	int startY;
	int endX;
	int endY;
	int deltaX;
	float ave_dir;    // ave gradient direction
	float line_lngth; // length of line rather than number of pixels (account for greater distance between diagonal pixels)     
	float smoothness; // ave delta X
	float curvature;  // ave rate of change of gradient direction
	float delta_curvature; //ave change in curvature (should be small for smooth line)
	float delta2_curvature; //2nd derivative curvature (large if curvature direction keeps changing) 
	float sum_curvature;
	float sum_delta_curvature;
	float sum_delta2_curvature;
} EDGE_STRING;

typedef struct 
{
	bool part_of_string;
	bool branch_node;
	int strng_no;
} EDGE_POINT_STATUS;

typedef struct 
{
	int no_of_strngs;
	//int *strng;
	int strng[MAX_STRINGS_PER_EDGE];
	int lftmost;
	int rtmost;
	int top;
	int btm;
	int startX;
	int startY;
	int endX;
	int endY;
	float ave_dir;
	float curv;
	float smoothness;
} LFT_EDGE;

typedef struct 
{
	int no_of_strngs;
	//int *strng;
	int strng[MAX_STRINGS_PER_EDGE];
	int lftmost;
	int rtmost;
	int top;
	int btm;
	int startX;
	int startY;
	int endX;
	int endY;
	float ave_dir;
	float curv;
	float smoothness;
} RT_EDGE;

typedef struct 
{
	int lft_edge;
	int rt_edge;
	float max_width;
	float max_height;
	bool elps_fit;
	float elps_hc;		// center coordinate
	float elps_vc;
	float elps_maj_rad; // radius (major axis)
	float elps_min_rad;	// radius (minor axis)
	float elps_angle;	// ellipse angle
	float elps_aspect;
	float elps_fit_err;
	float scaled_diam;
	float lft_edge_curv;
	float rt_edge_curv;
	float lft_smoothness;
	float rt_smoothness;
	float smoothness;
	float delta;
	int fuzzSmoothness;
	int fuzzAspect;
	int fuzzPD;
	int total_fuzz;
} OBJECT;

typedef struct
{
	int no_lft_edges;
	int no_rt_edges;
	//int *lft_edge;
	//int *rt_edge;
	int lft_edge[MAX_EDGES_PER_OBJECT];
	int rt_edge[MAX_EDGES_PER_OBJECT];
	bool elps_fit;
	float hc;			// horiz center
	float vc;
	float elps_maj_rad;
	float elps_min_rad;
	float elps_angle;
	float elps_aspect;
	float elps_fit_err;
	float scaled_diam;
	int fuzzAspect;
	int fuzzPD;
	int total_fuzz;
} CR_OBJECT;

typedef struct 
{
	PNT edgepoint[500];
} STRING_ARRAY;

typedef struct 
{
	int dir[500];
} STRNG_DIR_ARRAY;

typedef struct 
{
	int importance;          /* a number 0-15 signifying relative strength of result */
	int binSize;             /* size of factor for each bin (total range == binSize * FUZZ_BINS) */
	int factor[5];            /* collection of bin factors */
} tFuzz;

typedef struct
{
	tFuzz Smoothness;
	tFuzz Aspect;
	tFuzz PD;
}FuzzParam;

typedef struct 
{
	PNT pnt[1000]; //MAX_NO_PNTS_PER_AREA
	int no_of_pnts;
}LABELED_AREA;

typedef struct 
{
	bool label[100];
}EquivTable;

typedef struct 
{
	PNT edgepnt[100];
	int no_edge_pnts;
	PNT elpspnt[100];
	int no_elps_pnts;
	int area;//no pixels
	int lft;
	int rt;
	int top;
	int btm;
	double rough_aspect;
	bool elps_fit;
	float hc;			// center horiz
	float vc;			// center vert
	float elps_maj_rad;	// radius (major axis)
	float elps_min_rad;	// radius (minor axis)
	float elps_angle;	// ellipse angle
	float elps_aspect;
	float elps_fit_err;
	float scaled_diam;
	int fuzzAspect;
	int fuzzPD;
	int total_fuzz;
}AREA;

typedef struct 
{
	int DirBuf8[800000];
	int Ismooth_array[800000];
	int grad_array[800000];
	int dir_array_8[800000];
	int dir_array_400[800000];
	int thresh_grad_array[800000];
	int thin_array[800000];
	int dual_thresh_array[800000];
	int grayhist[256];
	int ghist[32767];

	DVBUF dvbuf[100000];
	DVBUF dvbuf2[100000];
} CANNY_EDGE_DETECT_DATA;

//////////////////////////////////////////////////////////////////////////

void DoEyeImageProcessing
( int *Isizeh, int *Isizev, int Iarray[], int *SkernSize, int Skern[5][5], int *GkernSize, int GkernH[5][5], int GkernV[5][5], 
 EDGE_DETECTION_PARAM_FPGA *edp_fpga, EDGE_DETECTION_PARAM_DSP *edp_dsp, OBJECT_RECOGNITION_PARAMETERS *orp_P, OBJECT_RECOGNITION_PARAMETERS *orp_CR, FuzzParam *Pfuzz,
 int *nstring, int *no_lft_edges, int *no_rt_edges, int *no_of_objects, int *nCRstring, int *no_lft_CR_edges, int *no_rt_CR_edges, int *no_of_CR_objects, 
 int *pupil_object_no, int *CR_object_no, bool *pupil_found, bool *CR_found,
 PNT Pbuf[], PNT CRbuf[], EDGE_POINT_STATUS edgepoint[], EDGE_POINT_STATUS CR_edgepoint[], int *Pbufsize, int *CRbufsize, 
 LFT_EDGE lft_edge[], LFT_EDGE lft_CR_edge[], RT_EDGE rt_edge[], RT_EDGE rt_CR_edge[], OBJECT objt[], CR_OBJECT CR_objt[], EDGE_STRING es[], EDGE_STRING CR_es[],
 int PDirBuf[], int CrDirBuf[], int grad_H_array[], int grad_V_array[], int CRarray[], PNT pupil_ellipse_points[], int *no_p_elps_pnts, PNT CR_ellipse_points[], int *no_CR_elps_pnts,
 LBLD_AREA_BUF CRthreshBuf[], int *cr_threshbuf_size, AREA crArea[], int *no_cr_areas, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int labld_array[], int validAreaCRs[], int *noValidCRAreas,
 CANNY_EDGE_DETECT_DATA* data );

void GetParameters
(int *SkernSize, int *GkernSize, int Skern[5][5], int GkernH[5][5], int GkernV[5][5], EDGE_DETECTION_PARAM_FPGA *edp_fpga, EDGE_DETECTION_PARAM_DSP *edp_dsp, 
 OBJECT_RECOGNITION_PARAMETERS *orp_P, OBJECT_RECOGNITION_PARAMETERS *orp_CR, FuzzParam *Pfuzz);

// top level function
typedef struct
{
	bool found;			// object has been identified succesfully
	float center_horz;	// center coordinate horiz
	float center_vert;	// center coordinate vert
	float major_radius;	// radius (major axis)
	float minor_radius;	// radius (minor axis)
	float elps_angle;	// ellipse angle
	PNT ellipse_points[MAX_NO_ELLIPSE_PLOT_PNTS];
	int no_ellipse_points;
} PUPIL_CR_OBJECT;

#define MAX_CR_OBJECTS  5
// note: declare PUPIL_CR_OBJECT CR[MAX_CR_OBJECTS] in the calling routine
void FindPupilCR (int Isizeh, int Isizev, int Iarray[], PUPIL_CR_OBJECT* pupil,	PUPIL_CR_OBJECT CR[], int* number_CR_objects);

#endif

#pragma once

class ImageProcessingUtilities_DSP
{
public:
	ImageProcessingUtilities_DSP(void);
	~ImageProcessingUtilities_DSP(void);
	static void ComputeThreshFromHist
		(  int hist[], int *hist_size, float *thresh_hist_fraction, int *min_pup_pix, int *max_mag, int *cmx, int *xpeak, int *ypeak, int *brt_pup_bndry, int *thresh);
	static void DoThresholdHysteresis2(DVBUF dvbuf[], int *dvbuf_size);
	static void MakeEdgeBuf_from_DVBUF( DVBUF dvbuf[], int *dvbuf_size, bool *use_Thresh_Hyst, int *maxsize, int *EdgeBuf_size, PNT EdgeBuf[] );
	static void GrayAreaEdgePoints4
			( PNT EdgeBuf[], int *EdgeBufSize, int dir8[], int dir400[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB, LBLD_AREA_BUF thresh_buf[], int *thresh_buff_size, 
			  int *dist, int *featuretype, int *FeatureBuf_size, PNT FeatureBuf[], int FeatureDirBuf[] );
	static void FindStrings
			( OBJECT_RECOGNITION_PARAMETERS *p, PNT EdgeBuf[], EDGE_POINT_STATUS edge_point[], int *edgebuf_size, int DirBuf[], int *nstring, EDGE_STRING *es);
	static void MakeEdges
			( OBJECT_RECOGNITION_PARAMETERS *p, EDGE_STRING strng[], int *no_of_strngs, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges, int *no_rt_edges );
	static void MakeObjects
            ( OBJECT_RECOGNITION_PARAMETERS *p, OBJECT objt[], int *no_of_objects, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges, 
			  int *no_rt_edges, EDGE_STRING strng[], int *nstring, PNT EdgeBuf[], int *EdgeBufSize );
	static void MakeCRObjects
            (OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT objt[], int *no_of_objects, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges, int *no_rt_edges);
	static bool FitPupilEllipse
			(OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p);
	static bool FitCREllipse(CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p);
	static bool FindBestPupilObject(OBJECT objt[], int *no_of_objects, int *best_objt, FuzzParam *Pfuzz);
	static void plot_ellipse(PNT ellipse_point[], int *no_list_elements, float *xc, float *yc, float *a, float *b, float *angle);
	static float AbsFloat(float f);
	static void FindCRsWithLabeledRegionAnalysis( OBJECT_RECOGNITION_PARAMETERS *orp_CR, int Array[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB,
		                                          LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int lbld_array[], AREA crArea[], int *no_cr_areas, int validAreaCRs[], 
												  int *no_validAreaCRs );
	static void FindCRsWithLabeledRegionAnalysis2( OBJECT_RECOGNITION_PARAMETERS *orp_CR, LBLD_AREA_BUF threshbuf[], int *threshbuf_size,  
                                                   LABELED_AREA labeled_area[], EquivTable EquivToLabel[], 
												   AREA crArea[], int *no_cr_areas, int validAreaCRs[], int *no_validAreaCRs  );
	static void FitElipseToLabeledAreaCRs( OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, AREA crArea[], int validAreaCRs[], int *no_validAreaCRs );
	static void FitElipseToEdgeObjectCRs( CR_OBJECT CR_objt[], OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, LFT_EDGE lft_CR_edge[], RT_EDGE rt_CR_edge[],
		                                                         EDGE_STRING CR_es[], PNT CRbuf[], int *CRbuf_size, int *no_of_CR_objects, int validCRs[], int *noValidCRs );



private:
	static void FollowEdge2(DVBUF dvbuf[], int *i, int *dvbuf_size, int *di, int *dj);
	static bool SearchBuff2(DVBUF dvbuf[], int *dvbuf_size, int *i, int *col, int *row, int *dir, int *j);
	static bool SearchBuff(LBLD_AREA_BUF buff[], int *buffsize, int *colno, int *rowno);
	static void StartString
            (int *strt_pnt_index, int *adjacent_pnt_index, int *last_strng_no, EDGE_STRING *es,
             EDGE_POINT_STATUS edge_point[], PNT EdgeBuf[], int *edgebuf_size, int dir[], int *pupiltype);
	static void FollowString
			(OBJECT_RECOGNITION_PARAMETERS *p, PNT EdgeBuf[], int *edgebuf_size, int *strng_no, EDGE_STRING *es, EDGE_POINT_STATUS edge_point[], int dir[]);
	static void AddToString(int *edge_pnt_index, int *strng_no, PNT EdgeBuf[], EDGE_STRING es[], EDGE_POINT_STATUS edge_point[], int dir[], int *pupiltype);
	static void AbortCurrentString(int *last_strng_no, EDGE_STRING es[], EDGE_POINT_STATUS edge_point[]);
	static void EndString
            (int *end_pnt_index, int *strng_no, EDGE_STRING es[], PNT EdgeBuf[], EDGE_POINT_STATUS edge_point[], int *pupiltype);
	static bool CloseOppositeEdge( int *ipnt, PNT EdgeBuf[], int *edgebuf_size, int DirBuf[], int *HairWidth );
	static int Convert800LevelDirTo8LevelDir(int *I800LevelDir);
	static bool fits_CR_object(OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *edge_no, int side);
	static float closest_end_points(int *Xstrt1, int *Ystrt1, int *Xend1, int *Yend1, int *Xstrt2, int *Ystrt2, int *Xend2, int *Yend2);
	static void MakePupilObjectEdgePointList
            ( OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p,
				PNT objt_point[], int *no_points_in_list );
              //PNT **objt_point, int *no_points_in_list );
	static void MakeCRObjectEdgePointList
            ( CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p,
				PNT objt_point[], int *no_points_in_list );
              //PNT **objt_point, int *no_points_in_list );
	static bool ComputeEllipse(PNT edge_point[], int *no_of_points, float *xc, float *yc, float *major, float *minor, float *alpha_major);
	static int pupil_fuzz(OBJECT *objt, FuzzParam *Pfuzz);
	static int calc_fuzz( int val, tFuzz *f);
	static void startCRobject(OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT *objt, int *i, int side);
	static void add_to_CR_object(CR_OBJECT *objt, int *edge_no, int side);
	static void LabeledRegions_from_array(int array[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB, int *max_labels, 
		                                  int *max_no_pnts_per_area, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int lbld_array[], int *no_of_regions, int *no_labels );
	static void LabeledRegions_from_thresh_buf(LBLD_AREA_BUF threshbuf[], int *threshbuf_size, int *max_labels, 
		                                       int *max_no_pnts_per_area, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int *no_of_regions, int *no_labels );

	static void MakeEdgePointListFromLabeledArea( LABELED_AREA *L, PNT list[], int *max_no_pnts_per_cr_edge, int *no_pnts_in_list, int *leftmost, int *rtmost, int *top, int *btm );
	static void FindValidCRsFromAreas(OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, AREA crArea[], int *no_CR_areas, int *max_labels, int validAreaCRs[], int *no_validAreaCRs);
	static bool AcceptableCR( CR_OBJECT *CR_obj, OBJECT_RECOGNITION_PARAMETERS *CR_rec_param );
};
#define LEFT 0
#define RIGHT 1

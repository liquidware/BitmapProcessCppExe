#include <stdio.h>
#include <string.h>
#include "Clock.h"
#include "ProcessImage.h"
#include "ImageProcessingUtilities_FPGA.h"
#include "ImageProcessingUtilities_DSP.h"
#define _USE_MATH_DEFINES // for C++
#include <math.h>
ImageProcessingUtilities_DSP::ImageProcessingUtilities_DSP(void)
{
}

ImageProcessingUtilities_DSP::~ImageProcessingUtilities_DSP(void)
{
}
    void ImageProcessingUtilities_DSP::ComputeThreshFromHist
	(  int hist[], int *hist_size, float *thresh_hist_fraction, int *min_pup_pix, int *max_mag, int *cmx, int *xpeak, int *ypeak, int *brt_pup_bndry, int *thresh)
    {
        //Outputs are max_mag, cmx, xpeak, ypeak, brt_pup_bndry, and thresh
        // count pixels with non zero values
        int npixels = 0;
        int sum_xy = 0;
        int sum_y = 0;
		bool sstop = false;//debug
        *xpeak = 0;
        *ypeak = 0;
        *brt_pup_bndry = 0;
        
        *max_mag = 0;
        clock_start("ImageProcessingUtilities_DSP::ComputeThreshFromHist");
        for (int i = 1; i < *hist_size; i++)
        {
            if ( i == 1030 )//debug
				sstop = true;//debug
			if (hist[i] != 0) 
				*max_mag = i;//find brightest pixel value
            npixels +=  hist[i];//count non zero pixels
            sum_xy += hist[i] * i;//collect sums for "center of mass"
            sum_y += hist[i];//collect sums for "center of mass"
            if (hist[i] >= *ypeak)
            {
                *ypeak = hist[i]; //number of pixels at peak
                *xpeak = i; //hist peak (use brightest peak in case of tie)
            }
            if ( hist[i] >= *min_pup_pix ) *brt_pup_bndry = i; //bndry will be lowest intensity with at least min_pup_pix pixels
        }                                                  //bright pupil should have at least min_pup_pix at its brightest intensity 
                                                          //and most of bright pupil is probably just below this brightness boundary

        //for floating point exception, but could be expensive
        if (sum_y > 0) {
        	*cmx = sum_xy / sum_y; //x axis center of mass
        }
        //compute specified fraction of non zero pixels
        int maxpixels = (int)((double)npixels * *thresh_hist_fraction + 0.5);

        // find histogram bin (pixel value) where sum of pixels first exceed maxpixels
        int r;
        npixels = 0;
        for (r = 1; r < *max_mag; r++)
        {
            npixels += hist[r];
            if (npixels >= maxpixels) break;
        }
         *thresh = r;
        clock_end();
    }

    void  ImageProcessingUtilities_DSP::DoThresholdHysteresis2(DVBUF dvbuf[], int *dvbuf_size)
    {
        int di[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
        int dj[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
        clock_start("ImageProcessingUtilities_DSP::DoThresholdHysteresis2");
		for (int i = 0; i < *dvbuf_size; i++)
        {
            if (dvbuf[i].value == 2)
                FollowEdge2(dvbuf, &i, dvbuf_size, di, dj );
        }
		clock_end();
    }

		void ImageProcessingUtilities_DSP::FollowEdge2(DVBUF dvbuf[], int *i, int *dvbuf_size, int di[], int dj[])
        {//This is a re-entrant routine 
            if (*i < *dvbuf_size)
            {
                for (int dir = 0; dir < 8; dir++)//check 8 neighboring pixels
                {
                    int col = dvbuf[*i].X + di[dir];
                    int row = dvbuf[*i].Y + dj[dir];
                    int j;
                    if (SearchBuff2( dvbuf, dvbuf_size, i, &col, &row, &dir, &j ))
                    {
                        if (dvbuf[j].value == 1)
                        {
                            dvbuf[j].value = 2;
                            FollowEdge2( dvbuf, &j, dvbuf_size, di, dj );
                        }
                    }
                }
            }
        }

		bool ImageProcessingUtilities_DSP::SearchBuff2(DVBUF dvbuf[], int *dvbuf_size, int *i, int *col, int *row, int *dir, int *j)
        {//find buff entry with specified col and row; set index number (j) if found; return false if not found
            *j = 0;
            switch (*dir)
            {
                case 0: //find neighbor to right of pixel i (fwd search)
                    for (*j = *i + 1; *j < *dvbuf_size; (*j)++)
                    {
                        if ((dvbuf[*j].X == *col) && (dvbuf[*j].Y == *row)) return (true);
                        if ((dvbuf[*j].X > *col) || (dvbuf[*j].Y > *row)) return (false);
                    }
                    break;
                case 1://find neighbor to lower right of pixel i (fwd search)
                case 2://find neighbor below pixel i (fwd search)
                case 3://find neighbor to lower left of pixel i (fwd search)
                    for (*j = *i + 1; *j < *dvbuf_size; (*j)++)
                    {
                        if ((dvbuf[*j].X == *col) && (dvbuf[*j].Y == *row)) return (true);
                        if ((dvbuf[*j].Y > *row)) return (false);
                        if ((dvbuf[*j].Y == *row) && (dvbuf[*j].X > *col)) return (false);
                    }
                    break;
                case 4:
                    for (*j = *i - 1; *j >= 0; *j--)
                    {
                        if ((dvbuf[*j].X == *col) && (dvbuf[*j].Y == *row)) return (true);
                        if ((dvbuf[*j].X < *col) || (dvbuf[*j].Y < *row)) return (false);
                    }
                    break;
                case 5:
                case 6:
                case 7:
                    for (*j = *i - 1; *j >= 0; *j--)
                    {
                        if ((dvbuf[*j].X == *col) && (dvbuf[*j].Y == *row)) return (true);
                        if ((dvbuf[*j].Y < *row)) return (false);
                        if ((dvbuf[*j].Y == *row) && (dvbuf[*j].X < *col)) return (false);
                    }
                    break;
            }
            return (false);
        }

        void ImageProcessingUtilities_DSP::MakeEdgeBuf_from_DVBUF( DVBUF dvbuf[], int *dvbuf_size, bool *use_Thresh_Hyst, int *maxsize, int *EdgeBuf_size, PNT EdgeBuf[] )
        {
            //Output EdgeBuf_size and EdgeBuf
			int j = 0;
			clock_start("ImageProcessingUtilities_DSP::MakeEdgeBuf_from_DVBUF");
            for (int i = 0; i < *dvbuf_size; i++)
            {
                if ( i >= *maxsize ) break; //don't overrun EdgBuf[]
				if( (dvbuf[i].value == 2) || !*use_Thresh_Hyst )
                {
                    EdgeBuf[j].X = dvbuf[i].X;
                    EdgeBuf[j].Y = dvbuf[i].Y;
                    j++;
                }
            }
            *EdgeBuf_size = j;
            clock_end();
        }

        void ImageProcessingUtilities_DSP::GrayAreaEdgePoints4
			( PNT EdgeBuf[], int *EdgeBufSize, int dir8[], int dir400[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB, LBLD_AREA_BUF thresh_buf[], int *thresh_buff_size, 
			  int *dist, int *featuretype, int *FeatureBuf_size, PNT FeatureBuf[], int FeatureDirBuf[] )
        {
            //Output FeatureBuf and FeatureBuf_size
			//Add pixels from the edge buffer list to the feature buff list only if the pixel a specified distance in the gradient direction is in thresh_buff (was found to be above a threshold gray level)
            //If looking for dark object (!bright), "dist" is opposite to gradient direction and must be in the thresh_buff
            //Zero all other gradient values.
            //Direction map is assumed to be 8 point map starting with 1 = gradient increasing to right, and proceeding in clockwise rotation
            int colno, rowno;
            int diarray[9] = { 0, 1, 1, 0, -1, -1, -1, 0, 1 };
            int djarray[9] = { 0, 0, 1, 1, 1, 0, -1, -1, -1 };
            int j = 0;
            FeatureBuf[0].X = 0;
            FeatureBuf[0].Y = 0;
            clock_start("ImageProcessingUtilities_DSP::GrayAreaEdgePoints4");
            for (int i = 0; i < *EdgeBufSize; i++)
             {
                if (EdgeBuf[i].X < *dist || EdgeBuf[i].X > (*Isizeh - *dist - 1) || EdgeBuf[i].Y < *dist || EdgeBuf[i].Y > (*Isizev - *dist - 1)) continue;
                else
                {
                    int di = diarray[ dir8[ i ] ];
                    int dj = djarray[ dir8[ i ] ];
                    if( *featuretype = BRIGHT )
                    {
                        colno = EdgeBuf[i].X + di * *dist;
                        rowno = EdgeBuf[i].Y + dj * *dist;
                    }
                    else //dark
                    {
                        colno = EdgeBuf[i].X - di * *dist;
                        rowno = EdgeBuf[i].Y - dj * *dist;
                    }
					if( (colno >= *wL) && (colno <= *wR) && (rowno >= *wT) && (rowno <= *wB) )
						if (SearchBuff(thresh_buf, thresh_buff_size, &colno, &rowno)) 
						{
							FeatureBuf[j].X = EdgeBuf[i].X;
							FeatureBuf[j].Y = EdgeBuf[i].Y;
							FeatureDirBuf[j] = dir400[i];
							++j;
						}
                }
            }
            *FeatureBuf_size = j;
            clock_end();
        }

		bool ImageProcessingUtilities_DSP::SearchBuff(LBLD_AREA_BUF buff[], int *buffsize, int *colno, int *rowno)
        {
            bool point_found = false;
            for (int n = 0; n < *buffsize; n++)
            {
                if (buff[n].Y < *rowno) continue;
                else 
                    if (buff[n].Y > *rowno) 
                    {
                        point_found = false;
                        break;
                    }
                if (buff[n].Y == *rowno)
                {
                    if (buff[n].X == *colno)
                    {
                        point_found = true;
                        break;
                    }
               }
            }
            return (point_found);
        }

		void ImageProcessingUtilities_DSP::FindStrings
			(OBJECT_RECOGNITION_PARAMETERS *p, PNT EdgeBuf[], EDGE_POINT_STATUS edge_point[], int *edgebuf_size, int DirBuf[], int *nstring, EDGE_STRING *es)
        {
            //Output nstring, es, and modify edge_point as appropriate.
			// es[] is an array of edge strings, each consisting of a set of adjacent edge points
            //      specified by a list of EdgeBuf indecies. There are nstring elements in the array.
            //      A string will end at any branch, and new strings will start for each branch.
            // edge_point[] is an array corresponding to each point in EdgeBuf, and containing information 
            //      about the status of that point, such as whether it has been assigned to a string, etc. 
            // Adjustable parameters for the FindStrings function are passed to it in OBJECT_RECOGNITION_PARAMETERS structure.
			clock_start("ImageProcessingUtilities_DSP::FindStrings");
            *nstring = 0; // number of strings found
            int strng_no = *nstring - 1;
            /* cladden: unused variable removed
            int istring = 0;
			*/

			//LARGE_INTEGER t1, t2, t3;
			//double dt1, dt2, dt3, deltat;

			//QueryPerformanceFrequency(&t1);
			//QueryPerformanceCounter(&t2);

            clock_start("ImageProcessingUtilities_DSP::FindStrings:init es");
            /* cladden: loop unroll optimizations */
            //for (int i = 0; i < p->string_max_no_strings; i++)// initialize es array
            for (int i = 0; i < (p->string_max_no_strings/100); i++)// initialize es array
            {
                //es[i].edge_buf_index = new int[p->string_max_points_per_string];

                es[i].pixel_lngth = 0;
                es[i+1].pixel_lngth = 0;
                es[i+2].pixel_lngth = 0;
                es[i+3].pixel_lngth = 0;
                es[i+4].pixel_lngth = 0;
                es[i+5].pixel_lngth = 0;
                es[i+6].pixel_lngth = 0;
                es[i+7].pixel_lngth = 0;
                es[i+8].pixel_lngth = 0;
                es[i+9].pixel_lngth = 0;

                es[i+10].pixel_lngth = 0;
                es[i+11].pixel_lngth = 0;
                es[i+12].pixel_lngth = 0;
                es[i+13].pixel_lngth = 0;
                es[i+14].pixel_lngth = 0;
                es[i+15].pixel_lngth = 0;
                es[i+16].pixel_lngth = 0;
                es[i+17].pixel_lngth = 0;
                es[i+18].pixel_lngth = 0;
                es[i+19].pixel_lngth = 0;

                es[i+20].pixel_lngth = 0;
                es[i+21].pixel_lngth = 0;
                es[i+22].pixel_lngth = 0;
                es[i+23].pixel_lngth = 0;
                es[i+24].pixel_lngth = 0;
                es[i+25].pixel_lngth = 0;
                es[i+26].pixel_lngth = 0;
                es[i+27].pixel_lngth = 0;
                es[i+28].pixel_lngth = 0;
                es[i+29].pixel_lngth = 0;

                es[i+30].pixel_lngth = 0;
                es[i+31].pixel_lngth = 0;
                es[i+32].pixel_lngth = 0;
                es[i+33].pixel_lngth = 0;
                es[i+34].pixel_lngth = 0;
                es[i+35].pixel_lngth = 0;
                es[i+36].pixel_lngth = 0;
                es[i+37].pixel_lngth = 0;
                es[i+38].pixel_lngth = 0;
                es[i+39].pixel_lngth = 0;

                es[i+40].pixel_lngth = 0;
                es[i+41].pixel_lngth = 0;
                es[i+42].pixel_lngth = 0;
                es[i+43].pixel_lngth = 0;
                es[i+44].pixel_lngth = 0;
                es[i+45].pixel_lngth = 0;
                es[i+46].pixel_lngth = 0;
                es[i+47].pixel_lngth = 0;
                es[i+48].pixel_lngth = 0;
                es[i+49].pixel_lngth = 0;

                es[i+50].pixel_lngth = 0;
                es[i+51].pixel_lngth = 0;
                es[i+52].pixel_lngth = 0;
                es[i+53].pixel_lngth = 0;
                es[i+54].pixel_lngth = 0;
                es[i+55].pixel_lngth = 0;
                es[i+56].pixel_lngth = 0;
                es[i+57].pixel_lngth = 0;
                es[i+58].pixel_lngth = 0;
                es[i+59].pixel_lngth = 0;

                es[i+60].pixel_lngth = 0;
                es[i+61].pixel_lngth = 0;
                es[i+62].pixel_lngth = 0;
                es[i+63].pixel_lngth = 0;
                es[i+64].pixel_lngth = 0;
                es[i+65].pixel_lngth = 0;
                es[i+66].pixel_lngth = 0;
                es[i+67].pixel_lngth = 0;
                es[i+68].pixel_lngth = 0;
                es[i+69].pixel_lngth = 0;

                es[i+70].pixel_lngth = 0;
                es[i+71].pixel_lngth = 0;
                es[i+72].pixel_lngth = 0;
                es[i+73].pixel_lngth = 0;
                es[i+74].pixel_lngth = 0;
                es[i+75].pixel_lngth = 0;
                es[i+76].pixel_lngth = 0;
                es[i+77].pixel_lngth = 0;
                es[i+78].pixel_lngth = 0;
                es[i+79].pixel_lngth = 0;

                es[i+80].pixel_lngth = 0;
                es[i+81].pixel_lngth = 0;
                es[i+82].pixel_lngth = 0;
                es[i+83].pixel_lngth = 0;
                es[i+84].pixel_lngth = 0;
                es[i+85].pixel_lngth = 0;
                es[i+86].pixel_lngth = 0;
                es[i+87].pixel_lngth = 0;
                es[i+88].pixel_lngth = 0;
                es[i+89].pixel_lngth = 0;

                es[i+90].pixel_lngth = 0;
                es[i+91].pixel_lngth = 0;
                es[i+92].pixel_lngth = 0;
                es[i+93].pixel_lngth = 0;
                es[i+94].pixel_lngth = 0;
                es[i+95].pixel_lngth = 0;
                es[i+96].pixel_lngth = 0;
                es[i+97].pixel_lngth = 0;
                es[i+98].pixel_lngth = 0;
                es[i+99].pixel_lngth = 0;
            }
            clock_end();
            for (int i = 0; i < *edgebuf_size; i++) // initialize edge_point array
            {
                edge_point[i].part_of_string = false;
                edge_point[i].branch_node = false;
            }
			//QueryPerformanceCounter(&t2);
            for (int ipnt = 0; ipnt < *edgebuf_size; ipnt++) //scroll through all edge points
            {
                if (*nstring >= p->string_max_no_strings) break;
				//if (CloseOppositeEdge(&ipnt, EdgeBuf, edgebuf_size, DirBuf, &p->string_HairWidth))
				//{
				//	continue;//if this seems likely to be a point on one edge of a thin hair (thin dark line), don't start a string with it.
				//}            //we don't want to have to sort through eye brow or eyelash hairs when we look for the edges that define the pupil and CR
				//if( (EdgeBuf[ipnt].Y > 60)&&(EdgeBuf[ipnt].Y < 160)&&(EdgeBuf[ipnt].X > 80)&&(EdgeBuf[ipnt].X < 240)) //debug
                if (!edge_point[ipnt].part_of_string) //if edge point already part of a string, go on to next point
                {
                    //does point have at least minstring other points attached to it
                    //int[] tempStrng = new int[1000];
                    //int n_adjacent = 0;
                    for (int jpnt = ipnt + 1; jpnt < *edgebuf_size; jpnt++)
                    {
                        if ((EdgeBuf[jpnt].Y - EdgeBuf[ipnt].Y) <= p->string_minVbreak)
                        {
							if ((ImageProcessingUtilities_FPGA::AbsInt(EdgeBuf[jpnt].X - EdgeBuf[ipnt].X)) <= p->string_minHbreak)
                            {
                                if (!edge_point[jpnt].part_of_string) //if adjacent pt already part of string, don't use it again
                                {
                                    //if (!CloseOppositeEdge(&jpnt, EdgeBuf, edgebuf_size, DirBuf, &p->string_HairWidth))//if this seems likely to be a point on one edge of a thin hair, don't make it part of string.
                                    //{
                                        StartString(&ipnt, &jpnt, &strng_no, es, edge_point, EdgeBuf, edgebuf_size, DirBuf, &p->pupil_type);
                                        FollowString(p, EdgeBuf, edgebuf_size, &strng_no, es, edge_point, DirBuf);
                                        break; //don't start more than one string from same point.
                                   // }
                                }
                            }
                        }
                        else break; // we have reached a buff entry with vert coord more than 1 row beyond EdgeBuf[ipnt]; 
                    }               // no subsequent buf entries will be neigbors.
                }
            }
			//QueryPerformanceCounter(&t3);
			//dt1 = t1.QuadPart;
			//dt2 = t2.QuadPart;
			//dt3 = t3.QuadPart;
			//deltat = (dt3 - dt2) / dt1;
			//WCHAR str[80];
			//swprintf_s(str, 80, L"deltat=%f", deltat);
			//OutputDebugString(str);

            *nstring = strng_no + 1;
            clock_end();
        }

		void ImageProcessingUtilities_DSP::StartString
            (int *strt_pnt_index, int *adjacent_pnt_index, int *last_strng_no, EDGE_STRING * es,
             EDGE_POINT_STATUS edge_point[], PNT EdgeBuf[], int *edgebuf_size, int dir[], int *pupiltype)
        {
            int strng_no = *last_strng_no + 1;
            double delta;
            //int[] lookupHdir = new int[9] { -10, -2, -1, 0, 1, 2, 1, 0, -1 };//zero at top & btm, increasing to 2 at rt side & -2 at left side.
            //int[] lookupCurv = new int[9] { -10, -2, -1, 0, 1, 2, 3, 4, -3 };//zero at top & btm, inceasing from top to 3 at btm rt & -3 at btm lft.
            es[strng_no].edge_buf_index[0] = *strt_pnt_index;
            es[strng_no].edge_buf_index[1] = *adjacent_pnt_index;
            es[strng_no].pixel_lngth = 2;
            es[strng_no].opened = true;
            es[strng_no].smoothness = 0;
            es[strng_no].startX = EdgeBuf[*strt_pnt_index].X;
            es[strng_no].startY = EdgeBuf[*strt_pnt_index].Y;
            edge_point[*strt_pnt_index].part_of_string = edge_point[*adjacent_pnt_index].part_of_string = true;
            es[strng_no].deltaX = EdgeBuf[*adjacent_pnt_index].X - EdgeBuf[*strt_pnt_index].X;
			if ((ImageProcessingUtilities_FPGA::AbsInt(es[strng_no].deltaX) > 0) && (EdgeBuf[*adjacent_pnt_index].Y > EdgeBuf[*strt_pnt_index].Y)) //new point is diagonal to previous point
                es[strng_no].line_lngth = 1.4142f;
            else es[strng_no].line_lngth = 1.0f;
            if ( *pupiltype == DARK )
            {
                es[strng_no].ave_dir = (float)( dir[*strt_pnt_index] - dir[*adjacent_pnt_index] );
                delta = -dir[*adjacent_pnt_index] + dir[*strt_pnt_index];
            }
            else
            {
                es[strng_no].ave_dir = (float)( dir[*strt_pnt_index] + dir[*adjacent_pnt_index] );
                delta = dir[*adjacent_pnt_index] - dir[*strt_pnt_index];
            }
            if (EdgeBuf[*strt_pnt_index].X < EdgeBuf[*adjacent_pnt_index].X)
            {
                es[strng_no].hmin = EdgeBuf[*strt_pnt_index].X;
                es[strng_no].hmax = EdgeBuf[*adjacent_pnt_index].X;
            }
            else
            {
                es[strng_no].hmax = EdgeBuf[*strt_pnt_index].X;
                es[strng_no].hmin = EdgeBuf[*adjacent_pnt_index].X;
            }
            if (EdgeBuf[*strt_pnt_index].Y < EdgeBuf[*adjacent_pnt_index].Y)
            {
                es[strng_no].top = EdgeBuf[*strt_pnt_index].Y;
                es[strng_no].btm = EdgeBuf[*adjacent_pnt_index].Y;
            }
            else
            {
                es[strng_no].btm = EdgeBuf[*strt_pnt_index].Y;
                es[strng_no].top = EdgeBuf[*adjacent_pnt_index].Y;
            }
            es[strng_no].delta2_curvature = 0;//current delta curv = last delta curv
            es[strng_no].delta_curvature = 0;// current curvature - last curvature
            es[strng_no].curvature = (float)(delta);
            es[strng_no].sum_curvature = (float)(delta);
            es[strng_no].sum_delta_curvature = 0;
            es[strng_no].sum_delta2_curvature = 0;

            *last_strng_no = strng_no;
        }

		void ImageProcessingUtilities_DSP::FollowString
			(OBJECT_RECOGNITION_PARAMETERS *p, PNT EdgeBuf[], int *edgebuf_size, int *strng_no, EDGE_STRING *es, EDGE_POINT_STATUS edge_point[], int dir[])
        {
            const int MULTIPLE_ADJACENT_POINTS = 2;
 
            int strng_indx = es[*strng_no].pixel_lngth - 1;
            int ieb = es[*strng_no].edge_buf_index[strng_indx]; //edge_buf index of last point in string
            //int eblngth = EdgeBuf.GetLength(0);
            int temp[5];
            int dist[5];;
            int tielist[5]; 
            int jeb = 0;
            int Ydist = 0;
            int Xdist = 0;
            int keb;
            int point_with_closest_dir, closest_point, closest_dist, tieindex, dirdif, mindirdif, n;
            bool continue_strng = true;

            while (continue_strng)
            {
				if( es[*strng_no].pixel_lngth > MAX_PTS_PER_STRING )
					break;
                // find first edgebuf entry with same vert position as edgebuf[ieb]; 
                // then search forward through edge buf (remembering to skip edgebuf[ieb])
                //if (ieb > 0) //if we are not already at he beginning of the edgebuf, look  backwards
                //{
                //    for (keb = ieb - 1; keb > 0; keb--)
                //    {
                //        if (EdgeBuf[keb].X <= EdgeBuf[ieb].X)  //we found a left neighbor (if y value is same)
                //        {
                //            if (EdgeBuf[keb].Y < EdgeBuf[ieb].Y) keb++;//increment to get back to earliest entry with same Y value.
                //            break;
                //        }
                //        if ((EdgeBuf[keb].Y < EdgeBuf[ieb].Y)) // we found the last entry with smaller Y value
                //        {
                //            keb++;//increment to get back to earliest entry with same Y value.
                //            break;
                //        }
                //    }
                //}
                keb = ieb - 1; //If there is an adjacent point to the left, it must the previous buffer entry.
                //EdgeBuf[keb] is the point in the edge buf just before EdgeBuf[ieb].

                int nadj = 0;
                for (jeb = keb; jeb < *edgebuf_size; jeb++)
                { // make list of points that are adjacent to end of string (should be max of 5)
                    if (jeb == ieb) continue; //skip point that is current end of string 
                    Ydist = EdgeBuf[jeb].Y - EdgeBuf[ieb].Y;
                    Xdist = EdgeBuf[jeb].X - EdgeBuf[ieb].X;
                    if ((Ydist >= 1) && (Xdist > 1)) break;
                    if (Ydist <= 1)
                    {
						if (ImageProcessingUtilities_FPGA::AbsInt(Xdist) <= 1) //if this point is adjacent to current end-of-string, consider adding it to the string
							//if( (Ydist > 0) && (Xdist > 1) ) break;
                            if (!edge_point[jeb].part_of_string) //if the point is part of another string, don't add it to this one
								if (!(ImageProcessingUtilities_FPGA::AbsInt(dir[jeb] - dir[ieb]) > p->string_max_delta_dir)) //if this point would make too big a "kink" in the string, don't use it
                                    //if (!CloseOppositeEdge(&jeb, EdgeBuf, edgebuf_size, dir, &p->string_HairWidth))//if this seems likely to be a point on one edge of a thin hair (thin dark line), don't make it part of string.
                                    {// add this point to list of adjacent points
                                        temp[nadj] = jeb;
										dist[nadj] = ImageProcessingUtilities_FPGA::AbsInt(Xdist) + ImageProcessingUtilities_FPGA::AbsInt(Ydist);
                                        nadj++;
                                        if (nadj >= 4) break;
                                    }
                    }
                    else break; //once we have a vertical break we might as well stop (Y values only increase as we go down the EdgeBuf )
                } // Leave "for loop" when we have found all the adjacent points we can to current ieb

                int num_of_adjacent_points = nadj;
                if (nadj > 1) num_of_adjacent_points = MULTIPLE_ADJACENT_POINTS;
                switch (num_of_adjacent_points)
                {
                    case 0: // No adjacent points
                        if (es[*strng_no].pixel_lngth < p->string_minstringlength)
                            AbortCurrentString(strng_no, es, edge_point); //the existing string is too short to be used as a string -- throw it out!
                        else //keep existing string, but end it here
                            EndString(&ieb, strng_no, es, EdgeBuf, edge_point, &p->pupil_type);
                        continue_strng = false; //don't look for any more points to add to this string
                        break;

                    case 1: // we have one adjacent point, with index in temp[0], that is not already part of another string
                        AddToString(&temp[0], strng_no, EdgeBuf, es, edge_point, dir, &p->pupil_type); //add the adjacent point to string
                        ieb = temp[0]; // reset end of string index 
                        break;

                    case MULTIPLE_ADJACENT_POINTS:
                        // there is more than one adjacent point -- find closest point or point with most consistent gradient direction 
                        //first try to find the closest point
                        closest_point = temp[0]; //initialize "winner" to first point in list
                        closest_dist = dist[0];
                        tielist[0] = temp[0];
                        tieindex = 0;
                        for (n = 1; n < nadj; ++n) //check the other points in the list
                        {
                            if (dist[n] < closest_dist)
                            {
                                closest_dist = dist[n];
                                closest_point = temp[n];
                            }
                            if (dist[n] == closest_dist) // if we have any ties for the "winner", make a list of these 
                            {
                                ++tieindex;
                                tielist[tieindex] = temp[n];
                            }
                            else tieindex = 0; //if new single "winner", zero tie list.
                        }
                        if (tieindex == 0) // if we have a single closest point, add it to string 
                        {
                            AddToString(&closest_point, strng_no, EdgeBuf, es, edge_point, dir, &p->pupil_type); //add closest point to string
                            ieb = closest_point; // reset end of string index 
                        }
                        else // If two or more points are equi-distant form string end, compare gradient direction values (find smallest change from current end of string)
                        {    // It two or more points compare equally in direction difference from current end of string, just use the first one we come to.
							mindirdif = ImageProcessingUtilities_FPGA::AbsInt(dir[ieb] - dir[tielist[0]]); // direction change for first point in list
                            point_with_closest_dir = tielist[0]; // inialize "winner" to first point in list
                            for (n = 1; n < tieindex; ++n) // check the other points in the list
                            {
                                dirdif = ImageProcessingUtilities_FPGA::AbsInt(dir[ieb]) - dir[tielist[n + 1]];
                                if (dirdif < mindirdif) point_with_closest_dir = tielist[n + 1];
                            }
                            AddToString(&point_with_closest_dir, strng_no, EdgeBuf, es, edge_point, dir, &p->pupil_type); //add point with most consistant direction to string
                            ieb = point_with_closest_dir; // reset end of string index
                        }
                        break;
                }//end switch
            }//end while
        }

		void ImageProcessingUtilities_DSP::AddToString(int *edge_pnt_index, int *strng_no, PNT EdgeBuf[], EDGE_STRING es[], EDGE_POINT_STATUS edge_point[], int dir[], int *pupiltype)
        {
            int smoothness_value, edge_lngth, previous_pnt_index, newDeltaX;
            float delta, delta_curv, delta_line_lngth;
			
			smoothness_value = 0;
            edge_lngth = es[*strng_no].pixel_lngth;
            previous_pnt_index = es[*strng_no].edge_buf_index[edge_lngth - 1];

            es[*strng_no].edge_buf_index[edge_lngth] = *edge_pnt_index; // add new EdgeBuf index to string
            edge_point[*edge_pnt_index].part_of_string = true;
            es[*strng_no].pixel_lngth = edge_lngth + 1;
            if (es[*strng_no].hmin > EdgeBuf[*edge_pnt_index].X) es[*strng_no].hmin = EdgeBuf[*edge_pnt_index].X;//keep track of leftmost point
            if (es[*strng_no].hmax < EdgeBuf[*edge_pnt_index].X) es[*strng_no].hmax = EdgeBuf[*edge_pnt_index].X;//keep track of rightmost point
            if (es[*strng_no].top > EdgeBuf[*edge_pnt_index].Y) es[*strng_no].top = EdgeBuf[*edge_pnt_index].Y;//keep track of top point
            if (es[*strng_no].btm < EdgeBuf[*edge_pnt_index].Y) es[*strng_no].btm = EdgeBuf[*edge_pnt_index].Y;//keep track of btm point
            newDeltaX = EdgeBuf[*edge_pnt_index].X - EdgeBuf[previous_pnt_index].X;
            if ((ImageProcessingUtilities_FPGA::AbsInt(newDeltaX) > 0) && (EdgeBuf[*edge_pnt_index].Y > EdgeBuf[previous_pnt_index].Y)) //new point is diagonal to previous point
                delta_line_lngth = 1.4142f;
            else delta_line_lngth = 1.0f; // new point is directly below or alongside previous point
            smoothness_value = ImageProcessingUtilities_FPGA::AbsInt(newDeltaX - es[*strng_no].deltaX);
            es[*strng_no].deltaX = newDeltaX;

            //Accumulate sum for ave. Sign of ave will indicate whether edge is a probable left edge (pos) or rt edge (neg).
            if (*pupiltype == DARK)
            {
                es[*strng_no].ave_dir -= dir[*edge_pnt_index];
                delta = (float)(-dir[*edge_pnt_index] + dir[previous_pnt_index] );
            }
            else
            {
                es[*strng_no].ave_dir += dir[*edge_pnt_index];
                delta = (float)( dir[*edge_pnt_index] - dir[previous_pnt_index] );
            }
            es[*strng_no].line_lngth += delta_line_lngth;
            es[*strng_no].smoothness += smoothness_value;
            delta_curv = AbsFloat( delta - es[*strng_no].curvature );
            if (edge_lngth > 2) es[*strng_no].delta2_curvature = delta_curv - es[*strng_no].delta_curvature;//current delta curv = last delta curv
            es[*strng_no].delta_curvature = delta_curv;// current curvature - last curvature
            es[*strng_no].curvature = delta;
            es[*strng_no].sum_curvature += delta;
            es[*strng_no].sum_delta_curvature += delta_curv;
            es[*strng_no].sum_delta2_curvature += es[*strng_no].delta2_curvature;
        }


		void ImageProcessingUtilities_DSP::AbortCurrentString( int *last_strng_no, EDGE_STRING es[], EDGE_POINT_STATUS edge_point[] )
        {
            for (int i = 0; i < es[*last_strng_no].pixel_lngth; i++)
            {
                int ieb = es[*last_strng_no].edge_buf_index[i];
                edge_point[ieb].part_of_string = false;
            }
            es[*last_strng_no].opened = false;
            --(*last_strng_no);
        }

	    void ImageProcessingUtilities_DSP::EndString
            (int *end_pnt_index, int *strng_no, EDGE_STRING es[], PNT EdgeBuf[], EDGE_POINT_STATUS edge_point[], int *pupiltype)
        {
            float sum, ssum;
			es[*strng_no].opened = false;
            edge_point[*end_pnt_index].part_of_string = true;
            edge_point[*end_pnt_index].strng_no = *strng_no;
            sum = es[*strng_no].ave_dir;
            ssum = es[*strng_no].smoothness;
            es[*strng_no].endX = EdgeBuf[*end_pnt_index].X;
            es[*strng_no].endY = EdgeBuf[*end_pnt_index].Y;
            es[*strng_no].ave_dir = sum / (float)es[*strng_no].pixel_lngth; //positive means probable left edge & visa versa;
            es[*strng_no].smoothness = ssum / (float)es[*strng_no].pixel_lngth;
            es[*strng_no].curvature = es[*strng_no].sum_curvature / es[*strng_no].line_lngth; //if bright pupil, ave curvature should be pos for rt edge, neg for lft edge.
            es[*strng_no].delta_curvature = es[*strng_no].sum_delta_curvature / es[*strng_no].line_lngth;
            es[*strng_no].delta2_curvature = es[*strng_no].sum_delta2_curvature / es[*strng_no].line_lngth;
            //check for special case; if gradient points to left (right edge) and all points at same vertical pos, 
            //change sign of curvature.  (We really should have evaluted this one from left to right). 
            if (*pupiltype == DARK)
            {
                if ((es[*strng_no].startY == es[*strng_no].endY) && (es[*strng_no].ave_dir > 0))// all points are at same vertical position & this is a right edge
                    es[*strng_no].curvature = -es[*strng_no].curvature; //we should have evaluated this from left to right instead of right to left
            }
            else
                if ((es[*strng_no].startY == es[*strng_no].endY) && (es[*strng_no].ave_dir < 0))// all points are at same vertical position & this is a right edge
                    es[*strng_no].curvature = -es[*strng_no].curvature; //we should have evaluated this from left to right instead of right to left
        }

		bool ImageProcessingUtilities_DSP::CloseOppositeEdge( int *ipnt, PNT EdgeBuf[], int *edgebuf_size, int DirBuf[], int *HairWidth )
            //This function is used to avoid considering edges that are probably edges of eye brow or eye lash hairs.  
            //A hair usually appears in the image as a thin dark line. After edge detection, it will usually appear as two, closely spaced, parallel edges,
            //(one edge at the light to dark transition, and a parallel edge at the dark to bright transition). 
            //If two very closely spaced points have light to dark gradients pointing towards each other, they are likely to be points on the two parallel edges of a thin dark line. 
            //"HairWidth" is currently measured in pixels.  In future should be real measurement unit, scaled to pixels by magnification and resolution factor for the image
            //(although that scaling could be done prior to this function call). 
        {
            bool closeOpposite = false;
			int direction;
            //values in DirBuf vary from 0 to 400 one goes from top to bottom, around the left edge of a bright circle. 
            //The DirBuf values vary from 0 to -400 from top to btoom, around the right edge.
            direction = Convert800LevelDirTo8LevelDir(&DirBuf[*ipnt]); //Convert DifBuf value to one of 8 cardinal directions 
            switch (direction)                                            
            {
                case 1: //look in up dir for a drk to brt edge [direction 5](look for top edge of horizontal hair)
                    for (int jpnt = *ipnt-1; jpnt > 0; --jpnt) //only need to look backward in buf to find point above current point
                    {
                        int deltaY = EdgeBuf[*ipnt].Y - EdgeBuf[jpnt].Y;
                        if( (deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if (EdgeBuf[jpnt].X == EdgeBuf[*ipnt].X)//we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 5 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == *HairWidth) && (EdgeBuf[jpnt].X < EdgeBuf[jpnt].X) ) break; //once we back up to left of current point, on "HairWidth" line above, we won't find another 
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 2: //look diagonally up to left for drk to brt edge [direction 6](look for top edge of hair oriented from lower left to upper rt)
                    for (int jpnt = *ipnt-1; jpnt > 0; --jpnt) //only need to look backward in buf to find point above current point
                    {
                        int deltaY = EdgeBuf[*ipnt].Y - EdgeBuf[jpnt].Y;
                        if( (deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X - deltaY)) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 6 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == *HairWidth) && (EdgeBuf[jpnt].X < (EdgeBuf[*ipnt].X - deltaY)) ) break; //we have backed up to the left of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 3: //look to left for drk to brt edge (look for left edge of vert hair)
                    for (int jpnt = *ipnt-1; jpnt > 0; --jpnt) //only need to look backward in buf to find point to left of current point
                    {
                        int deltaY = EdgeBuf[*ipnt].Y - EdgeBuf[jpnt].Y;
                        if( (deltaY == 0) ) // we are lookingto left along on same line
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X - deltaY)) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 7 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == 0) && ( (EdgeBuf[jpnt].X < (EdgeBuf[*ipnt].X) - *HairWidth) ) ) break; //we have backed up to the left of the last possible point
                        if (deltaY > 0) break; //we only want to look on same line
                    }
                    break;
                case 4: //look dwn to left for drk to brt edge (look for left edge of hair oriented from upper lft to lwr right)
                    for (int jpnt = *ipnt + 1; jpnt < *edgebuf_size; ++jpnt) //only need to look forward in buf to find point below current point
                    {
                        int deltaY = EdgeBuf[jpnt].Y - EdgeBuf[*ipnt].Y;
                        if( (deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X - deltaY)) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 8 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == *HairWidth) && (EdgeBuf[jpnt].X > (EdgeBuf[*ipnt].X - deltaY)) ) break; //we have advanced to the right of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 5: //look dwn [direction 1] for drk to brt edge
                    for (int jpnt = *ipnt + 1; jpnt < *edgebuf_size; ++jpnt) //only need to look forward in buf to find point below current point
                    {
                        int deltaY = EdgeBuf[jpnt].Y - EdgeBuf[*ipnt].Y;
                        if( (deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if ( EdgeBuf[jpnt].X == EdgeBuf[*ipnt].X ) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 1 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == *HairWidth) && (EdgeBuf[jpnt].X > (EdgeBuf[*ipnt].X + deltaY)) ) break; //we have advanced to the right of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 6: //look dwn to rt [direction 2] for drk to brt edge
                    for (int jpnt = *ipnt+1; jpnt < *edgebuf_size; ++jpnt) //only need to look forward in buf to find point below current point
                    {
                        int deltaY = EdgeBuf[jpnt].Y - EdgeBuf[*ipnt].Y;
                        if( (deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X + deltaY)) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 2 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == *HairWidth) && (EdgeBuf[jpnt].X > (EdgeBuf[*ipnt].X + deltaY)) ) break; //we have advanced to the right of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 7: //look right [direction 3] for drk to brt edge
                    for (int jpnt = *ipnt + 1; jpnt < *edgebuf_size; ++jpnt) //only need to look forward in buf to find point below current point
                    {
                        int deltaY = EdgeBuf[jpnt].Y - EdgeBuf[*ipnt].Y;
                        if( (deltaY == 0) )
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X + deltaY)) //we found another edge point a hair's width away in grad direction
                                if ( Convert800LevelDirTo8LevelDir(&DirBuf[jpnt]) == 3 ) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ( (deltaY == 0) && (EdgeBuf[jpnt].X > (EdgeBuf[*ipnt].X + deltaY)) ) break; //we have advanced to the right of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break;
                case 8: //look up to right [direction 4] for drk to brt edge 
                    for (int jpnt = *ipnt - 1; jpnt > 0; --jpnt) //only need to look backward in buf to find point above current point
                    {
                        int deltaY = EdgeBuf[*ipnt].Y - EdgeBuf[jpnt].Y;
                        if ((deltaY > 0) && (deltaY <= *HairWidth))
                        {
                            if (EdgeBuf[jpnt].X == (EdgeBuf[*ipnt].X + deltaY)) //we found another edge point a hair's width away in grad direction
                                if (Convert800LevelDirTo8LevelDir( &DirBuf[jpnt] ) == 4) //the edge point has grad direction opposite to original point (could be opposite edge of thin dark obj)
                                {
                                    closeOpposite = true;
                                    break;
                                }
                        }
                        if ((deltaY == *HairWidth) && (EdgeBuf[jpnt].X < (EdgeBuf[*ipnt].X + deltaY))) break; //we have backed up to the left of the last possible point
                        if (deltaY > *HairWidth) break; //if we are already more than HairWidth lines away, stop looking
                    }
                    break; 
            }
            return (closeOpposite);
        }

	    int ImageProcessingUtilities_DSP::Convert800LevelDirTo8LevelDir(int *I800LevelDir)
        //I800LevelDir are gradient directdion values that vary from 0 to 400 as one goes from top to bottom, around the left edge of a bright circle, 
        //and from 0 to -400 from top to btoom, around the right edge.
        //This routine Converts the diredtion value to one of 8 cardinal directions 
        {
            int direction = 0;
            if ((*I800LevelDir <= 50) && (*I800LevelDir >= -50)) direction = 1;  //drk to brt grad points down; brt to drk grad points up
            else if ((*I800LevelDir > 50) && (*I800LevelDir <= 150)) direction = 2;  //drk to brt grad points diag down to right
				else if ((*I800LevelDir > 150) && (*I800LevelDir <= 250)) direction = 3; //drk to brt grad points to right
					else if ((*I800LevelDir > 250) && (*I800LevelDir <= 350)) direction = 4;  //drk to brt grad points up to right
						else if ((*I800LevelDir > 350) && (*I800LevelDir <= 400)) direction = 5;  //drk to brt grad points up
							else if ((*I800LevelDir < -350) && (*I800LevelDir >= -400)) direction = 5; //drk to brt grad points up
								else if ((*I800LevelDir < -250) && (*I800LevelDir >= -350)) direction = 6;  //drk to brt grad points up to left
									else if ((*I800LevelDir < -150) && (*I800LevelDir >= -250)) direction = 7;    //drk to brt grad points left
										else if ((*I800LevelDir < -50) && (*I800LevelDir >= -150)) direction = 8;    //drk to brt grad down to left
            return (direction);
        }

		void ImageProcessingUtilities_DSP::MakeEdges
			(OBJECT_RECOGNITION_PARAMETERS *p, EDGE_STRING strng[], int *no_of_strngs, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges, int *no_rt_edges )
        {
			clock_start("ImageProcessingUtilities_DSP::MakeEdges");
            //Output no_lft_edges, no_rt_edges, lft_edge, and rt_edge. lft_edge and rt_edge are arrays of LFT_EDGE and RT_EDGE structures. Other two outputs specify number of elements in each array.
			// One string per edge  (no string combination logic yet)
            //const double mindir = 100.0; //left & rt edges should have ave gradients pointing more or less left or right,  
            //const double maxdir = 300.0; //not up or down.
            //const int minedgelength = 15; // minstring length to use as an edge
            //const double max_delta_curvature = 15.0; //we want consistant curvature (rate of change of direction); ave changes in curvature (delta curvature) should be small
            //const int mincurv = 2;//previous .01
            //const int maxcurv = 8;//previous .05
            *no_lft_edges = 0;
            *no_rt_edges = 0;
            //initialize edge arrays
            for (int i = 0; i < p->edge_max_no_edges; i++)
            {
                lft_edge[i].no_of_strngs = 0;
                rt_edge[i].no_of_strngs = 0;
            }
            for (int istrng = 0; istrng < *no_of_strngs; istrng++)
            {
                if (strng[istrng].pixel_lngth >= p->edge_minedgelength)//don't bother making an edge out of a very short string
                    if (strng[istrng].pixel_lngth <= p->edge_maxedgelength) //don't make an edge out of a ridiculously long string
                        if (strng[istrng].delta_curvature < p->edge_max_delta_curvature) // don't use strings that don't have reasonably consistent curvature
                        {//Note: if ave_dir close to 0, this is probably an eyelid boundary & not a string that we want to form a pupil "edge" with.
                            if ((strng[istrng].ave_dir > p->edge_mindir) && (strng[istrng].ave_dir < p->edge_maxdir)) //Left Edge
                                if ((strng[istrng].curvature > p->edge_mincurv) && (strng[istrng].curvature < p->edge_maxcurv))
                                {
                                    //lft_edge[*no_lft_edges].strng = new int[100];
                                    lft_edge[*no_lft_edges].strng[lft_edge[*no_lft_edges].no_of_strngs] = istrng;
                                    //lft_edge[no_lft_edges].lftmost = strng[lft_edge[no_lft_edges].no_of_strngs].hmin;
                                    lft_edge[*no_lft_edges].lftmost = strng[istrng].hmin;
                                    lft_edge[*no_lft_edges].rtmost = strng[istrng].hmax;
                                    lft_edge[*no_lft_edges].top = strng[istrng].top;
                                    lft_edge[*no_lft_edges].btm = strng[istrng].btm;
                                    lft_edge[*no_lft_edges].startX = strng[istrng].startX;
                                    lft_edge[*no_lft_edges].startY = strng[istrng].startY;
                                    lft_edge[*no_lft_edges].endX = strng[istrng].endX;
                                    lft_edge[*no_lft_edges].endY = strng[istrng].endY;
                                    lft_edge[*no_lft_edges].ave_dir = strng[istrng].ave_dir;
                                    ++( lft_edge[*no_lft_edges].no_of_strngs );
                                    ++( *no_lft_edges );
                                }
                            if ((strng[istrng].ave_dir < -p->edge_mindir) && (strng[istrng].ave_dir > -p->edge_maxdir)) //Right Edge
                                if ((strng[istrng].curvature < -p->edge_mincurv) && (strng[istrng].curvature > -p->edge_maxcurv))
                                {
                                    //rt_edge[*no_rt_edges].strng = new int[100];
                                    rt_edge[*no_rt_edges].strng[rt_edge[*no_rt_edges].no_of_strngs] = istrng;
                                    //rt_edge[no_rt_edges].rtmost = strng[rt_edge[no_rt_edges].no_of_strngs].hmax;
                                    rt_edge[*no_rt_edges].lftmost = strng[istrng].hmin;
                                    rt_edge[*no_rt_edges].rtmost = strng[istrng].hmax;
                                    rt_edge[*no_rt_edges].top = strng[istrng].top;
                                    rt_edge[*no_rt_edges].btm = strng[istrng].btm;
                                    rt_edge[*no_rt_edges].startX = strng[istrng].startX;
                                    rt_edge[*no_rt_edges].startY = strng[istrng].startY;
                                    rt_edge[*no_rt_edges].endX = strng[istrng].endX;
                                    rt_edge[*no_rt_edges].endY = strng[istrng].endY;
                                    rt_edge[*no_rt_edges].ave_dir = strng[istrng].ave_dir;
                                    ++( rt_edge[*no_rt_edges].no_of_strngs );
                                    ++( *no_rt_edges );
                                }
                        }
            }
            clock_end();
        }

		void ImageProcessingUtilities_DSP::MakeObjects
              (OBJECT_RECOGNITION_PARAMETERS *p, OBJECT objt[], int *no_of_objects, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges,
              int *no_rt_edges, EDGE_STRING strng[], int *nstring, PNT EdgeBuf[], int *EdgeBufSize )
        {
            //Output no_of_objects and objt.  objt is an array of OBJECT strutures; no_of_objects is the number of elements in the array.

			//For each left edge, see if any right edges overlap (have points on same vertical line).
            //If leftmost point of left edge is within reasonable distance of rightmost point of other, 
            //then these two edges form an object. 
            //Note that the same edge can be in multiple objects.
			// CAUTION:  Currently niether EdgeBufSize nor nstring is used.  Although we never write to EdgeBuf or strng in this routine, 
			//           EdgeBuf[] and string] index should really always be checked against EdgeBufSize or nstring to insure that we are not trying to read past the array boundary. 
			int toplftstringindex, toplftpntindex, btmlftstringindex, btmlftpntindex, toplft, btmlft;
			int toprtstringindex, toprtpntindex, btmrtstringindex, btmrtpntindex, toprt, btmrt, objht, objwdth;
			bool possible_object;
			//PNT BXtoplft, BXbtmlft, BXtoprt, BXbtmrt, top, btm;
            int n = 0;
			int ll, rr, ii;
			clock_start("ImageProcessingUtilities_DSP::MakeObjects");
            for ( ll = 0; ll < *no_lft_edges; ll++)
                for ( rr = 0; rr < *no_rt_edges; rr++)
                {
                    toplftstringindex = lft_edge[ll].strng[0];
                    toplftpntindex = strng[toplftstringindex].edge_buf_index[0];
                    btmlftstringindex = lft_edge[ll].strng[lft_edge[ll].no_of_strngs - 1];
                    btmlftpntindex = strng[btmlftstringindex].edge_buf_index[strng[btmlftstringindex].pixel_lngth - 1];
                    toplft = EdgeBuf[toplftpntindex].Y;
                    btmlft = EdgeBuf[btmlftpntindex].Y;

                    toprtstringindex = rt_edge[rr].strng[0];
                    toprtpntindex = strng[toprtstringindex].edge_buf_index[0];
                    btmrtstringindex = rt_edge[rr].strng[rt_edge[rr].no_of_strngs - 1];
                    btmrtpntindex = strng[btmrtstringindex].edge_buf_index[strng[btmrtstringindex].pixel_lngth - 1];
                    toprt = EdgeBuf[toprtpntindex].Y;
                    btmrt = EdgeBuf[btmrtpntindex].Y;
                    objht = 0;
                    objwdth = 0;
                    possible_object = false;

                    //We want the left edge top to be left of the right edge top, and the same for the bottoms.
                    //In other words we don't want left and rt edges to cross.
                    if ((EdgeBuf[toplftpntindex].X < EdgeBuf[toprtpntindex].X) && (EdgeBuf[btmlftpntindex].X < EdgeBuf[btmrtpntindex].X))
                    {
                        //If there are any lft & rt pnts on the same horizontal line (top rt is between top & btm left ),
                        if ((toplft <= toprt) && (btmlft >= toprt))
                        {
                            //find left edge point on same line as top rt -- make sure it is to the left of toprt
                            for ( ii = 0; ii < strng[lft_edge[ll].strng[0]].pixel_lngth; ii++)
                                if (EdgeBuf[strng[lft_edge[ll].strng[0]].edge_buf_index[ii]].Y >= toprt)//found the line
                                {
                                    if (EdgeBuf[strng[lft_edge[ll].strng[0]].edge_buf_index[ii]].X < EdgeBuf[toprtpntindex].X)
                                        possible_object = true;
                                    break;
                                }
                        }
                        //If there are any lft & rt pnts on the same horizontal line (top lft is between top & btm rt ),
                        else if ((toprt < toplft) && (btmrt > toplft))
                        {
                            //find right edge point on same line as top left -- make sure it is to the right of toplft
                            for (int ii = 0; ii < strng[rt_edge[rr].strng[0]].pixel_lngth; ii++)
                                if (EdgeBuf[strng[rt_edge[rr].strng[0]].edge_buf_index[ii]].Y >= toplft) //found the line
                                {
                                    if (EdgeBuf[strng[rt_edge[rr].strng[0]].edge_buf_index[ii]].X > EdgeBuf[toplftpntindex].X)
                                        possible_object = true;
                                    break;
                                }
                        }
                    }
                    if (possible_object)
                    // If we are satisied that left and right edges overlap & the left edge is really to the left of the right edge
                    // then check for edges in between them, and check for sensible size and curvature values
                    {
                        //if there are other edges in between (in the box formed by the end points of the pair), 
                        //then don't clasify this pair as an object
                        //BXtoplft.X = lft_edge[ll].startX;
                        //BXtoplft.Y = lft_edge[ll].startY;
                        //BXbtmlft.X = lft_edge[ll].endX;
                        //BXbtmlft.Y = lft_edge[ll].endY;
                        //BXtoprt.X = rt_edge[rr].startX;
                        //BXtoprt.Y = rt_edge[rr].startY;
                        //BXbtmrt.X = rt_edge[rr].endX;
                        //BXbtmrt.Y = rt_edge[rr].endY;

                        //for (int ii = 0; ii < no_lft_edges; ii++)
                        //{
                        //    if (ii == ll) continue;//skip the current left edge
                        //    top.X = lft_edge[ii].startX;
                        //    top.Y = lft_edge[ii].startY;
                        //    btm.X = lft_edge[ii].endX;
                        //    btm.Y = lft_edge[ii].endY;
                        //    if (!LineSegmentIsInsideBox(BXtoplft, BXtoprt, BXbtmlft, BXbtmrt, top, btm))
                        //        continue;
                        //    else reject_pair = true;
                        //}
                        //if ( !reject_pair )
                        //    for (int ii = 0; ii < no_rt_edges; ii++)
                        //    {
                        //        if (ii == rr) continue;//skip the current rt edge
                        //        top.X = rt_edge[ii].startX;
                        //        top.Y = rt_edge[ii].startY;
                        //        btm.X = rt_edge[ii].endX;
                        //        btm.Y = rt_edge[ii].endY;
                        //        if (!LineSegmentIsInsideBox(BXtoplft, BXtoprt, BXbtmlft, BXbtmrt, top, btm))
                        //            continue;
                        //        else reject_pair = true;
                        //    }
                        //if (!reject_pair)
                        {
                            objwdth = rt_edge[rr].rtmost - lft_edge[ll].lftmost;
                            if (btmlft > btmrt) objht = btmlft - toplft;
                            else objht = btmrt - toplft;
                            // if width and height are within reasonable range,
                            if ((objwdth < p->obj_maxwdth) && (objwdth > p->obj_minwdth))
                                if ((objht < p->obj_maxht) && (objht > p->obj_minht))
                                    //and if left and right edges have roughly equal and opposite curvatures
                                    if (((strng[rt_edge[rr].strng[0]].curvature < 0.0) && (strng[lft_edge[ll].strng[0]].curvature > 0.0)) ||
                                        ((strng[rt_edge[rr].strng[0]].curvature > 0.0) && (strng[lft_edge[ll].strng[0]].curvature < 0.0)))
                                        if (AbsFloat(strng[rt_edge[rr].strng[0]].curvature + strng[lft_edge[ll].strng[0]].curvature) < p->obj_max_curv_dif)
                                        {
                                            double ave_curv = (strng[lft_edge[ll].strng[0]].curvature - strng[rt_edge[rr].strng[0]].curvature) / 2.0;
                                            double PDest = (2.0 / 3.14) * (400.0 / ave_curv); //circum = pi*d; 1/2 circum = 400/curv; d = (2/pi) * (400/curv);
                                            double PDest_to_objwdth_ratio = PDest / objwdth;
                                            if ((PDest_to_objwdth_ratio > p->obj_minPD_to_width) && (PDest_to_objwdth_ratio < p->obj_maxPD_to_width))
                                            //and if curvature is roughly consistant with object width 
                                            { //this edge pair is an object.
                                                objt[n].lft_edge = ll;
                                                objt[n].rt_edge = rr;
                                                n++;
                                            }
                                        }
                        }
                    }
                }
            *no_of_objects = n;
            clock_end();
        }

		void ImageProcessingUtilities_DSP::MakeCRObjects
            (OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT objt[], int *no_of_objects, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *no_lft_edges, int *no_rt_edges)

			//Output objt, and no_of_objects.
        {
			bool rt_edge_part_of_obj [1000];
			bool lft_edge_part_of_obj [1000];
            //initialize
            int objt_no = 0;
            int iledge, jredge, iredge, jledge, ii;
            clock_start("ImageProcessingUtilities_DSP::MakeCRObjects");
            *no_of_objects = 0;
            for (ii = 0; ii < *no_lft_edges; ii++)
                lft_edge_part_of_obj[ii] = false;
            for (ii = 0; ii < *no_rt_edges; ii++)
                rt_edge_part_of_obj[ii] = false;

            //For each edge, see if other edges are close, & have roughly same center 
            for (iledge = 0; iledge < *no_lft_edges; iledge++) //first try each left edge as first edge in an object
            {
                if (!lft_edge_part_of_obj[iledge]) //don't start an object with an edge that is already part of another object
                {
                    startCRobject(p, &objt[objt_no], &iledge, LEFT); //start an object with first left edge, and see if there are rt edges that can be part of the object
                    lft_edge_part_of_obj[iledge] = true;
                    for (jredge = 0; jredge < *no_rt_edges; jredge++) //find at least one right edge to form obj with left edge
                    {
                        if (!rt_edge_part_of_obj[jredge])//if not already part of a CR object
                            if (fits_CR_object(p, &objt[objt_no], lft_edge, rt_edge, &jredge, RIGHT)) //add it to the object if it meets size, pos, and shape criteria
                            {
                                add_to_CR_object(&objt[objt_no], &jredge, RIGHT);
                                rt_edge_part_of_obj[jredge] = true;
                            }
                    }
                    if (objt[objt_no].no_rt_edges > 0) //don't bother adding left edges unless we have included at least 1 rt edge
                    {//if we have at least 1 right edge, see if there are aother lft edges that can be added
                        for (jledge = 0; jledge < *no_lft_edges; jledge++)
                        {
                            if (!(jledge == iledge))
                            {
                                if (!lft_edge_part_of_obj[jledge])
                                    if (fits_CR_object(p, &objt[objt_no], lft_edge, rt_edge, &jledge, LEFT))
                                    {
                                        add_to_CR_object(&objt[objt_no], &jledge, LEFT);
                                        lft_edge_part_of_obj[jledge] = true;
                                    }
                            }
                        }
                    }
                }
                if (objt[objt_no].no_rt_edges < 1)//don't include this as an object without at least 1 lft and 1 rt edge
                {
                    lft_edge_part_of_obj[iledge] = false;//remove part of obj flag
                    objt[objt_no].no_lft_edges = 0;//reset to no edges in object
                }
                else //This is a CR object
                {
                    ++objt_no;
                    ++( *no_of_objects );
                }
            }
            for (iredge = 0; iredge < *no_rt_edges; iredge++) //Now try each right edge as first in an object
            {
                if (!rt_edge_part_of_obj[iredge]) //don't start an object with an edge that is already part of another object
                {
                    startCRobject(p, &objt[objt_no], &iredge, RIGHT); //start an object with first rt edge, and see if there are lft edges that can be part of the object
                    for (jledge = 0; jledge < *no_lft_edges; jledge++) //find at least one left edge to form obj with rt edge
                    {
                        if (!lft_edge_part_of_obj[jledge])
                            if (fits_CR_object(p, &objt[objt_no], lft_edge, rt_edge, &jledge, LEFT))
                            {
                                add_to_CR_object(&objt[objt_no], &jledge, LEFT);
                                lft_edge_part_of_obj[jledge] = true;
                            }
                    }
                    if (objt[objt_no].no_lft_edges > 0) //don't bother adding rt edges unless we have included at least 1 lft edge
                    {//if we have at least 1 lft edge, see if there are aother rt edges that can be added
                        for (jredge = 0; jredge < *no_rt_edges; jredge++)
                        {
                            if (!(jredge == iredge)) //skip the edge that we used to start the object
                            {
                                if (!rt_edge_part_of_obj[jredge])
                                    if (fits_CR_object(p, &objt[objt_no], lft_edge, rt_edge, &jredge, RIGHT))
                                    {
                                        add_to_CR_object(&objt[objt_no], &jredge, RIGHT);
                                        rt_edge_part_of_obj[jredge] = true;
                                    }
                            }
                        }
                    }
                }
                if (objt[objt_no].no_lft_edges < 1)//don't include this as an object without at least 1 lft and 1 rt edge
                {
                    rt_edge_part_of_obj[iredge] = false; //remove part of obj flag  
                    objt[objt_no].no_rt_edges = 0;
                }
                else //This is a CR object
                {
                    ++objt_no;
                    ++( *no_of_objects );
                }
            }
            clock_end();
        }

		void ImageProcessingUtilities_DSP::startCRobject(OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT *objt, int *i, int side)
        {

            //objt->lft_edge = new int[p->obj_max_edges_per_object];
            //objt->rt_edge = new int[p->obj_max_edges_per_object];
            if (side == LEFT)
            {
                objt->no_lft_edges = 1;
                objt->no_rt_edges = 0;
                objt->lft_edge[0] = *i;
            }
            if (side == RIGHT)
            {
                objt->no_rt_edges = 1;
                objt->no_lft_edges = 0;
                objt->rt_edge[0] = *i;
            }
        }

		void ImageProcessingUtilities_DSP::add_to_CR_object(CR_OBJECT *objt, int *edge_no, int side)
        {
			int n;
            if (side == LEFT)
            {
                n = objt->no_lft_edges;
                objt->lft_edge[n] = *edge_no;
                (objt->no_lft_edges)++;
            }
            if (side == RIGHT)
            {
                n = objt->no_rt_edges;
                objt->rt_edge[n] = *edge_no;
                (objt->no_rt_edges)++;
            }
        }

		bool ImageProcessingUtilities_DSP::fits_CR_object(OBJECT_RECOGNITION_PARAMETERS *p, CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], int *edge_no, int side)
        {
            bool accept = true;
            float proximity, delta_dir;
            float minproximity = 10000; // initialize to impossibly large value
            int twice_mean_hpos_rt_edg, twice_mean_hpos_lft_edg, twice_mean_vpos_rt_edg, twice_mean_vpos_lft_edg, hdiff, vdiff;

            if (side == LEFT)
            {
                for (int ll = 0; ll < objt->no_lft_edges; ll++)//compare to other same-side edges in object
                {
                    //if candidate edge overlaps same side edge, reject
                    if (((lft_edge[objt->lft_edge[ll]].top > lft_edge[*edge_no].top) && (lft_edge[objt->lft_edge[ll]].top < lft_edge[*edge_no].btm))
                        || ((lft_edge[objt->lft_edge[ll]].btm > lft_edge[*edge_no].top) && (lft_edge[objt->lft_edge[ll]].btm < lft_edge[*edge_no].btm))
                       )
                    {
                        accept = false;
                        break;
                    }

                    if (lft_edge[objt->lft_edge[ll]].btm <= lft_edge[*edge_no].top) //if candidate is below same side edge it should have larger dir value
                    {
                        if (lft_edge[*edge_no].ave_dir < lft_edge[objt->lft_edge[ll]].ave_dir)
                        {
                            accept = false;
                            break;
                        }
                    }
                    if (lft_edge[objt->lft_edge[ll]].top >= lft_edge[*edge_no].btm) //if candidate is above same side edge it should have smaller dir value
                    {
                        if (lft_edge[*edge_no].ave_dir > lft_edge[objt->lft_edge[ll]].ave_dir)
                        {
                            accept = false;
                            break;
                        }
                    }
                    proximity = closest_end_points
                        (&lft_edge[*edge_no].startX, &lft_edge[*edge_no].startY, &lft_edge[*edge_no].endX, &lft_edge[*edge_no].endY,
                          &lft_edge[objt->lft_edge[ll]].startX, &lft_edge[objt->lft_edge[ll]].startY, &lft_edge[objt->lft_edge[ll]].endX, &lft_edge[objt->lft_edge[ll]].endY);

                    if (ll == 0) minproximity = proximity;
                    else
                        if (proximity < minproximity) minproximity = proximity;
                }
                for (int rr = 0; rr < objt->no_rt_edges; rr++)//compare to other opposite side edges in object
                {
                    twice_mean_hpos_rt_edg = rt_edge[objt->rt_edge[rr]].lftmost + rt_edge[objt->rt_edge[rr]].rtmost;
                    twice_mean_hpos_lft_edg = lft_edge[*edge_no].lftmost + lft_edge[*edge_no].rtmost;
                    hdiff = twice_mean_hpos_rt_edg - twice_mean_hpos_lft_edg;
                    if (hdiff < 0 || (hdiff > (2 * p->obj_maxwdth)))
                    {//if lft edge is further to right than right edge, or they are separated by > maximum CR diam, reject
                        accept = false;
                        break;
                    }
                    delta_dir = (float)( rt_edge[objt->rt_edge[rr]].ave_dir + lft_edge[*edge_no].ave_dir );
                    twice_mean_vpos_rt_edg = rt_edge[objt->rt_edge[rr]].top + rt_edge[objt->rt_edge[rr]].btm;
                    twice_mean_vpos_lft_edg = lft_edge[*edge_no].top + lft_edge[*edge_no].btm;
                    vdiff = twice_mean_vpos_rt_edg - twice_mean_vpos_lft_edg;
					if (ImageProcessingUtilities_FPGA::AbsInt(vdiff) > 2 * p->obj_maxht)  //if vertical center of candidate is farther than max CR diam from opposite edge, reject
                    {
                        accept = false;
                        break;
                    }
                    if (AbsFloat(delta_dir) > p->obj_neutral_delta_dir)  //neutral delta dir means they both point more or less equally up or equally down
                    { // If ave_dir for one edge points mostly up and the other points mostly down, 
                        // the one pointing up should have a vertical center below the one pointing down.
                        if (delta_dir < 0) //right edge vert center should be lower than left edge
                        {
                            if ((vdiff < 0))
                            {
                                accept = false;
                                break;
                            }
                        }
                        else //right edge vert center should be higher than left edge
                        {
                            if (vdiff > 0)
                            {
                                accept = false;
                                break;
                            }
                        }
                    }//if (delta_dir > neutral_delta_dir)
                    proximity = closest_end_points
                        (&lft_edge[*edge_no].startX, &lft_edge[*edge_no].startY, &lft_edge[*edge_no].endX, &lft_edge[*edge_no].endY,
                          &rt_edge[objt->rt_edge[rr]].startX, &rt_edge[objt->rt_edge[rr]].startY, &rt_edge[objt->rt_edge[rr]].endX, &rt_edge[objt->rt_edge[rr]].endY);

                    if ((rr == 0) && (objt->no_lft_edges == 0)) minproximity = proximity;
                    else
                        if (proximity < minproximity) minproximity = proximity;
                }//for loop -- compare to other opposite side edges
                if (minproximity < p->obj_max_edge_end_point_separation_sqrd) accept = false;
                return (accept);
            }//if candidate is a left edge

            else //if (side == Side.right)
            {
                for (int rr = 0; rr < objt->no_rt_edges; rr++)//compare to other same-side edges in object
                {
                    if (rr == *edge_no) continue; //don't compare edge to itself
                    //if candidate edge overlaps same side edge, reject
                    if (((rt_edge[objt->rt_edge[rr]].top > rt_edge[*edge_no].top) && (rt_edge[objt->rt_edge[rr]].top < rt_edge[*edge_no].btm))
                        || ((rt_edge[objt->rt_edge[rr]].btm > rt_edge[*edge_no].top) && (rt_edge[objt->rt_edge[rr]].btm < rt_edge[*edge_no].btm))
                       )
                    {
                        accept = false;
                        break;
                    }

                    if (rt_edge[objt->rt_edge[rr]].btm <= rt_edge[*edge_no].top) //if candidate is below same side edge it should have more negative dir value
                    {
                        if (rt_edge[*edge_no].ave_dir > rt_edge[objt->rt_edge[rr]].ave_dir)
                        {
                            accept = false;
                            break;
                        }
                    }
                    if (rt_edge[objt->rt_edge[rr]].top >= rt_edge[*edge_no].btm) //if candidate is above same side edge it should have a less negative dir value
                    {
                        if (rt_edge[*edge_no].ave_dir < rt_edge[objt->rt_edge[rr]].ave_dir)
                        {
                            accept = false;
                            break;
                        }
                    }
                    proximity = closest_end_points
                        (&rt_edge[*edge_no].startX, &rt_edge[*edge_no].startY, &rt_edge[*edge_no].endX, &rt_edge[*edge_no].endY,
                          &rt_edge[objt->rt_edge[rr]].startX, &rt_edge[objt->rt_edge[rr]].startY, &rt_edge[objt->rt_edge[rr]].endX, &rt_edge[objt->rt_edge[rr]].endY);

                    if (rr == 0) minproximity = proximity;
                    else
                        if (proximity < minproximity) minproximity = proximity;
                }
                for (int ll = 0; ll < objt->no_lft_edges; ll++)//compare to other opposite side edges in object
                {

                    twice_mean_hpos_lft_edg = lft_edge[objt->lft_edge[ll]].lftmost + lft_edge[objt->lft_edge[ll]].rtmost;
                    twice_mean_hpos_rt_edg = rt_edge[*edge_no].lftmost + rt_edge[*edge_no].rtmost;
                    hdiff = twice_mean_hpos_rt_edg - twice_mean_hpos_lft_edg;
                    if (hdiff < 0 || (hdiff > (2 * p->obj_maxwdth)))
                    {//if lft edge is further to right than right edge, or they are separated by > maximum CR diam, reject
                        accept = false;
                        break;
                    }
                    delta_dir = lft_edge[objt->lft_edge[ll]].ave_dir + rt_edge[*edge_no].ave_dir;
                    twice_mean_vpos_lft_edg = lft_edge[objt->lft_edge[ll]].top + lft_edge[objt->lft_edge[ll]].btm;
                    twice_mean_vpos_rt_edg = rt_edge[*edge_no].top + rt_edge[*edge_no].btm;
                    vdiff = twice_mean_vpos_lft_edg - twice_mean_vpos_rt_edg;
					if (ImageProcessingUtilities_FPGA::AbsInt(vdiff) > 2 * p->obj_maxht) //if vertical center of candidate is farther than max CR diam from opposite edge, reject
                    {
                        accept = false;
                        break;
                    }
                    if ( AbsFloat(delta_dir)  > (float)( p->obj_neutral_delta_dir ) )  //neutral delta dir means they both point more or less equally up or equally down
                    { // If ave_dir for one edge points mostly up and the other points mostly down, 
                        // the one pointing up should have a vertical center below the one pointing down.
                        if (delta_dir < 0.0f) //right edge vert center should be lower than left edge
                        {
                            if ((vdiff > 0))
                            {
                                accept = false;
                                break;
                            }
                        }
                        else //right edge vert center should be higher than left edge
                        {
                            if (vdiff < 0)
                            {
                                accept = false;
                                break;
                            }
                        }
                    }//if (delta_dir > neutral_delta_dir)
                    proximity = closest_end_points
                        (&rt_edge[*edge_no].startX, &rt_edge[*edge_no].startY, &rt_edge[*edge_no].endX, &rt_edge[*edge_no].endY,
                          &lft_edge[objt->lft_edge[ll]].startX, &lft_edge[objt->lft_edge[ll]].startY, &lft_edge[objt->lft_edge[ll]].endX, &lft_edge[objt->lft_edge[ll]].endY);
                    if ((ll == 0) && (objt->no_rt_edges == 0)) minproximity = proximity;
                    else
                        if (proximity < minproximity) minproximity = proximity;

                }//for loop -- compare to other opposite side edges
                if (minproximity > p->obj_max_edge_end_point_separation_sqrd) accept = false;
                return (accept);
            }//if candidate is a left edge
        }

		float ImageProcessingUtilities_DSP::closest_end_points(int *Xstrt1, int *Ystrt1, int *Xend1, int *Yend1, int *Xstrt2, int *Ystrt2, int *Xend2, int *Yend2)
        {
            //find the shortest distance from either end of string 1 to either end of string 2 and return that distance value
            // (combinations are strt1 to strt2, strt1 to end2, end1 to strt2, end1 to end2) 
            float dist_sqrd;
            float mindist_sqrd = 10000.0; //start at impossibly high value
			float x1[2], x2[2],y1[2],y2[2];
            int i, j;
            x1[0] = (float)( *Xstrt1 );
            y1[0] = (float)( *Ystrt1 );
            x1[1] = (float)( *Xend1 );
            y1[1] = (float)( *Yend1 );
            x2[0] = (float)( *Xstrt2 );
            y2[0] = (float)( *Ystrt2 );
            x2[1] = (float)( *Xend2 );
            y2[1] = (float)( *Yend2 );
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    dist_sqrd = ((x1[i] - x2[j]) * (x1[i] - x2[j])) + ((y1[i] - y2[j]) * (y1[i] - y2[j]));
                    if ((i == 0) && (j == 0)) mindist_sqrd = dist_sqrd;
                    else
                        if (dist_sqrd < mindist_sqrd) mindist_sqrd = dist_sqrd;
                }
            return (mindist_sqrd);
        }

		bool  ImageProcessingUtilities_DSP::FitPupilEllipse
			(OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p)
        {
            //Make list of edge points in object
            //PNT *obj_point;
			PNT obj_point[100000];
            int no_points_in_list;
            float xc, yc, major, minor, alpha_major;
            clock_start("ImageProcessingUtilities_DSP::FitPupilEllipse");
            MakePupilObjectEdgePointList(objt, lft_edge, rt_edge, es, EdgeBuf, edgebuf_size, p, obj_point, &no_points_in_list);
            if ( ComputeEllipse(obj_point, &no_points_in_list, &xc, &yc, &major, &minor, &alpha_major) )
            {
                objt->elps_hc = xc;
                objt->elps_vc = yc;
                objt->elps_maj_rad = major;
                objt->elps_min_rad = minor;
                objt->elps_angle = alpha_major;
                objt->elps_aspect = minor / major;
                objt->elps_fit = true;
				//delete[] obj_point; //free memory allocated in MakePupilObjectEdgePointList
                clock_end();
                return (true);
            }
            else
            {
                objt->elps_fit = false;
				//delete[] obj_point; //free memory allocated in MakePupilObjectEdgePointList
                clock_end();
                return (false); // error
            }
            clock_end();
        }

	    bool ImageProcessingUtilities_DSP::FitCREllipse(CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], 
			                                            int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p )
        {
            //Make list of edge points in object
            //PNT *obj_point;
			PNT obj_point[100000];
            int no_points_in_list;
            float xc, yc, major, minor, alpha_major;
            MakeCRObjectEdgePointList(objt, lft_edge, rt_edge, es, EdgeBuf, edgebuf_size, p, obj_point, &no_points_in_list);
            if ( ComputeEllipse(obj_point, &no_points_in_list, &xc, &yc, &major, &minor, &alpha_major) )
            {
                objt->hc = xc;
                objt->vc = yc;
                objt->elps_maj_rad = major;
                objt->elps_min_rad = minor;
                objt->elps_angle = alpha_major;
                objt->elps_aspect = major / minor;
				objt->elps_fit = true;
				//delete[] obj_point; //free memory allocated in MakeCRObjectEdgePointList
                return (true);
            }
            else 
			{
                objt->elps_fit = false;
				//delete[] obj_point; //free memory allocated in MakeCRObjectEdgePointList
				return (false); //error
			}
        }


		void ImageProcessingUtilities_DSP::MakePupilObjectEdgePointList
            (OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p,
			PNT objt_point[], int *no_points_in_list)
              //PNT **objt_point, int *no_points_in_list)
        {
			//Make list of edge points to be used by ellipse fit routine.
			//NOTE: OBJECT_RECOGNITION_PARAMETERS struct is passed to this routine so that p->feature_max_no_of_points_in_elps can be used
			//      in future to limit no of points used for curve fit.  It is not currently used. 
            int ii = 0; //points in objt_point array
            int maxdim = 0;
            *no_points_in_list = 0;
            for (int nl = 0; nl < lft_edge[objt->lft_edge].no_of_strngs; nl++)
                maxdim += es[lft_edge[objt->lft_edge].strng[nl]].pixel_lngth;
            for (int nr = 0; nr < rt_edge[objt->rt_edge].no_of_strngs; nr++)
                maxdim += es[rt_edge[objt->rt_edge].strng[nr]].pixel_lngth;
            //*objt_point = new PNT[maxdim];
            for (int lls = 0; lls < lft_edge[objt->lft_edge].no_of_strngs; lls++) //for each string in left edge
			{
                for (int ps = 0; ps < es[lft_edge[objt->lft_edge].strng[lls]].pixel_lngth; ps++) //for each pixel in string
                {
                    objt_point[ii].X = EdgeBuf[es[lft_edge[objt->lft_edge].strng[lls]].edge_buf_index[ps]].X;//get col value from edge buffer & add to objt point array
                    objt_point[ii].Y = EdgeBuf[es[lft_edge[objt->lft_edge].strng[lls]].edge_buf_index[ps]].Y;//get row value from edge buffer & add to objt point array
                    ++ii;
					if ( ii >= 100000) break;
                }
				if ( ii >= 100000) break;
			}
				if( ii < 100000 ) 
				{
					for (int rrs = 0; rrs < rt_edge[objt->rt_edge].no_of_strngs; rrs++) //for each string in right edge
					{
						for (int ps = 0; ps < es[rt_edge[objt->rt_edge].strng[rrs]].pixel_lngth; ps++)//for each pixel in string
						{
							objt_point[ii].X = EdgeBuf[es[rt_edge[objt->rt_edge].strng[rrs]].edge_buf_index[ps]].X;//get col value from edge buffer & add to objt point array
							objt_point[ii].Y = EdgeBuf[es[rt_edge[objt->rt_edge].strng[rrs]].edge_buf_index[ps]].Y;//get row value from edge buffer & add to objt point array
							++ii;
							if ( ii >= 100000) break;
						}
						if ( ii >= 100000) break;
					}
				}
            *no_points_in_list = ii;
        }

		void ImageProcessingUtilities_DSP::MakeCRObjectEdgePointList
            (CR_OBJECT *objt, LFT_EDGE lft_edge[], RT_EDGE rt_edge[], EDGE_STRING es[], PNT EdgeBuf[], int *edgebuf_size, OBJECT_RECOGNITION_PARAMETERS *p,
				PNT objt_point[], int *no_points_in_list)
              //PNT **objt_point, int *no_points_in_list)
        {
            int ii = 0; //points in objt_point array
            *no_points_in_list = 0;
            int maxdim = 0;
            //compute no of points in list (add no of points in each specified string)
            for (int nnle = 0; nnle < objt->no_lft_edges; nnle++)
                for (int nls = 0; nls < lft_edge[objt->lft_edge[nnle]].no_of_strngs; nls++)
                    maxdim += es[lft_edge[objt->lft_edge[nnle]].strng[nls]].pixel_lngth;
            for (int nnre = 0; nnre < objt->no_rt_edges; nnre++)
                for (int nrs = 0; nrs < rt_edge[objt->rt_edge[nnre]].no_of_strngs; nrs++)
                    maxdim += es[rt_edge[objt->rt_edge[nnre]].strng[nrs]].pixel_lngth;

            //*objt_point = new PNT[maxdim];
            for (int nle = 0; nle < objt->no_lft_edges; nle++) //for each left edge in object
                for (int lls = 0; lls < lft_edge[objt->lft_edge[nle]].no_of_strngs; lls++) //for each string in left edge
                    for (int ps = 0; ps < es[lft_edge[objt->rt_edge[nle]].strng[lls]].pixel_lngth; ps++) //for each pixel in string
                    {
                        objt_point[ii].X = EdgeBuf[es[lft_edge[objt->lft_edge[nle]].strng[lls]].edge_buf_index[ps]].X;//get col value from edge buffer & add to objt point array
                        objt_point[ii].Y = EdgeBuf[es[lft_edge[objt->lft_edge[nle]].strng[lls]].edge_buf_index[ps]].Y;//get row value from edge buffer & add to objt point array
                        ++ii;
                    }
            for (int nre = 0; nre < objt->no_rt_edges; nre++)
                for (int rrs = 0; rrs < rt_edge[objt->rt_edge[nre]].no_of_strngs; rrs++) //for each string in right edge
                    for (int ps = 0; ps < es[rt_edge[objt->rt_edge[nre]].strng[rrs]].pixel_lngth; ps++)//for each pixel in string
                    {
                        objt_point[ii].X = EdgeBuf[es[rt_edge[objt->rt_edge[nre]].strng[rrs]].edge_buf_index[ps]].X;//get col value from edge buffer & add to objt point array
                        objt_point[ii].Y = EdgeBuf[es[rt_edge[objt->rt_edge[nre]].strng[rrs]].edge_buf_index[ps]].Y;//get row value from edge buffer & add to objt point array
                        ++ii;
                    }
            *no_points_in_list = ii;
        }



		bool ImageProcessingUtilities_DSP::ComputeEllipse(PNT edge_point[], int *no_of_points, float *xc, float *yc, float *major, float *minor, float *alpha_major)
        {
            //compute sums
            double a[5][6];
            double y, y2, y3, x, x2, x3;
            double tmp, aa, bb, dd, ee, ff, det, aa1, bb1, tg2alpha, sin2alpha, f, a2, b2, ra, rb, alpha;
            int i, j, k;
			for (i = 0; i<5; i++)//initialize A matrix
				for (j=0; j<6; j++)
					a[i][j] = 0.0;
            for (int ii = 0; ii < *no_of_points; ii++)
            {
                x = edge_point[ii].X;
                y = edge_point[ii].Y;
                y2 = y * y;
                y3 = y * y2;

                x2 = x * x;
                x3 = x * x2;

				a[0][0] += x2 * x2;
                a[0][1] += x3 * y;
                a[0][2] += x3;
                a[0][3] += x2 * y;
                a[0][4] += x2;
                a[0][5] += x2 * y2;

                a[1][3] += x * y2;
                a[1][4] += x * y;
                a[1][5] += x * y3;

                a[2][4] += x;

                a[3][3] += y2;
                a[3][4] += y;
                a[3][5] += y3;

                a[4][4] += 1;
            }
            //compute ellipse
            // complete the matrix
            a[1][0] = a[0][1];
            a[1][1] = a[0][5];
            a[1][2] = a[0][3];

            a[2][0] = a[0][2];
            a[2][1] = a[0][3];
            a[2][2] = a[0][4];
            a[2][3] = a[1][4];
            a[2][5] = a[1][3];

            a[3][0] = a[0][3];
            a[3][1] = a[1][3];
            a[3][2] = a[1][4];

            a[4][0] = a[0][4];
            a[4][1] = a[1][4];
            a[4][2] = a[2][4];
            a[4][3] = a[3][4];
            a[4][5] = a[3][3];


            // make triangular matrix
            for (k = 4; k > 0; k--)
            {
                for (i = 5 - k; i < 5; i++) // rows 1-4
                {
                    tmp = a[i][k] / a[4 - k][k];
                    for (j = 0; j < 6; j++)
                        a[i][j] -= a[4 - k][j] * tmp;
                }
            }
            // solve simulteneous equations (now when matrix is triangular it becomes trivial)
            aa = a[4][5] / a[4][0];
            bb = (a[3][5] - a[3][0] * aa) / a[3][1];
            dd = (a[2][5] - a[2][0] * aa - a[2][1] * bb) / a[2][2];
            ee = (a[1][5] - a[1][0] * aa - a[1][1] * bb - a[1][2] * dd) / a[1][3];
            ff = (a[0][5] - a[0][0] * aa - a[0][1] * bb - a[0][2] * dd - a[0][3] * ee) / a[0][4];

            // calculate ellipse data
            // verify that equation returns ellipse
            det = bb * bb + 4 * aa;
            if (det >= 0)
            {
                // equation returned by MLS does not represent an ellipse
                *xc = 0.0;
                *yc = 0.0;
                ra = 0;
                rb = 0;
                alpha = 0;
                *major = 0.0;
                *minor = 0.0;
                *alpha_major = 0.0;
                return (false);//error -- no valid elipse
            }
            // ellipse center
            *xc = (float)( -(2 * dd + bb * ee) / det );
            *yc = (float)( -(bb * dd - 2 * aa * ee) / det );

            // Trigonometry
            aa1 = aa + 1.0;
            if (aa1 == 0.0) aa1 = 0.00001;//do this to prevent divide by zero
            tg2alpha = bb / aa1;
            sin2alpha = tg2alpha / sqrt(1 + tg2alpha * tg2alpha);
            alpha = sin2alpha; //*alpha = 0.5 * asin(sin2alpha);

            // ellipse axis
            f = 2 * (ff - aa * *xc * *xc - bb * *xc * *yc + *yc * *yc);
            bb1 = bb / sin2alpha;
            a2 = f / (1 - aa + bb1);
            b2 = f / (1 - aa - bb1);
            rb = sqrt(b2);
            ra = sqrt(a2);

            // find major axis
            if (rb >= ra)
            {
                *major = (float)(rb);
                *minor = (float)(ra);
                *alpha_major = (float)(alpha);
            }
            else
            {
                *major = (float)(ra);
                *minor = (float)(rb);
                *alpha_major = (float)( alpha - M_PI / 2.0 );
            }
            if ((*xc <= 0.0) || (*xc > 300) || (*yc <= 0.0) || (*yc > 300) || (*major > 500) || (*major <= 0.0)
                || (*minor > 500) || (*minor <= 0.0) || (alpha > 6.5) || (alpha < -6.5))
                return (false); //invalid ellipse -- can't be an valid image object
            return (true); //valid 
        }

		bool ImageProcessingUtilities_DSP::FindBestPupilObject(OBJECT objt[], int *no_of_objects, int *best_objt, FuzzParam *Pfuzz)
        {
            bool possible_pupil = false;
            int max_fuzz = 0;
            *best_objt = 0;
            int ii;
            clock_start("ImageProcessingUtilities_DSP::FindBestPupilObject");
            for (ii = 0; ii < *no_of_objects; ii++)
            {
                if (objt[ii].elps_fit)
                {
                    possible_pupil = true;
                    objt[ii].total_fuzz = pupil_fuzz(&objt[ii], Pfuzz);
                    if (objt[ii].total_fuzz > max_fuzz)
                    {
                        max_fuzz = objt[ii].total_fuzz;
                        *best_objt = ii;
                    }
                }
            }
            clock_end();
            return (possible_pupil);
        }

		int ImageProcessingUtilities_DSP::pupil_fuzz(OBJECT *objt, FuzzParam *Pfuzz)
        {

            int aspect, PD;
			//int smoothness;
			int total_fuzz = 0;
            aspect = (int)(objt->elps_aspect * 100 + 0.5);
            objt->fuzzAspect = calc_fuzz(aspect, &Pfuzz->Aspect);
            total_fuzz += objt->fuzzAspect;

            //smoothness = (int)(objt->smoothness * 1000 + 0.5); 
            //objt->fuzzSmoothness = calc_fuzz(smoothness, &Pfuzz->Smoothness);
            //total_fuzz += objt->fuzzSmoothness;

            PD = (int)(objt->scaled_diam * 100);
            objt->fuzzPD = calc_fuzz(PD, &Pfuzz->PD);
            total_fuzz += objt->fuzzPD;

            return (total_fuzz);
        }

		int ImageProcessingUtilities_DSP::calc_fuzz( int val, tFuzz *f)
        {
			int f_bin, f_remain, f_adjust, f_total;

            /* saturate the reading if value is larger than bins */
            if (ImageProcessingUtilities_FPGA::AbsInt(val) > f->binSize * (FUZZ_BINS - 1))
                val = (f->binSize * (FUZZ_BINS - 1)) - 1;
            {
                f_bin = val / f->binSize;
                f_remain = (int)((INT_MAX / FUZZ_MAX) * (int)(val % f->binSize)) / f->binSize;
                f_adjust = (f_remain * ((int)f->factor[f_bin + 1] - (int)f->factor[f_bin])) / (INT_MAX / FUZZ_MAX);
                f_total = f->factor[f_bin] + f_adjust;
                return (f_total * f->importance);
            }
        }

		void ImageProcessingUtilities_DSP::plot_ellipse(PNT ellipse_point[], int *no_list_elements, float *xc, float *yc, float *a, float *b, float *angle)
        {
            //a = major axis radius (in pixel units)
            //b = minor axis radius (in pixel units)
            //angle = major axis angle wrt horizontal
			bool repeat;
            float x, y, sin_angle, cos_angle, sin_t, cos_t, step_angle_increment, t;
            int ii, no_steps;
            clock_start("ImageProcessingUtilities_DSP::plot_ellipse");
            no_steps = (int)((6.0 * *a) + 0.5); // number of steps will be very roughly the number of pixels in the circumference
			if( no_steps > MAX_NO_ELLIPSE_PLOT_PNTS ) no_steps = MAX_NO_ELLIPSE_PLOT_PNTS; //don't exceed the fixed size of the ellpse_point array
            step_angle_increment = (float)( (2 * M_PI)/no_steps );
            t = 0.0; //radians
            sin_angle = sin(*angle);
            cos_angle = cos(*angle);
            for (ii = 0; ii < no_steps; ii++)
            {
                repeat = true;
                while (repeat)
                {
                    sin_t = sin(t);
                    cos_t = cos(t);
                    x = *xc + *a * cos_t * cos_angle - *b * sin_t * sin_angle;
                    y = *yc + *a * cos_t * sin_angle + *b * sin_t * cos_angle;
                    ellipse_point[ii].X = (int)(x + 0.5); //round off to nearest pixel
                    ellipse_point[ii].Y = (int)(y + 0.5);
                    t = t + step_angle_increment;
                    if (ii > 0)
                        if ((ellipse_point[ii].X == ellipse_point[ii - 1].X) && (ellipse_point[ii].Y == ellipse_point[ii - 1].Y))
                            repeat = true;
                        else repeat = false;
                    else repeat = false;
                    if (t > (2.0 * M_PI)) break;
                }
                if (t > (2.0 * M_PI)) break;
            }
            *no_list_elements = ii;
            clock_end();
        }

		float ImageProcessingUtilities_DSP::AbsFloat(float f)
		{
			if(f<0) return(-f);
			else return(f);
		}

    void ImageProcessingUtilities_DSP::LabeledRegions_from_array(int Array[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB, int *max_labels, 
		                                                          int *max_no_pnts_per_area, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int lbld_array[], int *no_of_regions, int *no_labels )
    {
        //array is a binary image array
        //labeled_area index is binary image array label minus 1 (array label zero is background, not a labeled region, and has no labeled area buffer list).
        //equivalence table indecies (column numbers and row numbers) are also binary image array label minus 1.
        int w = *Isizeh;
        int h = *Isizev;

        int diarray[4] = { -1, -1, 0, 1 };
        int djarray[4] = { 0, -1, -1, -1 };
        int /*label,*/ lowest_label, this_pixel_label, ij, id, jd, idjd;
        int n_label = 0;
        //bool sstop = false; //debug
        //int stpflg = 0;//debug
        //Cnct cnctns[];
        //EquivTable EquivToLabel[100];//equivalence table indecies (column numbers and row numbers) are also binary image array label minus 1.
        //LABELED_AREA labeled_area[100];
        //for( int n = 0; n < MAX_LABELS; n++)
        //{
        //    EquivToLabel[n].label = new bool[MAX_LABELS];
        //    labeled_area[n].pnt = new Point[MAX_NO_PNTS_PER_AREA];
        //}
        int cncted_labels[4];
        int n_cncted_labels, cl;

		//initialize
        EquivToLabel[0].label[0] = false;
        labeled_area[0].no_of_pnts = 0;
        
        for (int j = *wT; j <= *wB; j++)
            for (int i = *wL; i <= *wR; i++)
            {
                ij = i + (w * j);
				//if ((i == 175) && (j == 110))//debug
                //{                            //debug
                //    sstop = true;            //debug
                //}                            //debug
                if (Array[ij] > 0)
                {
                    //are any previously examained adjacent pixels labeled? (previously examained adjacent pixels are to left, upper left, above, or upper right of current point)
                    int array_label = 0;
                    lowest_label = *max_labels;
                    n_cncted_labels = 0;
                    cl = 0;
                    for (int n = 0; n < 4; n++)
                    {
                        id = i + diarray[n];
						jd = j + djarray[n];
						idjd = id + (w * jd);
						if( ( (i + diarray[n]) >= 0 ) && ((i+diarray[n]) < w) && ((j + djarray[n]) >= 0 ) && ((j+diarray[n]) < h) )
                            if (lbld_array[idjd] > 0) // if so, keep track of all adjacent labels in cncted_labels[] list. 
                            {
                                array_label = lbld_array[idjd];
                                if (array_label < lowest_label) lowest_label = array_label;
                                cncted_labels[cl] = array_label; // make list (cncted_labels[]) of labeled neigbors, and also keep track of lowest label no.  
                                ++cl;
                            }
                    }
                    if (cl > 0)// if there were any adjacent labels
                    {
                        lbld_array[ij] = lowest_label;
                        this_pixel_label = lowest_label;
                        for (int nn = 0; nn < cl; nn++)
                        {
                            if (cncted_labels[nn] != lowest_label)
                            {
                                EquivToLabel[lowest_label].label[cncted_labels[nn]] = true; // add other adjacent labels to equivalence table for this pixel's label.
                                EquivToLabel[cncted_labels[nn]].label[lowest_label] = true; //keep table symmetric
                            }
                        }
                    }
                    else //if no adjacent labels, assign next label to this pixel, increment label counter, and initialize equivalence table for the mew label. 
                    {
                        lbld_array[ij] = n_label + 1;
                        ++n_label;
                        this_pixel_label = n_label;
                        for (int nn = 0; nn < n_label; nn++) // initialize new column of equivalence table
                            EquivToLabel[n_label].label[nn] = false;
                        for (int cn = 0; cn < (n_label ); cn++) // initialize a new row on the equivalence table
                            EquivToLabel[cn].label[n_label ] = false;
                        labeled_area[n_label ].no_of_pnts = 0;//start new labeled area buffer index. 
						//++n_label;
                    }
                    // add this pixel to labeled area buffer
                    int pn = 0;
                    if( labeled_area[this_pixel_label].no_of_pnts > 0) 
                        pn = labeled_area[this_pixel_label].no_of_pnts;
                    if (pn < *max_no_pnts_per_area - 1)
                    {
                        labeled_area[this_pixel_label].pnt[pn].X = i;
                        labeled_area[this_pixel_label].pnt[pn].Y = j;
                        labeled_area[this_pixel_label].no_of_pnts = labeled_area[this_pixel_label].no_of_pnts + 1;
                    }
                }
                else lbld_array[ij] = 0; //the pixel is below threshold.  Set it to zero.  Make no addition to label buffer or equivalence table.
            } //end of first pass

        // Take a pass through equivalence table to be sure every column has the lowest label that it is equivalent to.  
        for( int nl = n_label-1; nl > 0; nl--)// check each column, starting with highest label column
        {
            for (int rn = 0; rn < n_label; rn++) //for each row in this column
            {
                if( EquivToLabel[nl].label[rn] ) //if another label (rn) is equiv to this col lablel (nl), find lowest equiv label in col rn
                    for( int nll = 0; nll <nl; nll++ )// start from the lowest row in column rn, and check only up to current column number (nl)
                    {
                        if( EquivToLabel[rn].label[nll] ) //if we come to an equivalence
                        {
                            EquivToLabel[nl].label[nll] = true; //be sure this label is set true in col nl
                            break; //leave for loop (first equivalence we come to is the lowest)
                        }
                    }
            }
        }
        //Now copy each label buffer into the lowest equivalent label buffer. If a lable is copied into a lower equivalent declare it empty by zeroing "no_of_points" for that label.
        int np; 
        *no_of_regions = n_label;
        for ( int lb = 1; lb < n_label; lb++) // start with 1 (no need to check lowest label). 
            for( int rn = 0; rn < lb; rn++ )
                if( EquivToLabel[lb].label[rn] )//copy label buf lb to end of label buf rn, and declare label buf lb empty
                {
                    for ( int nn = 0; nn < labeled_area[lb].no_of_pnts; nn++)
                    {
                        np = labeled_area[rn]. no_of_pnts;
                        if(np <  *max_no_pnts_per_area - 1)
                        {
                            labeled_area[rn].pnt[np].X = labeled_area[lb].pnt[nn].X;
                            labeled_area[rn].pnt[np].Y = labeled_area[lb].pnt[nn].Y;
                            labeled_area[rn]. no_of_pnts = labeled_area[rn]. no_of_pnts +1;
                        }
                    }
                    labeled_area[lb].no_of_pnts = 0;
                    -- *no_of_regions; // keep track of number of unique regions (not counting labels that are now empty).
                }
        *no_labels = n_label;
        //if (sstop)//debug
            //stpflg = 1;//debug
    }

    void ImageProcessingUtilities_DSP::LabeledRegions_from_thresh_buf(LBLD_AREA_BUF threshbuf[], int *threshbuf_size, int *max_labels, 
		                                                          int *max_no_pnts_per_area, LABELED_AREA labeled_area[], EquivTable EquivToLabel[], int *no_of_regions, int *no_labels )
    {
        //threshbuf is a list of pixels with gray scale values that were above a CR threshold.  Each array element has an X (column) and Y (row) address, and a label value. Label values are initially 0. 
        //labeled_area index is the pixel label.
        //equivalence table indecies (column numbers and row numbers) are also labels.

        int /*label,*/ lowest_label, this_pixel_label, cn, rn, cnn, rnn, ii;
        int n_label = 0;
        //bool sstop = false; //debug
        //int stpflg = 0;//debug
        int cncted_labels[4];
        int /*n_cncted_labels,*/ cl;

		//initialize
        EquivToLabel[0].label[0] = false;
        labeled_area[0].no_of_pnts = 0;

		if( *threshbuf_size > 0 ) 
		{
			threshbuf[0].label = 1;
			labeled_area[1].no_of_pnts = 1;
		}
        
        for (int i = 0; i < *threshbuf_size; i++)
        {
			//are any previously examained adjacent pixels labeled? (previously examained adjacent pixels are to left, upper left, above, or upper right of current point)
			cn = threshbuf[i].X;
			rn = threshbuf[i].Y;
			cl = 0;
			lowest_label = *max_labels;
			for( ii = i-1; ii >= 0; ii-- )
			{
				cnn = threshbuf[ii].X;
				rnn = threshbuf[ii].Y;
				if (rnn == rn)
				{
					if( ii == i - 1 )
						if( cnn == cn-1 )
						{
							lowest_label = threshbuf[ii].label;
							cncted_labels[cl] = threshbuf[ii].label;
							++cl;
						}
				}
				if (rnn < (rn -1) ) break;
				if (rnn == rn - 1 ) //if we have found one neigbor on same row, any other neighbors will be on previous row
				{
					if( ( cnn <= (cn + 1) ) && ( cnn >= (cn - 1) ) ) 
					{
						if( threshbuf[ii].label  < lowest_label ) 
							lowest_label = threshbuf[ii].label;
						cncted_labels[cl] = threshbuf[ii].label;
						++cl;
					}
				}
				if( cnn <= cn -1 ) break; 
			}
            if (cl > 0)// if there were any adjacent labels, use lowest and add others to equivalence table.
            {
                threshbuf[i].label = lowest_label;
                this_pixel_label = lowest_label;
				if( cl > 1)
					for (int nn = 0; nn < cl; nn++)
					{
						if (cncted_labels[nn] != lowest_label)
						{
							EquivToLabel[lowest_label].label[cncted_labels[nn]] = true; // add other adjacent labels to equivalence table for this pixel's label.
							EquivToLabel[cncted_labels[nn]].label[lowest_label] = true; //keep table symmetric
						}
					}
            }
            else //if no adjacent labels, assign next label to this pixel, increment label counter, and initialize equivalence table for the mew label. 
            {
                ++n_label;
				this_pixel_label = n_label;
				threshbuf[i].label = n_label;
                for (int nn = 0; nn < n_label; nn++) // initialize new column of equivalence table
                    EquivToLabel[n_label].label[nn] = false;
                for (int cn = 0; cn < (n_label ); cn++) // initialize a new row on the equivalence table
                    EquivToLabel[cn].label[n_label ] = false;
                labeled_area[n_label ].no_of_pnts = 0;//start new labeled area buffer index. 
            }
            // add this pixel to labeled area buffer
            int pn = 0;
            if( labeled_area[n_label].no_of_pnts > 0) 
                pn = labeled_area[this_pixel_label].no_of_pnts;
            if (pn < (*max_no_pnts_per_area - 1) )
			{
                labeled_area[this_pixel_label].pnt[pn].X = cn;
                labeled_area[this_pixel_label].pnt[pn].Y = rn;
				labeled_area[this_pixel_label].no_of_pnts = labeled_area[this_pixel_label].no_of_pnts + 1;
			}
        } //end of first pass

        // Take a pass through equivalence table to be sure every column has the lowest label that it is equivalent to.  
        for( int nl = n_label-1; nl > 0; nl--)// check each column, starting with highest label column
        {
            for (int rn = 0; rn < n_label; rn++) //for each row in this column
            {
                if( EquivToLabel[nl].label[rn] ) //if another label (rn) is equiv to this col lablel (nl), find lowest equiv label in col rn
                    for( int nll = 0; nll <nl; nll++ )// start from the lowest row in column rn, and check only up to current column number (nl)
                    {
                        if( EquivToLabel[rn].label[nll] ) //if we come to an equivalence
                        {
                            EquivToLabel[nl].label[nll] = true; //be sure this label is set true in col nl
                            break; //leave for loop (first equivalence we come to is the lowest)
                        }
                    }
            }
        }
        //Now copy each label buffer into the lowest equivalent label buffer. If a lable is copied into a lower equivalent declare it empty by zeroing "no_of_points" for that label.
        int np; 
        *no_of_regions = n_label;
        for ( int lb = 1; lb < n_label; lb++) // start with 1 (no need to check lowest label). 
            for( int rn = 0; rn < lb; rn++ )
                if( EquivToLabel[lb].label[rn] )//copy label buf lb to end of label buf rn, and declare label buf lb empty
                {
                    for ( int nn = 0; nn < labeled_area[lb].no_of_pnts; nn++)
                    {
                        np = labeled_area[rn].no_of_pnts;
                        if(np <  *max_no_pnts_per_area - 1)
                        {
                            labeled_area[rn].pnt[np].X = labeled_area[lb].pnt[nn].X;
                            labeled_area[rn].pnt[np].Y = labeled_area[lb].pnt[nn].Y;
                            labeled_area[rn]. no_of_pnts = labeled_area[rn]. no_of_pnts +1;
                        }
                    }
                    labeled_area[lb].no_of_pnts = 0;
                    -- *no_of_regions; // keep track of number of unique regions (not counting labels that are now empty).
                }
        *no_labels = n_label;
    }


	void ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis( OBJECT_RECOGNITION_PARAMETERS *orp_CR, int Array[], int *Isizeh, int *Isizev, int *wL, int *wR, int *wT, int *wB, 
		                                                                 LABELED_AREA labeled_area[], EquivTable EquivToLabel[], 
																		 int lbld_array[], AREA crArea[], int *no_cr_areas, int validAreaCRs[], int *no_validAreaCRs  )
	{
	 //****** Do Labeled Region analysis to find potential CRs ***********************//
        int no_unique_regions, no_labels, n;
        clock_start("ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis");
        LabeledRegions_from_array(Array, Isizeh, Isizev, wL, wR, wT, wB, &orp_CR->max_labels, &orp_CR->max_no_pnts_per_labeled_area, labeled_area, EquivToLabel, lbld_array, &no_unique_regions, &no_labels);
        //MakeLabeledArrayfromLabeledAreas(labeled_area, no_labels, OrigBMPw, OrigBMPh, out lbld_array2);//for debugging only
        int nn = 0;
        *no_cr_areas = 0;
       // int no_pnts_in_list, leftmost, rtmost, top, btm;
        for ( n = 1; n <= no_labels; n++ ) //note: there is no lable 0 (0 is background), so  max label index = no_labels (rather than no_labels - 1).
        {
            if (labeled_area[n].no_of_pnts > 0)
            {
                MakeEdgePointListFromLabeledArea(&labeled_area[n], crArea[nn].edgepnt, &orp_CR->max_no_pnts_per_labeled_area_edge,
                                                 &crArea[nn].no_edge_pnts, &crArea[nn].lft, &crArea[nn].rt, &crArea[nn].top, &crArea[nn].btm);
                crArea[nn].area = labeled_area[n].no_of_pnts;
				crArea[nn].no_elps_pnts = 1; //initialize this to 1 instead of 0 because the value will be used to allocate memory for managed memory calling program (can't allocate 0).
                nn++;
            }
            *no_cr_areas = nn;
        }
        FindValidCRsFromAreas(orp_CR, crArea, no_cr_areas, &orp_CR->max_labels, validAreaCRs, no_validAreaCRs);
        //FindBestCRTriadFromAreas();
        //FindBestCRArea();
        clock_end();
    }

	void ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis2( OBJECT_RECOGNITION_PARAMETERS *orp_CR, LBLD_AREA_BUF threshbuf[], int *threshbuf_size,  
		                                                                 LABELED_AREA labeled_area[], EquivTable EquivToLabel[], 
																		 AREA crArea[], int *no_cr_areas, int validAreaCRs[], int *no_validAreaCRs  )
	{
	 //****** Do Labeled Region analysis to find potential CRs ***********************//
        int no_unique_regions, no_labels, n;
        clock_start("ImageProcessingUtilities_DSP::FindCRsWithLabeledRegionAnalysis2");
		LabeledRegions_from_thresh_buf(threshbuf, threshbuf_size, &orp_CR->max_labels, &orp_CR->max_no_pnts_per_labeled_area, labeled_area, 
		                                                                  EquivToLabel, &no_unique_regions, &no_labels );
        //MakeLabeledArrayfromLabeledAreas(labeled_area, no_labels, OrigBMPw, OrigBMPh, out lbld_array2);//for debugging only
        int nn = 0;
        *no_cr_areas = 0;
       // int no_pnts_in_list, leftmost, rtmost, top, btm;
        for ( n = 1; n <= no_labels; n++ ) //note: there is no lable 0 (0 is background), so  max label index = no_labels (rather than no_labels - 1).
        {
            if (labeled_area[n].no_of_pnts > 0)
            {
                MakeEdgePointListFromLabeledArea(&labeled_area[n], crArea[nn].edgepnt, &orp_CR->max_no_pnts_per_labeled_area_edge,
                                                 &crArea[nn].no_edge_pnts, &crArea[nn].lft, &crArea[nn].rt, &crArea[nn].top, &crArea[nn].btm);
                crArea[nn].area = labeled_area[n].no_of_pnts;
				crArea[nn].no_elps_pnts = 1; //initialize this to 1 instead of 0 because the value will be used to allocate memory for managed memory calling program (can't allocate 0).
                nn++;
            }
            *no_cr_areas = nn;
        }
        FindValidCRsFromAreas(orp_CR, crArea, no_cr_areas, &orp_CR->max_labels, validAreaCRs, no_validAreaCRs);
        //FindBestCRTriadFromAreas();
        //FindBestCRArea();
        clock_end();
    }


    void ImageProcessingUtilities_DSP::MakeEdgePointListFromLabeledArea( LABELED_AREA *L, PNT list[], int *max_no_pnts_per_cr_edge, int *no_pnts_in_list, int *leftmost, int *rtmost, int *top, int *btm )
    {
        int ln, cp, n;
        bool start_looking_for_rt_edge, out_of_order_element;
        ln = 0;
        cp = 0;
        start_looking_for_rt_edge = true; //0 when first point on new line; 1 when we have already entered a left edge point on this line
        list[0].X = L->pnt[0].X;
        list[0].Y = L->pnt[0].Y;
        *leftmost = 10000;
        *rtmost = 0;
        for ( n = 0; n < (L->no_of_pnts - 1); n++ )
        {

            if (L->pnt[n + 1].Y > L->pnt[n].Y)//if n+1 is a new line, advance the list pointer
            {
                ln++;
                start_looking_for_rt_edge = true;
                out_of_order_element = false;
            }
            else
                if (L->pnt[n + 1].Y == L->pnt[n].Y)
                {
                    if (start_looking_for_rt_edge)
                    //if we just entered a left edge point and are starting to look for a rt edge, advance the list pointer.
                    //If we are replacing a rt edge point with a point further to the right (on the same line), do not advance the list pointer.
                    {
                        ln++;
                        start_looking_for_rt_edge = false;
                    }
                    out_of_order_element = false;
                }
                else //next buffer entry has a lower line number (buffer element is "out of order")
                {
                    out_of_order_element = true;
                    for(int lnn = ln; lnn > 0; lnn--)//change >= to >; jb 10-29-10
                        if (list[lnn].Y == L->pnt[n + 1].Y) //found edge list entry from same line as next buffer entry
                        {
                            if( list[lnn-1].Y == L->pnt[n+1].Y) //if previous edge list entry also on same line, then list[lnn] must be a rt edge
                                if (L->pnt[n + 1].X > list[lnn].X) //if buffer entry [n+1] is to rt of list{lnn], replace list[lnn] with buffer entry [n+1]. 
                                {
                                    list[lnn].X = L->pnt[n + 1].X;
                                    list[lnn].Y = L->pnt[n + 1].Y;
                                }
                                else
                                {
                                    if (L->pnt[n + 1].X < list[lnn - 1].X) // if buf entry [n+1] is to left of list[lnn-1], replace list[lnn-1] with buffer entry [n+1].
                                    {
                                        list[lnn-1].X = L->pnt[n + 1].X;
                                        list[lnn-1].Y = L->pnt[n + 1].Y;
                                    }
                                }
                        }
                }

            if (!out_of_order_element)
            {
                list[ln].X = L->pnt[n + 1].X;
                list[ln].Y = L->pnt[n + 1].Y;
            }
            if (L->pnt[n + 1].X < *leftmost) *leftmost = L->pnt[n + 1].X;
            if (L->pnt[n + 1].X > *rtmost) *rtmost = L->pnt[n + 1].X;
            if (ln >= (*max_no_pnts_per_cr_edge-1) ) break;
        }
        *no_pnts_in_list = ln + 1;
        *top = list[0].Y;
        *btm = list[ln].Y;
    }

	void ImageProcessingUtilities_DSP::FindValidCRsFromAreas(OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, AREA crArea[], int *no_CR_areas, int *max_labels, int validAreaCRs[], int *no_validAreaCRs)
    {
        double height, width;
        int max_pixels = (int)(CR_rec_param->area_max_area_size * (CR_rec_param->EyeImageScaleFactor * CR_rec_param->EyeImageScaleFactor) );
        int min_pixels = (int)(CR_rec_param->area_min_area_size * (CR_rec_param->EyeImageScaleFactor * CR_rec_param->EyeImageScaleFactor) );
        if (max_pixels < 1) max_pixels = 1;
        if (min_pixels < 1) min_pixels = 1;
        int n = 0;
        for (int i = 0; i < *no_CR_areas; i++)
        {
            if ((crArea[i].area < max_pixels) && (crArea[i].area > min_pixels))
            {
                height = (double)(crArea[i].btm - crArea[i].top);
                width = (double)(crArea[i].rt - crArea[i].lft);
                if (height > width) crArea[i].rough_aspect = height / width;
                else crArea[i].rough_aspect = width / height;
                if (crArea[i].rough_aspect <= CR_rec_param->area_rough_aspect_max)
                {
                    validAreaCRs[n] = i;
                    n++;
                }
            }
        }
        *no_validAreaCRs = n;
    }


	void ImageProcessingUtilities_DSP::FitElipseToLabeledAreaCRs( OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, AREA crArea[], int validAreaCRs[], int *no_validAreaCRs )
	{
        //***** fit elipse to "valid" CR areas found from labeled region algorithm.  For the time being, just show all of them. ********************
		int CRemax;
		clock_start("ImageProcessingUtilities_DSP::FitElipseToLabeledAreaCRs");
        if (*no_validAreaCRs > CR_rec_param->feature_max_no_of_elps_computations) CRemax = CR_rec_param->feature_max_no_of_elps_computations;
        else CRemax = *no_validAreaCRs;

        if (CRemax > 0)
        {
            for (int n = 0; n < CRemax; n++)
            {
                if (ComputeEllipse(crArea[validAreaCRs[n]].edgepnt, &crArea[validAreaCRs[n]].no_edge_pnts, &crArea[validAreaCRs[n]].hc, &crArea[validAreaCRs[n]].vc,
                    &crArea[validAreaCRs[n]].elps_maj_rad, &crArea[validAreaCRs[n]].elps_min_rad, &crArea[validAreaCRs[n]].elps_angle))
                {
                    crArea[validAreaCRs[n]].elps_aspect = crArea[validAreaCRs[n]].elps_min_rad / crArea[validAreaCRs[n]].elps_maj_rad;
                    crArea[validAreaCRs[n]].elps_fit = true;
                    //no_of_CRelps_found++;
                }
                else //too few points for valid ellipse fit, or ellipse fit failed for some other reason. Just find simple center (half way between top & btm, and left & rt).
                {

                    crArea[validAreaCRs[n]].vc = crArea[validAreaCRs[n]].top + (crArea[validAreaCRs[n]].btm - crArea[validAreaCRs[n]].top) / 2.0f;
                    crArea[validAreaCRs[n]].hc = crArea[validAreaCRs[n]].lft + (crArea[validAreaCRs[n]].rt - crArea[validAreaCRs[n]].lft) / 2.0f;
                    if ((crArea[validAreaCRs[n]].btm - crArea[validAreaCRs[n]].top) > (crArea[validAreaCRs[n]].rt - crArea[validAreaCRs[n]].lft))
                    {
                        crArea[validAreaCRs[n]].elps_maj_rad = (float)(crArea[validAreaCRs[n]].btm - crArea[validAreaCRs[n]].top);
                        crArea[validAreaCRs[n]].elps_angle = 6.28319f; //2.0 * 3.14159;
                    }
                    else
                    {
                        crArea[validAreaCRs[n]].elps_maj_rad = (float)(crArea[validAreaCRs[n]].rt - crArea[validAreaCRs[n]].lft);
                        crArea[validAreaCRs[n]].elps_angle = 0.0;
                    }
                }
            }
        }
        clock_end();
	}

	
	void ImageProcessingUtilities_DSP::FitElipseToEdgeObjectCRs( CR_OBJECT CR_objt[], OBJECT_RECOGNITION_PARAMETERS *CR_rec_param, LFT_EDGE lft_CR_edge[], RT_EDGE rt_CR_edge[],
		                                                         EDGE_STRING CR_es[], PNT CRbuf[], int *CRbuf_size, int *no_of_CR_objects, int validCRs[], int *noValidCRs )
	{
		int CRemax, i, no_of_elps_found;
		int list[100];
		clock_start("ImageProcessingUtilities_DSP::FitElipseToEdgeObjectCRs");
		if (*no_of_CR_objects > CR_rec_param->feature_max_no_of_elps_computations) CRemax = CR_rec_param->feature_max_no_of_elps_computations;
		if (CRemax > 100) CRemax = 100;
		else CRemax = *no_of_CR_objects;
		no_of_elps_found = 0;
		*noValidCRs = 0;
		if (CRemax > 0)
		{
			for (int n = 0; n < CRemax; n++)
			{
				if( ImageProcessingUtilities_DSP::FitCREllipse( &CR_objt[n], lft_CR_edge, rt_CR_edge, CR_es, CRbuf, CRbuf_size, CR_rec_param) )
				{
					CR_objt[n].scaled_diam = (2 * CR_objt[n].elps_maj_rad) / CR_rec_param->EyeImageScaleFactor;
					list[no_of_elps_found] = n;
					++no_of_elps_found;
				}
			}
		   //CR_object_no = ImageProcessingUtilities_DSP::FindBestCRObject(CR_objt, &no_of_CR_objects, &CR_object_no, &CRfuzz);
			for( i = 0; i < no_of_elps_found; i++ )
			{
				if (CR_objt[list[i]].elps_fit)
				{
					if( ImageProcessingUtilities_DSP::AcceptableCR( &CR_objt[list[i]], CR_rec_param ) )
					{
						validCRs[i] = list[i];
						++(*noValidCRs);
					}
				}
			}
		}
		clock_end();
	}

    bool ImageProcessingUtilities_DSP::AcceptableCR( CR_OBJECT *CR_objt, OBJECT_RECOGNITION_PARAMETERS *CR_rec_param )
	{
		bool acceptable = false;
		if( CR_objt->scaled_diam <= CR_rec_param->feature_max_scaled_diam )
			if( CR_objt->elps_aspect <= CR_rec_param->elps_aspect_max )
				acceptable = true;
		return( acceptable );
	}

























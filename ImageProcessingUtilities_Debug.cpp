#include "ImageProcessingUtilities_Debug.h"

ImageProcessingUtilities_Debug::ImageProcessingUtilities_Debug(void)
{
}

ImageProcessingUtilities_Debug::~ImageProcessingUtilities_Debug(void)
{
}

        void ImageProcessingUtilities_Debug::MakeReadableStringArrayForDebugging(EDGE_STRING es[], int *nstring, PNT PupilEdgeBuf[], int *pbufsize, STRING_ARRAY S[])
        {
            int iimax, jjmax;
			if (*nstring > 500) iimax = 500;
			else iimax = *nstring;
			for (int ii = 0; ii < iimax; ii++)
            {
				if( es[ii].pixel_lngth > 500 ) jjmax = 500;
				else jjmax = es[ii].pixel_lngth;
				for (int jj = 0; jj < es[ii].pixel_lngth; jj++)
                {
                    S[ii].edgepoint[jj].X = PupilEdgeBuf[es[ii].edge_buf_index[jj]].X;
                    S[ii].edgepoint[jj].Y = PupilEdgeBuf[es[ii].edge_buf_index[jj]].Y;
                }
            }
        }

        void ImageProcessingUtilities_Debug::MakeReadableDirArrayForDebugging(EDGE_STRING es[], int *nstring, int DirArray[], int *pbufsize, STRNG_DIR_ARRAY D[])
        {
			int iimax, jjmax;
			if( *nstring > 500 ) iimax = 500;
			else iimax = *nstring;
            for (int ii = 0; ii < iimax; ii++)
            {
				if( es[ii].pixel_lngth > 500 ) jjmax = 500;
				else jjmax = es[ii].pixel_lngth;
                for (int jj = 0; jj < jjmax; jj++)
                {
                    D[ii].dir[jj] = DirArray[es[ii].edge_buf_index[jj]];
                }
            }
        }

        void ImageProcessingUtilities_Debug::MakeReadableHgradVgradStringForDebugging(EDGE_STRING es[], int *nstring, PNT PupilEdgeBuf[], int grad_H_array[], int grad_V_array[], 
			                                                                          int *Isizeh, int *Isizev, STRING_ARRAY Ghv[])
        {
			int ij, iimax, jjmax;
			if(*nstring > 500) iimax = 500;
			else iimax = *nstring;
            for (int ii = 0; ii < *nstring; ii++)
            {
				if( es[ii].pixel_lngth > 500 ) jjmax = 500;
				else jjmax = es[ii].pixel_lngth;
                for (int jj = 0; jj < es[ii].pixel_lngth; jj++)
                {
					ij = PupilEdgeBuf[es[ii].edge_buf_index[jj]].X + ( *Isizeh *  PupilEdgeBuf[es[ii].edge_buf_index[jj]].Y );
                    Ghv[ii].edgepoint[jj].X = grad_H_array[ij];
                    Ghv[ii].edgepoint[jj].Y = grad_V_array[ij];
                }
            }
        }


#ifndef ImageProcessingUtilities_Debug_h
#define ImageProcessingUtilities_Debug_h

#pragma once

#include "ProcessImage.h"

class ImageProcessingUtilities_Debug
{
public:
	ImageProcessingUtilities_Debug(void);
	~ImageProcessingUtilities_Debug(void);
	 static void MakeReadableStringArrayForDebugging(EDGE_STRING es[], int *nstring, PNT PupilEdgeBuf[], int *pbufsize, STRING_ARRAY S[] );
	 static void MakeReadableDirArrayForDebugging( EDGE_STRING es[], int *nstring, int DirArray[], int *pbufsize, STRNG_DIR_ARRAY D[] );
	 static void MakeReadableHgradVgradStringForDebugging
		 ( EDGE_STRING es[], int *nstring, PNT PupilEdgeBuf[], int grad_H_array[], int grad_V_array[], int *Isizeh, int *Isizev, STRING_ARRAY Ghv[] );

};

#endif

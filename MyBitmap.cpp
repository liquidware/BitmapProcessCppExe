#include "MyBitmap.h"

CMyBitmap::CMyBitmap(void)
{
}

CMyBitmap::~CMyBitmap(void)
{
}

BOOL CMyBitmap::LoadBitmap(LPCTSTR szFilename) 
{ 
	//ASSERT(szFilename);
	//DeleteObject();

	//HBITMAP hBitmap = NULL;
	//hBitmap = (HBITMAP)LoadImage(NULL, szFilename, IMAGE_BITMAP, 0, 0,
	//	LR_LOADFROMFILE | LR_CREATEDIBSECTION | LR_DEFAULTSIZE);
	//return Attach(hBitmap);
	return true;
}

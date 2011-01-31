#ifndef MyBitmap_h
#define MyBitmap_h

#pragma once
typedef int BOOL;
#define LPCTSTR char*

class CMyBitmap
{
public:
	CMyBitmap(void);
	~CMyBitmap(void);

	BOOL LoadBitmap(LPCTSTR szFilename); 
};

#endif

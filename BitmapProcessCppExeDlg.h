// BitmapProcessCppExeDlg.h : header file
//

#pragma once
#include "afxwin.h"

// definition of PNT structure
#include "ProcessImage.h"

// CBitmapProcessCppExeDlg dialog
class CBitmapProcessCppExeDlg : public CDialog
{
// Construction
public:
	CBitmapProcessCppExeDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_BITMAPPROCESSCPPEXE_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	CString m_imageFilename;
	afx_msg void OnBnClickedBtnBrowse();
	afx_msg void OnBnClickedBtnDoProcess();
	CStatic m_picture1;
	CStatic m_picture2;

private:
	void BrightPupilProcessing(CBitmap* pBitmap);
	int CreatePixelArray(BYTE* in, int** out, int width, int height, int bytesPerPixel);
	void SavePixelArray(int arr[], int arrSize, int width, int height, int bytesPerPixel);
	int CopyBytes(BYTE** pb, int size);
	void DrawCross(BYTE* bits, int width, int height, int crossHorz, int crossVert, int crossSize, int R, int G, int B);
	void DrawBorder(BYTE* bits, int width, int height, PNT border[], int borderCount, int R, int G, int B);
};

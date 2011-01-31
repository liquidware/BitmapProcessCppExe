// BitmapProcessCppExeDlg.cpp : implementation file
//

#include "stdafx.h"
#include "BitmapProcessCppExe.h"
#include "BitmapProcessCppExeDlg.h"
#include "MyBitmap.h"
#include "ProcessImage.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CBitmapProcessCppExeDlg dialog




CBitmapProcessCppExeDlg::CBitmapProcessCppExeDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CBitmapProcessCppExeDlg::IDD, pParent)
	, m_imageFilename(_T(""))
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CBitmapProcessCppExeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_IMAGE_FILE, m_imageFilename);
	DDX_Control(pDX, IDC_PICTURE_1, m_picture1);
	DDX_Control(pDX, IDC_PICTURE_2, m_picture2);
}

BEGIN_MESSAGE_MAP(CBitmapProcessCppExeDlg, CDialog)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BTN_BROWSE, &CBitmapProcessCppExeDlg::OnBnClickedBtnBrowse)
	ON_BN_CLICKED(IDC_BTN_DO_PROCESS, &CBitmapProcessCppExeDlg::OnBnClickedBtnDoProcess)
END_MESSAGE_MAP()


// CBitmapProcessCppExeDlg message handlers

BOOL CBitmapProcessCppExeDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	//m_imageFilename = "C:\\Users\\gkhrapunovich.ASL\\Documents\\Visual Studio 2008\\Projects\\ASL_bmp_processing_program\\bright_pupil_test_image\\Katie Good.bmp";
	//m_imageFilename = "C:\\Temp\\KatieGood.bmp";
	//m_imageFilename = "C:\\Temp\\Test.bmp";
	m_imageFilename = "C:\\asl\\Misc\\ASL_bmp_processing_program\\bright_pupil_test_images\\VAS_1.BMP";
	UpdateData(false);

	// Browse
	PostMessage(WM_COMMAND, MAKEWPARAM(IDC_BTN_BROWSE, BN_CLICKED));

	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CBitmapProcessCppExeDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CBitmapProcessCppExeDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CBitmapProcessCppExeDlg::OnBnClickedBtnBrowse()
{
	UpdateData(true);

	CFileDialog dlg(
		TRUE, // open
		NULL, // default extension
		m_imageFilename
		);
	if (dlg.DoModal() == IDOK)
	{
		m_imageFilename = dlg.GetPathName();
		UpdateData(false);

		// Display image
		CMyBitmap bmp;
		BOOL res = bmp.LoadBitmap(m_imageFilename);
		if (res == FALSE)
		{
			MessageBox("Error opening image file", "", MB_OK | MB_ICONERROR);
			return;
		}
		m_picture1.SetBitmap(bmp);
		m_picture2.SetBitmap(NULL);
	}
}

void CBitmapProcessCppExeDlg::OnBnClickedBtnDoProcess()
{
	// Load original image
	CMyBitmap bmp;
	BOOL res = bmp.LoadBitmap(m_imageFilename);
	if (res == FALSE)
	{
		MessageBox("Error opening image file", "", MB_OK | MB_ICONERROR);
		return;
	}

	// Process image
	BrightPupilProcessing(&bmp);
	m_picture2.SetBitmap(bmp);
}

void CBitmapProcessCppExeDlg::BrightPupilProcessing(CBitmap* pBitmap)
{
	// Create out bitmap with same dimetions, only 4-byte (int) per pixel
	BITMAP bmp;
	pBitmap->GetBitmap(&bmp);
	//ASSERT (bmp.bmBitsPixel == 24); // 3 bytes per pixel in old bitmap
	int bytesPerPixel = bmp.bmBitsPixel / 8;
	int width = bmp.bmWidth;
	int height = bmp.bmHeight;

	// Create pixel array
	int* intArr;
	int arrSize = CreatePixelArray((BYTE*)bmp.bmBits, &intArr, width, height, bytesPerPixel);

	//** temp
	SavePixelArray(intArr,arrSize, width, height, bytesPerPixel);

	// Find pupil and CR
	PUPIL_CR_OBJECT pupil;
	PUPIL_CR_OBJECT CR[MAX_CR_OBJECTS];
	int number_CR_objects;
	FindPupilCR(width, height, intArr, &pupil, CR, &number_CR_objects);

	// Draw crosses
	BYTE* inpBytes = (BYTE*)bmp.bmBits;
	if (pupil.found)
	{
		DrawCross(inpBytes, width, height, (int)pupil.center_horz, (int)pupil.center_vert, 10, 255, 0, 0);
		DrawBorder(inpBytes, width, height, pupil.ellipse_points, pupil.no_ellipse_points, 255, 0, 0);
	}

	if (number_CR_objects > 0)
	{
		for (int i=0; i<number_CR_objects; i++)
		{
			DrawCross(inpBytes, width, height, (int)CR[i].center_horz, (int)CR[i].center_vert, 10, 100, 255, 0);
		}
	}

	delete[] intArr;
}

int CBitmapProcessCppExeDlg::CreatePixelArray(BYTE* in, int** out, int width, int height, int bytesPerPixel)
{
	int pixelSize = width * height;
	*out = new int[pixelSize];
	for (int row=0; row<height; row++)
	{
		for (int col=0; col<width; col++)
		{
			int index1 = ((height - 1 - row) * width + col) * bytesPerPixel;
			BYTE* pIn = in + index1;
			int val = CopyBytes(&pIn, bytesPerPixel);
			int index2 = row * width + col;
			(*out)[index2] = val;
		}
	}
	return pixelSize;
}


int CBitmapProcessCppExeDlg::CopyBytes(BYTE** pb, int size)
{
	int b = **pb;
	(*pb)++;
	int g = **pb;
	(*pb)++;
	int r = **pb;
	(*pb)++;

	if (size == 4)
	{
		int a = **pb;
		(*pb)++;
	}
	// Build greyscale image
	return (int)(r * 0.299 + g * 0.587 + b * 0.114);
}

void CBitmapProcessCppExeDlg::DrawCross(BYTE* bits, int width, int height,
		int crossHorz, int crossVert, int crossSize, int R, int G, int B)
{
	// Invert verical coordinate to acount for different row order
	int crossV2 = height - crossVert - 1;
	
	// Horiz line
	for (int i = crossHorz - crossSize; i < crossHorz + crossSize; i++)
	{
		int index = 3 * (i + crossV2 * width);
		bits[index] = B;		// Blue
		bits[index + 1] = G;	// Green 
		bits[index + 2] = R;	// Red
	}
	// Vert line
	for (int i = crossV2 - crossSize; i < crossV2 + crossSize; i++)
	{
		int index = 3 * (crossHorz + i * width);
		bits[index] = B;
		bits[index + 1] = G;
		bits[index + 2] = R;
	}
}

void CBitmapProcessCppExeDlg::DrawBorder(BYTE* bits, int width, int height,
		PNT border[], int borderCount, int R, int G, int B)
{
	for (int i=0; i<borderCount; i++)
	{
		PNT* pt = &(border[i]);
		int index = 3 * (pt->X + (height - pt->Y - 1) * width);
		bits[index] = B;
		bits[index + 1] = G;
		bits[index + 2] = R;
	}
}

void CBitmapProcessCppExeDlg::SavePixelArray(int arr[], int arrSize, int width, int height, int bytesPerPixel)
{
	// Obtain filename from orig bitmap file
	CString filename = m_imageFilename + ".pxl";
	CFile file(filename, CFile::modeCreate | CFile::modeWrite);
	for (int i=0; i<arrSize; i++)
	{
		unsigned short val = (arr[i] & 0xFF) << 8;
		file.Write(&val, 2);
	}
	file.Close();
}

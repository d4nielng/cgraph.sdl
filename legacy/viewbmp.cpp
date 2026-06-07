// ======================================================================
// VIEWBMP.CPP
// ----------------------------------------------------------------------
// Simple Windows Bitmap Viewer - v0.4
// Support reading of 
//	. uncompressed monochrome bitmap
//	. uncompressed 16 colors bitmap
//	. uncompressed 256 colors bitmap
//	. uncompressed 24-bit and 32-bit bitmap
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 viewbmp.cpp -o viewbmp -s -lgdi32
// ======================================================================

#include "include/cgraph.h"

using namespace std;
using namespace CGraph;

ConsoleGraphics g;
CBitmap bmp;
bool stretch = false;

string	gauge_text[] = {VT_COLOR(220)+string(">")+VT_COLOR(240)+"---------", 
						VT_COLOR(220)+string(">>")+VT_COLOR(240)+"--------",
						VT_COLOR(220)+string(">>>")+VT_COLOR(240)+"-------",
						VT_COLOR(220)+string(">>>>")+VT_COLOR(240)+"------",
						VT_COLOR(220)+string(">>>>>")+VT_COLOR(240)+"-----",
						VT_COLOR(220)+string(">>>>>>")+VT_COLOR(240)+"----",
						VT_COLOR(220)+string(">>>>>>>")+VT_COLOR(240)+"---",
						VT_COLOR(220)+string(">>>>>>>>")+VT_COLOR(240)+"--",
						VT_COLOR(220)+string(">>>>>>>>>")+VT_COLOR(240)+"-",
						VT_COLOR(154)+string(">>>>>>>>>>")+VT_COLOR(240)};

void callback(int current, int total) {
	if (current % 10 == 0) {
		cout << "- loading " << gauge_text[current*10/total] << " " << VT_COLOR(196) << (current+1)*100/total << "%" << VT_DEFAULTATTR << "\r";	
	}
}

void refresh() {
	system("cls");
	system("dir images /w");
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    switch (uMsg) {
    case WM_KEYDOWN:
		if (wParam == VK_ESCAPE || wParam == VK_RETURN) {
			bmp.Clear();
		    DestroyWindow(hwnd);
			PostQuitMessage(0);
		}
		if (wParam == 0x53) {	// S key to toggle stretching
			stretch = !stretch;
			PostMessage(hwnd, WM_PAINT, 0, 0);
		}
		break;
    case WM_SIZE:
	    bmp.Render(hwnd, stretch);
        break;
	case WM_PAINT:
		bmp.Render(hwnd, stretch);
	    break;
	case WM_CLOSE:
		bmp.Clear();
	    DestroyWindow(hwnd);
		PostQuitMessage(0);
	    return 0;			    
    default:
    	return DefWindowProc(hwnd, uMsg, wParam, lParam);
	}    
	return 1;
}

void ShowBitmap(const CBitmap & img) {
	const char CLASS_NAME[]  = "Bitmap Viewer";
	const char TITLE[]  = "Bitmap Viewer";
	WNDCLASS wc = { };
	HINSTANCE hInstance = GetModuleHandle(NULL);
	HWND hwnd;
	int th, cw;
		
	wc.style		 = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc   = WindowProc;
	wc.hInstance     = hInstance;
	wc.lpszClassName = CLASS_NAME;
	wc.hCursor		 = LoadCursor(NULL, IDC_ARROW);
	RegisterClass(&wc);

	hwnd = CreateWindowEx(0, CLASS_NAME, TITLE,
	    WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 
		CW_USEDEFAULT, CW_USEDEFAULT,
	    NULL, NULL, hInstance, NULL);

	if (hwnd == NULL) return;

	cw = GetSystemMetrics(SM_CXSIZEFRAME) * 2;
	th = GetSystemMetrics(SM_CYCAPTION);
	ShowWindow(hwnd, SW_SHOW);
 	SetWindowPos(hwnd, 0, 0, 0, 
	 			 img.Width() + cw, 
				 img.Height() + th + cw, SWP_NOMOVE);
	SetForegroundWindow(hwnd);
	
	MSG msg = { };
    while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }		
}

int main() {
	string filename;
	bool quit = false;
		
	while (!quit) {
		refresh();
		cout << endl << "stretch = " << (stretch ? (VT_COLOR(222)+"yes") : (VT_COLOR(154) + "no")) << VT_DEFAULTATTR << endl;
		cout << VT_COLOR(220) << ":q" << VT_DEFAULTATTR << " = quit, " 
		     << VT_COLOR(220) << ":r" << VT_DEFAULTATTR << " = refresh file list " 
			 << VT_COLOR(220) << ":?" << VT_DEFAULTATTR << " = short help" << endl;
		cout << ">";	// otherwise, tell the users to input the desired file name
		cin >> filename;
		cin.clear();			// clear input buffers
		cin.ignore(INT_MAX,'\n'); 
	
		if (filename == ":q") {
			quit = true;
			break;
		}

		if (filename == ":r") {
			continue;
		}

		if (filename == ":?") {
			MessageBox(	NULL, 
						"Type in the bitmap image on the list to view\n"
						"Commands:\n"
						":r - refresh file list\n"
						":? - display this help message\n"
						":q - quit this app", "Help", MB_OK);
			continue;
		}
		
		if (filename.find(".bmp") == string::npos) 
			filename += ".bmp";
			
		// attempt to read the bitmap file
		if (bmp.Load("images/"+filename, callback)) {
			cout << endl;
			bmp.Info();
			ShowBitmap(bmp);	// and show it if success
		}
		else {
			// otherwise, report error
			cout << "Cannot load image file from \"" << VT_COLOR(196) << filename << VT_DEFAULTATTR << "\"" << endl;
			ConsoleGraphics::wait();
		}
	}

	return 0;
}

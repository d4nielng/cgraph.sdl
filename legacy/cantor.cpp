// ======================================================================
// CANTOR.CPP
// ----------------------------------------------------------------------
// Cantor Set Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 cantor.cpp -o cantor -s -lgdi32
// ======================================================================

#include "include/cgraph.h"

using namespace CGraph;

#define	WIDTH		640
#define	HEIGHT		400
#define	THICKNESS	(HEIGHT/(DEPTH*2))
#define	SIZE		(WIDTH - 10)
#define	GAP			(THICKNESS * 2)
#define	DEPTH		7
#define	CX			(WIDTH - SIZE) / 2
#define	CY			(HEIGHT - DEPTH * GAP) / 2
#define HUE_STEP	(360/(DEPTH+1))

ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);

void cantor(CBitmap & bmp, int x, int y, int length, int depth) {
	if (depth == DEPTH)
		return;
	
	double v = depth;
	
	// map color to depth
	CHSV hsv(CMath::map(v, 0.0, (double) DEPTH, 0.0, 360.0), 0.6f, 1.0f);
	bmp.SetColor(hsv.ToRGB());
	
	// new length is one third of previous length
	int newlength = length / 3;
	
	// draw the bar
	bmp.RectangleD(x, y, length, THICKNESS);
	
	// progress recursively down
	cantor(bmp, x              , y + GAP, newlength, depth + 1);
	cantor(bmp, x + 2*newlength, y + GAP, newlength, depth + 1);
}

void render() {
	// render the initial cantor set using pastel-like color
	cantor(bmp, CX, CY, SIZE, 0);
	Sleep(100);
	bmp.Render(g);
//	bmp.Write("cantor.bmp");
}

int main() {
	render();
	ConsoleGraphics::wait();
	return 0;
}

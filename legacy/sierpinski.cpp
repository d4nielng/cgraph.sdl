// ======================================================================
// SIERPINSKI.CPP
// ----------------------------------------------------------------------
// Sierpinksi Gasket Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 sierpinski.cpp -o sierpinksi -s -lgdi32
// ======================================================================
#include "include/cgraph.h"

using namespace CGraph;

#define	WIDTH		640
#define	HEIGHT		400
#define MIN_HEIGHT	8

const float SQ3 = sqrt(3);

void triangle(CBitmap & bmp, float x, float y, float h) {
    bmp.DrawLine(x - h / SQ3, y - h / 3, x + h / SQ3, y - h / 3);
    bmp.DrawLine(x - h / SQ3, y - h / 3, x, y + 2 * h / 3);
    bmp.DrawLine(x, y + 2 * h / 3, x + h / SQ3, y - h / 3);
}
  
void sierpinski(CBitmap & bmp, float x, float y, float h, CGraph::CRGB color) { 
    if (h < MIN_HEIGHT) return;

	CHSV hsv(CMath::map(y, 0, HEIGHT, 0, 360), 1.0, 1.0);
    bmp.SetColor(hsv.ToRGB());

    triangle(bmp, x, y, h);
    
    sierpinski(bmp, x          , y - 2 * h / 3, h / 2, color);
    sierpinski(bmp, x - h / SQ3, y + h / 3    , h / 2, color);
    sierpinski(bmp, x + h / SQ3, y + h / 3    , h / 2, color);
}

ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);

void render() {
	sierpinski(bmp,
			   bmp.Width() / 2, 
			   2 * bmp.Height() / 3 - 4, 
			   bmp.Height() / 2, 
			   CGraph::CRGB(1.0, 0.5, 1.0));
	Sleep(100);
	bmp.Render(g);
//	bmp.Write("sierpinski.bmp");
}

int main() {	
	render();
	ConsoleGraphics::wait();
	return 0;
}

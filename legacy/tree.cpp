// ======================================================================
// TREE.CPP
// ----------------------------------------------------------------------
// Fractal Tree Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 tree.cpp -o tree -s -lgdi32
// ======================================================================
#include "include/cgraph.h"

#define	WIDTH		640
#define	HEIGHT		400
#define SIZE        WIDTH
#define SCALE		5
#define BRANCHES    10
#define ROTSCALE	0.85
#define INITLENGTH	HEIGHT/7

using namespace CGraph;
ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);
 
void tree(	double ofsx, double ofsy,  double dirx, double diry, 
			double size, double rot, int depth) {
	
	int v = depth * 3;
	CHSV hsv(CMath::map(v, 0, 48, 0, 360), 0.8, 1.0);
	bmp.SetColor(hsv.ToRGB());
	
	bmp.MoveTo((int)ofsx, (int)ofsy);
	bmp.LineTo((int)(ofsx + dirx * size), (int)(ofsy + diry * size));
	if (depth > 0) {
	    // left branch
	    tree(	ofsx + dirx * size, ofsy + diry * size, 
				dirx * cos(rot) + diry * sin(rot), 
				dirx * -sin(rot) + diry * cos(rot), 
				size * 0.5/ SCALE + size * (SCALE - 1) / SCALE,
		        rot * ROTSCALE,  depth - 1);
	    // right branch
	    tree(	ofsx + dirx * size, ofsy + diry * size, 
				dirx * cos(-rot) + diry * sin(-rot), 
				dirx * -sin(-rot) + diry * cos(-rot), 
				size * 0.5/ SCALE + size * (SCALE - 1) / SCALE, 
				rot * ROTSCALE, depth - 1);
	}
}
 
void render() { 
	tree(WIDTH / 2.0, HEIGHT - 5.0, 0.0, -1.0, INITLENGTH, M_PI / 5, BRANCHES);
	Sleep(100);
	bmp.Render(g);
//	bmp.Write("tree.bmp");
}

int main() {
	render();
	ConsoleGraphics::wait();
	return 0;
}

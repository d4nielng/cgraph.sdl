// ======================================================================
// MANDELBROT.CPP
// ----------------------------------------------------------------------
// Mandelbrot Set Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 mandelbrot.cpp -o mandelbrot -s -lgdi32
// ======================================================================

#include <iostream>
#include "include/cgraph.h"

using namespace CGraph;
using namespace std;

#define	MIN_ITER		100
#define	MAX_ITER		2000

ConsoleGraphics g;
CBitmap bmp(640, 400);

int iterations = MAX_ITER;

static string gauge_text[] = {"#---------",
					  	 	  "##--------",
				      	 	  "###-------",
					  	 	  "####------",
			    	  	 	  "#####-----",
					  	 	  "######----",
				      	 	  "#######---",
					  	 	  "########--",
					  	 	  "#########-",
					  	 	  "##########"};

void mandelbrot(CBitmap & bmp, double left, double top, double xside, double yside) {
    double xscale, yscale, zx, zy, cx, tempx, cy;
    int x, y, i, j, maxx, maxy, count;
  
    maxx = bmp.Width();
    maxy = bmp.Height();  
    xscale = xside / maxx;
    yscale = yside / maxy;
   
    cout << fixed;
    cout.precision(0);
    for (y = 1; y <= maxy - 1; y++) {
		cout << "Rendering " << (int) y * 100.0 / (maxy-1) << "% " << gauge_text[y * 10 / maxy] << "\r";
        for (x = 1; x <= maxx - 1; x++) {
            cx = x * xscale + left;	// c_real
            cy = y * yscale + top;	// c_imaginary
            zx = 0;	// z_real
            zy = 0;	// z_imaginary
            count = 0;
  			
            while ((zx * zx + zy * zy < 4) && (count < iterations)) {
                tempx = zx * zx - zy * zy + cx;
                zy = 2 * zx * zy + cy;               
				zx = tempx;
                count++;
            }
  
            CHSV hsv(CMath::map(count*32, 0.0, iterations, 0.0, 360.0), 1.0, 1.0);
            if (count == iterations)
            	bmp.SetPixel(x, y, 0, 0, 0);
            else
            	bmp.SetPixel(x, y, hsv.ToRGB());
        }
    }
}

void menu() {
	bool valid;
	do {
		cout << "Number of iterations (" << MIN_ITER << ".." << MAX_ITER << "): ";
		if (!(cin >> iterations)) {
			valid = false;
		}
		else {
			valid = (iterations >= MIN_ITER && iterations <= MAX_ITER);
		}
		if (!valid) 
			printf("Please select within the range\n");
		ConsoleGraphics::flush();
	} while (!valid);
		
}

void render() {
	menu();
	CRect rc(-2.5, -1.2, 1.5, 1.2);
	mandelbrot(bmp, rc.Left(), rc.Top(), rc.Width(), rc.Height());
    Sleep(100);
    bmp.Render(g);
//	bmp.Write("mandelbrot.bmp");
}  

int main() {
	render();	
	ConsoleGraphics::wait();
	return 0;
}

// This program generates a Mandelbrot set and saves it to a BMP file
// compile: g++ -O2 mandelbrot.cpp -o mandelbrot -lSDL2
// run: ./mandelbrot
#include "cgraph.h"

using namespace daniel;
using namespace std;

#define FACTOR          9/16
#define WIDTH           720
#define HEIGHT          WIDTH*FACTOR
#define	MIN_ITER		100
#define	MAX_ITER		2000

int iterations = MAX_ITER;

static string gauge_text[] = {  ">---------",
                                ">>--------",
                                ">>>-------",
                                ">>>>------",
                                ">>>>>-----",
                                ">>>>>>----",
                                ">>>>>>>---",
                                ">>>>>>>>--",
                                ">>>>>>>>>-",
                                ">>>>>>>>>>"};
                                
class Mandelbrot: public CGraph {
private:
    bool redraw;

public:
    Mandelbrot() {
        CGraph::create(WIDTH, HEIGHT, "Mandelbrot Set");
        redraw = true;
    }

    virtual void render() {
        CRect rc(-2.5, -1.2, 1.5, 1.2);
        mandelbrot(rc.left(), rc.top(), rc.width(), rc.height());
        redraw = false;
    }

    void mandelbrot(double left, double top, double xside, double yside) {
        if (!redraw) return;

        double xscale, yscale, zx, zy, cx, tempx, cy;
        int x, y, i, j, maxx, maxy, count;
      
        maxx = width;
        maxy = height;
        xscale = xside / maxx;
        yscale = yside / maxy;
       
        cout << fixed;
        cout.precision(0);
        for (y = 1; y <= maxy - 1; y++) {
            cout << "Rendering " << (int) y * 100.0 / (maxy-1) << "% " << gauge_text[y * 10 / maxy] << "\r";
            cout << flush;
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
      
                CHSV hsv(CMath::map(count*32, 0.0, iterations >> 1, 0.0, 360.0), 1.0, 1.0);
                plotPixel(x, y, ((count == iterations) ? MAKERGB(0, 0, 0) : hsv.toRGB().getColor()));
            }
        }
        CBitmap * bmp = getBitmap();
        bmp->write("mandel.bmp");
        delete bmp;
    }
};

void menu() {
	bool valid;
	do {
		cout << "Number of iterations (" << MIN_ITER << ".." << MAX_ITER << "): ";
		cin >> iterations;
		valid = (iterations >= MIN_ITER && iterations <= MAX_ITER);
		if (!valid) 
			printf("Please select within the range\n");
	} while (!valid);
		
}

int main(int argc, char* argv[]) {
    menu();
    Mandelbrot app;
    app.loop();
    return 0;
}
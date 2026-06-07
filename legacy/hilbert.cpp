// ======================================================================
// HILBERT.CPP
// ----------------------------------------------------------------------
// Hilbert Curve Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 hilbert.cpp -o hilbert -s -lgdi32
// ----------------------------------------------------------------------
// Reference: https://youtu.be/dSK-MW-zuAc
// ======================================================================

#include "include/cgraph.h"

using namespace CGraph;
using namespace std;

#define	WIDTH		640
#define	HEIGHT		480
#define	DEFAULT		7		// default Hilbert Curve order

ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);
int order = DEFAULT;

CVector2D points[] = {
	CVector2D(0, 0), 
	CVector2D(0, 1),
	CVector2D(1, 1),
	CVector2D(1, 0)
};

CVector2D hilbert(int i, int order) {
	int index = i & 3;
	CVector2D v = points[index];
	for (int j = 1; j < order; j++) {
		i >>= 2;
		index = i & 3;
		double len = pow(2, j);
		if (index == 0)
			v.shuffle();
		else if (index == 1)
			v += CVector2D(0.0f, len);
		else if (index == 2)
			v += CVector2D(len, len);
		else if (index == 3) {
			CVector2D temp(len-1.0f-v.getY(), len-1.0f-v.getX());
			v = temp;
			v += CVector2D(len, 0.0f);
		}
	}
	return v;
}

void curve(int order = DEFAULT) {
	int N = pow(2, order);
	int total = N * N;
	CVector2D * path = new CVector2D[total];
	float len = (WIDTH / N) * 0.6f;
	int cw = (N + 1) * len;
	int sx = (WIDTH - cw) / 2;
	int sy = (HEIGHT - cw) / 2;
		
	for (int i = 0; i < total; i++) {
		path[i] = hilbert(i, order);
		path[i].scale(len);
		path[i] += CVector2D(len, len);
	}
	
	for (int i = 1; i < total; i++) {
		
		CHSV hsv(CMath::map(i, 0, total, 0, 360), 1.0f, 1.0f);
		bmp.SetColor(hsv.ToRGB());
		
		bmp.DrawLine(sx + path[i].getX(), sy + path[i].getY(),
				 	 sx + path[i-1].getX(), sy + path[i-1].getY());
	}
	delete[] path;	
}

void menu() {
	bool valid;
	do {
		cout <<	"Select "<< VT_COLOR(15) << "Hilbert Curve " << VT_DEFAULTATTR << "order:" << endl <<
				"[" << VT_COLOR(196) << "2"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR(202) << "3"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR(220) << "4"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR(190) << "5"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR( 48) << "6"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR( 33) << "7"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR( 21) << "8"<< VT_DEFAULTATTR << "] " << 
				"[" << VT_COLOR(165) << "9"<< VT_DEFAULTATTR << "] " << endl;
		cout << ">";
			if (!(cin >> order)) {
				valid = false;
			}
			else {
				valid = (order >= 2 && order <= 9);
			}
		ConsoleGraphics::flush();
		if (!valid) 
			printf("Please select within the range\n");
	} while (!valid);
}

void render() { 
	menu();
	curve(order);
	Sleep(100);
	bmp.Render(g);
//	bmp.Write("hilbert.bmp");
}

int main() {
	render();
	ConsoleGraphics::wait();
	return 0;
}

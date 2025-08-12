// This program demonstrates the Hilbert Space Filling Curve
// compile: g++ -O2 -o hilbert hilbert.cpp -lSDL2
// run: ./hilbert
#include "cgraph.h"

using namespace daniel;
using namespace std;

#define	WIDTH		720
#define	HEIGHT		WIDTH/2
#define	DEFAULT		7		// default Hilbert Curve order

int order = DEFAULT;

CVector2D points[] = {
	CVector2D(0, 0), 
	CVector2D(0, 1),
	CVector2D(1, 1),
	CVector2D(1, 0)
};

class Hilbert: public CGraph {
private:
	int order;
public:
	Hilbert() {
		CGraph::create(WIDTH, HEIGHT, "Hilber Space Filling Curve");
	}

	void setOrder(int ord) { order = ord; }
	int getOrder() const { return order; }
	virtual void render() { curve(order); }

	CVector2D hilbert(int i, int order) {
		long index = i & 3;
		CVector2D v = points[index];
		for (long j = 1; j < order; j++) {
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
		long N = pow(2, order);
		long total = N * N;
		CVector2D * path = new CVector2D[total];
		float len = (WIDTH / N) * .5f;
		long cw = (N + 1) * len;
		long sx = (WIDTH - cw) / 2;
		long sy = (HEIGHT - cw) / 2;
			
		for (long i = 0; i < total; i++) {
			path[i] = hilbert(i, order);
			path[i].scale(len);
			path[i] += CVector2D(len, len);
		}
		
		moveTo(sx + path[0].getX(), sy + path[0].getY());
		for (long i = 1; i < total; i++) {
			CHSV hsv(CMath::map(i, 0, total, 0, 360), 1.0f, 1.0f);
			setColor(hsv.toRGB());
			lineTo(sx + path[i].getX(), sy + path[i].getY());
		}
		delete[] path;	
	}

};

int menu() {
	bool valid;
	int order;
	do {
		cout <<	"Select "<< "Hilbert Curve " << "order:" << endl <<
				"[" << "2"<< "] " << 
				"[" << "3"<< "] " << 
				"[" << "4"<< "] " << 
				"[" << "5"<< "] " << 
				"[" << "6"<< "] " << 
				"[" << "7"<< "] " << 
				"[" << "8"<< "] " << 
				"[" << "9"<< "] " << endl;
		cout << ">";
		cin >> order;
		valid = (order >= 2 && order <= 9);
		if (!valid) 
			printf("Please select within the range\n");
	} while (!valid);
	return order;
}

int main() {
	int choice = menu();
	Hilbert app;
	app.setOrder(choice);
	app.loop();
	return 0;
}

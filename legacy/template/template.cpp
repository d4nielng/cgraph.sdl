#include "cgraph.h"

using namespace CGraph;

#define	WIDTH		640
#define	HEIGHT		400

// ===============================
// your graphics drawing code here
// ===============================

void paint(CBitmap & bmp) {

}

// =============================== 

ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);
 
void render() { 
	paint(bmp);
	bmp.Render(g);
//	bmp.Write("template.bmp");
}

int main() {
	render();
	ConsoleGraphics::wait();
	return 0;
}

// This program draws the Cantor set using recursion, a simple fractal.
// compile: g++ -O2 -o cantor cantor.cpp -lSDL2
// run: ./cantor
#include "cgraph.h"

using namespace daniel;

#define WIDTH       640
#define HEIGHT      480
#define	DEPTH		7
#define	THICKNESS	(HEIGHT/(DEPTH*2))
#define	SIZE		(WIDTH - 10)
#define	GAP			(THICKNESS * 2)
#define	CX			(WIDTH - SIZE) / 2
#define	CY			(HEIGHT - DEPTH * GAP) / 2
#define HUE_STEP	(360/(DEPTH+1))

class Cantor: public CGraph {
public:
    Cantor() {
        CGraph::create(WIDTH, HEIGHT, "Cantor Set");
    }

    virtual void render() {
        cantor(CX, CY, SIZE, 0);
    }
    
    void cantor(int x, int y, int length, int depth) {
        if (depth == DEPTH)
            return;
        
        double v = depth;
        
        // map color to depth
        CHSV hsv(CMath::map(v, 0.0, (double) DEPTH, 0.0, 360.0), 0.6f, 1.0f);
        setColor(hsv.toRGB());
        
        // new length is one third of previous length
        int newlength = length / 3;
        
        // draw the bar
        rectangle(x, y, length, THICKNESS);
        
        // progress recursively down
        cantor(x              , y + GAP, newlength, depth + 1);
        cantor(x + 2*newlength, y + GAP, newlength, depth + 1);
    }
};

int main(int argc, char* argv[]) {
    Cantor app;
    app.loop();
    return 0;
}
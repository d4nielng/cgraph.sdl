// This program draws a Sierpinski triangle using recursive calls.
// compile: g++ -O2 sierpinski.cpp -o sierpinski -lSDL2
// run: ./sierpinski
#include "cgraph.h"

using namespace daniel;

#define WIDTH       720
#define HEIGHT      WIDTH/2
#define MIN_HEIGHT  8
#define SQ3         sqrt(3)

class Sierpinski: public CGraph {
public:
    Sierpinski() {
        CGraph::create(WIDTH, HEIGHT, "Sierpinski Triangle");
    }

    virtual void render() {
        sierpinski(width / 2,  2 * height / 3 - 4,  height / 2);
    }

    void triangle(float x, float y, float h) {
        line(x - h / SQ3, y - h / 3, x + h / SQ3, y - h / 3);
        line(x - h / SQ3, y - h / 3, x, y + 2 * h / 3);
        line(x, y + 2 * h / 3, x + h / SQ3, y - h / 3);
    }
      
    void sierpinski(float x, float y, float h) { 
        if (h < MIN_HEIGHT) return;

        CHSV hsv(CMath::map(y, 0, height, 0, 360), 1.0, 1.0);
        setColor(hsv.toRGB());

        triangle(x, y, h);

        sierpinski(x          , y - 2 * h / 3, h / 2);
        sierpinski(x - h / SQ3, y + h / 3    , h / 2);
        sierpinski(x + h / SQ3, y + h / 3    , h / 2);
    }
};

int main(int argc, char* argv[]) {
    Sierpinski app;
    app.loop();
    return 0;
}
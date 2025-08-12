// This example demonstrates a fractal tree using the CGraph class.
// compile: g++ -O2 tree.cpp -o tree -lSDL2
// run: ./tree
#include <random>
#include "cgraph.h"

using namespace daniel;

#define WIDTH 800
#define HEIGHT 600
#define DEPTH 10

class FractalTree : public CGraph {
private:
    double left_factor;
    double right_factor;
public:
    FractalTree() {
        CGraph::create(WIDTH, HEIGHT, "Fractal Tree");

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
    
        left_factor = CMath::clampf(dis(gen), 0.0, 0.9);
        right_factor = CMath::clampf(dis(gen), 0.0, 0.9);
    }

    virtual void render() {
        draw_tree(WIDTH / 2, HEIGHT - 50, M_PI / 2, DEPTH);
    }

    void draw_tree(int x1, int y1, double angle, int depth) {
        if (depth == 0) return;
    
        CHSV hsv(CMath::map(depth * 3, 0, 48, 0, 360), 0.8, 1.0);
        setColor(hsv.toRGB());
    
        int x2 = x1 + (int)(cos(angle) * depth * 10.0);
        int y2 = y1 - (int)(sin(angle) * depth * 10.0);
    
        line(x1, y1, x2, y2);
    
        draw_tree(x2, y2, angle - left_factor, depth - 1);
        draw_tree(x2, y2, angle + right_factor, depth - 1);
    }
    
};

int main(int argc, char *argv[]) {
    FractalTree app;
    app.loop();
    return 0;
}

#include "cgraph.h"

using namespace daniel;

#define WIDTH 640
#define HEIGHT 480

class Demo: public CGraph {
public:
    Demo() {
        CGraph::create(WIDTH, HEIGHT, "CGraph Demo");
    }

    virtual void render() {
        clear();
        static int count = 16;
        int ch = height / count;
        for (int i = 0; i < count; i++) {
            CHSV hsv = CHSV(CMath::map(i, 0, count, 0, 360), 1, 1);
            setColor(hsv.toRGB());
            rectangle(0, i * ch, width, ch);
        }
    }
};

int main(int argc, char* argv[]) {
    Demo app;
    app.loop();
    return 0;
}
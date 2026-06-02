// This program reads a bitmap file and displays it on the screen, using the cgraph 2D engine.
// compile: g++ -O2 viewbmp.cpp -o viewbmp -lSDL2
// usage: viewbmp [filename]
//
// Controls:
//   Arrow keys  - Pan image
//   + / -       - Zoom in / out
//   R           - Reset view (center, 1:1 zoom)
//   ESC         - Quit
//
// TODO:
// - Add a single-column file picker with previews
//

#include "cgraph.h"

using namespace std;
using namespace daniel;

#ifdef _WIN32
    #define CMD_CLEAR   "cls"
    #define CMD_LS      "dir images /w"
#else
    #define CMD_CLEAR   "clear"
    #define CMD_LS      "ls images/"
#endif

// BitmapObject - wraps a CBitmap as a 2D engine game object.
// Position and zoom are managed through CTransform (position + scale).
class BitmapObject : public CGameObject {
private:
    CBitmap& bmp;
public:
    BitmapObject(CBitmap& b, int x = 0, int y = 0)
        : CGameObject(x, y), bmp(b) {}

    virtual void update(double deltaTime) {
        transform.update(deltaTime);
    }

    // Render the bitmap at the current transform position and scale.
    // Uses direct pixel buffer writes with nearest-neighbour scaling —
    // the same technique as CBitmap::render(), guaranteed to match the
    // SDL window surface format without any SDL_Blit format conversion.
    virtual void render(CGraph& gfx) {
        if (!bmp.width() || !bmp.height()) return;

        int dstX  = (int)transform.position.getX();
        int dstY  = (int)transform.position.getY();
        double s  = transform.scale;
        int srcW  = (int)bmp.width();
        int srcH  = (int)bmp.height();
        int dstW  = (int)(srcW * s);
        int dstH  = (int)(srcH * s);

        uint32_t* dst  = gfx.getPixels();
        int winW = (int)gfx.getWidth();
        int winH = (int)gfx.getHeight();
        uint32_t* src  = bmp.getPixels();   // {b,g,r,a} packed as uint32_t

        for (int dy = 0; dy < dstH; dy++) {
            int py = dstY + dy;
            if (py < 0 || py >= winH) continue;
            int sy = (int)(dy / s);
            if (sy >= srcH) continue;
            for (int dx = 0; dx < dstW; dx++) {
                int px = dstX + dx;
                if (px < 0 || px >= winW) continue;
                int sx = (int)(dx / s);
                if (sx >= srcW) continue;
                dst[py * winW + px] = src[sy * srcW + sx];
            }
        }
    }

    int imgWidth()  const { return bmp.width(); }
    int imgHeight() const { return bmp.height(); }
};

class ViewBMP : public CGraph {
private:
    CBitmap      img;
    CScene       scene;
    std::string  filename;
    bool         loaded;

    static const int    PAN_STEP  = 20;
    static const double ZOOM_STEP;  // factor per keypress
    static const double ZOOM_MIN;
    static const double ZOOM_MAX;

    void resetView() {
        if (!loaded || scene.getObjectCount() == 0) return;
        CTransform& t = scene.getObjects()[0]->getTransform();
        t.setScale(1.0);
        // Center on window
        t.setPosition((int)(getWidth() - img.width())  / 2,
                      (int)(getHeight() - img.height()) / 2);
    }

public:
    ViewBMP(const std::string& fname) : filename(fname), loaded(false) {
        // Load the bitmap before creating the window so we know the size
        loaded = img.load(fname);

        uint32_t winW = loaded ? max((uint32_t)img.width(),  (uint32_t)640) : 640;
        uint32_t winH = loaded ? max((uint32_t)img.height(), (uint32_t)480) : 480;
        CGraph::create(winW, winH, "Bitmap Viewer - " + fname);
        setClipping(true);

        if (loaded) {
            BitmapObject* obj = new BitmapObject(img);
            scene.addObject(obj);
            resetView();
        }
    }

    virtual void onKeyDown(SDL_Keysym keysym) {
        if (!loaded || scene.getObjectCount() == 0) return;
        CTransform& t = scene.getObjects()[0]->getTransform();
        double s = t.scale;

        switch (keysym.sym) {
        case SDLK_a:     t.translate(-PAN_STEP, 0);             break;
        case SDLK_d:    t.translate( PAN_STEP, 0);             break;
        case SDLK_w:       t.translate(0, -PAN_STEP);             break;
        case SDLK_s:     t.translate(0,  PAN_STEP);             break;
        case SDLK_EQUALS:
        case SDLK_KP_PLUS:  t.setScale(min(s / (1.0 - ZOOM_STEP), ZOOM_MAX)); break;
        case SDLK_MINUS:
        case SDLK_KP_MINUS: t.setScale(max(s * (1.0 - ZOOM_STEP), ZOOM_MIN)); break;
        case SDLK_r:        resetView();                            break;
        }
    }

    virtual void render() {
        // Dark grey background (use setColor so SDL_MapRGB maps the value correctly)
        setColor(0x20, 0x20, 0x20);
        clear();

        if (!loaded) {
            setColor(0xFF5555);
            drawText("Cannot load: " + filename, 10, 10);
            update();
            return;
        }

        // Let the engine render the bitmap game object (handles pan + zoom)
        scene.render(*this);

        // Info bar at bottom
        if (scene.getObjectCount() > 0) {
            CTransform& t = scene.getObjects()[0]->getTransform();
            char info[256];
            snprintf(info, sizeof(info),
                     "%s  %dx%d  ZOOM:%.0f%%  [WASD]PAN [+/-]ZOOM [R] RESET",
                     filename.c_str(), img.width(), img.height(), t.scale * 100.0);
            // Black backing strip
            setColor(0x000000);
            rectangle(0, getHeight() - 10, getWidth(), 10);
            setColor(0xFFFF00);
            drawText(info, 2, getHeight() - 9);
        }

        update();
    }
};

const double ViewBMP::ZOOM_STEP = 0.25;
const double ViewBMP::ZOOM_MIN  = 0.10;
const double ViewBMP::ZOOM_MAX  = 8.00;

// ────────────────────────────────────────────────────────────────────────────

string gauge_text[] = {string(">") +"---------",
                       string(">>")+"--------",
                       string(">>>")+  "-------",
                       string(">>>>")+  "------",
                       string(">>>>>")+  "-----",
                       string(">>>>>>")+  "----",
                       string(">>>>>>>")+  "---",
                       string(">>>>>>>>")+  "--",
                       string(">>>>>>>>>")+  "-",
                       string(">>>>>>>>>>")};

void callback(int current, int total) {
    cout << "- loading [" << gauge_text[current * 10 / total] << "] "
         << (current + 1) * 100 / total << "%" << "\r" << flush;
}

void refresh() {
    system(CMD_CLEAR);
    system(CMD_LS);
}

int main(int argc, char** argv) {
    string filename;

    if (argc > 1) {
        filename = argv[1];
    } else {
        refresh();
        cout << "Enter filename: ";
        cin >> filename;
        cin.clear();
        cin.ignore(INT_MAX, '\n');
    }

    if (filename.find(".bmp") == string::npos)
        filename += ".bmp";

    // Prepend images/ prefix only when no path separator is present
    string path = (filename.find('/') == string::npos && filename.find('\\') == string::npos)
                  ? "images/" + filename
                  : filename;

    // Show a progress gauge in the terminal while loading
    if (argc <= 1) {
        CBitmap probe;
        probe.load(path, callback);
        cout << endl;
        probe.info();
    }

    ViewBMP viewer(path);
    viewer.loop();
    return 0;
}


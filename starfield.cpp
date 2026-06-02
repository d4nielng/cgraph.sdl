// This program renders a 3D starfield with perspective projection.
// Stars move toward the viewer and their projected radius grows from 1 to 3 pixels.
// compile: g++ -O2 -o starfield starfield.cpp -lSDL2
// run: ./starfield
#include "cgraph.h"

#include <cstdlib>
#include <ctime>
#include <vector>

using namespace daniel;

#define WIDTH       800
#define HEIGHT      600
#define STAR_COUNT  256

class Starfield: public CGraph {
    struct Star {
        CVector3D position;
        double speed;
        double baseRadius;
    };

    std::vector<Star> stars;

    double centerX;
    double centerY;
    double focalLength;
    double nearPlane;
    double farPlane;
    double spreadX;
    double spreadY;

    bool blurEnabled;
    SDL_Texture * trailTexture;
    bool hwBlurAvailable;
    std::string hudMessage;
    double hudMessageUntil;

public:
    Starfield() {
        CGraphOptions options(WIDTH, HEIGHT, "Starfield (CVector3D Demo)");
        options.allowFullScreen = true;
        options.resizable = true;
        options.minWidth = 480;
        options.minHeight = 320;
        options.targetFPS = 120;
        options.backgroundColor = MAKERGB(2, 5, 14);
        CGraph::create(options);

        centerX = width * 0.5;
        centerY = height * 0.5;

        focalLength = 300.0;
        nearPlane = 1.0;
        farPlane = 82.0;

        spreadX = 56.0;
        spreadY = spreadX * (double)height / (double)width;

        blurEnabled = false;
        trailTexture = NULL;
        hwBlurAvailable = false;
        hudMessage = "SPACE blur  F12 screenshot  Cmd/Alt+Enter fullscreen";
        hudMessageUntil = 0.0;

        rebuildTrailTexture();

        std::srand((unsigned)std::time(NULL));
        stars.resize(STAR_COUNT);
        for (size_t i = 0; i < stars.size(); i++) {
            respawn(stars[i], false);
        }

        updateTitle();
    }

    ~Starfield() {
        destroyTexture(trailTexture);
    }

    virtual void onKeyDown(SDL_Keysym keysym) {
        if (keysym.sym == SDLK_SPACE) {
            blurEnabled = !blurEnabled;
            if (hwBlurAvailable && !blurEnabled) {
                setRenderTarget(trailTexture);
                clearTargetRGBA(2, 5, 14, 255);
                setRenderTarget(NULL);
            }
            updateTitle();
        } else if (keysym.sym == SDLK_F12) {
            if (saveScreenshot("starfield-screenshot.bmp")) {
                hudMessage = "Saved starfield-screenshot.bmp";
            } else {
                hudMessage = std::string("Screenshot failed: ") + getLastError();
            }
            hudMessageUntil = getElapsedTime() + 2.5;
        }
    }

    virtual void onResize(int w, int h) {
        updateTitle();
    }

    virtual void render() {
        double dt = getDeltaTime();

        if (dt > 0.05) dt = 0.05;

        // Slight blue-black background helps depth perception.
        if (hwBlurAvailable) {
            setRenderTarget(trailTexture);
            setBlendMode(SDL_BLENDMODE_BLEND);
            if (!blurEnabled) {
                clearTargetRGBA(2, 5, 14, 255);
            } else {
                // Hardware-accelerated persistence fade.
                fillRectRGBA(0, 0, (int)width, (int)height, 2, 5, 14, 32);
            }
        } else {
            clear(2, 5, 14);
        }

        for (size_t i = 0; i < stars.size(); i++) {
            Star & s = stars[i];

            // Move star toward the camera (viewer at z = 0).
            s.position.incZ(-s.speed * dt);

            // Recycle star once it reaches the near plane.
            if (s.position.getZ() <= nearPlane) {
                respawn(s, true);
            }

            CVector2D p = s.position.projectTo2D(focalLength, centerX, centerY, 0.0);
            int px = (int)p.getX();
            int py = (int)p.getY();

            if (px < -3 || px >= (int)width + 3 || py < -3 || py >= (int)height + 3) {
                respawn(s, true);
                continue;
            }

            double scale = focalLength / (s.position.getZ());
            int radius = (int)(s.baseRadius * scale + 0.5);
            radius = (int)CMath::clampl(radius, 1, 3);

            // Far stars are darker blue; near stars shift to medium blue.
            double depthFactor = 1.0 - CMath::clampf((s.position.getZ() - nearPlane) / (farPlane - nearPlane), 0.0, 1.0);
            int r = (int)CMath::map(depthFactor, 0.0, 1.0, 10.0, 70.0);
            int g = (int)CMath::map(depthFactor, 0.0, 1.0, 24.0, 122.0);
            int b = (int)CMath::map(depthFactor, 0.0, 1.0, 78.0, 205.0);

            // Extra brightness kick for near stars.
            double nearBoost = depthFactor * depthFactor;
            r = (int)CMath::clampl(r + (long)(nearBoost * 95.0), 0, 255);
            g = (int)CMath::clampl(g + (long)(nearBoost * 120.0), 0, 255);
            b = (int)CMath::clampl(b + (long)(nearBoost * 70.0), 0, 255);

            if (hwBlurAvailable) {
                setColor(r, g, b);
                filledCircle(px, py, radius);
            } else {
                setColor(r, g, b);
                filledCircle(px, py, radius);
            }

            // Tiny bloom for the closest stars only.
            if (radius >= 3) {
                int br = (int)CMath::clampl((long)(r + 45), 0, 255);
                int bg = (int)CMath::clampl((long)(g + 55), 0, 255);
                int bb = (int)CMath::clampl((long)(b + 35), 0, 255);
                if (hwBlurAvailable) {
                    setColor(br, bg, bb);
                    drawCircle(px, py, radius + 1);
                } else {
                    setColor(br, bg, bb);
                    drawCircle(px, py, radius + 1);
                }
            }
        }

        if (hwBlurAvailable) {
            drawRendererStatus(12, 12);
            setColor(180, 200, 230);
            drawTextRight("SPACE Blur  F12 Shot", (int)width - 12, (int)height - 16);
            if (getElapsedTime() < hudMessageUntil) {
                setColor(220, 235, 255);
                drawTextCentered(hudMessage, (int)(width * 0.5), 12);
            }
            setRenderTarget(NULL);
            clearBackbuffer(0, 0, 0, 255);
            copyTextureToBackbuffer(trailTexture);
        } else {
            drawRendererStatus(12, 12);
            setColor(180, 200, 230);
            drawTextRight("SPACE Blur  F12 Shot", (int)width - 12, (int)height - 16);
            if (getElapsedTime() < hudMessageUntil) {
                setColor(220, 235, 255);
                drawTextCentered(hudMessage, (int)(width * 0.5), 12);
            }
        }
    }

private:
    void rebuildTrailTexture() {
        destroyTexture(trailTexture);
        trailTexture = NULL;
        hwBlurAvailable = false;

        if (!hasRenderer())
            return;

        trailTexture = createRenderTargetTexture(width, height, true);
        if (!trailTexture)
            return;

        setRenderTarget(trailTexture);
        setBlendMode(SDL_BLENDMODE_BLEND);
        clearTargetRGBA(2, 5, 14, 255);
        setRenderTarget(NULL);
        hwBlurAvailable = true;
    }

    void updateTitle() {
        setTitle(std::string("Starfield (CVector3D Demo) - Blur ") + (blurEnabled ? "ON" : "OFF") +
                 (hwBlurAvailable ? " [HW]" : " [SW]") +
                 " - " + getHardwareAcceleratorName());
    }

    double randf(double minv, double maxv) {
        return minv + (maxv - minv) * ((double)std::rand() / (double)RAND_MAX);
    }

    void respawn(Star & s, bool farDepthOnly) {
        s.position.setX(randf(-spreadX, spreadX));
        s.position.setY(randf(-spreadY, spreadY));

        if (farDepthOnly) {
            s.position.setZ(randf(farPlane * 0.70, farPlane));
        } else {
            s.position.setZ(randf(nearPlane + 0.1, farPlane));
        }

        s.speed = randf(12.0, 28.0);
        s.baseRadius = randf(0.035, 0.070);
    }

};

int main(int argc, char* argv[]) {
    Starfield app;
    app.loop();
    return 0;
}

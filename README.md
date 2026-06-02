# CGRAPH.SDL
Simple work-in-progress single-header hardware accelerated 2D graphics library/engine to use in my classes.

Recent engine usability upgrades:
- `CGraphOptions` for self-documenting window/engine configuration
- `getLastError()` for engine and SDL setup failures
- Built-in frame timing via `getDeltaTime()`, `getElapsedTime()`, `setTargetFPS()`
- Frame-based input helpers: `wasKeyPressed()`, `wasKeyReleased()`, mouse position/button/wheel helpers
- Resize support with `onResize()` and window helpers such as `toggleFullscreen()`, `setResizable()`, `setMinSize()`, `centerWindow()`, `showCursor()`
- Camera helpers: `setCamera()`, `moveCamera()`, `setZoom()`, `worldToScreen()` and world-space drawing wrappers
- Bitmap helpers: `drawBitmap()`, `drawBitmapScaled()`, `drawBitmapRegion()`, `drawBitmapRotated()`
- Text layout helpers: `measureText()`, `drawTextCentered()`, `drawTextRight()`
- Screenshot helper: `saveScreenshot()`
- Lifecycle hooks: `onInit()`, `onUpdate()`, `onShutdown()`

Demos included:
- Cantor set (`cantor.cpp`)
- Mandelbrot set (`mandel.cpp`)
- Sierpinski triangle (`sierpinski.cpp`)
- Hilber curve (`hilbert.cpp`)
- Tree generator (`tree.cpp`)
- 3D starfield (`starfield.cpp`)
- N-queens solver (`nqueens.cpp`)
- Starfield (`starfield.cpp`)
- Simple rainbow demo (`demo.cpp`)

Application:
- Bitmap viewer (`viewbmp.cpp`)

Games:
- Simple Snake game (`snake.cpp`)
- Simple Tetris game (`tetris.cpp`)

# Prerequisite
- SDL2
- GCC/MinGW/Clang

To install SLD2 on macOS, you can use Homebrew (`brew.sh`), type `brew install sdl2`.

On Ubuntu, type `sudo apt install libsdl2-dev`.

I'm too lazy to test the library on Windows. You can do it yourself.

# Building
To build them all, type `make` and hit ENTER.

The makefile does not track header dependencies, so after changing `cgraph.h` use `make clean && make`.

# Example
Minimal setup using the new options object:

```cpp
daniel::CGraphOptions options(800, 600, "My App");
options.allowFullScreen = true;
options.resizable = true;
options.targetFPS = 60;

daniel::CGraph gfx(options);
gfx.run();
```

# WIP
Planned games
- Sudoku
- Tic-Tac-Toe

/*
    Simple Snake Game - Demonstrates cgraph.h capabilities
    
    Controls:
    - Arrow keys or WASD to move
    - R to restart after game over
    - + or ] to increase speed
    - - or [ to decrease speed
    - ESC to quit
*/

#include "cgraph.h"
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace daniel;

class SnakeGame : public CGraph {
private:
    // Game constants
    static const int GRID_WIDTH = 32;
    static const int GRID_HEIGHT = 24;
    static const int CELL_SIZE = 16;
    
    // Game state
    int snakeSpeed;  // frames between moves (lower = faster)
    struct Point {
        int x, y;
        Point() : x(0), y(0) {}
        Point(int x, int y) : x(x), y(y) {}
    };
    
    std::vector<Point> snake;
    Point food;
    int dirX, dirY;
    int nextDirX, nextDirY;
    int score;
    bool gameOver;
    int frameCount;
    
public:
    SnakeGame() : dirX(1), dirY(0), nextDirX(1), nextDirY(0), score(0), gameOver(false), frameCount(0), snakeSpeed(48) {
        create(GRID_WIDTH * CELL_SIZE, GRID_HEIGHT * CELL_SIZE, "Snake Game - Arrow Keys or WASD (+/- to adjust speed)");
        srand((unsigned)time(0));
        resetGame();
    }
    
    void spawnFood() {
        bool validSpot = false;
        while (!validSpot) {
            food.x = rand() % GRID_WIDTH;
            food.y = rand() % GRID_HEIGHT;
            validSpot = true;
            for (const auto& seg : snake) {
                if (seg.x == food.x && seg.y == food.y) {
                    validSpot = false;
                    break;
                }
            }
        }
    }
    
    void setSpeed(int speed) {
        snakeSpeed = (speed > 0 && speed <= 20) ? speed : snakeSpeed;
    }
    
    int getSpeed() const {
        return snakeSpeed;
    }
    
    virtual void onKeyDown(SDL_Keysym keysym) override {
        switch (keysym.sym) {
            case SDLK_UP:
            case SDLK_w:
                if (dirY == 0) { nextDirX = 0; nextDirY = -1; }
                break;
            case SDLK_DOWN:
            case SDLK_s:
                if (dirY == 0) { nextDirX = 0; nextDirY = 1; }
                break;
            case SDLK_LEFT:
            case SDLK_a:
                if (dirX == 0) { nextDirX = -1; nextDirY = 0; }
                break;
            case SDLK_RIGHT:
            case SDLK_d:
                if (dirX == 0) { nextDirX = 1; nextDirY = 0; }
                break;
            case SDLK_r:
                if (gameOver) resetGame();
                break;
            case SDLK_EQUALS:
            case SDLK_RIGHTBRACKET:
                snakeSpeed = (snakeSpeed > 1) ? snakeSpeed - 1 : 1;
                break;
            case SDLK_MINUS:
            case SDLK_LEFTBRACKET:
                snakeSpeed = (snakeSpeed < 20) ? snakeSpeed + 1 : 20;
                break;
            default:
                break;
        }
    }
    
    void resetGame() {
        snake.clear();
        snake.push_back(Point(GRID_WIDTH / 2, GRID_HEIGHT / 2));
        snake.push_back(Point(GRID_WIDTH / 2 - 1, GRID_HEIGHT / 2));
        snake.push_back(Point(GRID_WIDTH / 2 - 2, GRID_HEIGHT / 2));
        dirX = 1;
        dirY = 0;
        nextDirX = 1;
        nextDirY = 0;
        score = 0;
        gameOver = false;
        frameCount = 0;
        spawnFood();
    }
    
    void updateGame() {
        if (gameOver) return;
        
        frameCount++;
        if (frameCount < snakeSpeed) return;
        frameCount = 0;
        
        // Update direction
        dirX = nextDirX;
        dirY = nextDirY;
        
        // Calculate new head position with wrapping
        Point newHead = snake[0];
        newHead.x = (newHead.x + dirX + GRID_WIDTH) % GRID_WIDTH;
        newHead.y = (newHead.y + dirY + GRID_HEIGHT) % GRID_HEIGHT;
        
        // Check self collision
        for (const auto& seg : snake) {
            if (seg.x == newHead.x && seg.y == newHead.y) {
                gameOver = true;
                return;
            }
        }
        
        snake.insert(snake.begin(), newHead);
        
        // Check food collision
        if (newHead.x == food.x && newHead.y == food.y) {
            score += 10;
            spawnFood();
        } else {
            snake.pop_back();
        }
    }
    
    virtual void render() override {
        // Update game logic
        updateGame();
        
        // Clear screen with black background
        clear(0, 0, 0);
        
        // Draw grid (light gray)
        setColor(40, 40, 40);
        for (int i = 0; i <= GRID_WIDTH; i++) {
            line(i * CELL_SIZE, 0, i * CELL_SIZE, (int)height);
        }
        for (int i = 0; i <= GRID_HEIGHT; i++) {
            line(0, i * CELL_SIZE, (int)width, i * CELL_SIZE);
        }
        
        // Draw food (red square - fills entire cell)
        setColor(255, 100, 100);
        rectangle(food.x * CELL_SIZE, food.y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
        
        // Draw snake with gradient color
        for (size_t i = 0; i < snake.size(); i++) {
            if (i == 0) {
                // Head in bright green
                setColor(100, 255, 100);
            } else {
                // Body with HSV color gradient
                double hue = CMath::map((double)i, 0.0, (double)snake.size(), 120.0, 60.0);
                CHSV col(hue, 0.8, 0.7);
                CRGB rgb = col.toRGB();
                setColor(rgb.red(), rgb.green(), rgb.blue());
            }
            // Draw snake segment - fills entire cell
            rectangle(snake[i].x * CELL_SIZE, snake[i].y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
        }
        
        // Draw UI text at bottom
        setColor(200, 200, 200);
        
        // Draw border frame to contain score info
        int infoY = GRID_HEIGHT * CELL_SIZE + 5;
        
        if (gameOver) {
            setColor(255, 100, 100);
            line(10, infoY, (int)width - 10, infoY);
            line(10, infoY + 20, (int)width - 10, infoY + 20);
            setColor(255, 255, 255);
            // Game Over indicator via visual feedback
        }
    }
};

int main() {
    SnakeGame game;
    game.loop();
    return 0;
}

/*
    Snake Game - Refactored to use the cgraph 2D engine
    
    Controls:
    - Arrow keys or WASD to move
    - R to restart
    - +/- to adjust speed
    - ESC to quit
*/

#include "cgraph.h"
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace daniel;

// Food game object with score value
class Food : public CGameObject {
private:
    int scoreValue;  // Points this food is worth (1-5)
    
public:
    Food(double x, double y, int score, double hue) : CGameObject(x, y), scoreValue(score) {
        // Color based on hue
        CHSV hsv(hue, 0.9, 0.8);
        CRGB rgb = hsv.toRGB();
        color = MAKERGB(rgb.red(), rgb.green(), rgb.blue());
        setCollider(new CCollider(12, 12));
    }
    
    int getScore() const { return scoreValue; }
    
    void render(CGraph& gfx) override {
        gfx.setColor((color >> 16) & 0xFF, (color >> 8) & 0xFF, color & 0xFF);
        gfx.rectangle((int)transform.position.getX() - 6, 
                      (int)transform.position.getY() - 6, 12, 12);
    }
};

// Snake segment
class SnakeSegment : public CGameObject {
private:
    int segmentIndex;
    
public:
    SnakeSegment(double x, double y, int index) 
        : CGameObject(x, y, MAKERGB(100, 255, 100)), segmentIndex(index) {
        setCollider(new CCollider(12, 12));
        setLayer(index);  // Higher indices render on top
    }
    
    void render(CGraph& gfx) override {
        if (segmentIndex == 0) {
            // Head - bright green
            gfx.setColor(100, 255, 100);
        } else {
            // Body with gradient
            double hue = CMath::map((double)segmentIndex, 0.0, 20.0, 120.0, 60.0);
            CHSV col(hue, 0.8, 0.7);
            CRGB rgb = col.toRGB();
            gfx.setColor(rgb.red(), rgb.green(), rgb.blue());
        }
        gfx.rectangle((int)transform.position.getX() - 6, 
                      (int)transform.position.getY() - 6, 12, 12);
    }
    
    void setSegmentIndex(int idx) { segmentIndex = idx; }
};

class SnakeGame : public CGraph {
private:
    static const int GRID_WIDTH = 32;
    static const int GRID_HEIGHT = 24;
    static const int CELL_SIZE = 16;
    
    CScene scene;
    CTimer timer;
    
    std::vector<SnakeSegment*> snakeSegments;
    std::vector<Food*> fruits;  // Multiple fruits active at once
    static const int MAX_FRUITS = 5;
    
    int dirX, dirY;
    int nextDirX, nextDirY;
    int score;
    int lives;
    bool gameOver;
    double moveTimer;
    double moveInterval;  // Time between moves (in seconds)
    
public:
    SnakeGame() : dirX(1), dirY(0), nextDirX(1), nextDirY(0), score(0), 
                  lives(3), gameOver(false), moveTimer(0), moveInterval(0.1) {
        create(GRID_WIDTH * CELL_SIZE + 150, GRID_HEIGHT * CELL_SIZE, 
               "Snake Game - Arrow Keys or WASD (+/- to adjust speed)");
        srand((unsigned)time(nullptr));
        resetGame();
    }
    
    ~SnakeGame() {
        scene.clear();  // Scene handles all object deletion
        fruits.clear();
    }
    
    void resetGame() {
        scene.clear();  // Scene handles all object deletion
        snakeSegments.clear();
        fruits.clear();
        
        // Create snake body (3 segments) - head first at index 0, centered in grid cells
        for (int i = 0; i <= 2; i++) {
            SnakeSegment* seg = new SnakeSegment(
                (GRID_WIDTH / 2 - i) * CELL_SIZE + CELL_SIZE / 2,
                (GRID_HEIGHT / 2) * CELL_SIZE + CELL_SIZE / 2,
                i
            );
            snakeSegments.push_back(seg);
            scene.addObject(seg);
        }
        
        // Create initial fruits
        spawnMultipleFruits();
        
        dirX = 1;
        dirY = 0;
        nextDirX = 1;
        nextDirY = 0;
        score = 0;
        lives = 3;
        gameOver = false;
        moveTimer = 0;
    }
    
    void resetGameAfterDeath() {
        scene.clear();  // Scene handles all object deletion
        snakeSegments.clear();
        fruits.clear();
        
        // Create snake body (3 segments) - head first at index 0, centered in grid cells
        for (int i = 0; i <= 2; i++) {
            SnakeSegment* seg = new SnakeSegment(
                (GRID_WIDTH / 2 - i) * CELL_SIZE + CELL_SIZE / 2,
                (GRID_HEIGHT / 2) * CELL_SIZE + CELL_SIZE / 2,
                i
            );
            snakeSegments.push_back(seg);
            scene.addObject(seg);
        }
        
        // Create initial fruits
        spawnMultipleFruits();
        
        dirX = 1;
        dirY = 0;
        nextDirX = 1;
        nextDirY = 0;
        // Keep score, decrement lives
        lives--;
        if (lives <= 0) {
            gameOver = true;
        }
        moveTimer = 0;
    }
    
    void drawNumber(int num, int x, int y, int cellSize) {
        // Draw single digit using block pattern (3x5 grid)
        static const bool digits[10][5][3] = {
            // 0
            { {1,1,1}, {1,0,1}, {1,0,1}, {1,0,1}, {1,1,1} },
            // 1
            { {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1} },
            // 2
            { {1,1,1}, {0,0,1}, {1,1,1}, {1,0,0}, {1,1,1} },
            // 3
            { {1,1,1}, {0,0,1}, {1,1,1}, {0,0,1}, {1,1,1} },
            // 4
            { {1,0,1}, {1,0,1}, {1,1,1}, {0,0,1}, {0,0,1} },
            // 5
            { {1,1,1}, {1,0,0}, {1,1,1}, {0,0,1}, {1,1,1} },
            // 6
            { {1,1,1}, {1,0,0}, {1,1,1}, {1,0,1}, {1,1,1} },
            // 7
            { {1,1,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1} },
            // 8
            { {1,1,1}, {1,0,1}, {1,1,1}, {1,0,1}, {1,1,1} },
            // 9
            { {1,1,1}, {1,0,1}, {1,1,1}, {0,0,1}, {1,1,1} }
        };
        
        if (num >= 0 && num <= 9) {
            for (int row = 0; row < 5; row++) {
                for (int col = 0; col < 3; col++) {
                    if (digits[num][row][col]) {
                        rectangle(x + col * cellSize, y + row * cellSize, cellSize, cellSize);
                    }
                }
            }
        }
    }
    
    void spawnMultipleFruits() {
        // Spawn up to MAX_FRUITS fruits
        while (fruits.size() < MAX_FRUITS) {
            spawnSingleFruit();
        }
    }
    
    void spawnSingleFruit() {
        if (fruits.size() >= MAX_FRUITS) return;
        
        bool validSpot = false;
        double foodX, foodY;
        int scoreValue = 0;
        
        while (!validSpot) {
            int gx = rand() % GRID_WIDTH;
            int gy = rand() % GRID_HEIGHT;
            foodX = gx * CELL_SIZE + CELL_SIZE / 2;  // Center in grid cell
            foodY = gy * CELL_SIZE + CELL_SIZE / 2;  // Center in grid cell
            scoreValue = (rand() % 5) + 1;  // Score 1-5
            
            validSpot = true;
            
            // Check snake collision
            for (auto seg : snakeSegments) {
                if (seg->getTransform().position.getX() == foodX &&
                    seg->getTransform().position.getY() == foodY) {
                    validSpot = false;
                    break;
                }
            }
            
            // Check existing fruit collision
            if (validSpot) {
                for (auto fruit : fruits) {
                    if (fruit->getTransform().position.getX() == foodX &&
                        fruit->getTransform().position.getY() == foodY) {
                        validSpot = false;
                        break;
                    }
                }
            }
        }
        
        // Create fruit with hue based on score (blue to purple)
        double hue = 240.0 + (scoreValue - 1) * 15.0;  // 240 (blue) to 300 (purple)
        Food* newFruit = new Food(foodX, foodY, scoreValue, hue);
        fruits.push_back(newFruit);
        scene.addObject(newFruit);
    }
    
    void moveSnake() {
        if (gameOver || snakeSegments.empty()) return;
        
        dirX = nextDirX;
        dirY = nextDirY;
        
        // Calculate new head position with wrapping
        double headX = snakeSegments[0]->getTransform().position.getX() + dirX * CELL_SIZE;
        double headY = snakeSegments[0]->getTransform().position.getY() + dirY * CELL_SIZE;
        
        // Wrap around edges
        while (headX < 0) headX += GRID_WIDTH * CELL_SIZE;
        while (headX >= GRID_WIDTH * CELL_SIZE) headX -= GRID_WIDTH * CELL_SIZE;
        while (headY < 0) headY += GRID_HEIGHT * CELL_SIZE;
        while (headY >= GRID_HEIGHT * CELL_SIZE) headY -= GRID_HEIGHT * CELL_SIZE;
        
        // Check self collision
        for (auto seg : snakeSegments) {
            if (seg->getTransform().position.getX() == headX &&
                seg->getTransform().position.getY() == headY) {
                resetGameAfterDeath();
                return;
            }
        }
        
        // Create new head
        SnakeSegment* newHead = new SnakeSegment(headX, headY, 0);
        snakeSegments.insert(snakeSegments.begin(), newHead);
        scene.addObject(newHead);
        
        // Update segment indices
        for (size_t i = 0; i < snakeSegments.size(); i++) {
            snakeSegments[i]->setSegmentIndex(i);
            snakeSegments[i]->setLayer((int)i);
        }
        
        // Check collision with any fruit
        bool fruitEaten = false;
        for (size_t i = 0; i < fruits.size(); i++) {
            if (newHead->checkCollision(*fruits[i])) {
                score += fruits[i]->getScore();  // Add fruit's score value
                scene.removeObject(fruits[i]);  // removeObject already deletes
                fruits.erase(fruits.begin() + i);
                fruitEaten = true;
                break;
            }
        }
        
        if (fruitEaten) {
            // Spawn a new fruit to maintain MAX_FRUITS
            spawnSingleFruit();
        } else {
            // Remove tail if no food eaten
            SnakeSegment* tail = snakeSegments.back();
            scene.removeObject(tail);
            snakeSegments.pop_back();
        }
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
                resetGame();
                break;
            case SDLK_EQUALS:
            case SDLK_RIGHTBRACKET:
                moveInterval = CMath::maxf(0.02, moveInterval - 0.02);
                break;
            case SDLK_MINUS:
            case SDLK_LEFTBRACKET:
                moveInterval = CMath::minf(0.3, moveInterval + 0.02);
                break;
            default:
                break;
        }
    }
    
    virtual void render() override {
        timer.update();
        
        // Update snake movement timer
        moveTimer += timer.getDeltaTime();
        if (moveTimer >= moveInterval) {
            moveTimer -= moveInterval;
            moveSnake();
        }
        
        // Update scene
        scene.update(timer.getDeltaTime());
        
        // Clear screen
        clear(0, 0, 0);
        
        // Draw grid (only within game field)
        setColor(40, 40, 40);
        int gameFieldHeight = GRID_HEIGHT * CELL_SIZE;
        int gameFieldWidth = GRID_WIDTH * CELL_SIZE;
        for (int i = 0; i <= GRID_WIDTH; i++) {
            line(i * CELL_SIZE, 0, i * CELL_SIZE, gameFieldHeight);
        }
        for (int i = 0; i <= GRID_HEIGHT; i++) {
            line(0, i * CELL_SIZE, gameFieldWidth, i * CELL_SIZE);
        }
        
        // Render all scene objects
        scene.render(*this);
        
        // Draw score panel on the right side
        int panelX = GRID_WIDTH * CELL_SIZE + 5;
        int panelY = 5;
        int panelWidth = 140;
        int panelHeight = GRID_HEIGHT * CELL_SIZE - 10;
        
        // Draw border only (no gray fill)
        setColor(100, 100, 100);
        line(panelX, panelY, panelX + panelWidth, panelY);
        line(panelX, panelY + panelHeight, panelX + panelWidth, panelY + panelHeight);
        line(panelX, panelY, panelX, panelY + panelHeight);
        line(panelX + panelWidth, panelY, panelX + panelWidth, panelY + panelHeight);
        
        // Draw score label and value
        
        setColor(200, 200, 200);
        
        // Display score as 5-digit number with leading zeros if needed
        int scoreDisplayX = panelX + 15;
        int scoreDisplayY = panelY + 20;
        int cellSize = 5;
        
        // Draw each digit of the score
        int scoreValue = score;
        int digits[5];
        digits[0] = (scoreValue / 10000) % 10;
        digits[1] = (scoreValue / 1000) % 10;
        digits[2] = (scoreValue / 100) % 10;
        digits[3] = (scoreValue / 10) % 10;
        digits[4] = scoreValue % 10;
        
        for (int i = 0; i < 5; i++) {
            drawNumber(digits[i], scoreDisplayX + i * 18, scoreDisplayY, cellSize);
        }
        
        // Draw lives label and count
        setColor(200, 200, 200);
        int livesDisplayX = panelX + 35;
        int livesDisplayY = panelY + 60;
        
        // Display lives as single digit
        drawNumber(lives, livesDisplayX, livesDisplayY, cellSize);
        
        if (gameOver) {
            setColor(255, 100, 100);
        }
    }
};

int main() {
    SnakeGame game;
    game.loop();
    return 0;
}

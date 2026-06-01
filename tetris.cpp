/*
    Simple Tetris Game - Demonstrates cgraph.h capabilities
    
    Controls:
    - Arrow keys or AD to move left/right
    - W or UP to rotate
    - S or DOWN to move down faster
    - SPACE to hard drop
    - R to restart after game over
    - ESC to quit
*/

#include "cgraph.h"
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace daniel;

class TetrisGame : public CGraph {
private:
    // Game constants
    static const int BOARD_WIDTH = 10;
    static const int BOARD_HEIGHT = 20;
    static const int CELL_SIZE = 24;
    static const int PREVIEW_SIZE = 4;
    
    // Tetris piece definitions (4x4 grid)
    struct Piece {
        int grid[4][4];
        int color;
    };
    
    // Standard Tetris pieces
    static const Piece PIECES[7];
    
    // Game state
    int board[BOARD_HEIGHT][BOARD_WIDTH];
    int currentPieceType;
    int nextPieceType;
    int pieceX, pieceY;
    int pieceRotation;
    int score;
    int lines;
    bool gameOver;
    int frameCount;
    int dropSpeed;
    int lockDelayFrames;  // Delay before piece is locked in place
    bool pieceLanded;     // Flag to prevent movement after piece lands
    const int LOCK_DELAY = 1;  // Frames to wait before locking, 1 is the best value
    
public:
    TetrisGame() : currentPieceType(0), nextPieceType(0), pieceX(3), pieceY(0), 
                   pieceRotation(0), score(0), lines(0), gameOver(false), 
                   frameCount(0), dropSpeed(90), lockDelayFrames(0), pieceLanded(false) {
        srand((unsigned)time(nullptr));  // Seed random once at startup
        create(CELL_SIZE * BOARD_WIDTH + 300, CELL_SIZE * BOARD_HEIGHT + 40, "Tetris - Arrow Keys to Move, W to Rotate, SPACE to Drop");
        resetGame();
    }
    
    void resetGame() {
        memset(board, 0, sizeof(board));
        score = 0;
        lines = 0;
        gameOver = false;
        frameCount = 0;
        lockDelayFrames = 0;
        pieceLanded = false;
        
        nextPieceType = rand() % 7;
        spawnNewPiece();
    }
    
    void rotatePieceGrid(int grid[4][4], int rotated[4][4]) {
        // Rotate 90 degrees clockwise
        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++) {
                rotated[x][3-y] = grid[y][x];
            }
        }
    }
    
    void spawnNewPiece() {
        currentPieceType = nextPieceType;
        nextPieceType = rand() % 7;
        pieceX = 3;
        pieceY = -3;  // Spawn above the field, piece drops into play area
        pieceRotation = 0;
        lockDelayFrames = 0;  // Reset lock delay for new piece
        pieceLanded = false;  // New piece is not landed
        
        if (checkCollision(pieceX, pieceY)) {
            gameOver = true;
        }
    }
    
    bool checkCollision(int x, int y, int rotation = -1) {
        if (rotation == -1) rotation = pieceRotation;
        
        const Piece& piece = PIECES[currentPieceType];
        
        // Create rotated grid
        int rotated[4][4];
        memcpy(rotated, piece.grid, sizeof(piece.grid));
        
        // Apply rotation
        for (int r = 0; r < rotation; r++) {
            int temp[4][4];
            rotatePieceGrid(rotated, temp);
            memcpy(rotated, temp, sizeof(temp));
        }
        
        // Check collision with rotated grid
        for (int py = 0; py < 4; py++) {
            for (int px = 0; px < 4; px++) {
                if (rotated[py][px] == 0) continue;
                
                int bx = x + px;
                int by = y + py;
                
                if (bx < 0 || bx >= BOARD_WIDTH || by >= BOARD_HEIGHT) {
                    if (by >= 0) return true;
                }
                if (by >= 0 && board[by][bx] != 0) {
                    return true;
                }
            }
        }
        return false;
    }
    
    void lockPiece() {
        const Piece& piece = PIECES[currentPieceType];
        
        // Create rotated grid
        int rotated[4][4];
        memcpy(rotated, piece.grid, sizeof(piece.grid));
        for (int r = 0; r < pieceRotation; r++) {
            int temp[4][4];
            rotatePieceGrid(rotated, temp);
            memcpy(rotated, temp, sizeof(temp));
        }
        
        for (int py = 0; py < 4; py++) {
            for (int px = 0; px < 4; px++) {
                if (rotated[py][px] == 0) continue;
                
                int bx = pieceX + px;
                int by = pieceY + py;
                
                if (by >= 0 && by < BOARD_HEIGHT && bx >= 0 && bx < BOARD_WIDTH) {
                    board[by][bx] = piece.color;
                }
            }
        }
        
        clearLines();
        spawnNewPiece();
    }
    
    void clearLines() {
        for (int y = BOARD_HEIGHT - 1; y >= 0; y--) {
            bool full = true;
            for (int x = 0; x < BOARD_WIDTH; x++) {
                if (board[y][x] == 0) {
                    full = false;
                    break;
                }
            }
            
            if (full) {
                // Move lines down
                for (int yy = y; yy > 0; yy--) {
                    for (int x = 0; x < BOARD_WIDTH; x++) {
                        board[yy][x] = board[yy-1][x];
                    }
                }
                for (int x = 0; x < BOARD_WIDTH; x++) {
                    board[0][x] = 0;
                }
                
                score += 100;
                lines++;
                y++; // Check this line again
            }
        }
    }
    
    virtual void onKeyDown(SDL_Keysym keysym) override {
        // R key restarts game anytime
        if (keysym.sym == SDLK_r) {
            resetGame();
            return;
        }
        
        // Handle hard drop before checking if piece is landed
        if (keysym.sym == SDLK_SPACE) {
            while (!checkCollision(pieceX, pieceY + 1)) {
                pieceY++;
                score += 2;
            }
            pieceLanded = true;  // Force piece to lock immediately
            lockDelayFrames = LOCK_DELAY;
            return;
        }
        
        // Don't allow other movement if piece has landed
        if (pieceLanded) return;
        
        switch (keysym.sym) {
            case SDLK_LEFT:
            case SDLK_a:
                if (!checkCollision(pieceX - 1, pieceY)) {
                    pieceX--;
                }
                break;
            case SDLK_RIGHT:
            case SDLK_d:
                if (!checkCollision(pieceX + 1, pieceY)) {
                    pieceX++;
                }
                break;
            case SDLK_UP:
            case SDLK_w:
                if (!checkCollision(pieceX, pieceY, (pieceRotation + 1) % 4)) {
                    pieceRotation = (pieceRotation + 1) % 4;
                }
                break;
            case SDLK_DOWN:
            case SDLK_s:
                // Soft drop - move down one cell
                if (!checkCollision(pieceX, pieceY + 1)) {
                    pieceY++;
                    score += 1;
                }
                break;
            default:
                break;
        }
    }
    
    void updateGame() {
        if (gameOver) return;
        
        frameCount++;
        if (frameCount < dropSpeed) return;
        frameCount = 0;
        
        if (checkCollision(pieceX, pieceY + 1)) {
            // Piece has hit something, start lock delay
            pieceLanded = true;  // Mark piece as landed
            lockDelayFrames++;
            if (lockDelayFrames >= LOCK_DELAY) {
                lockPiece();
            }
        } else {
            pieceY++;
            lockDelayFrames = 0;  // Reset delay if piece can still move
        }
    }
    
    void drawNumber(int num, int x, int y, int cellSize) {
        // Draw single digit using block pattern
        // Simple 3x5 grid representation
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
    
    virtual void render() override {
        updateGame();
        clear(0, 0, 0);
        
        // Center the board horizontally
        int windowWidth = CELL_SIZE * BOARD_WIDTH + 300;
        int boardWidth = BOARD_WIDTH * CELL_SIZE;
        int boardX = (windowWidth - boardWidth) / 2;
        int boardY = 20;
        
        // Draw preview box (left side) - enlarged for full piece display
        setColor(100, 100, 100);
        rectangle(5, 5, 115, 115);
        
        // Draw preview piece (centered in preview box)
        setColor(100, 200, 255);
        const Piece& nextPiece = PIECES[nextPieceType];
        int pieceSize = 4 * CELL_SIZE;  // 96 pixels
        int previewBoxSize = 115;
        int previewPieceOffsetX = 5 + (previewBoxSize - pieceSize) / 2;
        int previewPieceOffsetY = 5 + (previewBoxSize - pieceSize) / 2;
        for (int py = 0; py < 4; py++) {
            for (int px = 0; px < 4; px++) {
                if (nextPiece.grid[py][px] != 0) {
                    rectangle(previewPieceOffsetX + px * CELL_SIZE, previewPieceOffsetY + py * CELL_SIZE, CELL_SIZE, CELL_SIZE);
                }
            }
        }
        
        // Draw board background
        setColor(30, 30, 30);
        rectangle(boardX, boardY, BOARD_WIDTH * CELL_SIZE, BOARD_HEIGHT * CELL_SIZE);
        
        // Draw grid
        setColor(50, 50, 50);
        for (int x = 0; x <= BOARD_WIDTH; x++) {
            line(boardX + x * CELL_SIZE, boardY, boardX + x * CELL_SIZE, boardY + BOARD_HEIGHT * CELL_SIZE);
        }
        for (int y = 0; y <= BOARD_HEIGHT; y++) {
            line(boardX, boardY + y * CELL_SIZE, boardX + BOARD_WIDTH * CELL_SIZE, boardY + y * CELL_SIZE);
        }
        
        // Draw placed blocks
        for (int y = 0; y < BOARD_HEIGHT; y++) {
            for (int x = 0; x < BOARD_WIDTH; x++) {
                if (board[y][x] != 0) {
                    setColor((board[y][x] >> 16) & 0xFF, (board[y][x] >> 8) & 0xFF, board[y][x] & 0xFF);
                    rectangle(boardX + x * CELL_SIZE, boardY + y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
                }
            }
        }
        
        // Draw current piece
        const Piece& piece = PIECES[currentPieceType];
        
        // Create rotated grid
        int rotated[4][4];
        memcpy(rotated, piece.grid, sizeof(piece.grid));
        for (int r = 0; r < pieceRotation; r++) {
            int temp[4][4];
            rotatePieceGrid(rotated, temp);
            memcpy(rotated, temp, sizeof(temp));
        }
        
        setColor((piece.color >> 16) & 0xFF, (piece.color >> 8) & 0xFF, piece.color & 0xFF);
        for (int py = 0; py < 4; py++) {
            for (int px = 0; px < 4; px++) {
                if (rotated[py][px] != 0) {
                    int bx = boardX + (pieceX + px) * CELL_SIZE;
                    int by = boardY + (pieceY + py) * CELL_SIZE;
                    if (pieceY + py >= 0) {
                        rectangle(bx, by, CELL_SIZE, CELL_SIZE);
                    }
                }
            }
        }
        
        // Draw score on right side
        setColor(255, 200, 100);
        int scoreX = boardX + BOARD_WIDTH * CELL_SIZE + 30;
        int scoreY = 30;
        
        // Draw "Score" label using blocks
        drawNumber(score / 1000, scoreX, scoreY, 4);
        drawNumber((score / 100) % 10, scoreX + 20, scoreY, 4);
        drawNumber((score / 10) % 10, scoreX + 40, scoreY, 4);
        drawNumber(score % 10, scoreX + 60, scoreY, 4);
        
        // Draw game over message
        if (gameOver) {
            setColor(255, 100, 100);
            for (int i = 0; i < 5; i++) {
                line(boardX - 5, boardY - 5 - i, boardX + BOARD_WIDTH * CELL_SIZE + 5, boardY - 5 - i);
            }
        }
    }
};

// Initialize pieces
const TetrisGame::Piece TetrisGame::PIECES[7] = {
    // I - Cyan
    { { {0,0,0,0}, {1,1,1,1}, {0,0,0,0}, {0,0,0,0} }, 0x00FFFF },
    // O - Yellow
    { { {0,0,0,0}, {0,1,1,0}, {0,1,1,0}, {0,0,0,0} }, 0xFFFF00 },
    // T - Purple
    { { {0,0,0,0}, {0,1,1,1}, {0,0,1,0}, {0,0,0,0} }, 0xFF00FF },
    // S - Green
    { { {0,0,0,0}, {0,0,1,1}, {0,1,1,0}, {0,0,0,0} }, 0x00FF00 },
    // Z - Red
    { { {0,0,0,0}, {0,1,1,0}, {0,0,1,1}, {0,0,0,0} }, 0xFF0000 },
    // J - Blue
    { { {0,0,0,0}, {0,1,0,0}, {0,1,1,1}, {0,0,0,0} }, 0x0000FF },
    // L - Orange
    { { {0,0,0,0}, {0,0,0,1}, {0,1,1,1}, {0,0,0,0} }, 0xFF8800 }
};

int main() {
    TetrisGame game;
    game.loop();
    return 0;
}

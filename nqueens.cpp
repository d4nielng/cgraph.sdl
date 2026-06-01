// This program solves the N-Queens problem using backtracking with animation.
// compile: g++ -o nqueens nqueens.cpp -lSDL2
// run: ./nqueens
#include "cgraph.h"
#include <cstdlib>
#include <ctime>

using namespace daniel;

#define	WIDTH	480
#define	HEIGHT	492

#define	N	8
#define	BOARD_SIZE	480
#define	CW	BOARD_SIZE / N
#define	CH	BOARD_SIZE / N
#define BX	0
#define	BY	0

struct point {
    double x;
    double y;
} queen[8] =   {{1.0-0.20, 1.0-0.20}, 
                {1.0-0.10, 1.0-0.60},
                {1.0-0.30, 1.0-0.50},
                {1.0-0.45, 1.0-0.90},
                {1.0-0.60, 1.0-0.50},
                {1.0-0.80, 1.0-0.60}, 
                {1.0-0.70, 1.0-0.20},
                {1.0-0.20, 1.0-0.20}};

// Attempt tracker for visualization
struct Attempt {
    int row, col;
    bool succeeded;
};

// Queen game object for rendering on board
class Queen: public CGameObject {
private:
    int boardX, boardY;
public:
    Queen(int x, int y) : boardX(x), boardY(y) {
        transform.setPosition(BX + boardX * CW + CW * 0.5f, 
                             BY + boardY * CH + CW * 0.5f);
    }
    
    virtual void update(double deltaTime) {}
    
    virtual void render(CGraph& g) {
        g.setColor(0x00FF00);  // Green for placed queens
        g.moveTo(BX + boardX * CW + CW * queen[0].x, 
                 BY + boardY * CH + CW * queen[0].y);
        for (int i = 1; i < 8; i++)
            g.lineTo(BX + boardX * CW + CW * queen[i].x, 
                    BY + boardY * CH + CW * queen[i].y);
    }
};

class Queens: public CGraph {
private:
    bool board[N][N];
    bool solved;
    CScene scene;
    
    // Animation state
    std::vector<Attempt> attempts;
    int currentAttempt;
    double stepDelay;
    double stepTimer;
    int curCol;
    bool isRunning;
    bool isComplete;
    int stepCount;
    int firstRow;   // row where the first queen is placed
    
public:
    Queens() : currentAttempt(0), stepDelay(0.002), stepTimer(0), 
               curCol(0), isRunning(false), isComplete(false), stepCount(0), firstRow(0) {
        CGraph::create(WIDTH, HEIGHT, "N-Queens Solver - Solving...");
        memset(board, 0, sizeof(board));
        
        // Load the VGA 8x8 font
        loadFont("vga8x8.bmp");
        
        // Pick a random starting row for the first queen
        std::srand((unsigned)std::time(nullptr));
        firstRow = std::rand() % N;
        
        // Place first queen at the randomly chosen row in column 0
        board[firstRow][0] = true;
        scene.addObject(new Queen(0, firstRow));
        Attempt first;
        first.row = firstRow;
        first.col = 0;
        first.succeeded = true;
        attempts.push_back(first);
        
        // Solve from column 1 onward
        char titleBuf[128];
        if (animatedExplore(board, 1)) {
            snprintf(titleBuf, sizeof(titleBuf),
                     "N-Queens Solver - Solution Found! (first queen at row %d)", firstRow + 1);
        } else {
            snprintf(titleBuf, sizeof(titleBuf),
                     "N-Queens Solver - No solution from row %d", firstRow + 1);
        }
        isComplete = true;
        setTitle(titleBuf);
    }

    virtual void update(double deltaTime) {
        // Animation already complete during construction
    }

    virtual void render() {
        drawBoard();
        
        // Render attempted positions
        for (size_t i = 0; i < attempts.size(); i++) {
            int row = attempts[i].row;
            int col = attempts[i].col;
            if (attempts[i].succeeded) {
                setColor(0xFF6600);  // Orange for positions where queen was placed
            } else {
                setColor(0xFF0000);  // Red for failed attempts
            }
            int markerX = BX + col * CW + CW / 2;
            int markerY = BY + row * CH + CH / 2;
            filledCircle(markerX, markerY, 4);
        }
        
        // Render final queens in green
        scene.render(*this);
        
        // Fill stats area with black
        setColor(0x000000);
        rectangle(0, BOARD_SIZE, WIDTH, 12);
        
        // Draw statistics text in stats area
        setColor(0xFFFF00);
        char buf[128];
        snprintf(buf, sizeof(buf), "START ROW: %d | ATTEMPTS: %zu | STATES: %d", 
                 firstRow + 1, attempts.size(), stepCount);
        drawText(buf, 10, BOARD_SIZE + 2);
    }

    void drawCell(int x, int y, bool white) {
        if (white)
            setColor(0xFFFFFF);
        else
            setColor(0x202020);
        rectangle(BX + x * CW, BY + y * CH, CW, CH);
    }
    
    void drawBoard() {
        bool white = true;
        for (int i = 0; i < N; i++)	{
            white = !white;
            for (int j = 0; j < N; j++) {
                white = !white;
                drawCell(i, j, white);
            }
        }
    }  
    
    bool isSafe(bool board[N][N], int row, int col) {  
        // this row on left side
        for (int i = 0; i < col; i++)
            if (board[row][i]) return false; 
        // upper diagonal on left side
        for (int i = row, j = col; i >= 0 && j >= 0; i--, j--)  
            if (board[i][j]) return false;  
        // lower diagonal on left side
        for (int i = row, j = col; j >= 0 && i < N; i++, j--)  
            if (board[i][j]) return false;  
        return true;  
    }
    
    // Simple recursive solver that records attempts for animation
    bool animatedExplore(bool board[N][N], int col) {
        // Base case: all queens placed
        if (col >= N) {
            return true;
        }
        
        // Try each row in this column
        for (int row = 0; row < N; row++) {
            // Record attempt
            Attempt att;
            att.row = row;
            att.col = col;
            att.succeeded = false;
            
            if (isSafe(board, row, col)) {
                // Place queen
                board[row][col] = true;
                scene.addObject(new Queen(col, row));
                att.succeeded = true;
                attempts.push_back(att);
                curCol = col + 1;
                
                // Recursively solve for next column
                if (animatedExplore(board, col + 1)) {
                    return true;
                }
                
                // Backtrack
                board[row][col] = false;
                auto& objects = scene.getObjects();
                if (!objects.empty()) {
                    scene.removeObject(objects.back());
                }
            } else {
                // Failed attempt
                attempts.push_back(att);
            }
            
            stepCount++;  // Count each row attempt
        }
        
        curCol = col;
        return false;
    }
};


int main() {
    Queens app;
    app.loop();
    return 0;
}

// This program solves the N-Queens problem using backtracking.
// compile: g++ -o nqueens nqueens.cpp -lSDL2
// run: ./nqueens
#include "cgraph.h"

using namespace daniel;

#define	WIDTH	480
#define	HEIGHT	480

#define	N	8
#define	CW	HEIGHT / N
#define	CH	HEIGHT / N
#define BX	(WIDTH - N * CW)/2
#define	BY	(HEIGHT - N * CW)/2

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

class Queens: public CGraph {
public:
    Queens() {
        CGraph::create(WIDTH, HEIGHT, "N-Queens Solver");
    }

    virtual void render() {
        placeQueens();	
    }

    void drawCell(int x, int y, bool white) {
        if (white)
            setColor(0xFFFFFF);
        else
            setColor(0x202020);
        rectangle(BX + x * CW, BY + y * CH, CW, CH);
    }
    
    void drawQueen(int x, int y) {
        setColor(0x7F7F7F);
        moveTo(BX + x * CW + CW * queen[0].x, 
               BY + y * CH + CW * queen[0].y);
        for (int i = 1; i < 8; i++)
            lineTo(BX + x * CW + CW * queen[i].x, BY + y * CH + CW * queen[i].y);
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
    
    void showSolution(bool board[N][N]) {
        drawBoard();
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (board[i][j])  
                    drawQueen(i, j);
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
    
    bool explore(bool board[N][N], int col) { 
        // base case: if all queens are placed then return true 
        if (col >= N) return true;  
        // Consider this column and try placing this queen in all rows one by one
        for (int i = 0; i < N; i++) {  
            // Check if queen can be placed on board[i][col]
            if (isSafe(board, i, col)) {  
                // Place this queen in board[i][col]
                board[i][col] = true;  
                // recurse to place rest of the queens
                if (explore(board, col + 1)) return true;  
                // If failed to place queen in board[i][col], then remove her from board[i][col]
                board[i][col] = false; // BACKTRACK  
            }  
        }  
        return false;  // Cannot be placed in any row in this colum col then return false
    }  
    
    void placeQueens() {  
        bool board[N][N];
        for(int i = 0; i < N; i++)
               for(int j = 0; j < N; j++)
                   board[i][j] = false;        
        if (explore(board, 0) == false)
            printf("Solution does not exist\n");
        else
            showSolution(board);  
    }  
};

int main() {
    Queens app;
    app.loop();
	return 0;
}

// ======================================================================
// NQUEENS.CPP
// ----------------------------------------------------------------------
// N-Queens Placement Visualization
// Coded by Trinh D.D. Nguyen
// ----------------------------------------------------------------------
// Demonstration for CGRAPH library
// ----------------------------------------------------------------------
// Compilers supported: Dev-C++, TDM-GCC, LLVM MinGW64
// ----------------------------------------------------------------------
// Compile: 
//     g++ -std=c++11 -O2 nqueens.cpp -o nqueens -s -lgdi32
// ======================================================================
#include "include/cgraph.h"

using namespace CGraph;

#define	WIDTH	640
#define	HEIGHT	400

#define	N	8
#define	CW	HEIGHT / N
#define	CH	HEIGHT / N
#define BX	(WIDTH - N * CW)/2
#define	BY	(HEIGHT - N * CW)/2

typedef struct point {
	double x;
	double y;
} CPOINT;

CPOINT queen[8] = {{1.0-0.20, 1.0-0.20}, 
	        	   {1.0-0.10, 1.0-0.60},
			  	   {1.0-0.30, 1.0-0.50},
			       {1.0-0.45, 1.0-0.90},
			       {1.0-0.60, 1.0-0.50},
			       {1.0-0.80, 1.0-0.60}, 
			       {1.0-0.70, 1.0-0.20},
			       {1.0-0.20, 1.0-0.20}};

ConsoleGraphics g;
CBitmap bmp(WIDTH, HEIGHT);

void drawcell(int x, int y, bool white) {
	RECT rc = {BX + x*CW, BY + y*CH, BX + x*CW+CW, BY + y*CH+CH};
	if (white)
		bmp.SetColor(255, 255, 255);
	else
		bmp.SetColor(32, 32, 32);
	bmp.Rectangle(rc.left, rc.top, rc.right, rc.bottom);
}

void drawqueen(int x, int y) {
	bmp.SetColor(127, 127, 127);
	bmp.MoveTo(BX + x*CW + CW*queen[0].x, BY + y*CH+CW*queen[0].y);
	for (int i = 1; i < 8; i++)
		bmp.LineTo(BX + x*CW + CW*queen[i].x, BY + y*CH+CW*queen[i].y);
}

void drawboard() {
	bool white = true;
	for (int i = 0; i < N; i++)	{
		white = !white;
		for (int j = 0; j < N; j++) {
			white = !white;
			drawcell(i, j, white);
		}
	}
} 

void show_solution(bool board[N][N]) {
	drawboard();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (board[i][j])  
            	drawqueen(i, j);
}  

bool check_safe(bool board[N][N], int row, int col) {  
    int i, j;  
    // this row on left side
    for (i = 0; i < col; i++)
        if (board[row][i]) return false; 
    // upper diagonal on left side
    for (i = row, j = col; i >= 0 && j >= 0; i--, j--)  
        if (board[i][j]) return false;  
    // lower diagonal on left side
    for (i = row, j = col; j >= 0 && i < N; i++, j--)  
        if (board[i][j]) return false;  
    return true;  
}  

bool explore(bool board[N][N], int col) { 
	// base case: if all queens are placed then return true 
    if (col >= N) return true;  
    // Consider this column and try placing this queen in all rows one by one
    for (int i = 0; i < N; i++) {  
        // Check if queen can be placed on board[i][col]
        if (check_safe(board, i, col)) {  
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

void place_queens() {  
    bool board[N][N];
    
    for(int i = 0; i < N; i++)
   		for(int j = 0; j < N; j++)
   			board[i][j] = false;
    
	if (explore(board, 0) == false)
        printf("Solution does not exist\n");
	else
    	show_solution(board);  
}  

void render() {
	place_queens();	
    Sleep(100);
	bmp.Render(g);
//	bmp.Write("nqueens.bmp");	
}

int main() {
	render();
	ConsoleGraphics::wait();
	return 0;
}

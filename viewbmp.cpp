// This program reads a bitmap file and displays it on the screen, using SDL2 library.
// compile: g++ -O2 viewbmp.cpp -o viewbmp -lSDL2
// usage: viewbmp [filename]
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

class ViewBMP: public CGraph {
private:
	CBitmap img;
public:
	ViewBMP(const CBitmap & bmp): img(bmp) { 
		CGraph::create(bmp.width(), bmp.height(), "Bitmap Viewer"); 
	}
    virtual void render() { img.render(surface, true); }
};

string	gauge_text[] = {string(">")+"---------", 
						string(">>")+"--------", 
						string(">>>")+"-------",
						string(">>>>")+"------",
						string(">>>>>")+"-----",
						string(">>>>>>")+"----", 
						string(">>>>>>>")+"---", 
						string(">>>>>>>>")+"--", 
						string(">>>>>>>>>")+"-", 
						string(">>>>>>>>>>")};

void callback(int current, int total) {
	cout << "- loading [" << gauge_text[current*10/total] << "] " << (current+1)*100/total << "%" << "\r";	
	cout << flush;
}

void refresh() {
	system(CMD_CLEAR);
	system(CMD_LS);
}

int main(int argc, char ** argv) {
	string filename;
	CBitmap bmp;
	bool result;
	
	if (argc > 1)
		filename = argv[1];
	else {
		refresh();
		cout << ">";			// otherwise, tell the users to input the desired file name
		cin >> filename;
		cin.clear();			// clear input buffers
		cin.ignore(INT_MAX,'\n'); 
	}

	if (filename.find(".bmp") == string::npos) 
		filename += ".bmp";
	
	result = bmp.load((argc > 1) ? filename : "images/"+filename, callback);

	// attempt to read the bitmap file
	if (result) {
		cout << endl;
		bmp.info();

		ViewBMP viewbmp(bmp);
		viewbmp.loop();
	}
	else
		cout << "Cannot load image file from \"" << filename << "\"" << endl;

	return 0;
}

/*
    Single header for Simple Console Graphics
     _______  _______  ______    _______  _______  __   __        __   __ 
    |       ||       ||    _ |  |   _   ||       ||  | |  |      |  | |  |
    |       ||    ___||   | ||  |  |_|  ||    _  ||  |_|  |      |  |_|  |
    |       ||   | __ |   |_||_ |       ||   |_| ||       |      |       |
    |      _||   ||  ||    __  ||       ||    ___||       | ___  |       |
    |     |_ |   |_| ||   |  | ||   _   ||   |    |   _   ||   | |   _   |
    |_______||_______||___|  |_||__| |__||___|    |__| |__||___| |__| |__|
                                                Coded by Trinh D.D. Nguyen
                                                       For MS-Windows only

    Instructions:
    -------------
    - Create a ConsoleGraphics object, attatched to the current CMD window
      See included examples for more information.
    - Create a CBitmap and draw on this bitmap.
    - call CBitmap::Render() to display the content of the bitmap.
    - call CBitmap::Write() to save the content of the bitmap into a 32-bit Windows Bitmap file.
    - call CBitmap::Load() to load an uncompressed Windows Bitmap file.
    
    Requirements:
    - Windows 7 or later.
    - C++11
    - WinGDI

	To do:
	- Convert this into a simple 2D engine
	- Support rendering into a separated window
	- Port to Linux and macOS (pretty vague)
	
    Lastest update: Mar, 2024
*/

#ifndef __CONSOLE_GRAPPHICS__
#define __CONSOLE_GRAPPHICS__

#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <windows.h>

namespace CGraph {

#ifndef M_PI
	#define M_PI	3.14159265358979323846
#endif

// minimal ANSI stuffs
#define	VT_COLOR(v)		(std::string("\x1B[38;5;"#v"m"))
#define	VT_BKGD(v)		(std::string("\x1B[48;5;"#v"m"))
#define	VT_ATTR(v)		(std::string("\x1B["#v"m"))
#define	VT_CLEAR		(std::string("\x1B[H\x1B[2J"))
#define	VT_CURSORSHOW	(std::string("\x1B[?25h"))
#define	VT_CURSORHIDE	(std::string("\x1B[?25l"))
#define	VT_DEFAULTATTR	(std::string("\x1B[0m"))
#define	VT_CURSORHOME	(std::string("\x1B[H"))

typedef void CGRAPH_CALLBACK(int cur, int total);

// Mathematics utility class
class CMath {
public:
	// linear interpolation
	static double lerp(double a, double b, double f) {
		if (f <= 0.5f)
        	return a + (b - a) * f;
    	else
        	return b - (b - a) * (1.0 - f);
	}	
	
	// map a value from [imin, imax] to [omin, omax]
	static double map(double v, double imin, double imax, double omin, double omax	) {
		if (imax == imin) return omin;
		return (v - imin) * (omax - omin) / (imax - imin) + omin;
	}
	
	// clamping
	static double clampf(double fv, double minf, double maxf) {
		if (fv < minf) fv = minf;
		if (fv > maxf) fv = maxf;
		return fv;
	}
	
	static long clampl(long lv, long minl, long maxl) {
		if (lv < minl) lv = minl;
		if (lv > maxl) lv = maxl;
		return lv;
	}
	
	// find max
	static double maxf(double fa, double fb) {
		return (fa >= fb) ? fa : fb;
	}

	static long maxl(long la, long lb) {
		return (la >= lb) ? la : lb;
	}

	// find min
	static double minf(double fa, double fb) {
		return (fa < fb) ? fa : fb;
	}
	
	static long minl(long la, long lb) {
		return (la < lb) ? la : lb;
	}
	
	// absolute
	static double absf(double v) {
		return v < 0.0f ? -v : v;
	}
	
	static long absl(long l) {
		return l < 0 ? -l : l;
	}	
		
};

// CRect
// -----
// Represents a rectangular area
class CRect {
protected:
	double l, t, r, b;
	double width, height;
public:
	CRect() { l = t = r = b = width = height = 0.0f; }
	CRect(double l, double t, double r, double b) { SetRect(l, t, r, b); }
	CRect(const CRect & rhs) { SetRect(rhs.l, rhs.t, rhs.r, rhs.b); }
	void SetRect(double l, double t, double r, double b) {
		this->l = l; this->t = t; this->r = r; this->b = b;
		this->width  = r-l; this->height = b-t;
	}
	double Left() const { return this->l; }	
	double Top() const { return this->t; }	
	double Right() const { return this->r; }	
	double Bottom() const { return this->b; }	
	double Width() const { return width; }
	double Height() const { return height; }
};

// CPoint
// ------
// Represents a 2D point
class CPoint {
protected:
	double	x, y;
public:	
	CPoint() { x = y = 0.0f; }
	CPoint(double x, double y) { SetPosition(x, y); }
	CPoint(const CPoint & rhs) { SetPosition(rhs.x, rhs.y); }
	void SetPosition(double x, double y) { this->x = x; this->y = y; }
	double GetX() const { return x; }		
	double GetY() const { return y; }
	double Distance(const CPoint & p) {
		double dx = this->x - p.x, dy = this->y - p.y;
		return sqrt(dx*dx+dy*dy);
	}
};

// CVector2D
// ---------
// Represents a vector in 2D space
class CVector2D {
	union {
		struct {
			double x, y;
		};
		double v[2];
	};

public:

	CVector2D() : x(0.0f), y(0.0f){ }
	CVector2D(double x, double y) { this->x = x; this->y = y; }
	CVector2D(const CVector2D & rhs) {
		*this = rhs;
	}
	
	double getX() const { return this->x; }	
	double getY() const { return this->y; }	
	CVector2D & setX(double x) {
		this->x = x;
		return * this;
	}
	CVector2D & setY(double y) {
		this->y = y;
		return * this;
	}
	
	CVector2D & incX(double dx) {
		this->x += dx;
		return * this;
	}
	
	CVector2D & incY(double dy) {
		this->y += dy;
		return * this;
	}

	double magnitude() const { 
		return sqrt(x * x + y * y); 
	}

	CVector2D & normalize() { 
		double mag = magnitude();
		if (mag != 0) {
			x /= mag;
			y /= mag;
		}
		return *this;		
	}

	CVector2D & shuffle() {
		double t = x; x = y; y = t;
		return * this;
	}
	
	double dot(const CVector2D v) const {
		return this->x * v.x + this->y * v.y;
	}	
	
	// cross-product is not defined in 2D, however 
	// we return the magnitude of the z component in 3D space
	double cross(const CVector2D & v) const {
		return this->x * v.y - this->y * v.x;
	}	
	
	double operator[] (int index) const {
		return v[index];
	}

	CVector2D & operator = (const CVector2D & rhs) {
		if (&rhs != this) {
			this->x = rhs.x;
			this->y = rhs.y;
		}		
		return * this;
	}
	
	CVector2D & operator += (CVector2D v) {
		this->x += v.x;
		this->y += v.y;
		return * this;
	}
	
	CVector2D & operator += (double d) {
		this->x += d;
		this->y += d;
		return * this;	
	}
	
	CVector2D & move(double dx, double dy) {
		this->x += dx;
		this->y += dy;
		return * this;
	}
	
	CVector2D & scale(double factor) {
		this->x *= factor;
		this->y *= factor;
		return * this;
	}
	
	CVector2D & rotate(double theta) {
		double ox = x;
		double oy = y;
		this->x = ox * cos(theta) - oy * sin(theta);
		this->y = ox * sin(theta) + oy * cos(theta);
		return * this;	
	}
	
	double distance(const CVector2D & v) const {
		double dx = v.x - this->x;
		double dy = v.y - this->y;
		return sqrt(dx * dx + dy * dy);	
	}
	
	friend std::ostream & operator << (std::ostream & os, const CVector2D & v) {
		os << "(" << v.x << ", " << v.y << ")";
		return os;
	}
		
	friend std::istream & operator >> (std::istream & is, CVector2D & v) {
		is >> v.x >> v.y;
		return is;
	}
	
	friend double distance(const CVector2D & u, const CVector2D & v) {
		double dx = v.x - u.x;
		double dy = v.y - u.y;
		return sqrt(dx * dx + dy * dy);	
	}
	
	friend double dot(const CVector2D & u, const CVector2D & v) {
		return u.x * v.x + u.y * v.y;
	}

	friend CVector2D operator + (const CVector2D & u, const CVector2D & v) {
		return CVector2D(u.x + v.x, u.y + v.y);
	}	
	
	friend CVector2D operator - (const CVector2D & u, const CVector2D & v) {
		return CVector2D(u.x - v.x, u.y - v.y);
	}	
	
	friend bool operator == (const CVector2D & u, const CVector2D & v) {
		return u.magnitude() == v.magnitude();
	}

	friend bool operator != (const CVector2D & u, const CVector2D & v) {
		return u.magnitude() != v.magnitude();
	}
	
	static CVector2D zero;
};

CVector2D CVector2D::zero(0.0f, 0.0f);

class CRGB;	// forward declaration to be used inside CHSV

// CHSV
// ----
// Represents the HSV color space, allowing conversion to RGB.
// Requires class CRGB.
class CHSV {
protected:
	double h, s, v;	// h = [0, 359];  s, v in [0, 1]
public:
	CHSV() {  SetColor(0, 0, 0.5); }
	CHSV(double h, double s, double v) { SetColor(h, s, v); }
	CHSV(const CHSV & rhs) { SetColor(rhs); }
	void SetColor(double h, double s, double v) { this->h = h; this->s = s; this->v = v; }
	void SetColor(const CHSV & rhs) { this->h = rhs.h; this->s = rhs.s; this->v = rhs.v; }
	CRGB ToRGB();
	double Hue() const { return h; }
	double Saturation() const { return s; }
	double Value() const { return v; }
	CHSV & addHue(double dh) {
		this->h += dh;
		this->h = CMath::clampf(this->h, 0.0f, 359.0f);
		return *this;
	}
	CHSV & addSaturation(double ds) {
		s += ds;
		this->s = CMath::clampf(this->s, 0.0f, 1.0f);
		return *this;		
	}
	CHSV & addValue(double dv) {
		this->v += dv;
		this->v = CMath::clampf(this->v, 0.0f, 1.0f);
		return *this;	
	}
};

// CRGB
// ----
// Represents the sRGB color space, allowing conversion to HSV.
// Requires class CHSV.
class CRGB {
protected:
	unsigned char r, g, b;	// r, g, b in [0, 255]
public:
	CRGB() : r(0), g(0), b(0) { }
	CRGB(unsigned char r, unsigned char g, unsigned char b) { SetColor(r, g, b); }	
	CRGB(double r, double g, double b) { SetColorf(r, g, b); }
	CRGB(const CRGB & rhs) { SetColor(rhs); }
	void SetColor(unsigned char r, unsigned char g, unsigned char b) {
		this->r = r; this->g = g; this->b = b;
	}
	void SetColorf(double r, double g, double b) {
		this->r = r * 255.0f; this->g = g * 255.0f; this->b = b * 255.0f;
	}
	void SetColor(const CRGB & rhs) { SetColor(rhs.r, rhs.g, rhs.b); }
	COLORREF GetColor() const { return RGB(r, g, b); }
	CHSV ToHSV();
	double Redf() const { return r/255.0; }
	double Greenf() const { return g/255.0; }
	double Bluef() const { return b/255.0; }
	int Red() const { return r; }
	int Green() const { return g; }
	int Blue() const { return b; }
	
	CRGB & addRed(int dr) {
		this->r += dr;
		this->r = CMath::clampl(this->r, 0, 255);
		return *this;	
	}
	CRGB & addGreen(int dg) {
		this->g += dg;
		this->g = CMath::clampl(this->g, 0, 255);
		return *this;	
	}
	CRGB & addBlue(int db) {
		this->b += db;
		this->b = CMath::clampl(this->b, 0, 255);
		return *this;	
	}
};

CHSV CRGB::ToHSV() {
	double	min, max, delta, fh, fs, fv;
	double	fr = Redf(), fg = Greenf(), fb = Bluef();

    min = fr  < fg ? fr : fg;
    min = min < fb ? min: fb;

    max = fr  > fg ? fr : fg;
    max = max > fb ? max: fb;

    fv = max;
    delta = max - min;
    if (delta < 0.00001) {
        fs = 0;
        fh = 0; // undefined, maybe nan?
        return CHSV(fh, fs, fv);
    }
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        fs = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0              
        // s = 0, h is undefined
        fs = 0.0;
        fh = NAN;                            // its now undefined
        return CHSV(fh, fs, fv);
    }
    if( fr >= max )                           // > is bogus, just keeps compilor happy
        fh = ( fg - fb ) / delta;        // between yellow & magenta
    else
    if( fg >= max )
        fh = 2.0 + ( fb - fr ) / delta;  // between cyan & yellow
    else
        fh = 4.0 + ( fr - fg ) / delta;  // between magenta & cyan

    fh *= 60.0;                       // degrees

    if( fh < 0.0 )
        fh += 360.0;
    return CHSV(fh, fs, fv);
}

CRGB CHSV::ToRGB() {
	CRGB out;
	double hh, p, q, t, ff;
    long i;

    if (s <= 0.0)
    	return CRGB(v, v, v);
    hh = h;
    if (hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * ff));
    t = v * (1.0 - (s * (1.0 - ff)));

    switch(i) {
    case 0: out.SetColorf(v, t, p); break;
    case 1: out.SetColorf(q, v, p); break;
    case 2: out.SetColorf(p, v, t); break;
    case 3: out.SetColorf(p, q, v); break;
    case 4: out.SetColorf(t, p, v); break;
    case 5:
    default: out.SetColorf(v, p, q); break;
    }
	return out;
}

// ConsoleGraphics
// ---------------
// A simple trick to render graphics onto CMD window. 
// However, any changes to the CMD window would clear the rendered image.
// Only support CMD. Does not work in Windows Terminal or PowerShell.
class ConsoleGraphics {
protected:
	HWND wnd;
	HANDLE console;
	DWORD mode;
public:
	ConsoleGraphics() {
		this->wnd = GetConsoleWindow();
		this->console = GetStdHandle(STD_OUTPUT_HANDLE);
		if (this->console) {
			GetConsoleMode(this->console, &mode);
			SetConsoleMode(this->console, mode | 0x0004);
		}
		ShowCursor(false);
	}
	
	~ConsoleGraphics() {
		if (this->console)
			SetConsoleMode(this->console, mode);
		ShowCursor(true);
	}	
	
	HWND GetWindow() {
		return this->wnd;
	}
	
	void ShowCursor(bool show) {
   		CONSOLE_CURSOR_INFO info = {100, show};
   		SetConsoleCursorInfo(this->console, &info);
	}
	
	static void flush() {
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	
	static void wait() {
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
};

// CBitmap
// -------
// Core class for rendering stuffs into a bitmap object
class CBitmap {
protected:
	static const int signature = 0x4D42;	// 'MB', Windows Bitmap signature
	int width;
	int height;
	struct rgba { unsigned char b, g, r, a; };
	std::vector<struct rgba> data;
	BITMAPFILEHEADER bmhdr;
	BITMAPINFOHEADER bminf;
	bool loaded;

	// for plottings
	CRGB color;		// current plotting color
	CPoint cursor;	// virtual cursor

	bool InBounds(int x, int y) const {
		return x >= 0 && y >= 0 && x < width && y < height;
	}

public:
	// default constructor
	CBitmap() : width(0), height(0), loaded(false) {
		this->color.SetColor((unsigned char)255, (unsigned char)255, (unsigned char)255);
	}

	// parameterized constructor	
	CBitmap(int w, int h) : width(0), height(0), loaded(false) {
		Create(w, h);
		this->color.SetColor((unsigned char)255, (unsigned char)255, (unsigned char)255);
	}
	
	// constructor to load a Windows Bitmap from file
	CBitmap(std::string filename) : width(0), height(0), loaded(false) {
		Load(filename);	
		this->color.SetColor((unsigned char)255, (unsigned char)255, (unsigned char)255);		
	}
	
	// no copy constructor support at the moment
	
	// destructor
	~CBitmap() { Clear(); }
	
	// Creates a bitmap of size w and h
	bool Create(int w, int h) {	
		if (w <= 0 || h <= 0) {
			Clear();
			width = height = 0;
			loaded = false;
			return false;
		}
		size_t sw = static_cast<size_t>(w);
		size_t sh = static_cast<size_t>(h);
		if (sw > std::numeric_limits<size_t>::max() / sh)
			return false;
		width = w;
		height = h;
		data.resize(sw * sh);
		loaded = false;
		return true;
	}
	
	// Clears the content of the bitmap
	void Clear() {
		data.clear();
		loaded = false;
	}

	// Directly put pixel into the bitmap buffer
	void SetPixel(int x, int y, const CRGB & pix) {
		SetPixel(x, y, pix.Red(), pix.Green(), pix.Blue());
	}	
	
	// Same as previous methods, using r, g, b values
	void SetPixel(int x, int y, int r, int g, int b) {
		if (!InBounds(x, y)) return;
		int offset = (x + y * this->width);
		data[offset].r = r;
		data[offset].g = g;
		data[offset].b = b;
	}
	
	// Directly gets pixel from the bitmap buffer
	CRGB GetPixel(int x, int y) const {
		if (!InBounds(x, y)) return CRGB(0, 0, 0);
		int offset = (x + y * this->width);
		return CRGB(data[offset].r, data[offset].g, data[offset].b);
	}
	
	// Loads an uncompressed 1, 4, 8, 24 or 32-bit Windows Bitmap file
	bool Load(std::string filename, CGRAPH_CALLBACK * cb = NULL) {
		RGBQUAD	clut[256];	
		std::ifstream is;
		char * row = NULL;
		int linew, pixelw;
		loaded = false;
		
		// open the file for reading
		is.open(filename.c_str(), std::ios_base::binary);
		if (!is) return false;
	
		// read bitmap header
		is.read((char *)&bmhdr, sizeof(BITMAPFILEHEADER));
		if (!is) {
			is.close();
			return false;
		}
		if (bmhdr.bfType != signature) {
			is.close();
			return false;
		}
		
		// read bitmap info header
		is.read((char *)&bminf, sizeof(BITMAPINFOHEADER));
		if (!is) {
			is.close();
			return false;			
		}
		
		if (bminf.biCompression != BI_RGB || bminf.biPlanes != 1 || bminf.biSize != 40) {
			is.close();
			return false;
		}

		if (bminf.biWidth <= 0 || bminf.biHeight <= 0) {
			is.close();
			return false;
		}

		if (!(bminf.biBitCount == 1 || bminf.biBitCount == 4 || bminf.biBitCount == 8 ||
		      bminf.biBitCount == 24 || bminf.biBitCount == 32)) {
			is.close();
			return false;
		}
		
		// compute the width of a scanline, must be 32-bit aligned
		long long bitsPerLine = static_cast<long long>(bminf.biWidth) * bminf.biBitCount;
		if (bitsPerLine <= 0) {
			is.close();
			return false;
		}
		linew = static_cast<int>((((bitsPerLine + 31) & ~31LL) >> 3));
		pixelw = bminf.biBitCount >> 3;
		if (linew <= 0) {
			is.close();
			return false;
		}
		
		// allocate memory enough to hold a scanline
		row = new char[linew];

		// read up Color-Look-Up Table
		if (bminf.biBitCount <= 8) {
			is.read((char *)clut, sizeof(RGBQUAD) * (1 << bminf.biBitCount));
			if (!is) {
				delete []row;
				is.close();
				return false;
			}
		}
		
		// seek to bitmap data
		is.seekg (bmhdr.bfOffBits, std::ios::beg);
		if (!is) {
			delete []row;
			is.close();
			return false;
		}
				
		// initialize the bitmap
		Clear();
		if (!Create(bminf.biWidth, bminf.biHeight)) {
			delete []row;
			is.close();
			return false;
		}

		switch(bminf.biBitCount) {	
		case 1:		// monochrome bitmap
			for (int i = 0; i < bminf.biHeight; i++) {
				int offset = (bminf.biHeight-1-i)*bminf.biWidth;
				is.read(row, linew);
				if (!is) {
					delete []row;
					is.close();
					return false;
				}
				int j = 0; 
				while(j < bminf.biWidth) {
					unsigned char c = (unsigned char)row[j >> 3];
					int k = 0;
					while (k < 8 && j < bminf.biWidth) {
						if ((c & 0x80)) {
							data[offset].r = clut[1].rgbRed;
							data[offset].g = clut[1].rgbGreen;
							data[offset].b = clut[1].rgbBlue;							
						}
						else {
							data[offset].r = clut[0].rgbRed;
							data[offset].g = clut[0].rgbGreen;
							data[offset].b = clut[0].rgbBlue;												
						}
						c = c << 1;
						k++; j++; offset++;
					}		
				}
				if (cb) (*cb)(i, bminf.biHeight);	// invoke the callback function
			}
			break;
		case 4:		// indexed 4-bit
			for (int i = 0; i < bminf.biHeight; i++) {
				int offset = (bminf.biHeight-1-i)*bminf.biWidth;
				is.read(row, linew);
				if (!is) {
					delete []row;
					is.close();
					return false;
				}
				int j = 0; 
				while(j < bminf.biWidth) {
					unsigned char c = row[j >> 1];
					int k = 0;
					while (k < 2 && j < bminf.biWidth) {
						int v = (c & 0xF0) >> 4;
						data[offset].r = clut[v].rgbRed;
						data[offset].g = clut[v].rgbGreen;
						data[offset].b = clut[v].rgbBlue;
						c = c << 4;
						k++; j++; offset++;
					}		
				}
				if (cb) (*cb)(i, bminf.biHeight);	// invoke the callback function
			}
			break;
		case 8:		// indexed 8-bit 
			for (int i = 0; i < bminf.biHeight; i++) {
				int offset = (bminf.biHeight-1-i)*bminf.biWidth;
				is.read(row, linew);
				if (!is) {
					delete []row;
					is.close();
					return false;
				}
				for (int j = 0; j < bminf.biWidth; j++) {
					unsigned char c = row[j];
					data[offset].r = clut[c].rgbRed;
					data[offset].g = clut[c].rgbGreen;
					data[offset].b = clut[c].rgbBlue;
					offset++;
				}
				if (cb) (*cb)(i, bminf.biHeight);	// invoke the callback function
			}
			break;
		case 24:	// 24-bit RGB
		case 32:	// 32-bit RGBA
			for (int i = 0; i < bminf.biHeight; i++) {
				int offset = (bminf.biHeight-1-i) * bminf.biWidth;
				is.read(row, linew);
				if (!is) {
					delete []row;
					is.close();
					return false;
				}
				// swap BGR to RGB
				for (int j = 0; j < bminf.biWidth; j++) {
					data[offset].r = row[j*pixelw + 2];
					data[offset].g = row[j*pixelw + 1];
					data[offset].b = row[j*pixelw + 0];
					offset++;
				}
				if (cb) (*cb)(i, bminf.biHeight);	// invoke the callback function
			}
			break;
		default:
			delete []row;
			is.close();
			return false;
		}
		
		// finalizing
		delete []row;
		is.close();
		loaded = true;
		return true;
	}
	
	// Saves current buffer into a 32-bit uncompressed Windows Bitmap file
	bool Write(std::string filename, CGRAPH_CALLBACK * cb = NULL) {
		if (!data.size() || !width || !height)
			return false;

		std::ofstream os(filename.c_str(), std::ios::out | std::ios::binary);
		if (os.fail()) return false;

		// Prepares neccessaries informations, headers
		unsigned long size = data.size() * sizeof(struct rgba);
		unsigned long bits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
		unsigned long linew = width * sizeof(struct rgba);

		// setup Windows Bitmap File Header structure - refers to MSDN for more details
		BITMAPFILEHEADER bh = { signature,	// bfType
								bits + size,// bfSize
								0,			// bfReserved1
								0,			// bfReserved2
								bits		// bfOffBits
							};
							
		// setup Window Bitmap Infomation Header structure - refers to MSDN for more details
		BITMAPINFOHEADER bi = {	40,							// biSize
								width,						// biWidth
								height,						// biHeight
								1,							// biPlanes
								sizeof(struct rgba) << 3,	// biBitCount
								BI_RGB,						// biCompression
								size,						// biSizeImage
								2835,						// biXPelsPerMeter, 72dpi
								2835, 						// biYPelsPerMeter, 72dpi
								0, 							// biClrUsed
								0							// biClrImportant
							};
						
		// write the header
		os.write((char *)&bh, sizeof(BITMAPFILEHEADER));
		if (os.fail()) {
			os.close();
			return false;
		}
		
		// write the info header
		os.write((char *)&bi, sizeof(BITMAPINFOHEADER));
		if (os.fail()) {
			os.close();
			return false;
		}
		
		// write the image, upside down
		for (int row = 0; row < height; row++) {
			const struct rgba * line = &data[(height - row - 1) * width];
			// write down a bitmap row
			os.write((char *)line, linew);
			if (os.fail()) {
				os.close();
				return false;
			}
			if (cb) (*cb)(row, height);	// invoke the callback function
		}
		os.close();
		return true;
	}

	// Displays information about a loaded Windows Bitmap image
	void Info() {
		int depth;
		if (loaded) {
			std::cout << "[" << VT_COLOR(192) << bminf.biWidth << VT_DEFAULTATTR << " x " << VT_COLOR(192) << bminf.biHeight << VT_DEFAULTATTR << ", ";
			depth = bminf.biBitCount;
		}
		else {
			std::cout << "[" << VT_COLOR(192) << width << VT_DEFAULTATTR << " x " << VT_COLOR(192) << height << VT_DEFAULTATTR << ", ";
			depth = sizeof(struct rgba) * 8;
		}		
		std::cout << VT_COLOR(200) << depth << VT_DEFAULTATTR << " bbp (";
		switch(depth) {
		case  1: std::cout << VT_COLOR(231) << "Monochrome"; break;
		case  4: std::cout << VT_COLOR(104) << "16 colors"; break;
		case  8: std::cout << VT_COLOR(104) << "256 colors"; break;
		case 24: std::cout << VT_COLOR(111) << "16M colors"; break;
		case 32: std::cout << VT_COLOR(111) << "True-color"; break;
		default: std::cout << VT_COLOR(208) << "Unknown" ;break;
		}
		std::cout << VT_DEFAULTATTR << ")]" << std::endl;
	}
	
	// Renders the bitmap onto a console screen
	void Render(ConsoleGraphics & cg) {
		Render(cg.GetWindow());	
	}

	// Renders the bitmap onto a Window via its handle
	void Render(HWND wnd, bool stretch = false) {
		if (!wnd || width <= 0 || height <= 0 || data.empty())
			return;

		BITMAPINFO bmi = {
			.bmiHeader = {
				.biSize = sizeof(BITMAPINFOHEADER),
				.biWidth = width,
				.biHeight = -height,
				.biPlanes = 1,
				.biBitCount = 32,	
				.biCompression = BI_RGB,
				.biSizeImage = 0,
				.biXPelsPerMeter = 0,
				.biYPelsPerMeter = 0,
				.biClrUsed = 0,
				.biClrImportant = 0
			}
		};
		HDC dc = GetDC(wnd);
		RECT rc;
		struct rgba * bits = &data[0];		
		
		GetClientRect(wnd, &rc);
		FillRect(dc, &rc, (HBRUSH) GetStockObject(BLACK_BRUSH));
		if (stretch) {
			StretchDIBits(dc, 
						  0, 0, rc.right-rc.left+1, rc.bottom-rc.top-1, 
						  0, 0, width, height, 
						  bits, &bmi, DIB_RGB_COLORS, SRCCOPY);
		}
		else
			StretchDIBits(dc, 
						  0, 0, width, height, 
						  0, 0, width, height, 
						  bits, &bmi, DIB_RGB_COLORS, SRCCOPY);
		ReleaseDC(wnd, dc);
	}

	int Width() const { return width; }
	int Height() const { return height; }

	// Plots a pixel to position (x, y) in the bitmap buffer, using current color
	void PutPixel(int x, int y) { SetPixel(x, y, color); }

	// Set active color for plotting using R, G, B triplet
	void SetColor(unsigned char r, unsigned char g, unsigned char b) {
		this->color.SetColor(r, g, b);
	}
	
	// Set active color for plotting using a CRGB class
	void SetColor(const CRGB & pix) {
		this->color.SetColor(pix);
	}

	// Draws a horizontal line from x1 to x2 at row y
	void HLine(int x1, int x2, int y) {
		if (width <= 0 || height <= 0 || y < 0 || y >= height) return;
		if (x1 > x2) std::swap(x1, x2);
		x1 = CMath::maxl(x1, 0);
		x2 = CMath::minl(x2, width - 1);
		if (x1 > x2) return;
		int offset = (x1 + y * this->width);
		for (int i = x1; i <= x2; i++, offset++) {
			data[offset].r = color.Red();
			data[offset].g = color.Green();
			data[offset].b = color.Blue();
		}
	}
	
	// Draws a vertical line from y1 to y2 at column x
	void VLine(int x, int y1, int y2) {
		if (width <= 0 || height <= 0 || x < 0 || x >= width) return;
		if (y1 > y2) std::swap(y1, y2);
		y1 = CMath::maxl(y1, 0);
		y2 = CMath::minl(y2, height - 1);
		if (y1 > y2) return;
		int offset = (x + y1 * this->width);
		for (int i = y1; i <= y2; i++, offset += this->width) {
			data[offset].r = color.Red();
			data[offset].g = color.Green();
			data[offset].b = color.Blue();
		}
	}
	
	// Draws a rectangular frame at (l, r) to (r, b)
	void DrawRect(int l, int t, int r, int b) {
		DrawLine(l, t, r, t);
		DrawLine(l, b, r, b);
		DrawLine(l, t, l, b);
		DrawLine(r, t, r, b);
	}
	
	// Draws a rectangular frame at (l, r) of size (w, h)
	void DrawRectD(int l, int t, int w, int h) {
		DrawRect(l, t, l + w, t + h);
	}
	
	// Draws a rectangular box at (l, r) to (r, b)
	void Rectangle(int l, int t, int r, int b) {
		for (int i = t; i <= b; i++)
			HLine(l, r, i);
	}
	
	// Draw a rectangular box at (l, r) of size (w, h)
	void RectangleD(int l, int t, int w, int h) {
		Rectangle(l, t, l + w, t + h);
	}
	
	// Draws a generic line using MidPoint algorithm
	void DrawLine(int l, int t, int r, int b) {
    	int dx = r - l;
    	int dy = b - t;
    	int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);
		if (steps == 0) {
			PutPixel(l, t);
			return;
		}
    	double xinc = dx / (double)steps;
    	double yinc = dy / (double)steps;
    	double x = l;
    	double y = t;
    	for (int i = 0; i <= steps; i++) {
        	PutPixel(x, y);
        	x += xinc;
        	y += yinc;
    	}
	}

	// Moves current virtual cursor to (x, y)
	void MoveTo(int x, int y) {
		cursor.SetPosition(x, y);
	}

	// Draws a line from current position to (x, y) and then move the virtual cursor to (x, y)
	void LineTo(int x, int y) {
		DrawLine(this->cursor.GetX(), this->cursor.GetY(), x, y);
		MoveTo(x, y);
	}

	// Draws a circle at (xc, yc) with radius r
	void DrawCircle(int xc, int yc, int r) {
		int x = 0, y = r;
		int d = 3 - (r << 1);
		CirclePixels(xc, yc, x, y);
		while (y >= x) {
			x++;
			if (d > 0) {
				y--;
				d = d + ((x - y) << 2) + 10;
			}
			else
				d = d + (x << 2) + 6;
			CirclePixels(xc, yc, x, y);
		}
	}

	// Draws a filled circle at (xc, yc) with radius r
	void DrawFilledCircle(int xc, int yc, int r) {
		int x = 0, y = r;
		int d = 3 - (r << 1);
		CircleLines(xc, yc, x, y);
		while (y >= x) {
			x++;
			if (d > 0) {
				y--;
				d = d + ((x - y) << 2) + 10;
			}
			else
				d = d + (x << 2) + 6;
			CircleLines(xc, yc, x, y);
		}
	}

private:
	void CirclePixels(int xc, int yc, int x, int y)	{
		PutPixel(xc+x, yc+y);
		PutPixel(xc-x, yc+y);
		PutPixel(xc+x, yc-y);
		PutPixel(xc-x, yc-y);
		PutPixel(xc+y, yc+x);
		PutPixel(xc-y, yc+x);
		PutPixel(xc+y, yc-x);
		PutPixel(xc-y, yc-x);
	}

	void CircleLines(int xc, int yc, int x, int y)	{
		DrawLine(xc+x, yc+y, xc-x, yc+y);
		DrawLine(xc+x, yc-y, xc-x, yc-y);
		DrawLine(xc+y, yc+x, xc-y, yc+x);
		DrawLine(xc+y, yc-x, xc-y, yc-x);
	}
};


// WORK-IN-PROGRESS
class CApp {
	HWND	hWnd;
	CBitmap	surface;
	
	CApp(std::string title) {
	}
	
	~CApp() {
		
	}
	
	// create window
	virtual void create() {
		
	}
	
	// window paint
	virtual void paint() {
		
	}
	
	// window update
	virtual void update() {
		
	}
	
	// window timer
	virtual void timer() {
	}
};

}

#endif

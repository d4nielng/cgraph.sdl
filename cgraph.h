#if !(defined(__CGRAPH_H__))
#define __CGRAPH_H__

#define SDL_MAIN_HANDLED

/*
    Single header for simple graphics programming using SDL2
     _______  _______  ______    _______  _______  __   __        __   __ 
    |       ||       ||    _ |  |   _   ||       ||  | |  |      |  | |  |
    |       ||    ___||   | ||  |  |_|  ||    _  ||  |_|  |      |  |_|  |
    |       ||   | __ |   |_||_ |       ||   |_| ||       |      |       |
    |      _||   ||  ||    __  ||       ||    ___||       | ___  |       |
    |     |_ |   |_| ||   |  | ||   _   ||   |    |   _   ||   | |   _   |
    |_______||_______||___|  |_||__| |__||___|    |__| |__||___| |__| |__|
                                                Coded by Trinh D.D. Nguyen
    
    Requirements:
    - C++11
    - SDL2

    Changes:
    - v0.1:
        Initial development, based on Window GDI library
        Support Windows only
    - v0.2:
        Completely switch to SDL2 for cross-platform purpose
        Support Windows, macOS and Linux
	
    To do:
	- Convert this into a simple 2D engine
	
    Usages:
    - Include this header file in your project
    - Linking with SDL2 library (using the parameter -lSDL2)
    - Enjoy coding with simple graphics
    
    Permissions:
    - Feel free to modify the source code to fit your needs

    Lastest update: Mar, 2025
*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <SDL.h>

namespace daniel {

#ifndef M_PI
	#define M_PI	3.14159265358979323846
#endif

typedef uint32_t COLORREF;
typedef void CGRAPH_CALLBACK(int current, int total);

#ifndef RGBQUAD
    typedef struct __attribute__((packed)) rgbquad {
        uint8_t rgbBlue;
        uint8_t rgbGreen;
        uint8_t rgbRed;
        uint8_t rgbReserved;
    } RGBQUAD;
#endif

#ifndef BITMAPFILEHEADER
    typedef struct __attribute__((packed)) bitmapfileheader {
        uint16_t bfType;
        uint32_t bfSize;
        uint16_t bfReserved1;
        uint16_t bfReserved2;
        uint32_t bfOffBits;
    } BITMAPFILEHEADER;
#endif

#ifndef BITMAPINFOHEADER
    typedef struct __attribute__((packed)) tagBITMAPINFOHEADER {
        uint32_t biSize;
        uint32_t biWidth;
        uint32_t biHeight;
        uint16_t biPlanes;
        uint16_t biBitCount;
        uint32_t biCompression;
        uint32_t biSizeImage;
        uint32_t biXPelsPerMeter;
        uint32_t biYPelsPerMeter;
        uint32_t biClrUsed;
        uint32_t biClrImportant;
    } BITMAPINFOHEADER;
#endif

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
    double w, h;
public:
    CRect() { l = t = r = b = w = h = 0.0f; }
    CRect(double l, double t, double r, double b) { setRect(l, t, r, b); }
    CRect(const CRect & rhs) { setRect(rhs.l, rhs.t, rhs.r, rhs.b); }
    void setRect(double l, double t, double r, double b) {
        this->l = l; this->t = t; this->r = r; this->b = b;
        this->w  = r-l; this->h = b-t;
    }
    double left() const { return this->l; }	
    double top() const { return this->t; }	
    double right() const { return this->r; }	
    double bottom() const { return this->b; }	
    double width() const { return w; }
    double height() const { return h; }
};

// CPoint
// ------
// Represents a 2D point
class CPoint {
protected:
    double	x, y;
public:	
    CPoint() { x = y = 0.0f; }
    CPoint(double x, double y) { setPosition(x, y); }
    CPoint(const CPoint & rhs) { setPosition(rhs.x, rhs.y); }
    void setPosition(double x, double y) { this->x = x; this->y = y; }
    double getX() const { return x; }		
    double getY() const { return y; }
    double distance(const CPoint & p) {
        double dx = this->x - p.x;
        double dy = this->y - p.y;
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
	CVector2D(const CVector2D & rhs) { *this = rhs; }
	
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
		this->x = x * cos(theta) - y * sin(theta);
		this->y = x * sin(theta) + y * cos(theta);
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
    CHSV() {  setColor(0, 0, 0.5); }
    CHSV(double h, double s, double v) { setColor(h, s, v); }
    CHSV(const CHSV & rhs) { setColor(rhs); }
    void setColor(double h, double s, double v) { this->h = h; this->s = s; this->v = v; }
    void setColor(const CHSV & rhs) { this->h = rhs.h; this->s = rhs.s; this->v = rhs.v; }
    CRGB toRGB();
    double hue() const { return h; }
    double saturation() const { return s; }
    double value() const { return v; }
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

#define MAKERGB(r, g, b)    ((uint32_t)((r) << 16) | \
                            ((uint32_t)((g) << 8)) | \
                             (uint32_t)(b))
#define RGB_RED(color)      ((color >> 16) & 0xFF)
#define RGB_GREEN(color)    ((color >> 8) & 0xFF)
#define RGB_BLUE(color)     (color & 0xFF)
                
// CRGB
// ----
// Represents the sRGB color space, allowing conversion to HSV.
// Requires class CHSV.
class CRGB {
protected:
    uint8_t r, g, b;	// r, g, b in [0, 255]
public:
    CRGB() : r(0), g(0), b(0) { }
    CRGB(uint8_t r, uint8_t g, uint8_t b) { setColor(r, g, b); }	
    CRGB(double r, double g, double b) { setColorf(r, g, b); }
    CRGB(const CRGB & rhs) { setColor(rhs); }
    CRGB(COLORREF color) { setColor(RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color)); }
    void setColor(uint8_t r, uint8_t g, uint8_t b) {
        this->r = r; this->g = g; this->b = b;
    }
    void setColorf(double r, double g, double b) {
        this->r = r * 255.0f; this->g = g * 255.0f; this->b = b * 255.0f;
    }
    void setColor(const CRGB & rhs) { setColor(rhs.r, rhs.g, rhs.b); }
    uint32_t getColor() const { return MAKERGB(r, g, b); }
    CHSV toHSV();
    double redf() const { return r/255.0; }
    double greenf() const { return g/255.0; }
    double bluef() const { return b/255.0; }
    uint8_t red() const { return r; }
    uint8_t green() const { return g; }
    uint8_t blue() const { return b; }
    
    operator uint32_t() const { return MAKERGB(r, g, b); }

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

CHSV CRGB::toHSV() {
    double	min, max, delta, fh, fs, fv;
    double	fr = redf(), fg = greenf(), fb = bluef();

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

CRGB CHSV::toRGB() {
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
    case 0: out.setColorf(v, t, p); break;
    case 1: out.setColorf(q, v, p); break;
    case 2: out.setColorf(p, v, t); break;
    case 3: out.setColorf(p, q, v); break;
    case 4: out.setColorf(t, p, v); break;
    case 5:
    default: out.setColorf(v, p, q); break;
    }
    return out;
}

class CGraph;

// CBitmap
// -------
// Core class for rendering stuffs into a bitmap object
class CBitmap {
protected:
    static const uint32_t signature = 0x4D42;	// 'MB', Windows Bitmap signature
    uint32_t w;
    uint32_t h;
    struct rgba { uint8_t b, g, r, a; };
    std::vector<struct rgba> data;
    BITMAPFILEHEADER bmhdr;
    BITMAPINFOHEADER bminf;
    bool loaded;

public:
    // default constructor
    CBitmap() : w(0), h(0), loaded(false) {}

    // parameterized constructor	
    CBitmap(uint32_t w, uint32_t h) { create(w, h); }
    
    // constructor to load a Windows Bitmap from file
    CBitmap(std::string filename) { load(filename);	 }
    
    // no copy constructor support at the moment
    
    // destructor
    ~CBitmap() { clear(); }
    
    // Creates a 32-bit bitmap of size w and h
    bool create(uint32_t w, uint32_t h) {	
        this->w = w;
        this->h = h;
        data.resize(w * h * sizeof(struct rgba));
        return true;
    }
    
    // Clears the content of the bitmap
    void clear() { data.clear(); }

    // Directly put pixel into the bitmap buffer
    void setPixel(uint32_t x, uint32_t y, const CRGB & pix) {
        setPixel(x, y, pix.red(), pix.green(), pix.blue());
    }	
    
    // Same as previous methods, using r, g, b values
    void setPixel(uint32_t x, uint32_t y, uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255) {
        uint32_t offset = (x + y * this->w);
        data[offset].r = r;
        data[offset].g = g;
        data[offset].b = b;
        data[offset].a = a;
    }
    
    void setPixel(uint32_t x, uint32_t y, uint32_t color) {
        uint32_t offset = (x + y * this->w);
        memcpy(&data[offset], &color, sizeof(uint32_t));
    }

    // Directly gets pixel from the bitmap buffer
    CRGB getPixel(uint32_t x, uint32_t y) const {
        uint32_t offset = (x + y * this->w);
        return CRGB(data[offset].r, data[offset].g, data[offset].b);
    }   

    uint32_t * getPixels() {
        return (uint32_t*) data.data();
    }

    // Loads an uncompressed 1, 4, 8, 24 or 32-bit Windows Bitmap file
    bool load(std::string filename, CGRAPH_CALLBACK * cb = NULL) {
        RGBQUAD	clut[256];	
        std::ifstream is;
        char * row = NULL;
        uint32_t linew, pixelw;
        
        // open the file for reading
        is.open(filename.c_str(), std::ios_base::binary);
        if (!is) return false;

        // read bitmap header
        is.read((char *)&bmhdr, sizeof(BITMAPFILEHEADER));
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

        if (bminf.biCompression != 0 || bminf.biPlanes != 1 || bminf.biSize != 40) {
            is.close();
            return false;
        }
        
        // compute the width of a scanline, must be 32-bit aligned
        linew = ((((bminf.biWidth * bminf.biBitCount) + 31) & ~31) >> 3);
        pixelw = bminf.biBitCount >> 3;
        
        // allocate memory enough to hold a scanline
        row = new char[linew];
        if (!row) {
            is.close();
            return false;
        }

        // read up Color-Look-Up Table
        if (bminf.biBitCount <= 8) 
            is.read((char *)clut, sizeof(RGBQUAD) * (1 << bminf.biBitCount));
        
        // seek to bitmap data
        is.seekg (bmhdr.bfOffBits, std::ios::beg);
                
        // initialize the bitmap
        clear();
        create(bminf.biWidth, bminf.biHeight);

        switch(bminf.biBitCount) {	
        case 1:		// monochrome bitmap
            for (int i = 0; i < bminf.biHeight; i++) {
                uint32_t offset = (bminf.biHeight-1-i) * bminf.biWidth;
                is.read(row, linew);
                if (!is) {
                    delete []row;
                    is.close();
                    return false;
                }
                int j = 0; 
                while(j < bminf.biWidth) {
                    uint32_t c = row[j >> 3];
                    uint32_t k = 0;
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
                uint32_t offset = (bminf.biHeight-1-i)*bminf.biWidth;
                is.read(row, linew);
                if (!is) {
                    delete []row;
                    is.close();
                    return false;
                }
                int j = 0; 
                while(j < bminf.biWidth) {
                    uint8_t c = row[j >> 1];
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
                uint32_t offset = (bminf.biHeight-1-i)*bminf.biWidth;
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
        }
        
        // finalizing
        delete []row;
        is.close();
        loaded = true;
        return true;
    }
    
    SDL_Surface * createSurface() {
        if (!data.size() || !w || !h || !loaded)
            return NULL;
        SDL_Surface * surface = SDL_CreateRGBSurfaceFrom(data.data(), w, h, 32, w * sizeof(struct rgba), 0, 0, 0, 0);
        if (!surface) return NULL;
        return surface;
    }

    // Saves current buffer into a 32-bit uncompressed Windows Bitmap file
    bool write(std::string filename, CGRAPH_CALLBACK * cb = NULL) {
        if (!data.size() || !w || !h)
            return false;

        std::ofstream os(filename.c_str(), std::ios::out | std::ios::binary);
        if (os.fail()) return false;

        // Prepares neccessaries informations, headers
        uint32_t size = data.size() * sizeof(struct rgba);
        uint32_t bits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
        uint32_t linew = w * sizeof(struct rgba);

        // setup Windows Bitmap File Header structure - refers to MSDN for more details
        BITMAPFILEHEADER bh = { signature,	// bfType
                                bits + size,// bfSize
                                0,			// bfReserved1
                                0,			// bfReserved2
                                bits		// bfOffBits
                            };
                            
        // setup Window Bitmap Infomation Header structure - refers to MSDN for more details
        BITMAPINFOHEADER bi = {	40,							// biSize
                                w,						    // biWidth
                                h,						    // biHeight
                                1,							// biPlanes
                                sizeof(struct rgba) << 3,	// biBitCount
                                0,						    // biCompression
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
        for (int row = 0; row < h; row++) {
            const struct rgba * line = &data[(h - row - 1) * w];
            // write down a bitmap row
            os.write((char *)line, linew);
            if (os.fail()) {
                os.close();
                return false;
            }
            if (cb) (*cb)(row, bminf.biHeight);	// invoke the callback function
        }
        os.close();
        return false;
    }

    // Displays information about a loaded Windows Bitmap image
    void info() {
        uint32_t depth;
        if (loaded) {
            std::cout << "[" << bminf.biWidth << " x " << bminf.biHeight << ", ";
            depth = bminf.biBitCount;
        }
        else {
            std::cout << "[" << w << " x " << h << ", ";
            depth = sizeof(struct rgba) * 8;
        }		
        std::cout << depth << " bbp (";
        switch(depth) {
        case  1: std::cout << "Monochrome"; break;
        case  4: std::cout << "16 colors";  break;
        case  8: std::cout << "256 colors"; break;
        case 24: std::cout << "16M colors"; break;
        case 32: std::cout << "True-color"; break;
        default: std::cout << "Unknown" ;   break;
        }
        std::cout << ")]" << std::endl;
    }
    
    // Renders the bitmap onto a Window via its handle
    void render(SDL_Surface * surface, bool stretch = false) {
        uint32_t color;
        uint32_t * pixels = (uint32_t *) surface->pixels;
        for (int y = 0; y < h; y++) {
            uint32_t offset = y * surface->w;
            for (int x = 0; x < w; x++) {
                memcpy(&pixels[offset + x], 
                        &data  [offset + x], sizeof(uint32_t));
            }
        }
    }

    int width() const { return w; }
    int height() const { return h; }

};    
    

class CGraph {
protected:
    SDL_Window* window;
    SDL_Surface* surface;
    uint32_t* pixels;
    uint32_t width;
    uint32_t height;
    uint32_t pitch;
    std::string title;
    bool clipping;
    bool initialized;
    bool fullscreen;

    COLORREF color;
	CPoint cursor;	// virtual cursor

public:
    CGraph() {
        this->clipping = false;
        this->initialized = false;
        this->title = "SDL Window";
        this->width = 640;
        this->height = 480;
        this->color = 0;
    }

    CGraph(int w, int h, std::string title = "SDL Window") {
        create(w, h, title);
    }

    ~CGraph() {
        SDL_DestroyWindow(window);
        SDL_Quit();    
    }

    bool create(int w, int h, std::string title, 
                bool fullscreen = false,
                int flags = SDL_WINDOW_SHOWN | 
                            SDL_WINDOW_OPENGL) {
        if (SDL_Init(SDL_INIT_VIDEO) < 0)
            return false;
        
        this->width = w;
        this->height = h;
        this->title = title;
        this->clipping = true;

        if (fullscreen)
            flags |= SDL_WINDOW_FULLSCREEN;

        window = SDL_CreateWindow(
            this->title.c_str(),
            SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
            this->width, this->height, flags
        );
    
        if (window == NULL) 
            return false;
               
        this->surface = SDL_GetWindowSurface(window);
        this->pixels = (uint32_t*) surface->pixels;
        this->pitch = surface->pitch >> 2;
        this->color = SDL_MapRGB(surface->format, 255, 255, 255);
        initialized = true;
        return true;
    }

    SDL_Window * getWindow() {
        return window;
    }

    SDL_Surface * getSurface() {
        return surface;
    }

    void clear(uint32_t color) {
        SDL_FillRect(surface, NULL, color);
    }

    void clear() {
        clear(this->color);
    }

    void clear(int r, int g, int b) {
        clear(MAKERGB(r, g, b));
    }

    void update() {
        SDL_UpdateWindowSurface(window);
    }

    void setClipping(bool clip) {
        this->clipping = clip;
    }

    bool isClipping() {
        return this->clipping;
    }

    void setFullscreen(bool fs) {
        this->fullscreen = fs;
        if (fs)
            SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
        else
            SDL_SetWindowFullscreen(window, 0);
    }

    bool isFullscreen() {
        return this->fullscreen;
    }

    void setTitle(std::string title) {
        this->title = title;
        SDL_SetWindowTitle(window, title.c_str());
    }

    void setWidth(uint32_t w) {
        this->width = w;
        SDL_SetWindowSize(window, w, height);    
    }

    void setHeight(uint32_t h) {
        this->height = h;
        SDL_SetWindowSize(window, width, h);    
    }

    void setFlags(uint32_t flags) {
        SDL_SetWindowFullscreen(window, flags);
    }

    uint32_t getFlags() {
        return SDL_GetWindowFlags(window);
    }

    uint32_t getWidth() {
        return this->width;
    }

    uint32_t getHeight() {
        return this->height;
    }

    uint32_t getColor() {
        return this->color;
    }

    uint32_t* getPixels() {
        return pixels;
    }

    CBitmap * getBitmap() {
        CBitmap * bmp = new CBitmap(width, height);
        memcpy(bmp->getPixels(), pixels, width * height * sizeof(uint32_t));
        return bmp;
    }

    virtual void render() { }

    virtual void loop() {
		SDL_Event e;
		int quit = 0;
	
		while (!quit) {
			while (SDL_PollEvent(&e) != 0) {
				if (e.type == SDL_QUIT) {
					quit = 1;
				}
				if (e.type == SDL_KEYDOWN) {
					switch (e.key.keysym.sym) {
						case SDLK_ESCAPE:
							quit = 1;
							break;
					}
				}
			}
			render();
			update();
		}    
    }

    // rendering methods

    void setColor(int r, int g, int b) {
        color = MAKERGB(r, g, b);
        SDL_MapRGB(surface->format, r, g, b);
    }

    void setColor(COLORREF color) {
        setColor(CRGB(color));
    }

    void setColor(CRGB color) {
        setColor(color.red(), color.green(), color.blue());
    }

    void plotPixel(int x, int y, COLORREF color) {
        if (clipping)
            if (!(x >= 0 && x < surface->w && y >= 0 && y < surface->h)) return;
        pixels[y * surface->w + x] = color;
    }
    
    void plotPixel(uint32_t x, uint32_t y) {
        plotPixel(x, y, color);
    }

    void lineHorz(uint32_t x1, uint32_t x2, uint32_t y, COLORREF color) {
        if (clipping) {
            if ((y < 0 | y >= surface->h)) return;
    
            if (x1 > x2) {
                int temp = x1;
                x1 = x2;
                x2 = temp;
            }
    
            if (x1 < 0) x1 = 0;
            if (x2 >= surface->w) x2 = surface->w -1;
        }
    
        for (int x = x1; x <= x2; x++) {
            pixels[y * surface->w + x] = color;
        }
    }

    void lineHorz(uint32_t x1, uint32_t x2, uint32_t y) {
        lineHorz(x1, x2, y, color);
    }

    void lineVert(uint32_t x, uint32_t y1, uint32_t y2, COLORREF color) {
        if (clipping) {
            if (x < 0 || x >= surface->w) return;
            if (y1 > y2) {
                int temp = y1;
                y1 = y2;
                y2 = temp;
            }
            if (y1 < 0) y1 = 0;
            if (y2 >= surface->h) y2 = surface->h -1;
        }
    
        for (int y = y1; y <= y2; y++) {
            pixels[y * surface->w + x] = color;
        }
    }

    void lineVert(uint32_t x, uint32_t y1, uint32_t y2) {
        lineVert(x, y1, y2, color);
    }

    void rectangle(uint32_t x, uint32_t y, uint32_t width, uint32_t height, COLORREF color) {
        if (clipping) {
            if (width <= 0 || height <= 0) return;
            if (x < 0) x = 0;
            if (y < 0) y = 0;
            if (x + width > surface->w) width = surface->w - x;
            if (y + height > surface->h) height = surface->h - y;
        }
        
        for (int i = 0; i < height; i++) {
            uint32_t* row = &pixels[(y + i) * pitch + x];
            for(int j = 0; j < width; j++)
                row[j] = color;
        }
    }
    
    void rectangle(int x, int y, int width, int height) {
        rectangle(x, y, width, height, color);
    }

    void line(int x1, int y1, int x2, int y2, COLORREF color) {
        int dx = abs(x2 - x1);
        int dy = abs(y2 - y1);
        int sx = (x1 < x2) ? 1 : -1;
        int sy = (y1 < y2) ? 1 : -1;
        int err = dx - dy;
    
        while (1) {
            plotPixel(x1, y1, color);
    
            if (x1 == x2 && y1 == y2) break;
    
            int e2 = 2 * err;
            if (e2 > -dy) {
                err -= dy;
                x1 += sx;
            }
            if (e2 < dx) {
                err += dx;
                y1 += sy;
            }
        }
    }

    void line(int x1, int y1, int x2, int y2) {
        line(x1, y1, x2, y2, color);
    }

	// Moves current virtual cursor to (x, y)
	void moveTo(int x, int y) {
		cursor.setPosition(x, y);
	}

	// Draws a line from current position to (x, y) and then move the virtual cursor to (x, y)
	void lineTo(int x, int y, COLORREF color) {
		line(cursor.getX(), cursor.getY(), x, y, color);
		moveTo(x, y);
	}

    void lineTo(int x, int y) {
		lineTo(x, y, color);
		moveTo(x, y);
	}

	// Draws a circle at (xc, yc) with radius r
	void drawCircle(int xc, int yc, int r, COLORREF color) {
		int x = 0, y = r;
		int d = 3 - (r << 1);
		circlePixels(xc, yc, x, y, color);
		while (y >= x) {
			x++;
			if (d > 0) {
				y--;
				d = d + ((x - y) << 2) + 10;
			}
			else
				d = d + (x << 2) + 6;
			circlePixels(xc, yc, x, y, color);
		}
	}

    void drawCircle(int xc, int yc, int r) {
        drawCircle(xc, yc, r, color);
    }

	// Draws a filled circle at (xc, yc) with radius r
	void filledCircle(int xc, int yc, int r, COLORREF color) {
		int x = 0, y = r;
		int d = 3 - (r << 1);
		circleLines(xc, yc, x, y, color);
		while (y >= x) {
			x++;
			if (d > 0) {
				y--;
				d = d + ((x - y) << 2) + 10;
			}
			else
				d = d + (x << 2) + 6;
			circleLines(xc, yc, x, y, color);
		}
	}

    void filledCircle(int xc, int yc, int r) {
        filledCircle(xc, yc, r, color);
    }

private:
	void circlePixels(uint32_t xc, uint32_t yc, uint32_t x, uint32_t y, COLORREF color)	{
		plotPixel(xc+x, yc+y, color);
		plotPixel(xc-x, yc+y, color);
		plotPixel(xc+x, yc-y, color);
		plotPixel(xc-x, yc-y, color);
		plotPixel(xc+y, yc+x, color);
		plotPixel(xc-y, yc+x, color);
		plotPixel(xc+y, yc-x, color);
		plotPixel(xc-y, yc-x, color);
	}

	void circleLines(uint32_t xc, uint32_t yc, uint32_t x, uint32_t y, COLORREF color)	{
		line(xc+x, yc+y, xc-x, yc+y, color);
		line(xc+x, yc-y, xc-x, yc-y, color);
		line(xc+y, yc+x, xc-y, yc+x, color);
		line(xc+y, yc-x, xc-y, yc-x, color);
	}
};

} // namespace daniel

#endif // __CGRAPH_H__


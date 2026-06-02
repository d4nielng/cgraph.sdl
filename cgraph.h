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
    - v0.3:
        Upgrade cgraph into a simple 2D engine with CTransform, CCollider,
        CGameObject, CScene, CTimer, and CAnimator support
        Add bitmap font loading from vga8x8.bmp using a 1-bit LUT for colored text
        Add bitmap viewer and game demos that use the engine rendering pipeline
	
    To do:
	- Keep extending the engine with reusable objects and scenes
	
    Usages:
    - Include this header file in your project
    - Linking with SDL2 library (using the parameter -lSDL2)
    - Enjoy coding with simple graphics
    
    Permissions:
    - Feel free to modify the source code to fit your needs

    Latest update: June 1, 2026
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
	
	CVector2D & operator -= (CVector2D v) {
		this->x -= v.x;
		this->y -= v.y;
		return * this;
	}
	
	CVector2D & operator -= (double d) {
		this->x -= d;
		this->y -= d;
		return * this;	
	}
	
	CVector2D & operator *= (double factor) {
		this->x *= factor;
		this->y *= factor;
		return * this;
	}
	
	CVector2D & operator /= (double divisor) {
		if (divisor != 0.0) {
			this->x /= divisor;
			this->y /= divisor;
		}
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
	
	double angle() const {
		return atan2(y, x);
	}
	
	double angleDegrees() const {
		return atan2(y, x) * 180.0 / M_PI;
	}
	
	CVector2D & setMagnitude(double newMag) {
		double mag = magnitude();
		if (mag != 0.0) {
			double factor = newMag / mag;
			this->x *= factor;
			this->y *= factor;
		}
		return * this;
	}
	
	CVector2D & clampMagnitude(double maxMag) {
		double mag = magnitude();
		if (mag > maxMag) {
			double factor = maxMag / mag;
			this->x *= factor;
			this->y *= factor;
		}
		return * this;
	}
	
	CVector2D perpendicular() const {
		return CVector2D(-this->y, this->x);
	}
	
	CVector2D & reflect(const CVector2D & normal) {
		double d = 2.0 * dot(normal);
		this->x -= d * normal.x;
		this->y -= d * normal.y;
		return * this;
	}
	
	CVector2D & rotate(double theta) {
		double ox = x;
		this->x = ox * cos(theta) - y * sin(theta);
		this->y = ox * sin(theta) + y * cos(theta);
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
	
	friend CVector2D operator * (const CVector2D & v, double scalar) {
		return CVector2D(v.x * scalar, v.y * scalar);
	}
	
	friend CVector2D operator * (double scalar, const CVector2D & v) {
		return CVector2D(v.x * scalar, v.y * scalar);
	}
	
	friend CVector2D operator / (const CVector2D & v, double divisor) {
		if (divisor == 0.0) return CVector2D(v);
		return CVector2D(v.x / divisor, v.y / divisor);
	}
	
	friend CVector2D operator - (const CVector2D & v) {
		return CVector2D(-v.x, -v.y);
	}
	
	friend bool operator == (const CVector2D & u, const CVector2D & v) {
		return u.magnitude() == v.magnitude();
	}

	friend bool operator != (const CVector2D & u, const CVector2D & v) {
		return u.magnitude() != v.magnitude();
	}
	
	static CVector2D lerp(const CVector2D & u, const CVector2D & v, double t) {
		return CVector2D(u.x + (v.x - u.x) * t, u.y + (v.y - u.y) * t);
	}
	
	static double angle(const CVector2D & u, const CVector2D & v) {
		double d = u.dot(v);
		double mag = u.magnitude() * v.magnitude();
		if (mag == 0.0) return 0.0;
		return acos(CMath::clampf(d / mag, -1.0, 1.0));
	}
	
	static double angleDegrees(const CVector2D & u, const CVector2D & v) {
		return angle(u, v) * 180.0 / M_PI;
	}
	
	static CVector2D & zero() {
		static CVector2D z(0.0f, 0.0f);
		return z;
	}
};

typedef class CVector2D Vector2D;

// CVector3D
// ---------
// Represents a vector in 3D space
class CVector3D {
    union {
        struct {
            double x, y, z;
        };
        double v[3];
    };

public:

    CVector3D() : x(0.0f), y(0.0f), z(0.0f) { }
    CVector3D(double x, double y, double z) { this->x = x; this->y = y; this->z = z; }
    CVector3D(const CVector3D & rhs) { *this = rhs; }

    double getX() const { return this->x; }
    double getY() const { return this->y; }
    double getZ() const { return this->z; }

    CVector3D & setX(double x) {
        this->x = x;
        return * this;
    }

    CVector3D & setY(double y) {
        this->y = y;
        return * this;
    }

    CVector3D & setZ(double z) {
        this->z = z;
        return * this;
    }

    CVector3D & set(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
        return * this;
    }

    CVector3D & incX(double dx) {
        this->x += dx;
        return * this;
    }

    CVector3D & incY(double dy) {
        this->y += dy;
        return * this;
    }

    CVector3D & incZ(double dz) {
        this->z += dz;
        return * this;
    }

    double magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }

    CVector3D & normalize() {
        double mag = magnitude();
        if (mag != 0.0) {
            x /= mag;
            y /= mag;
            z /= mag;
        }
        return * this;
    }

    double dot(const CVector3D & v) const {
        return this->x * v.x + this->y * v.y + this->z * v.z;
    }

    CVector3D cross(const CVector3D & v) const {
        return CVector3D(
            this->y * v.z - this->z * v.y,
            this->z * v.x - this->x * v.z,
            this->x * v.y - this->y * v.x
        );
    }

    double distance(const CVector3D & v) const {
        double dx = v.x - this->x;
        double dy = v.y - this->y;
        double dz = v.z - this->z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    // Perspective projection from 3D to 2D screen space.
    // focalLength controls zoom, centerX/centerY define the screen origin,
    // and cameraZ is the camera position on the Z axis.
    CVector2D projectTo2D(double focalLength, double centerX = 0.0, double centerY = 0.0, double cameraZ = 0.0) const {
        double dz = this->z - cameraZ;
        if (dz == 0.0) dz = 0.000001;
        double scale = focalLength / dz;
        return CVector2D(centerX + this->x * scale, centerY - this->y * scale);
    }

    double operator[] (int index) const {
        return v[index];
    }

    CVector3D & operator = (const CVector3D & rhs) {
        if (&rhs != this) {
            this->x = rhs.x;
            this->y = rhs.y;
            this->z = rhs.z;
        }
        return * this;
    }

    CVector3D & operator += (const CVector3D & v) {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return * this;
    }

    CVector3D & operator += (double d) {
        this->x += d;
        this->y += d;
        this->z += d;
        return * this;
    }

    CVector3D & operator -= (const CVector3D & v) {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return * this;
    }

    CVector3D & operator -= (double d) {
        this->x -= d;
        this->y -= d;
        this->z -= d;
        return * this;
    }

    CVector3D & operator *= (double factor) {
        this->x *= factor;
        this->y *= factor;
        this->z *= factor;
        return * this;
    }

    CVector3D & operator /= (double divisor) {
        if (divisor != 0.0) {
            this->x /= divisor;
            this->y /= divisor;
            this->z /= divisor;
        }
        return * this;
    }

    friend std::ostream & operator << (std::ostream & os, const CVector3D & v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }

    friend std::istream & operator >> (std::istream & is, CVector3D & v) {
        is >> v.x >> v.y >> v.z;
        return is;
    }

    friend double distance(const CVector3D & u, const CVector3D & v) {
        double dx = v.x - u.x;
        double dy = v.y - u.y;
        double dz = v.z - u.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    friend double dot(const CVector3D & u, const CVector3D & v) {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }

    friend CVector3D cross(const CVector3D & u, const CVector3D & v) {
        return u.cross(v);
    }

    friend CVector3D operator + (const CVector3D & u, const CVector3D & v) {
        return CVector3D(u.x + v.x, u.y + v.y, u.z + v.z);
    }

    friend CVector3D operator - (const CVector3D & u, const CVector3D & v) {
        return CVector3D(u.x - v.x, u.y - v.y, u.z - v.z);
    }

    friend CVector3D operator * (const CVector3D & v, double scalar) {
        return CVector3D(v.x * scalar, v.y * scalar, v.z * scalar);
    }

    friend CVector3D operator * (double scalar, const CVector3D & v) {
        return CVector3D(v.x * scalar, v.y * scalar, v.z * scalar);
    }

    friend CVector3D operator / (const CVector3D & v, double divisor) {
        if (divisor == 0.0) return CVector3D(v);
        return CVector3D(v.x / divisor, v.y / divisor, v.z / divisor);
    }

    friend CVector3D operator - (const CVector3D & v) {
        return CVector3D(-v.x, -v.y, -v.z);
    }

    friend bool operator == (const CVector3D & u, const CVector3D & v) {
        return u.magnitude() == v.magnitude();
    }

    friend bool operator != (const CVector3D & u, const CVector3D & v) {
        return u.magnitude() != v.magnitude();
    }

    static CVector3D lerp(const CVector3D & u, const CVector3D & v, double t) {
        return CVector3D(
            u.x + (v.x - u.x) * t,
            u.y + (v.y - u.y) * t,
            u.z + (v.z - u.z) * t
        );
    }

    static double angle(const CVector3D & u, const CVector3D & v) {
        double d = u.dot(v);
        double mag = u.magnitude() * v.magnitude();
        if (mag == 0.0) return 0.0;
        return acos(CMath::clampf(d / mag, -1.0, 1.0));
    }

    static double angleDegrees(const CVector3D & u, const CVector3D & v) {
        return angle(u, v) * 180.0 / M_PI;
    }

    static CVector3D & zero() {
        static CVector3D z(0.0f, 0.0f, 0.0f);
        return z;
    }
};

typedef CVector3D Vector3D;

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
        while (this->h >= 360.0) this->h -= 360.0;
        while (this->h < 0.0)   this->h += 360.0;
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
        data.resize(w * h);
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
        return true;
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
    

static const uint8_t CGRAPH_DEFAULT_FONT_LUT[256][8] = {
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 0
    {0x81, 0x7E, 0x5A, 0x7E, 0x42, 0x66, 0x7E, 0x81}, // 1
    {0x81, 0x24, 0x00, 0x00, 0x3C, 0x18, 0x81, 0xFF}, // 2
    {0xC9, 0x80, 0x80, 0x80, 0xC1, 0xE3, 0xF7, 0xFF}, // 3
    {0xF7, 0xE3, 0xC1, 0x80, 0xC1, 0xE3, 0xF7, 0xFF}, // 4
    {0xE7, 0xC3, 0xE7, 0x99, 0x00, 0x99, 0xE7, 0xC3}, // 5
    {0xF7, 0xE3, 0xC1, 0x80, 0x80, 0xC1, 0xF7, 0xE3}, // 6
    {0xFF, 0xFF, 0xE7, 0xC3, 0xC3, 0xE7, 0xFF, 0xFF}, // 7
    {0x00, 0x00, 0x18, 0x3C, 0x3C, 0x18, 0x00, 0x00}, // 8
    {0xFF, 0xC3, 0x99, 0xBD, 0xBD, 0x99, 0xC3, 0xFF}, // 9
    {0x00, 0x3C, 0x66, 0x42, 0x42, 0x66, 0x3C, 0x00}, // 10
    {0xC0, 0xF2, 0xE3, 0xC1, 0x9C, 0x9C, 0xC1, 0xFF}, // 11
    {0xC1, 0x9C, 0x9C, 0xC1, 0xE3, 0x80, 0xE3, 0xFF}, // 12
    {0xF1, 0xF0, 0xF2, 0xF2, 0xE3, 0x80, 0xE3, 0xFF}, // 13
    {0xF0, 0xC4, 0xC8, 0xC4, 0xCC, 0xC8, 0x88, 0x8F}, // 14
    {0xE7, 0x24, 0xC3, 0x18, 0xC3, 0x24, 0xE7, 0xFF}, // 15
    {0x9F, 0x87, 0x81, 0x80, 0x81, 0x87, 0x9F, 0xFF}, // 16
    {0xFC, 0xF0, 0xC0, 0x80, 0xC0, 0xF0, 0xFC, 0xFF}, // 17
    {0xE7, 0xC3, 0x81, 0xE7, 0xE7, 0x81, 0xC3, 0xE7}, // 18
    {0x99, 0x99, 0x99, 0x99, 0x99, 0xFF, 0x99, 0xFF}, // 19
    {0xC0, 0x92, 0x92, 0xC2, 0xF2, 0xF2, 0xF2, 0xFF}, // 20
    {0xC0, 0x8F, 0xC1, 0x9C, 0x9C, 0xC1, 0xF8, 0x81}, // 21
    {0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0xFF}, // 22
    {0xC3, 0x81, 0xE7, 0xE7, 0x81, 0xC3, 0xE7, 0x81}, // 23
    {0xE7, 0xC3, 0x81, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF}, // 24
    {0xE7, 0xE7, 0xE7, 0xE7, 0x81, 0xC3, 0xE7, 0xFF}, // 25
    {0xE7, 0xF3, 0xF9, 0x80, 0xF9, 0xF3, 0xE7, 0xFF}, // 26
    {0xF3, 0xE7, 0xCF, 0x80, 0xCF, 0xE7, 0xF3, 0xFF}, // 27
    {0xFF, 0xFF, 0xFF, 0x9F, 0x9F, 0x80, 0xFF, 0xFF}, // 28
    {0xFF, 0xDB, 0x99, 0x00, 0x99, 0xDB, 0xFF, 0xFF}, // 29
    {0xF7, 0xE3, 0xE3, 0xC1, 0xC1, 0x80, 0x80, 0xFF}, // 30
    {0x80, 0x80, 0xC1, 0xC1, 0xE3, 0xE3, 0xF7, 0xFF}, // 31
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 32
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF, 0xE7, 0xFF}, // 33
    {0xCC, 0x99, 0x33, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 34
    {0xC9, 0xC9, 0x80, 0xC9, 0xC9, 0x80, 0xC9, 0xC9}, // 35
    {0xE7, 0x81, 0xA7, 0xC7, 0xE3, 0xE5, 0x81, 0xE7}, // 36
    {0x1C, 0x59, 0x13, 0xE7, 0xC8, 0x9A, 0x38, 0xFF}, // 37
    {0xC3, 0x99, 0xC3, 0xC7, 0x92, 0x99, 0xC2, 0xFF}, // 38
    {0xF9, 0xF3, 0xE7, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 39
    {0xE3, 0xCF, 0x9F, 0x9F, 0x9F, 0xCF, 0xE3, 0xFF}, // 40
    {0xC7, 0xF3, 0xF9, 0xF9, 0xF9, 0xF3, 0xC7, 0xFF}, // 41
    {0x9C, 0xC9, 0xE3, 0x80, 0xE3, 0xC9, 0x9C, 0xFF}, // 42
    {0xE7, 0xE7, 0xE7, 0x81, 0xE7, 0xE7, 0xE7, 0xFF}, // 43
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xE7, 0xE7, 0xCF}, // 44
    {0xFF, 0xFF, 0xFF, 0x81, 0xFF, 0xFF, 0xFF, 0xFF}, // 45
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xE7, 0xE7, 0xFF}, // 46
    {0xFD, 0xF9, 0xF3, 0xE7, 0xCF, 0x9F, 0xBF, 0xFF}, // 47
    {0xC1, 0x9C, 0x9C, 0x94, 0x9C, 0x9C, 0xC1, 0xFF}, // 48
    {0xE3, 0xC3, 0xF3, 0xF3, 0xF3, 0xF3, 0xF3, 0xFF}, // 49
    {0xC1, 0x9C, 0xFC, 0xF1, 0xC7, 0x9F, 0x80, 0xFF}, // 50
    {0x83, 0xF9, 0xF9, 0xC3, 0xF9, 0xF9, 0x83, 0xFF}, // 51
    {0xF1, 0xE1, 0xC9, 0x99, 0x80, 0xF9, 0xF9, 0xFF}, // 52
    {0x81, 0x9F, 0x83, 0xF9, 0xF9, 0xF9, 0x83, 0xFF}, // 53
    {0xE1, 0xCF, 0x9F, 0x81, 0x9C, 0x9C, 0xC1, 0xFF}, // 54
    {0x81, 0x99, 0xF3, 0xE7, 0xCF, 0xCF, 0xCF, 0xFF}, // 55
    {0xC1, 0x9C, 0xC9, 0xE3, 0xC9, 0x9C, 0xC1, 0xFF}, // 56
    {0xC1, 0x9C, 0x9C, 0xC0, 0xFC, 0xF9, 0x83, 0xFF}, // 57
    {0xFF, 0xFF, 0xE7, 0xE7, 0xFF, 0xE7, 0xE7, 0xFF}, // 58
    {0xFF, 0xFF, 0xE7, 0xE7, 0xFF, 0xE7, 0xE7, 0xCF}, // 59
    {0xF3, 0xE7, 0xCF, 0x9F, 0xCF, 0xE7, 0xF3, 0xFF}, // 60
    {0xFF, 0xFF, 0x81, 0xFF, 0x81, 0xFF, 0xFF, 0xFF}, // 61
    {0xCF, 0xE7, 0xF3, 0xF9, 0xF3, 0xE7, 0xCF, 0xFF}, // 62
    {0xC1, 0x9C, 0xF9, 0xF3, 0xF3, 0xFF, 0xF3, 0xFF}, // 63
    {0xC1, 0x9C, 0x90, 0x96, 0x90, 0x9F, 0xC0, 0xFF}, // 64
    {0xC3, 0x99, 0x99, 0x81, 0x99, 0x99, 0x99, 0xFF}, // 65
    {0x81, 0x9C, 0x9C, 0x81, 0x9C, 0x9C, 0x81, 0xFF}, // 66
    {0xE1, 0xCC, 0x9F, 0x9F, 0x9F, 0xCC, 0xE1, 0xFF}, // 67
    {0x83, 0x99, 0x9C, 0x9C, 0x9C, 0x99, 0x83, 0xFF}, // 68
    {0x81, 0x9F, 0x9F, 0x83, 0x9F, 0x9F, 0x81, 0xFF}, // 69
    {0x81, 0x9F, 0x9F, 0x83, 0x9F, 0x9F, 0x9F, 0xFF}, // 70
    {0xE1, 0xCC, 0x9C, 0x9F, 0x98, 0xCC, 0xE1, 0xFF}, // 71
    {0x9C, 0x9C, 0x9C, 0x80, 0x9C, 0x9C, 0x9C, 0xFF}, // 72
    {0xC3, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xC3, 0xFF}, // 73
    {0xF9, 0xF9, 0xF9, 0xF9, 0x99, 0x99, 0xC3, 0xFF}, // 74
    {0x9C, 0x99, 0x93, 0x87, 0x93, 0x99, 0x9C, 0xFF}, // 75
    {0xCF, 0xCF, 0xCF, 0xCF, 0xCF, 0xCF, 0xC0, 0xFF}, // 76
    {0x9C, 0x88, 0x80, 0x94, 0x9C, 0x9C, 0x9C, 0xFF}, // 77
    {0x9C, 0x8C, 0x84, 0x90, 0x98, 0x9C, 0x9C, 0xFF}, // 78
    {0xC1, 0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0xC1, 0xFF}, // 79
    {0x81, 0x9C, 0x9C, 0x81, 0x9F, 0x9F, 0x9F, 0xFF}, // 80
    {0xC1, 0x9C, 0x9C, 0x9C, 0x84, 0x91, 0xC3, 0xF8}, // 81
    {0x81, 0x9C, 0x9C, 0x81, 0x93, 0x99, 0x9C, 0xFF}, // 82
    {0xC1, 0x9C, 0xCF, 0xE3, 0xF9, 0x9C, 0xC1, 0xFF}, // 83
    {0x81, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF}, // 84
    {0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0xC1, 0xFF}, // 85
    {0x9C, 0x9C, 0x9C, 0x9C, 0xC9, 0xE3, 0xF7, 0xFF}, // 86
    {0x9C, 0x9C, 0x9C, 0x94, 0x94, 0x80, 0xC9, 0xFF}, // 87
    {0x9C, 0x9C, 0xC9, 0xE3, 0xC9, 0x9C, 0x9C, 0xFF}, // 88
    {0x9C, 0x9C, 0x9C, 0xC1, 0xF3, 0xF3, 0xF3, 0xFF}, // 89
    {0x80, 0xF9, 0xF3, 0xE7, 0xCF, 0x9F, 0x80, 0xFF}, // 90
    {0xC1, 0xCF, 0xCF, 0xCF, 0xCF, 0xCF, 0xC1, 0xFF}, // 91
    {0xBF, 0x9F, 0xCF, 0xE7, 0xF3, 0xF9, 0xFD, 0xFF}, // 92
    {0xC1, 0xF9, 0xF9, 0xF9, 0xF9, 0xF9, 0xC1, 0xFF}, // 93
    {0xE3, 0xC9, 0x9C, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 94
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x00}, // 95
    {0xCF, 0xE7, 0xF3, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 96
    {0xFF, 0xFF, 0xC3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 97
    {0x9F, 0x9F, 0x83, 0x99, 0x99, 0x99, 0x83, 0xFF}, // 98
    {0xFF, 0xFF, 0xC1, 0x9C, 0x9F, 0x9F, 0xC0, 0xFF}, // 99
    {0xF9, 0xF9, 0xC1, 0x99, 0x99, 0x99, 0xC1, 0xFF}, // 100
    {0xFF, 0xFF, 0xC3, 0x99, 0x83, 0x9F, 0xC1, 0xFF}, // 101
    {0xE1, 0xCF, 0xCF, 0x83, 0xCF, 0xCF, 0xCF, 0xFF}, // 102
    {0xFF, 0xFF, 0xC0, 0x9C, 0x9C, 0xC0, 0xFC, 0x81}, // 103
    {0x9F, 0x9F, 0x93, 0x89, 0x99, 0x99, 0x99, 0xFF}, // 104
    {0xE7, 0xFF, 0xC7, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF}, // 105
    {0xF9, 0xFF, 0xF9, 0xF9, 0xF9, 0xF9, 0x99, 0xC3}, // 106
    {0x9F, 0x9F, 0x99, 0x93, 0x87, 0x93, 0x99, 0xFF}, // 107
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xF3, 0xFF}, // 108
    {0xFF, 0xFF, 0x89, 0x80, 0x94, 0x94, 0x9C, 0xFF}, // 109
    {0xFF, 0xFF, 0x93, 0x89, 0x99, 0x99, 0x99, 0xFF}, // 110
    {0xFF, 0xFF, 0xC1, 0x9C, 0x9C, 0x9C, 0xC1, 0xFF}, // 111
    {0xFF, 0xFF, 0x83, 0x99, 0x99, 0x83, 0x9F, 0x9F}, // 112
    {0xFF, 0xFF, 0xC1, 0x99, 0x99, 0xC1, 0xF9, 0xF9}, // 113
    {0xFF, 0xFF, 0xC9, 0xC4, 0xCF, 0xCF, 0xCF, 0xFF}, // 114
    {0xFF, 0xFF, 0xC1, 0x8F, 0xC3, 0xF1, 0x83, 0xFF}, // 115
    {0xE7, 0xE7, 0x81, 0xE7, 0xE7, 0xE7, 0xF1, 0xFF}, // 116
    {0xFF, 0xFF, 0x99, 0x99, 0x99, 0x99, 0xC4, 0xFF}, // 117
    {0xFF, 0xFF, 0x99, 0x99, 0x99, 0xC3, 0xE7, 0xFF}, // 118
    {0xFF, 0xFF, 0x9C, 0x9C, 0x94, 0x80, 0xC9, 0xFF}, // 119
    {0xFF, 0xFF, 0x9C, 0xC9, 0xE3, 0xC9, 0x9C, 0xFF}, // 120
    {0xFF, 0xFF, 0x99, 0x99, 0x99, 0xC1, 0xF9, 0x83}, // 121
    {0xFF, 0xFF, 0x81, 0xF3, 0xE7, 0xCF, 0x81, 0xFF}, // 122
    {0xF1, 0xE7, 0xE7, 0x87, 0xE7, 0xE7, 0xF1, 0xFF}, // 123
    {0xE7, 0xE7, 0xE7, 0xFF, 0xE7, 0xE7, 0xE7, 0xFF}, // 124
    {0x8F, 0xE7, 0xE7, 0xF1, 0xE7, 0xE7, 0x8F, 0xFF}, // 125
    {0xC4, 0x91, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}, // 126
    {0xE7, 0xC3, 0x99, 0x3C, 0x3C, 0x00, 0xFF, 0xFF}, // 127
    {0xC1, 0x9C, 0x9F, 0x9F, 0x9C, 0xC1, 0xE7, 0xCF}, // 128
    {0x99, 0xFF, 0x99, 0x99, 0x99, 0x99, 0xC4, 0xFF}, // 129
    {0xF3, 0xEF, 0xC3, 0x99, 0x83, 0x9F, 0xC1, 0xFF}, // 130
    {0xC1, 0xBE, 0xE3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 131
    {0xC9, 0xFF, 0xC3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 132
    {0xE7, 0xF7, 0xC3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 133
    {0xC3, 0x99, 0xC3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 134
    {0xFF, 0xFF, 0xC1, 0x9F, 0x9F, 0xC1, 0xE7, 0xCF}, // 135
    {0x81, 0x3C, 0xC3, 0x99, 0x83, 0x9F, 0xC1, 0xFF}, // 136
    {0xC9, 0xFF, 0xC3, 0x99, 0x83, 0x9F, 0xC1, 0xFF}, // 137
    {0xE7, 0xF7, 0xC3, 0x99, 0x83, 0x9F, 0xC1, 0xFF}, // 138
    {0x93, 0xFF, 0xC7, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF}, // 139
    {0xC1, 0xBE, 0xE3, 0xF3, 0xF3, 0xF3, 0xF3, 0xFF}, // 140
    {0xF3, 0xFB, 0xC7, 0xE7, 0xE7, 0xE7, 0xE7, 0xFF}, // 141
    {0x99, 0xE7, 0xC3, 0x99, 0x81, 0x99, 0x99, 0xFF}, // 142
    {0xE7, 0xDB, 0xC3, 0x99, 0x81, 0x99, 0x99, 0xFF}, // 143
    {0xF3, 0x81, 0x9F, 0x83, 0x9F, 0x9F, 0x81, 0xFF}, // 144
    {0xFF, 0xFF, 0x81, 0xE6, 0xC0, 0xB3, 0x88, 0xFF}, // 145
    {0xE0, 0xC9, 0x99, 0x80, 0x99, 0x99, 0x98, 0xFF}, // 146
    {0x81, 0x7E, 0xC3, 0x99, 0x99, 0x99, 0xC3, 0xFF}, // 147
    {0x99, 0xFF, 0xC3, 0x99, 0x99, 0x99, 0xC3, 0xFF}, // 148
    {0xE7, 0xF7, 0xC3, 0x99, 0x99, 0x99, 0xC3, 0xFF}, // 149
    {0x81, 0x7E, 0x99, 0x99, 0x99, 0x99, 0xC4, 0xFF}, // 150
    {0xE7, 0xF7, 0x99, 0x99, 0x99, 0x99, 0xC4, 0xFF}, // 151
    {0xC9, 0xFF, 0x9C, 0x9C, 0x9C, 0xC0, 0xFC, 0x81}, // 152
    {0x9C, 0xFF, 0xC1, 0x9C, 0x9C, 0x9C, 0x9C, 0xC1}, // 153
    {0x9C, 0xFF, 0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0xC1}, // 154
    {0xF7, 0xF7, 0xC1, 0x9F, 0x9F, 0xC1, 0xF7, 0xF7}, // 155
    {0xF0, 0xE6, 0xE7, 0x81, 0xE7, 0xE7, 0x80, 0xFF}, // 156
    {0x99, 0x99, 0xC3, 0x81, 0xE7, 0x81, 0xE7, 0xE7}, // 157
    {0x83, 0x99, 0x99, 0x83, 0x9B, 0x90, 0x9B, 0x98}, // 158
    {0xF8, 0xF2, 0xF3, 0xC0, 0xF3, 0xF3, 0xB3, 0xC7}, // 159
    {0xF3, 0xEF, 0xC3, 0xF9, 0xC1, 0x99, 0xC0, 0xFF}, // 160
    {0xE7, 0xDF, 0xE3, 0xF3, 0xF3, 0xF3, 0xF3, 0xFF}, // 161
    {0xF3, 0xEF, 0xC1, 0x9C, 0x9C, 0x9C, 0xC1, 0xFF}, // 162
    {0xF3, 0xEF, 0x99, 0x99, 0x99, 0x99, 0xC4, 0xFF}, // 163
    {0x8E, 0x71, 0x93, 0x89, 0x99, 0x99, 0x99, 0xFF}, // 164
    {0x8C, 0x73, 0x9C, 0x8C, 0x94, 0x98, 0x9C, 0xFF}, // 165
    {0xC1, 0x99, 0xC0, 0xFF, 0x80, 0xFF, 0xFF, 0xFF}, // 166
    {0xC3, 0x99, 0xC3, 0xFF, 0x81, 0xFF, 0xFF, 0xFF}, // 167
    {0xE7, 0xFF, 0xE7, 0xE7, 0xCF, 0x9C, 0x9C, 0xC1}, // 168
    {0xFF, 0xFF, 0x80, 0x9F, 0x9F, 0x9F, 0x9F, 0xFF}, // 169
    {0xFF, 0xFF, 0x80, 0xFC, 0xFC, 0xFC, 0xFC, 0xFF}, // 170
    {0x9C, 0x99, 0x93, 0xE0, 0xCE, 0x98, 0xB3, 0xB0}, // 171
    {0x9C, 0x99, 0x93, 0xE4, 0xCA, 0x96, 0xB0, 0xBE}, // 172
    {0xE7, 0xFF, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7}, // 173
    {0xFF, 0xCC, 0x99, 0x33, 0x99, 0xCC, 0xFF, 0xFF}, // 174
    {0xFF, 0x33, 0x99, 0xCC, 0x99, 0x33, 0xFF, 0xFF}, // 175
    {0xEE, 0xBB, 0xEE, 0xBB, 0xEE, 0xBB, 0xEE, 0xBB}, // 176
    {0xAA, 0x55, 0xAA, 0x55, 0xAA, 0x55, 0xAA, 0x55}, // 177
    {0x22, 0x88, 0x22, 0x88, 0x22, 0x88, 0x22, 0x88}, // 178
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7}, // 179
    {0xE7, 0xE7, 0xE7, 0xE7, 0x07, 0xE7, 0xE7, 0xE7}, // 180
    {0xE7, 0xE7, 0x07, 0xE7, 0x07, 0xE7, 0xE7, 0xE7}, // 181
    {0xC9, 0xC9, 0xC9, 0xC9, 0x09, 0xC9, 0xC9, 0xC9}, // 182
    {0xFF, 0xFF, 0xFF, 0xFF, 0x01, 0xC9, 0xC9, 0xC9}, // 183
    {0xFF, 0xFF, 0x07, 0xE7, 0x07, 0xE7, 0xE7, 0xE7}, // 184
    {0xC9, 0xC9, 0x09, 0xF9, 0x09, 0xC9, 0xC9, 0xC9}, // 185
    {0xC9, 0xC9, 0xC9, 0xC9, 0xC9, 0xC9, 0xC9, 0xC9}, // 186
    {0xFF, 0xFF, 0x01, 0xF9, 0x09, 0xC9, 0xC9, 0xC9}, // 187
    {0xC9, 0xC9, 0x09, 0xF9, 0x01, 0xFF, 0xFF, 0xFF}, // 188
    {0xC9, 0xC9, 0xC9, 0xC9, 0x01, 0xFF, 0xFF, 0xFF}, // 189
    {0xE7, 0xE7, 0x07, 0xE7, 0x07, 0xFF, 0xFF, 0xFF}, // 190
    {0xFF, 0xFF, 0xFF, 0xFF, 0x07, 0xE7, 0xE7, 0xE7}, // 191
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE0, 0xFF, 0xFF, 0xFF}, // 192
    {0xE7, 0xE7, 0xE7, 0xE7, 0x00, 0xFF, 0xFF, 0xFF}, // 193
    {0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xE7, 0xE7, 0xE7}, // 194
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE0, 0xE7, 0xE7, 0xE7}, // 195
    {0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0xFF, 0xFF}, // 196
    {0xE7, 0xE7, 0xE7, 0xE7, 0x00, 0xE7, 0xE7, 0xE7}, // 197
    {0xE7, 0xE7, 0xE0, 0xE7, 0xE0, 0xE7, 0xE7, 0xE7}, // 198
    {0xC9, 0xC9, 0xC9, 0xC9, 0xC8, 0xC9, 0xC9, 0xC9}, // 199
    {0xC9, 0xC9, 0xC8, 0xCF, 0xC0, 0xFF, 0xFF, 0xFF}, // 200
    {0xFF, 0xFF, 0xC0, 0xCF, 0xC8, 0xC9, 0xC9, 0xC9}, // 201
    {0xC9, 0xC9, 0x08, 0xFF, 0x00, 0xFF, 0xFF, 0xFF}, // 202
    {0xFF, 0xFF, 0x00, 0xFF, 0x08, 0xC9, 0xC9, 0xC9}, // 203
    {0xC9, 0xC9, 0xC8, 0xCF, 0xC8, 0xC9, 0xC9, 0xC9}, // 204
    {0xFF, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF}, // 205
    {0xC9, 0xC9, 0x08, 0xFF, 0x08, 0xC9, 0xC9, 0xC9}, // 206
    {0xE7, 0xE7, 0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF}, // 207
    {0xC9, 0xC9, 0xC9, 0xC9, 0x00, 0xFF, 0xFF, 0xFF}, // 208
    {0xFF, 0xFF, 0x00, 0xFF, 0x00, 0xE7, 0xE7, 0xE7}, // 209
    {0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xC9, 0xC9, 0xC9}, // 210
    {0xC9, 0xC9, 0xC9, 0xC9, 0xC0, 0xFF, 0xFF, 0xFF}, // 211
    {0xE7, 0xE7, 0xE0, 0xE7, 0xE0, 0xFF, 0xFF, 0xFF}, // 212
    {0xFF, 0xFF, 0xE0, 0xE7, 0xE0, 0xE7, 0xE7, 0xE7}, // 213
    {0xFF, 0xFF, 0xFF, 0xFF, 0xC0, 0xC9, 0xC9, 0xC9}, // 214
    {0xC9, 0xC9, 0xC9, 0xC9, 0x00, 0xC9, 0xC9, 0xC9}, // 215
    {0xE7, 0xE7, 0x00, 0xE7, 0x00, 0xE7, 0xE7, 0xE7}, // 216
    {0xE7, 0xE7, 0xE7, 0xE7, 0x07, 0xFF, 0xFF, 0xFF}, // 217
    {0xFF, 0xFF, 0xFF, 0xFF, 0xE0, 0xE7, 0xE7, 0xE7}, // 218
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // 219
    {0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x00, 0x00}, // 220
    {0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F}, // 221
    {0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0}, // 222
    {0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF}, // 223
    {0xFF, 0xFF, 0x84, 0x31, 0x33, 0x31, 0x84, 0xFF}, // 224
    {0xC3, 0x99, 0x99, 0x91, 0x9C, 0x9C, 0x91, 0xFF}, // 225
    {0xC0, 0xCC, 0xCF, 0xCF, 0xCF, 0xCF, 0xCF, 0xFF}, // 226
    {0xFF, 0xFF, 0x80, 0x49, 0xC9, 0xC9, 0xC9, 0xFF}, // 227
    {0x80, 0x8E, 0xC7, 0xE3, 0xC7, 0x8E, 0x80, 0xFF}, // 228
    {0xFF, 0xFF, 0xC0, 0x99, 0x99, 0x99, 0xC3, 0xFF}, // 229
    {0xFF, 0xFF, 0x99, 0x99, 0x99, 0x84, 0x9F, 0x9F}, // 230
    {0xFF, 0xFF, 0xC0, 0xB3, 0xF3, 0xF3, 0xF3, 0xFF}, // 231
    {0xE7, 0x81, 0x24, 0x24, 0x24, 0x81, 0xE7, 0xFF}, // 232
    {0xC3, 0x99, 0x3C, 0x42, 0x3C, 0x99, 0xC3, 0xFF}, // 233
    {0xC3, 0x99, 0x3C, 0x3C, 0x3C, 0x99, 0x18, 0xFF}, // 234
    {0xC1, 0x8F, 0xC3, 0x99, 0x99, 0x99, 0xC3, 0xFF}, // 235
    {0xFF, 0xFF, 0x81, 0x24, 0x24, 0x81, 0xFF, 0xFF}, // 236
    {0xE7, 0xE7, 0x81, 0x24, 0x24, 0x81, 0xE7, 0xE7}, // 237
    {0xFF, 0xFF, 0xC1, 0x9F, 0xC3, 0x9F, 0xC1, 0xFF}, // 238
    {0xC1, 0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0x9C, 0xFF}, // 239
    {0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF, 0x80, 0xFF}, // 240
    {0xE7, 0xE7, 0x81, 0xE7, 0xE7, 0xFF, 0x81, 0xFF}, // 241
    {0xC7, 0xE3, 0xF1, 0xE3, 0xC7, 0xFF, 0x81, 0xFF}, // 242
    {0xE3, 0xC7, 0x8F, 0xC7, 0xE3, 0xFF, 0x81, 0xFF}, // 243
    {0xF1, 0xE4, 0xE4, 0xE7, 0xE7, 0xE7, 0xE7, 0xE7}, // 244
    {0xE7, 0xE7, 0xE7, 0xE7, 0xE7, 0x27, 0x27, 0x8F}, // 245
    {0xFF, 0xE7, 0xFF, 0x81, 0xFF, 0xE7, 0xFF, 0xFF}, // 246
    {0xFF, 0xFF, 0xC4, 0x91, 0xFF, 0xC4, 0x91, 0xFF}, // 247
    {0xFF, 0xC3, 0x99, 0xC3, 0xFF, 0xFF, 0xFF, 0xFF}, // 248
    {0xFF, 0xFF, 0xE7, 0xC3, 0xE7, 0xFF, 0xFF, 0xFF}, // 249
    {0xFF, 0xFF, 0xFF, 0xE7, 0xFF, 0xFF, 0xFF, 0xFF}, // 250
    {0xF0, 0xF3, 0xF3, 0x73, 0x33, 0x93, 0xC3, 0xE7}, // 251
    {0x93, 0xC9, 0xC9, 0xC9, 0xFF, 0xFF, 0xFF, 0xFF}, // 252
    {0x83, 0xF1, 0xC3, 0x9F, 0x81, 0xFF, 0xFF, 0xFF}, // 253
    {0xFF, 0xFF, 0x83, 0x83, 0x83, 0x83, 0x83, 0x83}, // 254
    {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF} // 255
};
class CGraph {
protected:
    enum RenderBackend {
        BackendSoftware,
        BackendHardwarePrimitives,
        BackendHardwarePixels
    };

    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* frameTexture;
    SDL_Surface* windowSurface;
    SDL_Surface* surface;
    uint32_t* pixels;
    uint32_t width;
    uint32_t height;
    uint32_t pitch;
    std::string title;
    bool clipping;
    bool initialized;
    bool fullscreen;
    bool hardwareAccelerated;
    bool ownsSurface;
    RenderBackend backend;

    COLORREF color;
	CPoint cursor;	// virtual cursor
	const uint8_t* keyboardState;  // SDL keyboard state
	uint8_t fontAtlasLUT[256][8];  // Font LUT: [character][row], each byte is 8 pixels (bits)
	bool fontLoaded;

public:
    CGraph()
        : window(NULL), renderer(NULL), frameTexture(NULL), windowSurface(NULL), surface(NULL), pixels(NULL),
          width(640), height(480), pitch(0),
          title("SDL Window"), clipping(false), initialized(false), fullscreen(false),
          hardwareAccelerated(false), ownsSurface(false), backend(BackendSoftware),
                    color(MAKERGB(0, 0, 0)), keyboardState(NULL), fontLoaded(true) {
                memcpy(fontAtlasLUT, CGRAPH_DEFAULT_FONT_LUT, sizeof(fontAtlasLUT));
    }

    CGraph(int w, int h, std::string title = "SDL Window") {
        create(w, h, title);
    }

    ~CGraph() {
        if (frameTexture)
            SDL_DestroyTexture(frameTexture);
        if (renderer)
            SDL_DestroyRenderer(renderer);
        if (ownsSurface && surface)
            SDL_FreeSurface(surface);
        if (window)
            SDL_DestroyWindow(window);
        if (initialized)
            SDL_Quit();
    }

    bool create(int w, int h, std::string title, 
                bool fullscreen = false,
                int flags = SDL_WINDOW_SHOWN,
                bool preferHardware = true) {
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

        // Always keep a CPU surface available for software backend and optional HW pixel mode.
        surface = SDL_CreateRGBSurfaceWithFormat(0,
                                                 this->width,
                                                 this->height,
                                                 32,
                                                 SDL_PIXELFORMAT_ARGB8888);
        if (!surface)
            return false;
        ownsSurface = true;

        if (preferHardware) {
            renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
            if (renderer) {
                // HW primitive backend draws directly to renderer backbuffer.
                backend = BackendHardwarePrimitives;
                hardwareAccelerated = true;
                frameTexture = NULL;
            }
        }

        // Software fallback when hardware initialization fails or is disabled.
        if (!renderer) {
            if (frameTexture) {
                SDL_DestroyTexture(frameTexture);
                frameTexture = NULL;
            }
            if (renderer) {
                SDL_DestroyRenderer(renderer);
                renderer = NULL;
            }

            if (surface && ownsSurface) {
                SDL_FreeSurface(surface);
                surface = NULL;
            }

            surface = SDL_GetWindowSurface(window);
            ownsSurface = false;
            hardwareAccelerated = false;
            backend = BackendSoftware;
        }

        if (!surface)
            return false;

        windowSurface = SDL_GetWindowSurface(window);
        this->pixels = (uint32_t*) surface->pixels;
        this->pitch = surface->pitch >> 2;
        this->color = SDL_MapRGB(surface->format, 0, 0, 0);
        initialized = true;
        
        // Initialize keyboard state
        int numKeys;
        keyboardState = SDL_GetKeyboardState(&numKeys);
        
        SDL_FillRect(surface, NULL, 0);
        update();
        return true;
    }

    SDL_Window * getWindow() {
        return window;
    }

    SDL_Surface * getSurface() {
        return surface;
    }

    SDL_Renderer * getRenderer() {
        return renderer;
    }

    bool hasRenderer() const {
        return renderer != NULL;
    }

    bool isHardwareAccelerated() const {
        return hardwareAccelerated;
    }

    std::string getHardwareAcceleratorName() const {
        if (!renderer || !hardwareAccelerated)
            return "Software (SDL_Surface)";

        SDL_RendererInfo info;
        if (SDL_GetRendererInfo(renderer, &info) == 0 && info.name)
            return std::string(info.name);

        return "Unknown SDL Renderer";
    }

    SDL_Texture * createRenderTargetTexture(uint32_t w, uint32_t h, bool blended = true) {
        if (!renderer) return NULL;
        SDL_Texture * tex = SDL_CreateTexture(renderer,
                                              SDL_PIXELFORMAT_RGBA8888,
                                              SDL_TEXTUREACCESS_TARGET,
                                              w,
                                              h);
        if (!tex) return NULL;
        if (blended)
            SDL_SetTextureBlendMode(tex, SDL_BLENDMODE_BLEND);
        return tex;
    }

    void destroyTexture(SDL_Texture * tex) {
        if (tex) SDL_DestroyTexture(tex);
    }

    bool setRenderTarget(SDL_Texture * target) {
        if (!renderer) return false;
        return SDL_SetRenderTarget(renderer, target) == 0;
    }

    void setBlendMode(SDL_BlendMode mode) {
        if (renderer)
            SDL_SetRenderDrawBlendMode(renderer, mode);
    }

    void setColorRGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
        if (renderer)
            SDL_SetRenderDrawColor(renderer, r, g, b, a);
    }

    void clearTargetRGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255) {
        if (!renderer) return;
        SDL_SetRenderDrawColor(renderer, r, g, b, a);
        SDL_RenderClear(renderer);
    }

    void fillRectRGBA(int x, int y, int w, int h, uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255) {
        if (!renderer) return;
        SDL_Rect rect = { x, y, w, h };
        SDL_SetRenderDrawColor(renderer, r, g, b, a);
        SDL_RenderFillRect(renderer, &rect);
    }

    void clearBackbuffer(uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255) {
        if (!renderer) return;
        SDL_SetRenderTarget(renderer, NULL);
        SDL_SetRenderDrawColor(renderer, r, g, b, a);
        SDL_RenderClear(renderer);
    }

    void copyTextureToBackbuffer(SDL_Texture * tex) {
        if (!renderer || !tex) return;
        SDL_RenderCopy(renderer, tex, NULL, NULL);
    }

    void clear(uint32_t color) {
        if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderClear(renderer);
        } else {
            SDL_FillRect(surface, NULL, color);
        }
    }

    void clear() {
        clear(this->color);
    }

    void clear(int r, int g, int b) {
        if (surface)
            clear(SDL_MapRGB(surface->format, r, g, b));
    }

    void update() {
        if (backend == BackendHardwarePrimitives) {
            SDL_RenderPresent(renderer);
        } else if (backend == BackendHardwarePixels) {
            SDL_UpdateTexture(frameTexture, NULL, pixels, surface->pitch);
            SDL_SetRenderTarget(renderer, NULL);
            SDL_RenderClear(renderer);
            SDL_RenderCopy(renderer, frameTexture, NULL, NULL);
            SDL_RenderPresent(renderer);
        } else {
            SDL_UpdateWindowSurface(window);
        }
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
        if (backend == BackendHardwarePrimitives) {
            // Raw pixel clients use the HW pixel backend (single upload during update).
            if (!frameTexture) {
                frameTexture = SDL_CreateTexture(renderer,
                                                 SDL_PIXELFORMAT_ARGB8888,
                                                 SDL_TEXTUREACCESS_STREAMING,
                                                 this->width,
                                                 this->height);
            }

            if (frameTexture) {
                backend = BackendHardwarePixels;
                SDL_FillRect(surface, NULL, 0);
                pixels = (uint32_t*) surface->pixels;
                pitch = surface->pitch >> 2;
            }
        }
        return pixels;
    }

    CBitmap * getBitmap() {
        CBitmap * bmp = new CBitmap(width, height);
        memcpy(bmp->getPixels(), pixels, width * height * sizeof(uint32_t));
        return bmp;
    }

    // Keyboard handling
    
    virtual void onKeyDown(SDL_Keysym keysym) {
        // Override in subclass to handle key down events
    }
    
    virtual void onKeyUp(SDL_Keysym keysym) {
        // Override in subclass to handle key up events
    }
    
    bool isKeyPressed(SDL_Keycode keycode) {
        if (!keyboardState) return false;
        SDL_Scancode scancode = SDL_GetScancodeFromKey(keycode);
        return keyboardState[scancode] != 0;
    }
    
    bool isKeyDown(SDL_Keycode keycode) {
        return isKeyPressed(keycode);
    }
    
    // Helper methods for common keys
    bool isKeyEscapePressed() { return isKeyPressed(SDLK_ESCAPE); }
    bool isKeySpacePressed() { return isKeyPressed(SDLK_SPACE); }
    bool isKeyEnterPressed() { return isKeyPressed(SDLK_RETURN); }
    bool isKeyUpPressed() { return isKeyPressed(SDLK_UP); }
    bool isKeyDownPressed() { return isKeyPressed(SDLK_DOWN); }
    bool isKeyLeftPressed() { return isKeyPressed(SDLK_LEFT); }
    bool isKeyRightPressed() { return isKeyPressed(SDLK_RIGHT); }
    bool isKeyWPressed() { return isKeyPressed(SDLK_w); }
    bool isKeyAPressed() { return isKeyPressed(SDLK_a); }
    bool isKeyDPressed() { return isKeyPressed(SDLK_d); }
    bool isKeySPressed() { return isKeyPressed(SDLK_s); }

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
					onKeyDown(e.key.keysym);
					switch (e.key.keysym.sym) {
						case SDLK_ESCAPE:
							quit = 1;
							break;
					}
				}
				if (e.type == SDL_KEYUP) {
					onKeyUp(e.key.keysym);
				}
			}
			render();
			update();
		}    
    }

    // rendering methods

    void setColor(int r, int g, int b) {
        if (surface)
            color = SDL_MapRGB(surface->format, r, g, b);
        else
            color = MAKERGB(r, g, b);
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
        if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderDrawPoint(renderer, x, y);
        } else {
            pixels[y * surface->w + x] = color;
        }
    }
    
    void plotPixel(uint32_t x, uint32_t y) {
        plotPixel(x, y, color);
    }

    void lineHorz(int x1, int x2, int y, COLORREF color) {
        if (clipping) {
            if (y < 0 || y >= surface->h) return;
    
            if (x1 > x2) {
                int temp = x1;
                x1 = x2;
                x2 = temp;
            }
    
            if (x1 < 0) x1 = 0;
            if (x2 >= surface->w) x2 = surface->w - 1;
            if (x1 > x2) return;
        }

        if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderDrawLine(renderer, x1, y, x2, y);
        } else {
            for (int x = x1; x <= x2; x++) {
                pixels[y * surface->w + x] = color;
            }
        }
    }

    void lineHorz(int x1, int x2, int y) {
        lineHorz(x1, x2, y, color);
    }

    void lineVert(int x, int y1, int y2, COLORREF color) {
        if (clipping) {
            if (x < 0 || x >= surface->w) return;
            if (y1 > y2) {
                int temp = y1;
                y1 = y2;
                y2 = temp;
            }
            if (y1 < 0) y1 = 0;
            if (y2 >= surface->h) y2 = surface->h - 1;
            if (y1 > y2) return;
        }

        if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderDrawLine(renderer, x, y1, x, y2);
        } else {
            for (int y = y1; y <= y2; y++) {
                pixels[y * surface->w + x] = color;
            }
        }
    }

    void lineVert(int x, int y1, int y2) {
        lineVert(x, y1, y2, color);
    }

    void rectangle(int x, int y, int width, int height, COLORREF color) {
        if (clipping) {
            if (width <= 0 || height <= 0) return;

            if (x >= surface->w || y >= surface->h) return;
            if (x + width <= 0 || y + height <= 0) return;

            if (x < 0) {
                width += x;
                x = 0;
            }
            if (y < 0) {
                height += y;
                y = 0;
            }
            if (x + width > surface->w) width = surface->w - x;
            if (y + height > surface->h) height = surface->h - y;
            if (width <= 0 || height <= 0) return;
        }

        if (backend == BackendHardwarePrimitives) {
            SDL_Rect rect = { x, y, width, height };
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderFillRect(renderer, &rect);
        } else {
            for (int i = 0; i < height; i++) {
                uint32_t* row = &pixels[(y + i) * pitch + x];
                for(int j = 0; j < width; j++)
                    row[j] = color;
            }
        }
    }
    
    void rectangle(int x, int y, int width, int height) {
        rectangle(x, y, width, height, color);
    }

    void line(int x1, int y1, int x2, int y2, COLORREF color) {
        if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
            return;
        }

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
	}

	// Draws a circle at (xc, yc) with radius r
	void drawCircle(int xc, int yc, int r, COLORREF color) {
    if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);
            int xh = 0, yh = r;
            int dh = 3 - (r << 1);
            while (yh >= xh) {
                SDL_RenderDrawPoint(renderer, xc + xh, yc + yh);
                SDL_RenderDrawPoint(renderer, xc - xh, yc + yh);
                SDL_RenderDrawPoint(renderer, xc + xh, yc - yh);
                SDL_RenderDrawPoint(renderer, xc - xh, yc - yh);
                SDL_RenderDrawPoint(renderer, xc + yh, yc + xh);
                SDL_RenderDrawPoint(renderer, xc - yh, yc + xh);
                SDL_RenderDrawPoint(renderer, xc + yh, yc - xh);
                SDL_RenderDrawPoint(renderer, xc - yh, yc - xh);
                xh++;
                if (dh > 0) {
                    yh--;
                    dh = dh + ((xh - yh) << 2) + 10;
                }
                else
                    dh = dh + (xh << 2) + 6;
            }
            return;
        }

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
    if (backend == BackendHardwarePrimitives) {
            SDL_SetRenderDrawColor(renderer, RGB_RED(color), RGB_GREEN(color), RGB_BLUE(color), 255);

            int xh = 0, yh = r;
            int dh = 3 - (r << 1);
            while (yh >= xh) {
                SDL_RenderDrawLine(renderer, xc - xh, yc + yh, xc + xh, yc + yh);
                SDL_RenderDrawLine(renderer, xc - xh, yc - yh, xc + xh, yc - yh);
                SDL_RenderDrawLine(renderer, xc - yh, yc + xh, xc + yh, yc + xh);
                SDL_RenderDrawLine(renderer, xc - yh, yc - xh, xc + yh, yc - xh);
                xh++;
                if (dh > 0) {
                    yh--;
                    dh = dh + ((xh - yh) << 2) + 10;
                }
                else
                    dh = dh + (xh << 2) + 6;
            }
            return;
        }

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

	// Load font atlas from 1-bit BMP file (256x64 pixels = 32x8 grid of 8x8 characters)
	// Reads raw bits directly to build 1-bit LUT, independent of CLUT/color interpretation.
	// Bit = 1 in LUT means "character pixel" (inverted BMP: character pixels have bit=0 in file).
	bool loadFont(const std::string& fontPath) {
		std::ifstream is(fontPath, std::ios::binary);
        if (!is) return false;

		// Read BMP file header and info header
		BITMAPFILEHEADER fhdr;
		BITMAPINFOHEADER ihdr;
		is.read((char*)&fhdr, sizeof(fhdr));
		is.read((char*)&ihdr, sizeof(ihdr));

        if (fhdr.bfType != 0x4D42 || ihdr.biBitCount != 1)
            return false;

		int bmpW = ihdr.biWidth;    // expected 256
		int bmpH = ihdr.biHeight;   // expected 64
		int gridCols = bmpW / 8;    // 32 characters per row
		int gridRows = bmpH / 8;    // 8 rows of characters

		// Each row of the BMP is padded to a 4-byte boundary
		int stride = ((bmpW + 31) / 32) * 4;

		// Read all raw pixel rows (BMP is stored bottom-to-top)
		std::vector<std::vector<uint8_t>> rows(bmpH, std::vector<uint8_t>(stride));
		is.seekg(fhdr.bfOffBits, std::ios::beg);
		for (int i = 0; i < bmpH; i++) {
			// BMP rows are bottom-to-top, so row 0 in file = last scanline
			is.read((char*)rows[bmpH - 1 - i].data(), stride);
		}
		is.close();

		// Determine polarity: sample a pixel inside character ' ' (ASCII 32)
		// Space should be blank. If the bit for space's first pixel is SET,
		// then bit=1 means background (uninverted file); invert our interpretation.
		int spaceCol = (32 % gridCols) * 8;
		int spaceRow = (32 / gridCols) * 8;
		uint8_t sampleByte = rows[spaceRow][spaceCol / 8];
		// bit 7 of that byte = leftmost pixel of space
		bool invertBits = (sampleByte >> 7) & 1;  // if set, bit=1 is background

		// Build 1-bit LUT: fontAtlasLUT[charCode][row]
		// Each byte: bit 7 = leftmost pixel, bit 0 = rightmost
		// A set bit means "draw this pixel with the text color"
		memset(fontAtlasLUT, 0, sizeof(fontAtlasLUT));
		for (int charCode = 0; charCode < 256; charCode++) {
			int charGridCol = charCode % gridCols;
			int charGridRow = charCode / gridCols;
			int startX = charGridCol * 8;
			int startY = charGridRow * 8;

			for (int row = 0; row < 8; row++) {
				// Extract the byte containing the 8 pixels for this row
				// Each group of 8 pixels is one byte in a 1-bit BMP
				int byteIndex = startX / 8;
				int bitShift  = 7 - (startX % 8);  // 0 when byte-aligned
				uint8_t rawByte = (bitShift == 7)
					? rows[startY + row][byteIndex]
					: (uint8_t)((rows[startY + row][byteIndex] << (7 - bitShift)) |
					            (rows[startY + row][byteIndex + 1] >> (bitShift + 1)));
				// For an inverted BMP, character bits are 0; flip to make them 1
				fontAtlasLUT[charCode][row] = invertBits ? rawByte : ~rawByte;
			}
		}

		fontLoaded = true;
		return true;
	}

	// Fast text rendering using 1-bit LUT with any specified color
	void drawChar(char ch, int x, int y, COLORREF textColor) {
		if (!fontLoaded) return;
		
		unsigned char charCode = (unsigned char)ch;
		
		// Use 1-bit LUT to render character with specified color
		for (int row = 0; row < 8; row++) {
			uint8_t rowByte = fontAtlasLUT[charCode][row];
			
			// Check each bit and draw pixel if set (1 = character pixel)
			for (int col = 0; col < 8; col++) {
				if (!(rowByte & (0x80 >> col))) {  // Bit 7 = left, Bit 0 = right
					plotPixel(x + col, y + row, textColor);
				}
			}
		}
	}
	
	void drawChar(char ch, int x, int y) {
		drawChar(ch, x, y, color);
	}
	
	void drawText(const std::string& text, int x, int y, COLORREF textColor) {
		int currentX = x;
		for (char ch : text) {
			drawChar(ch, currentX, y, textColor);
			currentX += 8;  // 8 pixels per character
		}
	}
	
	void drawText(const std::string& text, int x, int y) {
		drawText(text, x, y, color);
	}

    // Draws a tiny renderer status badge (HW/SW) at the given position.
    // If bitmap font is loaded, also prints the mode text.
    void drawRendererStatus(int x = 8, int y = 8) {
        COLORREF oldColor = color;

        if (hardwareAccelerated)
            setColor(64, 220, 120);
        else
            setColor(235, 170, 56);

        filledCircle(x, y, 4);

        setColor(18, 24, 36);
        drawCircle(x, y, 5);

        if (fontLoaded) {
            setColor(220, 230, 245);
            drawText(hardwareAccelerated ? "HW" : "SW", x + 9, y - 4);
        }

        setColor(oldColor);
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

// ============================================================================
// 2D ENGINE COMPONENTS
// ============================================================================

// CTransform - Manages position, rotation, and scale
class CTransform {
public:
	CVector2D position;
	CVector2D velocity;
	double rotation;  // radians
	double scale;
	
	CTransform() : position(0, 0), velocity(0, 0), rotation(0), scale(1.0) {}
	
	void update(double deltaTime) {
		position.incX(velocity.getX() * deltaTime);
		position.incY(velocity.getY() * deltaTime);
	}
	
	CVector2D getWorldPosition() const { return position; }
	void setPosition(double x, double y) { position.setX(x).setY(y); }
	void translate(double dx, double dy) { position.incX(dx).incY(dy); }
	
	void setVelocity(double vx, double vy) { velocity.setX(vx).setY(vy); }
	void setRotation(double rad) { rotation = rad; }
	void rotate(double dRad) { rotation += dRad; }
	void setScale(double s) { scale = s; }
};

// CCollider - Basic collision detection
class CCollider {
public:
	enum Type { AABB, Circle };
	
	Type type;
	double width, height;  // for AABB
	double radius;          // for circle
	
	CCollider(double w, double h) : type(AABB), width(w), height(h), radius(0) {}
	CCollider(double r) : type(Circle), width(0), height(0), radius(r) {}
	
	bool intersects(const CVector2D& pos1, const CCollider& other, const CVector2D& pos2) const {
		if (type == AABB && other.type == AABB) {
			// AABB-AABB collision
			return !(pos1.getX() + width < pos2.getX() ||
					 pos1.getX() > pos2.getX() + other.width ||
					 pos1.getY() + height < pos2.getY() ||
					 pos1.getY() > pos2.getY() + other.height);
		}
		else if (type == Circle && other.type == Circle) {
			// Circle-Circle collision
			double dx = pos1.getX() - pos2.getX();
			double dy = pos1.getY() - pos2.getY();
			double dist = sqrt(dx * dx + dy * dy);
			return dist < (radius + other.radius);
		}
		else if (type == AABB && other.type == Circle) {
			// AABB-Circle collision
			double closestX = CMath::clampf(pos2.getX(), pos1.getX(), pos1.getX() + width);
			double closestY = CMath::clampf(pos2.getY(), pos1.getY(), pos1.getY() + height);
			double dx = pos2.getX() - closestX;
			double dy = pos2.getY() - closestY;
			return (dx * dx + dy * dy) < (other.radius * other.radius);
		}
		else {
			// Circle-AABB collision (swap and recurse)
			return other.intersects(pos2, *this, pos1);
		}
	}
};

// CGameObject - Base class for all game objects
class CGameObject {
protected:
	CTransform transform;
	CCollider* collider;
	COLORREF color;
	bool active;
	bool visible;
	int layer;
	
public:
	CGameObject(double x = 0, double y = 0, COLORREF col = MAKERGB(255, 255, 255))
		: color(col), active(true), visible(true), layer(0), collider(nullptr) {
		transform.setPosition(x, y);
	}
	
	virtual ~CGameObject() {
		if (collider) delete collider;
	}
	
	virtual void update(double deltaTime) {
		transform.update(deltaTime);
	}
	
	virtual void render(CGraph& gfx) {
		if (!visible) return;
		// Draw a circle by default
		gfx.setColor((color >> 16) & 0xFF, (color >> 8) & 0xFF, color & 0xFF);
		gfx.filledCircle((int)transform.position.getX(), 
					     (int)transform.position.getY(), 5);
	}
	
	void setCollider(CCollider* col) {
		if (collider) delete collider;
		collider = col;
	}
	
	bool checkCollision(const CGameObject& other) const {
		if (!collider || !other.collider) return false;
		return collider->intersects(transform.position, *other.collider, other.transform.position);
	}
	
	CTransform& getTransform() { return transform; }
	const CTransform& getTransform() const { return transform; }
	
	void setActive(bool a) { active = a; }
	void setVisible(bool v) { visible = v; }
	void setLayer(int l) { layer = l; }
	
	bool isActive() const { return active; }
	bool isVisible() const { return visible; }
	int getLayer() const { return layer; }
};

// CAnimator - Simple sprite animation
class CAnimator {
private:
	std::vector<CBitmap> frames;
	int currentFrame;
	double frameTime;
	double currentTime;
	bool isPlaying;
	bool isLooping;
	
public:
	CAnimator() : currentFrame(0), frameTime(0.1), currentTime(0), 
				  isPlaying(false), isLooping(true) {}
	
	void addFrame(const CBitmap& frame) {
		frames.push_back(frame);
	}
	
	void play(bool loop = true) {
		isPlaying = true;
		isLooping = loop;
		currentFrame = 0;
		currentTime = 0;
	}
	
	void stop() {
		isPlaying = false;
		currentFrame = 0;
		currentTime = 0;
	}
	
	void update(double deltaTime) {
		if (!isPlaying || frames.empty()) return;
		
		currentTime += deltaTime;
		if (currentTime >= frameTime) {
			currentTime -= frameTime;
			currentFrame++;
			
			if (currentFrame >= (int)frames.size()) {
				if (isLooping) {
					currentFrame = 0;
				} else {
					isPlaying = false;
					currentFrame = frames.size() - 1;
				}
			}
		}
	}
	
	const CBitmap& getCurrentFrame() const {
		if (frames.empty()) {
			static CBitmap empty;
			return empty;
		}
		return frames[CMath::clampf(currentFrame, 0, frames.size() - 1)];
	}
	
	void setFrameTime(double t) { frameTime = t; }
};

// CScene - Manages collection of game objects
class CScene {
private:
	std::vector<CGameObject*> objects;
	int objectIdCounter;
	
public:
	CScene() : objectIdCounter(0) {}
	
	virtual ~CScene() {
		for (auto obj : objects) {
			delete obj;
		}
		objects.clear();
	}
	
	void addObject(CGameObject* obj) {
		if (obj) {
			objects.push_back(obj);
			std::sort(objects.begin(), objects.end(),
				[](CGameObject* a, CGameObject* b) {
					return a->getLayer() < b->getLayer();
				});
		}
	}
	
	void removeObject(CGameObject* obj) {
		auto it = std::find(objects.begin(), objects.end(), obj);
		if (it != objects.end()) {
			delete *it;
			objects.erase(it);
		}
	}
	
	void update(double deltaTime) {
		for (auto obj : objects) {
			if (obj->isActive()) {
				obj->update(deltaTime);
			}
		}
	}
	
	void render(CGraph& gfx) {
		for (auto obj : objects) {
			if (obj->isVisible()) {
				obj->render(gfx);
			}
		}
	}
	
	std::vector<CGameObject*>& getObjects() { return objects; }
	const std::vector<CGameObject*>& getObjects() const { return objects; }
	
	int getObjectCount() const { return objects.size(); }
	void clear() {
		for (auto obj : objects) delete obj;
		objects.clear();
	}
};

// CTimer - Frame timing utility
class CTimer {
private:
	uint32_t lastTicks;
	double deltaTime;
	double totalTime;
	
public:
	CTimer() : lastTicks(SDL_GetTicks()), deltaTime(0), totalTime(0) {}
	
	void update() {
		uint32_t currentTicks = SDL_GetTicks();
		deltaTime = (currentTicks - lastTicks) / 1000.0;
		totalTime += deltaTime;
		lastTicks = currentTicks;
		
		// Cap delta time to prevent large jumps
		if (deltaTime > 0.05) deltaTime = 0.05;
	}
	
	double getDeltaTime() const { return deltaTime; }
	double getTotalTime() const { return totalTime; }
};

} // namespace daniel

#endif // __CGRAPH_H__


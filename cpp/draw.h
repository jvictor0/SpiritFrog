#pragma once

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h> 

#define atomic_increment32(destination) __sync_add_and_fetch(destination, 1)
#define atomic_increment64(destination) __sync_add_and_fetch(destination, 1)

struct Point 
{
    Point(double x, double y) : m_x(x), m_y(y) { }
    double m_x, m_y;
};

struct Pixel 
{
    Pixel(int64_t x, int64_t y) : m_x(x), m_y(y) { }
    int64_t m_x, m_y;
};

int64_t Real2Pix(double x, double xmin, double xmax, int64_t width)
{
    return (x - xmin) / (xmax-xmin) * width;
}

double Pix2Real(int64_t pix, double xmin, double xmax, int64_t width)
{
    return (pix * (xmax - xmin))/width + xmin;
}

struct Rect
{
    Rect(double xmin, double xmax, double ymin, double ymax)
        : m_xmin(xmin)
        , m_xmax(xmax)
        , m_ymin(ymin)
        , m_ymax(ymax)
    {}
    Rect(double coord)
        : m_xmin(-coord)
        , m_xmax(coord)
        , m_ymin(-coord)
        , m_ymax(coord)
    {}
    
    Pixel Real2Pix(Point in, int64_t width, int64_t height)
    {
        return Pixel(::Real2Pix(in.m_x, m_xmin, m_xmax, width),
                     ::Real2Pix(in.m_y, m_ymin, m_ymax, height));
    }
    Point Pix2Real(Pixel in, int64_t width, int64_t height)
    {
        return Point(::Pix2Real(in.m_x, m_xmin, m_xmax, width),
                     ::Pix2Real(in.m_y, m_ymin, m_ymax, height));
    }

    Point Rand()
    {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        Point result =  Point(x * (m_xmax - m_xmin) + m_xmin, 
                              y * (m_ymax - m_ymin) + m_ymin);
        return result;
    }

    double m_xmin, m_xmax, m_ymin, m_ymax;
};

template<typename T>
struct Array3d
{
    Array3d(size_t x, size_t y, size_t z) 
        : m_x(x)
        , m_y(y)
        , m_z(z)
        , m_data(new T[x * y * z])
    { 
        memset(m_data, 0, sizeof(T) * x * y * z);
    }

    Array3d(const char* filename)
    {
        std::ifstream stream(filename);
        stream.read(reinterpret_cast<char*>(this), sizeof(size_t) * 3);
        m_data = new T[m_x * m_y * m_z];
        stream.read(reinterpret_cast<char*>(m_data), ByteLength());
        stream.close();
    }
    Array3d(const Array3d&) = delete; // fuck you
    
    ~Array3d()
    {
        delete m_data;
    }

    T& Get(size_t x, size_t y, size_t z)
    {
        assert (x >= 0 && x < m_x);
        assert (y >= 0 && y < m_y);
        assert (z >= 0 && z < m_z);
        return m_data[x * m_y * m_z + y * m_z + z];
    }
    
    void AtomicIncrement(size_t x, size_t y, size_t z)
    {
        assert (x >= 0 && x < m_x);
        assert (y >= 0 && y < m_y);
        assert (z >= 0 && z < m_z);
        if (sizeof(T) == 4)
        {
            atomic_increment32(&m_data[x * m_y * m_z + y * m_z + z]);
        }
        else if (sizeof(T) == 8)
        {
            atomic_increment64(&m_data[x * m_y * m_z + y * m_z + z]);
        }
        else
        {
            assert (false);
        }
    }

    void Serialize(const char* filename)
    {
        std::ofstream stream(filename);
        stream.write(reinterpret_cast<char*>(this), sizeof(size_t) * 3);
        stream.write(reinterpret_cast<char*>(m_data), ByteLength());
        stream.close();
    }

    size_t ByteLength()
    {
        return sizeof(T) * m_x * m_y * m_z;
    }

    size_t m_x, m_y, m_z;
    T* m_data;
};


#pragma once

#include "draw.h"
#include <vector>
#include <thread>
#include <time.h>
#include "CImg-1.6.8_pre110615/CImg.h"

struct BuddhabrotParams
{
    bool PreFilter(Point pt)
    {
        double x = pt.m_x;
        double y = pt.m_y;
        if ((x+1)*(x+1) + y * y < 1.0/16)
        {
            return true;
        }
        double q = (x - 0.25) * (x - 0.25) + y * y;
        if (q * (q + (x - 0.25)) < 0.25 * y * y)
        {
            return true;
        }
        return false;
    }

    Point NextPoint(Point z, Point c)
    {
        return Point(z.m_x * z.m_x - z.m_y * z.m_y + c.m_x,
                     2 * z.m_x * z.m_y + c.m_y);
    }

    bool Done(Point z)
    {
        return (z.m_x > 2.0 || z.m_x < -2.0 || z.m_y > 2.0 || z.m_y < -2.0);
    }
};

template<class Params>
void BuddhaAlgorithmThread(Array3d<uint32_t>& data, std::vector<uint64_t>& its, Params& params, Rect& rect, int64_t seconds)
{
    time_t t0 = time(0);
    while (difftime(time(0), t0) < seconds)
    {
        Point c = rect.Rand();
        if (params.PreFilter(c))
        {
            continue;
        }
        Point z = Point(0,0);
        std::vector<Point> boxes;
        for (uint64_t i = 0; i < its.back(); ++i)
        {
            z = params.NextPoint(z,c);
            if (params.Done(z))
            {
                uint64_t minj = 0;
                for (minj = 0; minj < its.size(); ++minj)
                {                    
                    if (i < its[minj])
                    {
                        break;
                    }
                }
                for (uint64_t k = 0; k < boxes.size(); ++k)
                {
                    Pixel p = rect.Real2Pix(boxes[k], data.m_x, data.m_y);                     
                    if (p.m_x >= 0 && p.m_y >= 0 && p.m_x < data.m_x && p.m_y < data.m_y)
                    {
                        for (uint64_t j = minj; j < its.size(); ++j)
                        {                    
                            data.AtomicIncrement(p.m_x, p.m_y, j);
                        }
                    }
                }
                break;
            }
            else
            {
                boxes.push_back(z);
            }
        }               
    }
}

template<class Params>
void BuddhaAlgorithm(int numProcs, Array3d<uint32_t>& data, std::vector<uint64_t>& its, Params& params, Rect& rect, int64_t seconds)
{
    std::vector<std::thread*> procs;

    for (int i = 0; i < numProcs; ++i)
    {
        procs.push_back(new std::thread([&]() { BuddhaAlgorithmThread(data, its, params, rect, seconds); }));
    }
    for (int i = 0; i < numProcs; ++i)
    {
        procs[i]->join();
    }    
    for (int i = 0; i < numProcs; ++i)
    {
        delete procs[i];
    }
}

void ChannelMaxs(Array3d<uint32_t>& in, std::vector<uint32_t>& maxes)
{
    for (size_t i = 0; i < in.m_z; ++i)
    {
        uint32_t mx = 0;
        for (size_t j = 0; j < in.m_x; ++j)
        {
            for (size_t k = 0; k < in.m_y; ++k)
            {
                mx = std::max(in.Get(j,k,i), mx);
            }
        }
        maxes.push_back(mx);
    }
}

struct RGB
{
    RGB(unsigned char r, unsigned char g, unsigned char b)
        : m_r(r), m_g(g), m_b(b) { }
    RGB(double r, double g, double b)
        : m_r(255 *r), m_g(255 * g), m_b(255 * b) { }
    
    unsigned char m_r, m_g, m_b;
};

struct HistomixParams
{
    std::vector<uint32_t> m_maxes;
    void Init(Array3d<uint32_t>& in)
    {
        ChannelMaxs(in, m_maxes);
        assert (m_maxes.size() == 1 || m_maxes.size() == 3);
    }
    RGB Get(uint32_t* channels)
    {
        return RGB(((double)channels[0])/m_maxes[0], 
                   ((double)channels[1])/m_maxes[1],
                   ((double)channels[2])/m_maxes[2]);
    }
};

template<typename Params>
void Draw(Array3d<uint32_t>& in, cimg_library::CImg<unsigned char>& out, Params& params)
{
    params.Init(in);
    std::cout << "init";
    for (size_t i = 0; i < in.m_x; ++i)
    {
        for (size_t j = 0; j < in.m_y; ++j)
        {
            RGB color = params.Get(&in.Get(i,j,0));
            const unsigned char colarr[] = { color.m_r, color.m_g, color.m_b };
            out.draw_point(i, j, colarr);
        }
    }    
}

Array3d<uint32_t>* GenerateHistogram(size_t width, int numProcs, std::vector<uint64_t>& its, double seconds)
{
    Array3d<uint32_t>* data = new Array3d<uint32_t>(width, width, its.size());
    BuddhabrotParams params;
    Rect rect(2.0);
    std::cout << "generating historgram" << std::endl;
    BuddhaAlgorithm(numProcs, *data, its, params, rect, seconds);
    std::cout << "done with histogram" << std::endl;
    return data;
}

void GenerateAndSaveHistogram(size_t width, int numProcs, std::vector<uint64_t>& its, double seconds, const char* filename)
{
    Array3d<uint32_t>* data = GenerateHistogram(width, numProcs, its, seconds);
    data->Serialize(filename);
    delete data;
}

void Buddhabrot1T()
{
    size_t width = 1024;
    std::vector<uint64_t> its;
    its.push_back(100);
    its.push_back(500);
    its.push_back(2000);
    Array3d<uint32_t> data(width, width, its.size());
    BuddhabrotParams params;
    Rect rect(2.0);
    std::cout << "generating historgram" << std::endl;
    BuddhaAlgorithm(8, data, its, params, rect, 5);
    std::cout << "done with histogram" << std::endl;

    cimg_library::CImg<unsigned char> img(width, width, 1, 3, 0);
    
    HistomixParams hmp;
    Draw(data, img, hmp);
    
    cimg_library::CImgDisplay disp(img, "buddha");
    while (!disp.is_closed())
    {
        sleep(1);
    }
    img.save_bmp("/vagrant/cppman.bmp");
}

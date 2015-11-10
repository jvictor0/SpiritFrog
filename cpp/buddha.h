#pragma once

#include "draw.h"
#include <vector>
#include <thread>
#include <time.h>
#include "CImg-1.6.8_pre110615/CImg.h"
#include <unistd.h>
#include <termios.h>


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
        : m_r(255 * std::max(0.0,std::min(1.0,r)))
        , m_g(255 * std::max(0.0,std::min(1.0,g)))
        , m_b(255 * std::max(0.0,std::min(1.0,b))) 
    { }
    
    unsigned char m_r, m_g, m_b;
};

struct HistomixParams
{
    HistomixParams(Array3d<uint32_t>& in) : m_mixmat(in.m_z), m_toggle(in.m_z)
    {
        m_maxes.clear();
        ChannelMaxs(in, m_maxes);
        m_xstart = 0;
        m_ystart = 0;
        m_xend = in.m_x;
        m_yend = in.m_y;
        
        for (size_t i = 0; i < m_mixmat.size(); ++i)
        {
            m_toggle[i] = false;
            for (size_t j = 0; j < 3; ++j)
            {
                m_mixmat[i].push_back(i==j?1:0);
            }
        }
    }
    void Init(Array3d<uint32_t>& in)
    {
    }
    RGB Get(uint32_t* channels)
    {
        double result[3];
        for (int i = 0; i < 3; ++i)
        {
            result[i] = 0;
        }
        for (size_t i = 0; i < m_mixmat.size(); ++i)
        {
            if (!m_toggle[i])
            {
                for (int j = 0; j < 3; ++j)
                {
                    result[j] += m_mixmat[i][j] * ((double)channels[i]) / m_maxes[i];
                }
            }
        }
        return RGB(result[0], result[1], result[2]);
    }

    void Zoom(int64_t amt)
    {
        if (static_cast<int64_t>(m_xend - m_xstart) > 2 * amt
            && static_cast<int64_t>(m_yend - m_ystart) > 2 * amt)
        {
            m_yend -= amt; m_xend -= amt; m_ystart += amt; m_xstart += amt;
        }
    }
    void Move(int64_t x, int64_t y)
    {
        if (m_xstart + x >= 0)
        {
            m_xend += x; m_xstart += x;
        }
        if (m_ystart + y >= 0)
        {
            m_yend += y; m_ystart += y;
        }
    }
    void UnZoom(int64_t amt)
    {
        if (static_cast<int64_t>(m_xstart) >= amt && static_cast<int64_t>(m_ystart) >= amt)
        {
            m_yend += amt; m_xend += amt; m_ystart -= amt; m_xstart -= amt;
        }
    }

    std::vector<uint32_t> m_maxes;
    std::vector<std::vector<double>> m_mixmat;
    std::vector<bool> m_toggle;
    size_t m_xstart, m_xend, m_ystart, m_yend;
};

template<typename Params>
void Draw(Array3d<uint32_t>& in, cimg_library::CImg<unsigned char>& out, Params& params)
{
    params.Init(in);
    out.fill(0);
    for (size_t i = params.m_xstart; i < params.m_xend; ++i)
    {
        for (size_t j = params.m_ystart; j < params.m_yend; ++j)
        {
            RGB color = params.Get(&in.Get(i,j,0));
            const unsigned char colarr[] = { color.m_r, color.m_g, color.m_b };
            out.draw_point(i - params.m_xstart, j - params.m_ystart, colarr);
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
    its.push_back(1000);
    its.push_back(20000);
    Array3d<uint32_t> data(width, width, its.size());
    BuddhabrotParams params;
    Rect rect(2.0);
    std::cout << "generating historgram" << std::endl;
    BuddhaAlgorithm(8, data, its, params, rect, 30);
    std::cout << "done with histogram" << std::endl;

    data.Serialize("buddy.a");

    std::cout << "done with serialiing" << std::endl;

    cimg_library::CImg<unsigned char> img(width, width, 1, 3, 0);    
    HistomixParams hmp(data);
    hmp.m_mixmat[0][0]=0;
    hmp.m_mixmat[1][1]=0;
    hmp.m_mixmat[2][2]=0;

    hmp.m_mixmat[0][0] = 207.0/255.0;
    hmp.m_mixmat[0][1] = 181.0/255.0;
    hmp.m_mixmat[0][2] = 59.0/255.0;

    hmp.m_mixmat[1][1] = 0.2;
    hmp.m_mixmat[1][2] = 0.2;

    hmp.m_mixmat[2][0] = 150.0/255.0;
    hmp.m_mixmat[2][2] = 260.0/255.0;
    
    hmp.m_xstart = 104;
    hmp.m_xend = 531;
    hmp.m_ystart = 420;
    hmp.m_yend = 764;
    Draw(data, img, hmp);
    cimg_library::CImgDisplay disp(img, "buddha");
    int mixer = 0;
    int color = 0;
    while (true)
    {
        std::string str;
        std::cin >> str;
        for (size_t i = 0; i < str.size(); ++i)
        {
            switch (str[i])
            {
                case '0': mixer = 0; break;
                case '1': mixer = 1; break;
                case '2': mixer = 2; break;
                case 'r': color = 0; break;
                case 'g': color = 1; break;
                case 'b': color = 2; break;
                case '+': hmp.m_mixmat[mixer][color] += 0.1; break;
                case '-': hmp.m_mixmat[mixer][color] -= 0.1; break;
                case 't': hmp.m_toggle[mixer] = !hmp.m_toggle[mixer]; break;
                case 'a': hmp.Move(-10,0); break;
                case 's': hmp.Move(0, -10); break;
                case 'd': hmp.Move(10, 0); break;
                case 'w': hmp.Move(0, 10); break;
                case 'z': hmp.Zoom(10); break;
                case 'x': hmp.UnZoom(10); break;
                default: break;
            }
        }
        std::cout << hmp.m_xstart << ", " << hmp.m_xend << ", " << hmp.m_ystart << ", " << hmp.m_yend << std::endl;
        std::cout << "current " << mixer << "," << color << std::endl;
        for (size_t i = 0; i < hmp.m_mixmat.size(); ++i)
        {
            std::cout << (hmp.m_toggle[i] ? "(OFF)  " : "(ON) ");
            for (size_t j = 0; j < 3; ++j)
            {
                std::cout << hmp.m_mixmat[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "drawing" << std::endl;        
        Draw(data, img, hmp);
        img.display(disp);
    }
    img.save_bmp("/vagrant/cppman.bmp");
}
///94, 521, 420, 764

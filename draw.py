import numpy as np
import scipy.misc.pilutil as smp
import time
import bigfloat
import random

class Rect:
    def __init__(self, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def real2pix(self, x,y, width, height):
        return real2pix(x, self.xmin, self.xmax, width), real2pix(y, self.ymin, self.ymax, height)

    def pix2real(self, x,y, width, height):
        return pix2real(x, self.xmin, self.xmax, width), pix2real(y, self.ymin, self.ymax, height)        

    def rand(self):
        return random.random() * (self.xmax - self.xmin) + self.xmin, random.random() * (self.ymax - self.ymin) + self.ymin


def FloatRect(size):
    return Rect(-size, size, -size, size)

def BigFloatRect(size):
    return Rect(bigfloat.BigFloat.exact(-size), bigfloat.BigFloat.exact(size), bigfloat.BigFloat.exact(-size),bigfloat.BigFloat.exact(size))

UnitRect = Rect(-1.0,1.0,-1.0,1.0)
    
def real2pix(x, xmin, xmax, width):
    return int ((x - xmin) / (xmax-xmin) * width)

def pix2real(pix, xmin, xmax, width):
    return (pix * (xmax - xmin))/width + xmin

class Draw:
    def __init__(self, width, height, rect):        
        self.data = np.zeros( (width, height,3), dtype=np.uint8 ) 
        self.width = width
        self.height = height
        self.rect = rect

    def real2pix(self, x,y):
        return self.rect.real2pix(x,y,self.width,self.height)

    def pix2real(self, x,y):
        return self.rect.read2pix(x,y,self.width,self.height)
                
    def drawreal(self, x, y, val):
        px, py = self.real2pix(x, y)
        self.data[px,py] = val

    def drawpix(self, px, py, val):
        self.data[px,py] = val

    def save(self, filename):
        smp.toimage(self.data).save(filename)

    def getiter(self):
        return DrawIterator(self)

    def drawfunc(self, fn):
        print "drawing"
        for px in xrange(self.width):
            print "%f percent done" % (float(px)/self.width)
            for py in xrange(self.height):
                self.data[px, py] = fn(*self.pix2real(px,py))
        return self

class DrawIterator:
    def __init__(self, draw):
        self.draw = draw
        self.x = -1
        self.y = 0

    def __iter__(self):
        return self

    def next(self):
        self.x += 1
        if self.x == self.draw.width:
            self.x = 0
            self.y += 1
            if self.y == self.draw.height:
                raise StopIteration()
        x,y = self.draw.pix2real(self.x,self.y)
        return self.x, self.y, x, y

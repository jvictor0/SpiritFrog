import draw
import cpl
import bigfloat
import numpy as np
import multiprocessing
import cPickle
import time

def InBulb(x,y):
    if (x+1)**2 + y**2 < 1.0/16:
        return True
    q = (x-0.25)**2 + y **2
    if q * (q + (x - 0.25)) < 0.25 * y**2:
        return True
    return False

class Mandlebrot:
    def __init__(self, its):
        self.its = its

    def __call__(self, x,y):
        z = (x * 0.0, y * 0.0)
        c = (x, y)
        if not InBulb(x,y):
            for i in xrange(self.its):
                z = cpl.c_add(cpl.c_mult(z,z), c)
                if abs(z[0]) >= 2.0 or abs(z[1]) >= 2.0:
                    return [ 0,0 , int (255 * (float(i) / self.its))]
        return [0,0,0]

class Buddhabrot:
    def __init__(self, secs, its, width, height, rect=draw.FloatRect(2)):
        self.secs = secs
        self.its = its
        self.height = height
        self.width = width
        self.data = np.zeros((width, height, len(its)), dtype=np.uint32)
        self.rect = rect

    def populate(self):
        t0 = time.time()
        i = 0
        while True:
            i = i + 1
            if i % 1000 == 0:
                if time.time() - t0 > self.secs:
                    return
            boxes = []
            c = self.rect.rand()
            if InBulb(*c):
                continue
            z = (c[0] * 0.0, c[1] * 0.0)
            for i in xrange(self.its[-1]):
                z = cpl.c_add(cpl.c_mult(z,z), c)
                if abs(z[0]) >= 2.0 or abs(z[1]) >= 2.0:
                    for px,py in boxes:
                        for j in xrange(len(self.its)):
                            if i < self.its[j]:
                                self.data[px,py,j] += 1
                    break
                else:
                    boxes.append(self.rect.real2pix(z[0], z[1], self.width, self.height))        

    def pickle(self, filename):
        cPickle.dump(self, open(filename, "wb"))
                    
    def draw(self):
        print "converting buddha bro to draw"
        d = draw.Draw(self.width, self.height, self.rect)
        scale = [0,0,0]
        for i in xrange(3):            
            scale[i] = 255.0 / max([self.data[x,y,i] for x in xrange(self.width) for y in xrange(self.height)])
        scale[0] = scale[0]/2
        scale[1] = scale[1]
        print "max"
        for i in xrange(self.width):
            for j in xrange(self.height):
                d.data[i, j] = [int(self.data[i,j,0] * scale[0]), int(self.data[i,j,1] * scale[1]), int(self.data[i,j,2] * scale[2])]
        return d

def MultiBuddhabrot(numprocs, secs, its, width, height, rect=draw.FloatRect(2)):
    manager = multiprocessing.Manager()
    rlist = manager.list()
    def fn(rlist):
        print "start"
        bud = Buddhabrot(secs, its, width, height, rect)
        bud.populate()
        rlist.append(bud)
    procs = [multiprocessing.Process(target=fn,args=(rlist,)) for i in xrange(numprocs)]
    [p.start() for p in procs]
    [p.join() for p in procs]
    result = Buddhabrot(secs, its, width, height, rect)
    for bud in rlist:
        result.data += bud.data
    return result
        

if __name__ == '__main__':
    bigfloat.setcontext(bigfloat.precision(52))
    #man = Mandlebrot(100)
    #d = draw.Draw(1024, 1024, draw.FloatRect(2.0))
    bud = MultiBuddhabrot(8, 10, [20,200,2000], 1024, 1024)
#    bud.pickle("man.pck")
#    bud = cPickle.load(open("man2.pck", "rb"))
    d = bud.draw()
    d.save("/vagrant/man.png")
    

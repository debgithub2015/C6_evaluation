import numpy as np
from math import ceil, floor, sqrt

class CUBE:
  def __init__(self, fname):
    f = open(fname, 'r')
    for i in range(2): f.readline() # echo comment
    tkns = f.readline().split() # number of atoms included in the file followed by the position of the origin of the volumetric data
    self.natoms = int(tkns[0])
    self.origin = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
# The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector.
    tkns = f.readline().split() #
    self.NX = int(tkns[0])
    self.X = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
    tkns = f.readline().split() # 
    self.NY = int(tkns[0])
    self.Y = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
    tkns = f.readline().split() #
    self.NZ = int(tkns[0])
    self.Z = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
# The last section in the header is one line for each atom consisting of 5 numbers, the first is the atom number, second (?), the last three are the x,y,z coordinates of the atom center.
    self.atoms = []
    for i in range(self.natoms):
      tkns = f.readline().split()
      self.atoms.append([tkns[0], tkns[2], tkns[3], tkns[4]])
# Volumetric data
    self.data = np.zeros((self.NX,self.NY,self.NZ))
    i=0
    for s in f:
      for v in s.split():
        self.data[i/(self.NY*self.NZ), (i/self.NZ)%self.NY, i%self.NZ] = float(v)
        i+=1
    if i != self.NX*self.NY*self.NZ: raise NameError, "FSCK!"

  def dump(self, f):
# output Gaussian cube into file descriptor "f".
# Usage pattern: f=open('filename.cube'); cube.dump(f); f.close()
    print >>f, "CUBE file\ngenerated by piton _at_ erg.biophys.msu.ru"
    print >>f, "%4d %.6f %.6f %.6f" % (self.natoms, self.origin[0], self.origin[1], self.origin[2])
    print >>f, "%4d %.6f %.6f %.6f"% (self.NX, self.X[0], self.X[1], self.X[2])
    print >>f, "%4d %.6f %.6f %.6f"% (self.NY, self.Y[0], self.Y[1], self.Y[2])
    print >>f, "%4d %.6f %.6f %.6f"% (self.NZ, self.Z[0], self.Z[1], self.Z[2])
    for atom in self.atoms:
      print >>f, "%s %d %s %s %s" % (atom[0], 0, atom[1], atom[2], atom[3])
    for ix in xrange(self.NX):
      for iy in xrange(self.NY):
         for iz in xrange(self.NZ):
            print >>f, "%.5e " % self.data[ix,iy,iz],
            if (iz % 6 == 5): print >>f, ''
         print >>f,  ""

  def mask_sphere(self, R, Cx,Cy,Cz):
# produce spheric volume mask with radius R and center @ [Cx,Cy,Cz]
# can be used for integration over spherical part of the volume
    m=0*self.data
    for ix in xrange( int(ceil((Cx-R)/self.X[0])), int(floor((Cx+R)/self.X[0])) ):
      ryz=sqrt(R**2-(ix*self.X[0]-Cx)**2)
      for iy in xrange( int(ceil((Cy-ryz)/self.Y[1])), int(floor((Cy+ryz)/self.Y[1])) ):
          rz=sqrt(ryz**2 - (iy*self.Y[1]-Cy)**2)
          for iz in xrange( int(ceil((Cz-rz)/self.Z[2])), int(floor((Cz+rz)/self.Z[2])) ):
              m[ix,iy,iz]=1
    return m


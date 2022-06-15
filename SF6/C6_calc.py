from numpy import *
import numpy as np
import Dens_tool as DT
import sys

def kF(n):
    return ((3*pi**2.*n)**(1./3.))


def omega_pl_n(n):
    return sqrt(4*pi*n)
    

#def omega_pl_rs(rs):
#  n= 3.0/(4.0*pi*rs**3)
#  return sqrt(4*pi*n)

cutoff = 1e-16


class C6:
    def __init__(self,filename,version=2):
        self.filename = filename
        self.nspin = 1  
        self.version = version

    def read_from_cube(self):
# Cube files are assumed to be stored in (Hartree) atomic units
        import load_cube
        file = self.filename
        cube = load_cube.CUBE(file)
        Nx = cube.NX
        Ny = cube.NY
        Nz = cube.NZ
        a1 = cube.X*Nx
        a2 = cube.Y*Ny
        a3 = cube.Z*Nz
        unitcell = zeros((3,3))
        unitcell[0,:] = a1
        unitcell[1,:] = a2
        unitcell[2,:] = a3
        self.unitcell = unitcell
        print unitcell
        self.old_dens = cube.data

    def compute_base(self):
        unitcell=self.unitcell

        a1=unitcell[0]
        a2=unitcell[1]
        a3=unitcell[2]
#        V_UnitCell=abs(dot(a1,cross(a2,a3)))
        V_UnitCell = np.linalg.det(unitcell)
        print "Volume of the unitcell=",V_UnitCell
        print "Bravais lattice: "
        print "a1: ",a1
        print "a2: ",a2
        print "a3: ",a3

        # Calculate the inverse unitcell
        Ma1=2*pi*cross(a2,a3)/V_UnitCell
        Ma2=2*pi*cross(a3,a1)/V_UnitCell
        Ma3=2*pi*cross(a1,a2)/V_UnitCell

        self.a1=a1
        self.a2=a2
        self.a3=a3

        I_UnitCell=np.zeros((3,3))

        I_UnitCell[0]=[Ma1[0], Ma1[1], Ma1[2]]
        I_UnitCell[1]=[Ma2[0], Ma2[1], Ma2[2]]
        I_UnitCell[2]=[Ma3[0], Ma3[1], Ma3[2]]
        print "Inverse bravais lattice: "
        print "b1: ",I_UnitCell[0]
        print "b2: ",I_UnitCell[1]
        print "b3: ",I_UnitCell[2]

        if self.nspin==1:
            self.chgR=self.old_dens  #/V_UnitCell
        else:
            self.chgR=self.old_dens  #/V_UnitCell
            self.chgR_up=self.old_dens_up  #/V_UnitCell
            self.chgR_down=self.old_dens_down  #/V_UnitCell

        self.nga1,self.nga2,self.nga3=self.chgR.shape

        self.bx=DT.DerivativeXYZ(self.chgR,I_UnitCell,("x"))
        self.by=DT.DerivativeXYZ(self.chgR,I_UnitCell,("y"))
        self.bz=DT.DerivativeXYZ(self.chgR,I_UnitCell,("z"))

        mag_grad = ((self.bx.real)**2.+ \
                                                (self.by.real)**2.+ \
                                                (self.bz.real)**2)
        n = abs(self.chgR) + cutoff
        self.n = n
        q0s = DT.evaluate_q0_nos(n, mag_grad,vers=self.version)
        self.grid_points = np.product(self.chgR.shape)
        self.dV = V_UnitCell/self.grid_points
        self.eGap = 9*q0s**2/(8*pi)

        print "ns",sum(sum(sum(self.chgR)))*V_UnitCell/float(self.grid_points)
        print "grid points", self.grid_points
        print "size new dens", self.chgR.shape


    def C6_num(self):
        def alpha2(u):
            chiArray = self.n/(u**2+self.eGap**2)*(self.n>cutoff)
            return (sum(chiArray)*self.dV)**2
        from scipy.integrate import  quad
        C6=quad(alpha2,0.0001,800)[0]*3/pi
        return C6

atoms = C6(sys.argv[1],version=2)
atoms.read_from_cube()
atoms.compute_base()
print 'C6',atoms.C6_num()

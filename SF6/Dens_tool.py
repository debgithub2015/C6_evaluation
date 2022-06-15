# Based on code by Jesper Kleis, Elisa Londero, K. Berland, P. Hyldgaard 
from numpy import *

cutoff = 1e-16

# evaluates the fermi velocity in a.u.
def vf(n):
    return((3.*pi**2.*n)**(1./3.) )

# evaluates the fermi wavevector in a.u. 
def kf(n):
    return ((3*pi**2.*n)**(1./3.))

a_bohr=0.529177

# returns the exchange part of the LDA energy density 
def eps_LDA_x(n):
    n=abs(n)
    return(-3./(4.*pi)*(3.*pi*pi*n)**(1./3.))

def q_VWN_LDA (A,p,c,d,x):
    Xp=p**2.+c*p+d
    Xx=x**2.+c*x+d
    Qcd=sqrt(4.*d-c**2.)    
    return(A*(log(x**2./Xx)+2.*c/Qcd*arctan(Qcd/(2.*x+c))-
        c*p/Xp*(log((x-p)**2./Xx)+2.*(c+2.*p)/Qcd*arctan(Qcd/(2.*x+c)))))

# Returns the LDA correlation energy in a.u. 
# Implementation: Vosko-Wilk-Nusair (1980)
# J. Phys. 58, 1200 (1980)

def eps_c_VWN_LDA (n_up,n_down):
    n=n_up+n_down
    x=(3./(4.*pi*n))**(1./6.)
    mzeta=(n_up-n_down)/n
#  print x
    mLambda=q_VWN_LDA(0.0310907,-0.10498,3.72744,12.9352,x)
#  print mLambda
    mlambda=q_VWN_LDA(0.01554535,-0.325,7.06042,18.0578,x)
    alpha=q_VWN_LDA(-1./(6.*pi)**2.,-0.0047584,1.13107,13.0045,x)
    y=9./8.*(1.+mzeta)**(4./3.)+9./8.*(1.-mzeta)**(4./3.)-9./4.
#  print y
    h=4.*(mlambda-mLambda)/(9.*(2.**(1./3.)-1.)*alpha)-1.
    
    return(mLambda+alpha*y*(1.+h*mzeta**4.))

def e_PW92_LDA (r,t,u,v,w,x,y,p):
    return(-2.*t*(1.+u*r)*log(1.+1./(2.*t*(v*sqrt(r)+w*r+x*r**(3./2.)+y*r**(p+1.)))))

# Returns the LDA correlation energy in a.u. 
# Implementation: Perdew and Wang (1992)
# Phys. Rev. B 45, 13244 (1992)

def eps_c_PW92_LDA (n_up,n_down):
    c=1.7099210
    n=abs(n_up+n_down)
    r=(3./(4.*pi*n))**(1./3.)
    zeta=(n_up-n_down)/n
    wz=((1.+zeta)**(4./3.)+(1.-zeta)**(4./3.)-2.)/(2.**(4./3.)-2.)
    res=e_PW92_LDA(r,0.031091,0.21370,7.5957,3.5876,1.6382,0.49294,1.)*(1.-wz*zeta**4.)
    res=res+e_PW92_LDA(r,0.015545,0.20548,14.1189,6.1977,3.3662,0.62517,1.)*wz*zeta**4.
    res=res-e_PW92_LDA(r,0.016887,0.11125,10.357,3.6231,0.88026,0.49671,1.)*wz*(1.-zeta**4.)/c
    return(res)

    
#GetSpinDownGrid ( self )
#GetUpDownGrid ( self )

# evaluates q0 in the generslized geometry formalism
# gradn2 is the squared gradient of the electron density
# n at a point r
# Zab=-0.8491        # vdW-DF1
# Zab=-1.887         # vdW-DF2
# iFx=rPW86 # vdW-sDF/vdW-DF3

mu_C09x=0.0617
kappa_C09x=1.245
alfa=0.0483
def F_C09x_x(s2):
    return 1+mu_C09x*s2*exp(-alfa*s2)+kappa_C09x*(1-exp(-alfa*s2/2.0))

a_PW86=0.0864
b_PW86=14
c_PW86=0.2
def F_PW86_x(s2):
    return (1 + 15*a_PW86*s2 + b_PW86*(s2)**2. +c_PW86*s2**3)**(1.0/15.0)

mu_B88=0.2743
c_B88=2**(4.0/3.0)*(3*pi*pi)**(1.0/3.0)
beta_B88=9*mu_B88*(6.0/pi)**(1./3)/(2*c_B88)
def F_B88_x(s):
    return 1+mu_B88*s**2/(1+beta_B88*s*arcsinh(c_B88*s))

a_PW86_refit=0.1234
b_PW86_refit=17.33
c_PW86_refit=0.163
def F_PW86_refit_x(s2):
    return (1 + 15*a_PW86_refit*s2 + b_PW86_refit*(s2)**2. +c_PW86_refit*(s2)**3)**(1.0/15.0)

Zab_v1=-0.8491      # vdW-DF1
Zab_v2=-1.887       # vdW-DF2
def evaluate_iFx (s2, vers):
    if vers == 1:
        val=1.-(Zab_v1/9.)*(s2)
    elif vers == 2:
       val=1.-(Zab_v2/9.)*(s2)
    return(val)

def evaluate_q0_old (n, gradn2,Zab):
    if n < 1e-16:
        val=100000000000
    else:
        eps_xc_LDA=eps_c_PW92_LDA(n/2.,n/2.)+eps_LDA_x(n)
        eps_0xc=eps_xc_LDA-eps_LDA_x(n)*(Zab/9.)*(gradn2/(2.*kf(n)*n)**2.)
        val=kf(n)*eps_0xc/eps_LDA_x(n)
    return(val)

def evaluate_q0_nos(n, gradn2, vers):
    svals=(gradn2/(2.*kf(n)*n)**2.)
    eps_0sxc=eps_c_PW92_LDA(n/2.,n/2.)
    eps_0sxc=eps_0sxc+eps_LDA_x(n)*evaluate_iFx(svals, vers)
    val=kf(n)*eps_0sxc/eps_LDA_x(n)
    return(val)

def evaluate_q0_spin (n_up, n_down, gradn_up2, gradn_down2, vers):
    n=abs(n_up+n_down)
    n_up = (n_up>cutoff)*n_up + (n_up<cutoff)*cutoff
    n_down = (n_down>cutoff)*n_down + (n_down<cutoff)*cutoff
 
    svals_up = ((2.)**(-2/3.))*(gradn_up2/(2.*kf(n_up)*n_up)**2.)
    svals_down = ((2.)**(-2/3.))*(gradn_down2/(2.*kf(n_down)*n_down)**2.)
    eps_0sxc=eps_c_PW92_LDA(n_up,n_down)
    eps_0sxc = eps_0sxc+0.5*eps_LDA_x(2.*n_up)*evaluate_iFx(svals_up, vers)*2*n_up/n
    delta = 0.5*eps_LDA_x(2.*n_down)*evaluate_iFx(svals_down, vers)*2*n_down/n
    eps_0sxc = eps_0sxc+ delta
    val=kf(n)*eps_0sxc/eps_LDA_x(n)
    return(val)
 
# evaluates the eps_uq relevant for the 2003 layered vdW-DF
def eps_q(n, q_perp, u):
# n is in elec/[Ang^3]
    n=abs(n) 

    vfn=(3.*pi**2.*n)**(1./3.) # atomic units
    omega_p2=4.*pi*n     # atomic units
    return (1.+omega_p2/(u**2.+vfn**2.*q_perp**2./3.+q_perp**4./4.))

def DerivativeXYZ(array,reciprocal,componentlist,kpoint=[0,0,0]):
    """Returns the derivative of an array

    This method returns the derivative of an array. Note that the method
    implicitly assumes the array to be periodically repeated. However,
    if the spatial variation of the array can be written on the form
    'exp(ikr)*array' the derivative may still be obtained by specifying a
    'kpoint' given in cartesiancoordinates.

    **An example:** 

    To obtain the x,y derivative of array:

    'DerivativeXYZ(array,reciprocal,("x","y"))'

    where 'reciprocal' is an array representing the reciprocal unit cell 
    of the space in which 'array' is defined. 
    """
    from numpy import fromfunction

    # Using the following procedure (with i=x,y,z)
    # di f = sum_G prefactor*c_G*exp(i*(k+G)*r), where
    # prefactor = i*G_i

    # Generates a list with numbers 0,1,...,range/2,-(range/2-1),...,-1
    indexfunction=lambda i,length:(i+(length/2-1))%length-(length/2-1)

    components={'x':0,'y':1,'z':2}
    rec_coordinates=CoordinateArrayFromUnitVectors(array.shape,reciprocal,kpoint,indexfunction)
    prefactor=1.0

    for component in componentlist:
        prefactor=-1.0j*rec_coordinates[components[component]]*prefactor
    # Finding c_G
    c_G=InverseFFT(array)
    return FFT3D(prefactor*c_G)

def CoordinateArrayFromUnitVectors(shape,gridunitvectors,origin=[0,0,0],indexfunction=lambda i,length:i):
    """
    This method can be used to obtain an array representing the coordinates
    of a space defined by 'gridunitvecors'. 'gridunitvectors' is in turn a
    list containing the vectors defining the cells of the grid, i.e. the
    vectors between neighboring grid points. These vectors are spanned
    according to the specified shape. 

    'origin' -- specifies the origin of the returned coordinate array. 
    
    'indexfunction' -- is a lambda expression that defines the indices 
    with which each of the specified gridunitvectors are to be multiplied. 
    'indexfunction' must take two arguments, 'i' and 'length' - default
    is 'lambda i,length:i'. During exection the input index 'i' will run 
    over the interval 0,1,..., 'length' -1.

    **An Example**

    To obtain a coordinate array of shape (10,10) with 
    'gridunitvectors' =[[2,0],[0,1]] and the origin at [10,0] use:

    'CoordinateArrayFromUnitVectors((10,10),[[2,0],[0,1],[10,0])'

    Note that the output array will be of shape 
    (< *dimension* > ,  < *spatialcoordinates* >).
    """
    from numpy import add,fromfunction,array,asarray
    coordinatelist=[]
    gridunitvectors=asarray(gridunitvectors)
    # Looping over the dimensionality of the vectors
    for dim in range(gridunitvectors.shape[1]):
        coordinates=origin[dim]
        # Contribution from each unitvector
        for nunitvector in range(gridunitvectors.shape[0]):
            # Finding the indices from which the coordinate grid
            # is spanned
            indices=map(lambda i,f=indexfunction,l=shape[nunitvector]:f(i,l),range(shape[nunitvector]))
            coordinatefunc=lambda i,v=gridunitvectors[nunitvector,dim]:i*v
            coordinates=add.outer(coordinates,map(coordinatefunc,indices))
        coordinatelist.append(coordinates)
    return array(coordinatelist)

def InverseFFT(array):
    """Returns the inverse FFT of an array

    This method can be used to obtain the inverse FFT of an array
    """
    from numpy import fft
    dim=array.shape
    for i in range(len(dim)):
        array=fft.ifft(array,dim[i],axis=i)
    return array

def FFT3D(array):
    """Returns the FFT of a three dimensional array
    
    This method can be used to obtain the FFT of a three dimensional array.
    """
    from numpy import fft
    N1,N2,N3=array.shape
    return fft.fft(fft.fft2(array,(N1,N2),axes=(0,1)),N3,axis=2)

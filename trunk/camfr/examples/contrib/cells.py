#!/usr/bin/env python

####################################################################
#
# Code to generate various PhC lattices.
# Contributed by K.C. Huang (kch23@mit.edu).
#
####################################################################

from Numeric import *
from camfr import *
from cmath import *

def squares_GM(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    ex = Expression()

    z_max = a*sqrt(2.).real*0.5
    x_max = a*sqrt(2.).real*0.5

    if abs(r*sqrt(2))>=z_max/2:
        for i in range(slices):
            z = z_max/slices*(i+0.5)

            lower = 0
            upper = 0

            if abs(z) <=abs(r)*sqrt(2.).real:
                x_lower = sqrt(2.).real*abs(r)-z
                lower = 1
                if abs(x_lower) <= 0.:
                    x_lower = 0.
                    lower = 0
               
            if abs(z_max-z)<=abs(r)*sqrt(2.).real:
                x_upper = sqrt(2.).real*r-(z_max-z)
                upper = 1
                if abs(x_upper) <=0.:
                    x_upper = 0.
                    upper = 0
               
            if lower and not upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))

            if lower and upper: 
                part = Slab(rod(x_lower)+ambient(x_max-x_lower-x_upper)+rod(x_upper))

            if not lower and not upper:
                part = Slab(ambient(x_max))

            if not lower and upper:
                part = Slab(ambient(x_max-x_upper)+rod(x_upper))

            waveguides.append(part)
            ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))

    else:
        R = abs(r*sqrt(2.))

        for i in range(slices):
            z = R/slices*(i+0.5)
            x = R-z
            part = Slab(rod(x)+ambient(x_max-x))            

            waveguides.append(part)
            ex.add(part(R/slices))

        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(z_max-2.*R))

        for i in range(slices):
            z = R/slices*(i+0.5)
            x = z
            part = Slab(ambient(x_max-x)+rod(x))

            waveguides.append(part)
            ex.add(part(R/slices))        

        for i in range(len(waveguides)-1,len(waveguides)-slices-1,-1):
            ex.add(waveguides[i](R/slices))

        ex.add(waveguides[len(waveguides)-slices-1](abs(z_max-2.*R)))

        for i in range(slices-1,0-1,-1):
            ex.add(waveguides[i](R/slices))
        
    xz[0] = abs(x_max)
    xz[1] = abs(z_max)

    return ex

def circles_GX(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    x_max = 0.5*a
    z_max = 0.5*a
    
    ex = Expression()

    if r>z_max:
        for i in range(slices):
            z = z_max/slices*(i+0.5)
            
            lower = 0
            
            if abs(z)<=abs(r):
                x_lower = sqrt(r*r-z*z)
                lower = 1
                if abs(x_lower)<=0.:
                    x_lower = 0.
                    lower = 0
        
            if lower:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))
            else:
                part = Slab(ambient(x_max))

            waveguides.append(part)
            ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))

    else:
        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-z*z))

            part = Slab(rod(x)+ambient(x_max-x))
            waveguides.append(part)
            ex.add(part(r/slices))
            
        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(2.*z_max-2.*r))

        for i in range(len(waveguides)-2,-1,-1):
            ex.add(waveguides[i](r/slices))

    xz[0] = x_max
    xz[1] = abs(z_max)

    return ex

def circles_GM(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    ex = Expression()

    z_max = a*sqrt(2.).real*0.5
    x_max = a*sqrt(2.).real*0.5

    if r>=z_max/2:
        for i in range(slices):
            z = z_max/slices*(i+0.5)

            lower = 0
            upper = 0
            
            if abs(z) <=abs(r):
                x_lower = sqrt(r*r-z*z)
                lower = 1
                if abs(x_lower)<=0.:
                    x_lower = 0.
                    lower = 0
            
            if abs(z_max-z)<=abs(r):
                x_upper = sqrt(r*r-(z_max-z)*(z_max-z))
                upper = 1
                if abs(x_upper)<=0.:
                    x_upper = 0.
                    upper = 0
                
            if lower and not upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))
                
            if lower and upper: 
                part = Slab(rod(x_lower)+ambient(x_max-x_lower-x_upper)+rod(x_upper))
           
            if not lower and not upper:
                part = Slab(ambient(x_max))

            if not lower and upper:
                part = Slab(ambient(x_max-x_upper)+rod(x_upper))
                
            waveguides.append(part)
            ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))

    else:
        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-z*z))
            part = Slab(rod(x)+ambient(x_max-x))            

            waveguides.append(part)
            ex.add(part(r/slices))

        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(z_max-2.*r))

        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r**2-(r-z)**2))
            
            part = Slab(ambient(x_max-x)+rod(x))

            waveguides.append(part)
            ex.add(part(r/slices))        

        for i in range(len(waveguides)-1,len(waveguides)-slices-1,-1):
            ex.add(waveguides[i](r/slices))

        ex.add(waveguides[len(waveguides)-slices-1](abs(z_max-2.*r)))

        for i in range(slices-1,0-1,-1):
            ex.add(waveguides[i](r/slices))
            
    xz[0] = abs(x_max)
    xz[1] = abs(z_max)

    return ex

def squares_GX(waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    x_max = 0.5*a
    z_max = 0.5*a
    
    ex = Expression()
    
    part1 = Slab(ambient(0.5*a - r) + rod(r))
    part2 = Slab(ambient(0.5*a))

    waveguides.append(part1)
    waveguides.append(part2)
 
    xz[0] = x_max
    xz[1] = abs(z_max)
   
    ex.add(part1(2*r))
    ex.add(part2(a-2*r))

    return ex

def triangular_lattice_circles_GM(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    z_max = abs(sqrt(3.)*a*0.5)
    x_max = a*0.5
    
    ex = Expression()

    if r>z_max*0.5:
        for i in range(slices):
            z = z_max/slices*(i+0.5)
            
            lower = 0
            upper = 0
            
            if abs(z)<=abs(r):
                x_lower = sqrt(r*r-z*z)
                lower = 1
                if abs(x_lower)<=0.:
                    x_lower = 0.
                    lower = 0

            dz = z_max-z

            if abs(dz)<=abs(r):
                x_upper = sqrt(r*r-dz*dz)
                upper = 1
                if abs(x_upper)<=0.:
                    x_upper = 0.
                    upper = 0

            if lower and not upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))

            if lower and upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower-x_upper)+rod(x_upper))
	 
            if not lower and not upper:
                part = Slab(ambient(x_max))
                
            if not lower and upper:
                part = Slab(ambient(x_max-x_upper)+rod(x_upper))

        waveguides.append(part)
        ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))
            
    else:
        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-z*z))
            part = Slab(rod(x)+ambient(x_max-x))            

            waveguides.append(part)
            ex.add(part(r/slices))

        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(z_max-2.*r))

        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-(r-z)**2))
            part = Slab(ambient(x_max-x)+rod(x))

            waveguides.append(part)
            ex.add(part(r/slices))        

        for i in range(len(waveguides)-1,len(waveguides)-slices-1,-1):
            ex.add(waveguides[i](r/slices))

        ex.add(waveguides[len(waveguides)-slices-1](abs(z_max-2.*r)))

        for i in range(slices-1,0-1,-1):
            ex.add(waveguides[i](r/slices))        
        
    xz[0] = x_max
    xz[1] = abs(z_max)

    return ex

def triangular_lattice_circles_GK(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    z_max = 0.5*a
    x_max = sqrt(3)*0.5*a
    
    ex = Expression()

    if r>z_max*0.5:
        for i in range(slices):
            z = z_max/slices*(i+0.5)

            lower = 0
            upper = 0
       
            if abs(z)<=abs(r):
                x_lower = sqrt(r*r-z*z)
                lower = 1
                if abs(x_lower)<=0.:
                    x_lower = 0
                    lower = 0

            dz = z_max-z

            if abs(dz)<=abs(r):
                x_upper = sqrt(r*r-dz*dz)
                upper = 1
                if abs(x_upper)<=0.:
                    x_upper = 0
                    upper = 0

            if lower and not upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))

            if lower and upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower-x_upper)+rod(x_upper))
	 
            if not lower and not upper:
                part = Slab(ambient(x_max))

            if not lower and upper:
                part = Slab(ambient(x_max-x_upper)+rod(x_upper))

            waveguides.append(part)
            ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))
            
    else:
        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-z*z))
            part = Slab(rod(x)+ambient(x_max-x))            

            waveguides.append(part)
            ex.add(part(r/slices))

        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(z_max-2.*r))

        for i in range(slices):
            z = r/slices*(i+0.5)
            x = abs(sqrt(r*r-(r-z)**2))
            part = Slab(ambient(x_max-x)+rod(x))

            waveguides.append(part)
            ex.add(part(r/slices))        

        for i in range(len(waveguides)-1,len(waveguides)-slices-1,-1):
            ex.add(waveguides[i](r/slices))

        ex.add(waveguides[len(waveguides)-slices-1](abs(z_max-2.*r)))

        for i in range(slices-1,0-1,-1):
            ex.add(waveguides[i](r/slices))
        
    xz[0] = abs(x_max)
    xz[1] = abs(z_max)

    return ex

def triangular_lattice_squares_GM(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    z_max = abs(sqrt(3.)*a*0.5)
    x_max = a*0.5
    
    ex = Expression()

    R = abs(sqrt(2.)*r)

    if R>=z_max*0.5:
        for i in range(slices):
            z = z_max/slices*(i+0.5)

            lower = 0
            upper = 0
       
            if abs(z) <=abs(r)*sqrt(2.).real:
                x_lower = sqrt(2.).real*r-z
                lower = 1
                if abs(x_lower) <= 0.:
                    x_lower = 0.
                    lower = 0
                
            if abs(z_max-z)<=abs(r)*sqrt(2.).real:
                x_upper = sqrt(2.).real*r-(z_max-z)
                upper = 1
                if abs(x_upper) <=0.:
                    x_upper = 0.
                    upper = 0

            if lower and not upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower))

            if lower and upper:
                part = Slab(rod(x_lower)+ambient(x_max-x_lower-x_upper)+rod(x_upper))
	 
            if not lower and not upper:
                part = Slab(ambient(x_max))

            if not lower and upper:
                part = Slab(ambient(x_max-x_upper)+rod(x_upper))
                
            waveguides.append(part)
            ex.add(part(z_max/slices))

        for i in range(len(waveguides)-1,0-1,-1):
            ex.add(waveguides[i](z_max/slices))

    else:
        for i in range(slices):
            z = R/slices*(i+0.5)
            x = R-z
            part = Slab(rod(x)+ambient(x_max-x))            

            waveguides.append(part)
            ex.add(part(R/slices))

        part = Slab(ambient(x_max))
        waveguides.append(part)
        ex.add(part(z_max-2.*R))

        for i in range(slices):
            z = R/slices*(i+0.5)
            x = z
            part = Slab(ambient(x_max-x)+rod(x))

            waveguides.append(part)
            ex.add(part(R/slices))        

        for i in range(len(waveguides)-1,len(waveguides)-slices-1,-1):
            ex.add(waveguides[i](R/slices))

        ex.add(waveguides[len(waveguides)-slices-1](abs(z_max-2.*R)))

        for i in range(slices-1,0-1,-1):
            ex.add(waveguides[i](R/slices))        

    xz[0] = x_max
    xz[1] = abs(z_max)

    return ex

def triangular_lattice_squares_GK(slices,waveguides,ambient,rod,r,a,xz):
    # define half the unit cell
    z_max = 0.5*a
    x_max = sqrt(3)*0.5*a

    ex = Expression()

    if (r<0.25*a):
        part1 = Slab(rod(r)+ambient(abs(x_max)-r))
        part2 = Slab(ambient(abs(x_max)))
        part3 = Slab(ambient(abs(x_max)-r)+rod(r))

        waveguides.append(part1)
        waveguides.append(part2)
        waveguides.append(part3)

        ex.add(part1(r))
        ex.add(part2(abs(z_max)-2*r))
        ex.add(part3(2*r))
        ex.add(part2(abs(z_max)-2*r))
        ex.add(part1(r))
        
    if (r==0.25*a):
        part1 = Slab(rod(r)+ambient(abs(x_max)-r))
        part2 = Slab(ambient(abs(x_max)-r)+rod(r))

        waveguides.append(part1)
        waveguides.append(part2)

        ex.add(part1(r))
        ex.add(part2(2*r))
        ex.add(part1(r))
        
    if (r>0.25*a)&(r<abs(x_max)*0.5):
        part1 = Slab(rod(r)+ambient(abs(x_max)-r))
        part2 = Slab(rod(r)+ambient(abs(x_max)-2*r)+rod(r))
        part3 = Slab(ambient(abs(x_max)-r)+rod(r))

        waveguides.append(part1)
        waveguides.append(part2)
        waveguides.append(part3)

        ex.add(part1(abs(z_max)-r))
        ex.add(part2(2*r-abs(z_max)))
        ex.add(part3(2*(abs(z_max)-r)))
        ex.add(part2(2*r-abs(z_max)))
        ex.add(part1(abs(z_max)-r))
        
    if (r>abs(x_max)*0.5):
        part1 = Slab(rod(r)+ambient(abs(x_max)-r))
        part2 = Slab(rod(abs(x_max)))
        part3 = Slab(ambient(abs(x_max)-r)+rod(r))
        
        waveguides.append(part1)
        waveguides.append(part2)
        waveguides.append(part3)

        ex.add(part1(abs(z_max)-r))
        ex.add(part2(2*r-abs(z_max)))
        ex.add(part3(2*(abs(z_max)-r)))
        ex.add(part2(2*r-abs(z_max)))
        ex.add(part1(abs(z_max)-r))

    xz[0] = abs(x_max)
    xz[1] = abs(z_max)
    
    return ex

def flip_k(k,xz,a,mode):
    
    h0 = mode.field(Coord(0.,0.,0.)).H2()
    h1 = mode.field(Coord(xz[0],0.,xz[1])).H2()
    
    x = xz[0]
    z = xz[1]
    
    hyp = sqrt(x**2+z**2).real

    #project k onto the lattice vector
    kproj = k*z/hyp
    ka = kproj*a
    expk = exp(-1j*ka)
    
    kp = pi/a-k
    kpa = kp*a
    expkp = exp(1j*kpa)

    h1d = abs(h1-expk*h0)
    h1dp = abs(h1-expkp*h0)

    ratio = 2.*xz[1]/hyp
    if h1dp<h1d:
        k = 2.*pi/a/ratio-k

    print "RATIO = ",ratio
        
    return k

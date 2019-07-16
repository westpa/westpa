#!/usr/bin/env python

import numpy
from numpy import pi
from numpy import sin, cos, exp
from numpy.random import normal as random_normal

SEG_STATUS_COMPLETE = 2

PI = pi
A, B, C, x0 = 2, 10, 0.5, 1
coord_len = 21

coord_dtype = numpy.float32
sigma = 0.001**(0.5)
reflect_at = 10.0

coords = numpy.empty(coord_len, dtype=coord_dtype)


print(repr('inside Propagate'))
print(repr(coord_len), repr(coord_dtype))
print(repr(coords))
print()

#coords[0] = 8.0

f = open('odld.crd', 'r')
coords[0] = f.readline()
f.close

print(repr('coords[0] = '), repr(coords[0]))

twopi_by_A = 2*PI/A
half_B = B/2
sigma = sigma
gradfactor = sigma*sigma/2
all_displacements = numpy.zeros(coord_len, dtype=coord_dtype)

print(repr('A B C x0'))
print(repr(A), repr(B), repr(C), repr(x0))
print(repr(twopi_by_A), repr(half_B), repr(sigma), repr(gradfactor), repr(coord_len), repr(reflect_at))
print(repr(coord_len))
print(repr(coords))


for istep in range(1,coord_len):
    x = coords[istep-1]
    
    xarg = twopi_by_A*(x - x0)
    
    #print repr(istep), repr(x), repr(x0)
    #print repr(C)
    #print repr(x)
    #print repr(xarg)

    eCx = numpy.exp(C*x)
    eCx_less_one = eCx - 1.0
   
    all_displacements[istep] = displacements = random_normal(scale=sigma, size=(1,))
    grad = half_B / (eCx_less_one*eCx_less_one)*(twopi_by_A*eCx_less_one*sin(xarg)+C*eCx*cos(xarg))
    
    newx = x - gradfactor*grad + displacements
    newy = x - gradfactor*grad + displacements #KFW

    if reflect_at is not None:
        # Anything that has moved beyond reflect_at must move back that much

        # boolean array of what to reflect
        to_reflect = newx > reflect_at

        # how far the things to reflect are beyond our boundary
        reflect_by = newx[to_reflect] - reflect_at

        # subtract twice how far they exceed the boundary by
        # puts them the same distance from the boundary, on the other side
        newx[to_reflect] -= 2*reflect_by

    coords[istep] = newx
    print(repr(istep), repr(x), repr(newx), repr(coords[istep]))
    
print(repr('newx'))
print(repr(newx))
print(repr(coords))

f = open('odld.rst', 'w')
f.write('{:12.8f}'.format(coords[coord_len-1])+'\n')
f.close

f = open('pcoords.dat', 'w')
for element in coords.flat:
    f.write('{:12.8f}'.format(element)+'\n')
f.close

f = open('displacements.dat', 'w')
for element in all_displacements.flat:
    f.write('{:12.8f}'.format(element)+'\n')
f.close


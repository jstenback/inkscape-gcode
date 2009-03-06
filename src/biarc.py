#!/usr/bin/python
# -*- Mode: Python; tab-width: 4; python-indent: 4; indent-tabs-mode: nil -*- */
#    planar spline to arc conversion
#    Copyright 2007 Jeff Epler <jepler@unpythonic.net>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import math

class P:
    def __init__(self, x, y = None):
        if (y != None):
            self.x, self.y = .0 + x, .0 + y
        else:
            self.x, self.y = .0 + x[0], .0 + x[1]
    def __getitem__(self,index):
        assert index < 2
        if index == 0:
            return self.x
        if index == 1:
            return self.y
    def __setitem__(self,index,value):
        assert index < 2
        if index == 0:
            self.x = value
        if index == 1:
            self.y = value
    def __add__(self, other):
        return P(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return P(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        if isinstance(other, P):
            return self.x * other.x + self.y * other.y
        return P(self.x * other, self.y * other)
    __rmul__ = __mul__
    def __div__(self, other):
        return P(self.x / other, self.y / other)

    def mag(self):
        return math.hypot(self.x, self.y)

    def unit(self):
        h = self.mag()
        if h: return self / h
        else: return P(0,0)

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def rot(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)

        return P(self.x * c - self.y * s,
                 self.x * s + self.y * c)

    def angle(self):
        return math.atan2(self.y, self.x)

    def __repr__(self): return "<%.6f %.6f>" % (self.x, self.y)

def biarc(P0, TS, P4, TE, r):
    TS = TS.unit()
    TE = TE.unit()

    v = P0 - P4

    c = v*v
    b = 2*v*(r*TS+TE)
    a = 2*r*(TS*TE-1)

    if a == 0:
        raise ValueError, (a,b,c)

    discr = b*b-4*a*c
    if discr < 0:
        raise ValueError, (a,b,c,discr)

    disq = discr**.5
    beta1 = (-b - disq) / 2 / a
    beta2 = (-b + disq) / 2 / a
    #print (beta1, beta2)
    if beta1 > 0 and beta2 > 0:
        raise ValueError, (a,b,c,disq,beta1,beta2)
    beta = max(beta1, beta2)

    if beta < 0:
        raise ValueError, (a,b,c,disq,beta)

    alpha = beta * r
    ab = alpha + beta
    P1 = P0 + alpha * TS
    P3 = P4 - beta * TE
    P2 = (beta / ab)  * P1 + (alpha / ab) * P3

    #print "(%f %f %f %f)" % (a, b, c, discr)
    #print "(%f %f %f %f)" % (P2.x, P2.y, alpha, beta)
    return P1, P2, P3, alpha, beta

small = 1e-10
def unit(a,b):
    h = math.hypot(a, b)
    if h: return a/h, b/h
    else: return a, b

class Arc:
    def __init__(self, PS, PE, TS):
        """\
        r(PS, PE, TS) -> TS

        Given the initial location (PS) the desired location (PE), and
        the direction of travel at the initial location (TS), find the
        line or arc satisfying these conditions

        In the case that the target point is exactly in the opposite
        direction of the initial direction of motion, the tangency
        condition cannot be maintained.  In this case, the points are
        joined with a line that violates the tangency condition.

        Returns the direction of travel at the final location
        """

        TS = TS.unit()
        #print("(r %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f)"
        #      % (PS.x, PS.y, PE.x, PE.y, TS.x, TS.y))

        self.start = PS
        self.end = PE
        self.isLine = False

        delta = PE - PS

        den = 2 * (delta.y * TS.x - delta.x * TS.y)

        if abs(den) < small:  # Line
            self.isLine = True
            return

        self.r = -(delta.x**2 + delta.y**2) / den

        if abs(self.r) < 0.0002:  # Line
            self.isLine = True
            return

        self.i = TS.y * self.r
        self.j = -TS.x * self.r

        # XXX
        if (abs((self.start - self.end).mag()) < precision and
            r < precision):
            self.isLine = True
            return

        self.center = PS + P(self.i, self.j)

    def printGCode(self, output):
        """
        Output a line of g-code to perform this arc
        """

        if self.isLine:
            output.write("G1 X%.6f Y%.6f\n" % (self.end.x, self.end.y))
            #return unit(delta.x, delta.y)
        elif self.r < 0:   # Counterclockwise
            output.write("G3 I%.6f J%.6f X%.6f Y%.6f\n"
                         % (self.i, self.j, self.end.x, self.end.y))
            #return unit(-(PE.y-c.y), PE.x-c.x)
        else:       # Clockwise
            output.write("G2 I%.6f J%.6f X%.6f Y%.6f\n"
                         % (self.i, self.j, self.end.x, self.end.y))
            #return unit(PE.y-c.y, -(PE.x-c.x))

    def radiusDifference(self, p):
        if self.isLine:
            return 0

        return (self.center - p).mag() - abs(self.r)


def giarc(PS, TS, PE, TE, r=1):
    TS = TS.unit()
    TE = TE.unit()
    P1, P2, P3, alpha, beta = biarc(PS, TS, PE, TE, r)

    if 0:
        print "G0 X%.6f Y%.6f" % (P1.x, P1.y)
        print "G0 X%.6f Y%.6f" % (P2.x, P2.y)
        print "G4 P0"
        print "G0 X%.6f Y%.6f" % (P3.x, P3.y)
        print "G0 X%.6f Y%.6f" % (PE.x, PE.y)
        print "G0 X%.6f Y%.6f" % (PS.x, PS.y)

    return Arc(PS, P2, TS), Arc(P2, PE, P3 - P2)

class Spline:
    def __init__(self, P0, P1, P2, P3, n=2):
        # P0 is the start of the spline
        self.P0 = P0

        # P1 is the first control point
        self.P1 = P1

        # P2 is the second control point
        self.P2 = P2

        # P3 is the end of the spline
        self.P3 = P3

        # n is the initial number of sub-divisions
        self.n = n

        self.Q0 = 3*(P1-P0)
        self.Q1 = 3*(P2-P1)
        self.Q2 = 3*(P3-P2)

        # The maximum allowable divergence from the curve, measured as the
        # difference betwee the arc radius and the distance from the
        # center of the arc to the middle of a spline section that the arc
        # represents.
        self.max_divergence = 0.0005

        self.arcs = []

    # Calculate the the location of the point on the curve at t
    def calculateSplinePoint(self, t):
        return (self.P0 * (1 - t) ** 3 +
                3 * self.P1 * t * (1 - t) ** 2 +
                3 * self.P2 * t ** 2 * (1-t) +
                self.P3 * t ** 3)

    # Calculate the the tangent on the curve at t
    def calculateTangent(self, t):
        return (self.Q0 * (1 - t) ** 2 +
                2 * self.Q1 * t * (1 - t) +
                self.Q2 * t ** 2)

    # Create biarcs for a piece of the curve. The first arc will start
    # on the curve at tS, the second one will end at tE. If the curve
    # distance from the center of the arcs to the midpoint of the half
    # of the curve segment the arc represents is larger than
    # max_divergence, split the curve and repeat until the arcs stay
    # within the limit. This function returns the location and tangent
    # at tE
    def createArcs(self, PS, TS, tS, tE):
        PE = self.calculateSplinePoint(tE)
        TE = self.calculateTangent(tE)

        tDelta = tE - tS

        a1, a2 = giarc(PS, TS, PE, TE)

        PM1 = self.calculateSplinePoint(tS + (tDelta * 0.25))
        PM2 = self.calculateSplinePoint(tS + (tDelta * 0.75))

        if ((abs(a1.radiusDifference(PM1) ) > self.max_divergence) or
            (abs(a2.radiusDifference(PM2) ) > self.max_divergence)):
            tM = tS + (tDelta * 0.5)

            PM, TM = self.createArcs(PS, TS, tS, tM)
            self.createArcs(PM, TM, tM, tE)
        else:
            self.arcs.append(a1)
            self.arcs.append(a2)

        return PE, TE

    # Output the G-code for the arcs describing this curve.
    def printGCode(self, output):
        if len(self.arcs) == 0:
            PS = self.P0
            TS = self.P1 - self.P0
            tS = 0

            for i in range(0, self.n):
                tE = (i + 1) * 1. / self.n

                PS, TS = self.createArcs(PS, TS, tS, tE)
                tS = tE

        for arc in self.arcs:
            arc.printGCode(output)

if __name__ == '__main__':
    import sys
    print "G0X0Y0Z0"
    print "F100"
    s = Spline(P(0,0), P(.25,1), P(1,1), P(2,0), 2)
    s.max_divergence = 2
    s.printGCode(sys.stdout)
    print "G0X0Y0Z0"
    s = Spline(P(0,0), P(.25,1), P(1,1), P(2,0), 4)
    s.max_divergence = 2
    s.printGCode(sys.stdout)
    print "G0X0Y0Z0"
    s = Spline(P(0,0), P(.25,1), P(1,1), P(2,0), 8)
    s.max_divergence = 2
    s.printGCode(sys.stdout)
    print "M2"

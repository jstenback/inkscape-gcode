#!/usr/bin/env python 

'''
Copyright (C) 2009 glen.harris@middlegable.org
based on gcode.py (C) 2007 jst... 
based on gcode.py (C) 2007 hugomatic... 
based on dots.py (C) 2005 Aaron Spike, aaron@ekips.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Note that this still needs to be extended in the following ways:
1) Convert path curves internally into straight lines and/or arcs.  This is what
    the machineTolerance parameter is designed for.
'''
import inkex, simplestyle, simplepath, simpletransform, biarc
import os, math, copy, re, sys
import gettext
_ = gettext.gettext
def nolog(text):
    pass

logDebug = nolog
#logDebug = inkex.debug
logError = inkex.debug
logMessage = inkex.debug
logStage = inkex.debug
logWarning = inkex.debug
g_toolPath = "ToolPath"
option_filenameBase = ""
option_machineTolerance = 0.005
option_traverseZ = 1.0
option_defaultCutZ = -1.0
option_defaultFeedRate  = 1.0
option_returnToXYOrigin  = True
option_useCutVariables = True
option_useScaleVariables = False
option_reverseYAxis = True
option_outputExtension = "ngc"
epsilon = 0.00001

#Math/Trig/Line helpers
def isPointStraightLineDistanceLessThanTolerance(straightLine, point, toleranceSquared):
    """Calculate the closest distance between the point and any other point
    on the straight line (noting that the straight line has a begin and end)."""
    # Using algorithm from http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
    deltax = straightLine.destination.x - straightLine.source.x
    deltay = straightLine.destination.y - straightLine.source.y
    denom = deltax**2 + deltay**2
    if denom < epsilon:
        # The line segment extremely short, so just do a point test
        return isPointPointDistanceLessThanTolerance( straightLine.source, point, toleranceSquared)
    numx = (point.x - straightLine.source.x)*deltax
    numy = (point.y - straightLine.source.y)*deltay
    u = (numx + numy ) / denom
    if u < 0.0:
        # Closest point is the start of the line
        return isPointPointDistanceLessThanTolerance(straightLine.source, point, toleranceSquared)
    if u > 1.0:
        # Closest point is the finish of the line
        return isPointPointDistanceLessThanTolerance(straightLine.destination, point, toleranceSquared)
    closex = straightLine.source.x + u * deltax
    closey = straightLine.source.y + u * deltay
    return isPointPointDistanceLessThanTolerance(Point(closex,closey), point, toleranceSquared)
                                                 
def isPointPointDistanceLessThanTolerance(point1, point2, toleranceSquared):
    deltax = point1.x - point2.x
    deltay = point1.y - point2.y
    distanceSquared = deltax**2 + deltay**2
    if distanceSquared < toleranceSquared: 
        return True
    return False
 
def interpolateBezier(segment, position):
    if position == 0.0:
        return segment.source
    if position == 1.0:
        return segment.destination
    # Source is P0, SourceTangent is P1, DestinationTangent is P2 and Destination is P3
    p0 = (1-position)**3
    p1 = 3 * position * (1-position)**2 
    p2 = 3 * position **2 * (1-position) 
    p3 = position **3
    x = segment.source.x * p0 + segment.sourceTangent.x * p1 + segment.destinationTangent.x * p2 + segment.destination.x * p3
    y = segment.source.y * p0 + segment.sourceTangent.y * p1 + segment.destinationTangent.y * p2 + segment.destination.y * p3
    return Point( x, y)

def calculateLineSegmentIntersectionPoint(line1,line2):
    """Returns a Point if the two lines intersect at a single point, else None.  If the two
    lines are coincident, then return None"""
    logDebug("Intersection of %s x %s" % (line1, line2))
    xslope1 = line1.destination.x-line1.source.x
    yslope1 = line1.destination.y-line1.source.y
    xslope2 = line2.destination.x-line2.source.x
    yslope2 = line2.destination.y-line2.source.y
    xstartdelta = line1.source.x-line2.source.x
    ystartdelta = line1.source.y-line2.source.y
    denominator = (-xslope2*yslope1+xslope1*yslope2)
    logDebug("Denominator: %.3f" % (denominator))
    if math.fabs(denominator) >= epsilon :
        alpha1 = (xslope2*ystartdelta - yslope2*xstartdelta)/denominator
        alpha2 = (-yslope1*xstartdelta + xslope1*ystartdelta)/denominator
        logDebug("Alpha: %.3f, %.3f" % (alpha1, alpha2))
        returnPoint = Point(line1.source.x + alpha1*xslope1, line1.source.y + alpha1*yslope1)
        #checkPoint = Point(line2.source.x + alpha2*xslope2, line2.source.y + alpha2*yslope2)
        logDebug(returnPoint)
        return returnPoint
    else:
        # lines are basically parallel, so there could be zero or an infinite
        # number of intersection points
        return None

def calculateMidpoint(line1,line2):
    logDebug("%s midpoint %s" % (line1,line2))
    returnPoint = Point()
    returnPoint.translate( line1.destination )
    returnPoint.translate( line2.source )
    returnPoint.x /= 2
    returnPoint.y /= 2
    logDebug(returnPoint)
    return returnPoint


class AbstractElementReader(object):
    """
    Abstract definition of the readers.  They should take an xml element, and
    attempt to return a list of CutPaths that represent that element.  Most readers
    will return one CutPath, but some (such as Inkscape paths) may have to return 
    multiple CutPaths because each CutPath is continuous.
    
    Note that readers should call coordinateTransformer.convertPosition() to create
    Points that are in the correct output units.
    """
    def convertElementToPaths(self, element, coordinateTransformer):
        pass

class AbstractSegmentConverter(object):
    """
    Abstract definition of the converters.  They should take a single Segment
    and return a list of 1 or more segments.  If no changes are to take place
    then return None
    """
    def convertSegment(self, segment):
        pass
    
class CoordinateTransformer(object):
    """
    Helper class to return Points in output units.  Internally, this class manages
    a stack of transforms, scales it to get to mm or inches, then handles a user-specified 
    origin (in case you dont like (0,0) at the top left of the page).
    
    This class is created during the parsing of the svg tree, and then used
    repeatedly by the conversion classes which call convertPosition() and convertVector().
    
    Uses separate scales for x and y to allow for coordinate flips (e.g. flip y axis)
    """
    def __init__(self):
        self.transforms = []
        self.scalex = 1.0
        self.scaley = 1.0
        self.origin = Point(0.0,0.0)
    def setOriginInOutputUnits(self,point):
        self.origin = point
    def setScale(self,scalex, scaley):
        self.scalex = scalex
        self.scaley = scaley
    def appendTransform(self,transform):
        self.transforms.append(transform)
    def convertPosition(self,x,y):
        #x and y are in document units
        if len(self.transforms) > 0:
            point = biarc.P(x,y)
            for transform in self.transforms:
                simpletransform.applyTransformToPoint(transform, point)
            # point is now in document units, so convert to output units
            return Point(point.x/self.scalex - self.origin.x, point.y/self.scaley -self.origin.y)
        else:
            return Point(x/self.scalex - self.origin.x, y/self.scaley -self.origin.y)
    
class Point(object):
    """
    Basic point class used by segments etc.  Normally should be immutable, but
    normalise() etc are present for convenience.  This class should probably
    be extended to 3D in the future (if we start doing 2.5D conversions).
    Most points will be in final output units (i.e. mm or inches).
    """
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y
    def normalise(self):
        length = self.magnitude()
        self.x /= length
        self.y /= length
    def magnitude(self):
        return math.sqrt(self.x*self.x + self.y*self.y)
    def translate(self, offset):
        self.x += offset.x
        self.y += offset.y
    def __str__(self):
        return "[%.3f,%.3f]" % (self.x, self.y)

class AbstractSegment(object):
    def __init__(self, source, destination):
        self.source = source
        self.destination = destination

class StraightLineSegment(AbstractSegment):
    def __init__(self, source, destination):
        self.source = source
        self.destination = destination

class ArcSegment(AbstractSegment):
    def __init__(self, source, destination, center, clockwise):
        self.source = source
        self.destination = destination
        self.center = center
        self.clockwise = clockwise

class BezierSegment(AbstractSegment):
    def __init__(self, source, destination, sourceTangent, destinationTangent):
        self.source = source
        self.destination = destination
        self.sourceTangent = sourceTangent
        self.destinationTangent = destinationTangent
    def __str__(self):
        return "S%s D%s ST%s DT%s" % (self.source, self.destination, self.sourceTangent, self.destinationTangent)
     
class CutAttributes(object):
    """
    Just a collection of name-value pairs, but abstracted from a regular tuple
    for readability.
    """
    def __init__(self):
        self.attributes = {}
    def setAttribute(self, name, value):
        logDebug("Setting attribute %s to %s" % (name, value) )
        assert not self.attributes.has_key(name)
        self.attributes[name] = value
    def hasAttribute(self, name):
        return self.attributes.has_key(name)
    def getAttribute(self, name):
        assert self.attributes.has_key(name)
        return self.attributes[name]
    def getAttributeWithDefault(self, name, default):
        if self.attributes.has_key(name):
            return self.attributes[name]
        return default
               
class CutList(object):
    """
    This class is really just a collection of CutPaths, where all CutPaths share
    a particular set of attributes.  Typical attributes that will be attached to
    a CutList may include 'FeedRate', 'CutDepth', 'ToolSelection' etc.
    """
    def __init__(self):
        self.cutPaths = []
        self.attributes = CutAttributes()
    def appendCutPath(self, cutPath):
        self.cutPaths.append(cutPath)
    def getCutPaths(self):
        return self.cutPaths
    def hasCutPaths(self):
        return len(self.cutPaths) > 0

class CutPath(object):
    """
    This class represents a continuous cut (or tool path).  Note that it is irrelevant
    whether the path is closed or open - closed is just the degenerate case where
    the final point is the same as the first point.
    
    Each segment of the CutPath should be an AbstractSegment and may be straight, arc, bezier, etc etc.  
    
    Not sure what attributes will be attached to a CutPath yet.
    """
    def __init__(self):
        self.segments = []
        self.attributes = CutAttributes()
    def appendSegment(self, segment):
        if len(self.segments) > 0:
            # Ensure that the path is contiguous - note that this uses object identity,
            # but could maybe be relaxed to use position identity.  I strongly
            # recommend we keep with object identity
            lastSegment = self.segments[-1]
            assert segment.source == lastSegment.destination
        self.segments.append(segment)
    def hasSegments(self):
        return len(self.segments) > 0
    def getSegments(self):
        return self.segments


class GcodeOperationList(object):        
    def __init__(self):
        self.lines = []
    def appendInstruction(self, code, **arguments):
        line = code + " "
        argnames = arguments.keys()
        argnames.sort()
        for name in argnames : 
            value = arguments[name]
            if isinstance(value,basestring):
                line +=  "%c%s " % (name[0].upper(),arguments[name])
            else:
                assert isinstance(value,float)
                line +=  "%c%.4f " % (name[0].upper(),arguments[name])
        self.lines.append( line )
        
    def appendToolChange(self, tool):
        self.lines.append("%s M6" % tool)

    def appendEnd(self):
        self.lines.append("M2")
        
    def appendRapidMove(self, **arguments):
        self.appendInstruction("G00", **arguments)
        
    def appendFeedMove(self,**arguments):
        self.appendInstruction("G01",**arguments)

    def appendFeedArc(self,clockwise,**arguments):
        if clockwise:
            self.appendInstruction("G02",**arguments)
        else:
            self.appendInstruction("G03",**arguments)

    def appendComment(self, comment):
        self.lines.append( "(" + comment + ")" )
        
    def setUnitsToMillimeters(self):
        self.lines.append("G21 (All units in mm)")
    def setUnitsToInches(self):
        self.lines.append("G20 (All units in inches)")
    def defineVariable(self,index,defaultValue,comment):
        self.lines.append( "#%i=%.4f (%s)" % (index, defaultValue, comment) )


def calculateBounds(cutPath):
    assert cutPath.hasSegments()
    point = cutPath.segments[0].source
    minX = maxX = point.x
    minY = maxY = point.y
    for segment in cutPath.segments:
        point = segment.destination
        minX = min(point.x, minX)
        minY = min(point.y, minY)
        maxX = max(point.x, maxX)
        maxY = max(point.y, maxY)
    return StraightLineSegment( Point(minX,minY), Point(maxX,maxY))

class PathElementReader(AbstractElementReader):
    def convertElementToPaths(self, element, coordinateTransformer):
        assert element.tag == inkex.addNS("path","svg")
        returnPaths = []
        d = element.attrib.get('d')
        if d is None:
            logWarning("Ignoring path with no d attribute")
            return returnPaths
            
        p = simplepath.parsePath(d)
        currentPath = None
        currentPoint = None
        # Loop through each element of the line, building up a path representation
        for cmd, params in p:
            if cmd == 'M':
                if currentPath is not None:
                    #Terminate the last path, and start a new one
                    returnPaths.append(currentPath)
                currentPath = CutPath()
                currentPoint = coordinateTransformer.convertPosition(params[0],params[1] )
            elif cmd == 'L':
                # Lineto - somewhere inside the path
                assert currentPath is not None
                assert currentPoint is not None
                newPoint = coordinateTransformer.convertPosition(params[0],params[1] )
                currentPath.appendSegment( StraightLineSegment( currentPoint, newPoint ) )
                currentPoint = newPoint
            elif cmd == 'C':
                # Cubicto - p1, p2, p3
                assert currentPath is not None
                assert currentPoint is not None
                currentTangent = coordinateTransformer.convertPosition(params[0],params[1] )
                newTangent = coordinateTransformer.convertPosition(params[2],params[3] )
                newPoint = coordinateTransformer.convertPosition(params[4],params[5] )
                currentPath.appendSegment( BezierSegment( currentPoint, newPoint, currentTangent, newTangent ) )
                currentPoint = newPoint
            elif cmd == 'Z':
                # Close the path
                assert currentPath is not None
                assert currentPath.hasSegments()
                #newPoint = currentPath.segments[0].source
                #currentPath.appendSegment( StraightLineSegment( currentPoint, newPoint ) )
                returnPaths.append( currentPath )
                currentPath = None
                currentPoint = None
            else:
                logWarning(_("Path contains invalid line segments (%s) and will be ignored" % cmd ))
                return []
            
        # clean  up in case the path did not finish with a Z
        if currentPath is not None:
            returnPaths.append(currentPath)
            currentPath = None
            currentPoint = None
        return returnPaths

class NullElementReader(AbstractElementReader):
    """Do nothing - useful for groups because the overall algorithm will
    still iterate through children."""
    def convertElementToPaths(self, element, coordinateTransformer):
        return []

class LineElementReader(AbstractElementReader):
    def convertElementToPaths(self, element, coordinateTransformer):
        assert element.tag == inkex.addNS("line","svg")
        startPoint = coordinateTransformer.convertPosition(node.get('x1'),node.get('y1') )
        endPoint = coordinateTransformer.convertPosition(node.get('x1'),node.get('y1') )
        currentPath = CutPath()
        currentPath.appendSegment( StraightLineSegment( startPoint, endPoint ) )
        returnPaths = []
        returnPaths.append( currentPath )
        return returnPaths


class RectElementReader(AbstractElementReader):
    def convertElementToPaths(self, element, coordinateTransformer):
        assert element.tag == inkex.addNS("rect","svg")
        rx = element.get('rx')
        ry = element.get('ry')
        # if only one is set, then it is a radius
        if rx is None:
            rx = ry
        if ry is None:
            ry = rx
        if rx is None:
            assert ry is None
            rx = 0.0
            ry = 0.0
        else:
            rx = float(rx)
            ry = float(ry)
        x = float(element.get("x"))
        y = float(element.get("y"))
        width = float(element.get("width"))
        height = float(element.get("height"))
        logDebug("Found rect x=%.4f y=%.4f width=%.4f height=%.4f rx=%.4f ry=%.4f" % (x,y,width,height,rx,ry))
        returnPaths = []
        if (width > 0.0) and (height >= 0.0):
            if ( rx < width / 2.0 ) and (ry < height / 2.0 ):
                currentPath = CutPath()
                if  rx + ry > 0.0:
                    p1 = coordinateTransformer.convertPosition(x+rx,y)
                    p2 = coordinateTransformer.convertPosition(x+width-rx,y)
                    p3 = coordinateTransformer.convertPosition(x+width,y+ry)
                    p4 = coordinateTransformer.convertPosition(x+width,y+height-ry)
                    p5 = coordinateTransformer.convertPosition(x+width-rx,y+height)
                    p6 = coordinateTransformer.convertPosition(x+rx,y+height)
                    p7 = coordinateTransformer.convertPosition(x,y+height-ry)
                    p8 = coordinateTransformer.convertPosition(x,y+ry)
                    if abs(rx - ry) < epsilon:
                        # Use arcs for corners
                        c23 = coordinateTransformer.convertPosition(x+width-rx,y+ry)
                        c45 = coordinateTransformer.convertPosition(x+width-rx,y+height-ry)
                        c67 = coordinateTransformer.convertPosition(x+rx,y+height-ry)
                        c81 = coordinateTransformer.convertPosition(x+rx,y+ry)
                        currentPath.appendSegment(StraightLineSegment(p1, p2))
                        currentPath.appendSegment(ArcSegment(p2, p3, c23, False))
                        currentPath.appendSegment(StraightLineSegment(p3, p4))
                        currentPath.appendSegment(ArcSegment(p4, p5, c45, False))
                        currentPath.appendSegment(StraightLineSegment(p5, p6))
                        currentPath.appendSegment(ArcSegment(p6, p7, c67, False))
                        currentPath.appendSegment(StraightLineSegment(p7, p8))
                        currentPath.appendSegment(ArcSegment(p8, p1, c81, False))
                    else:
                        # Use biarcs for corners
                        p2t = coordinateTransformer.convertPosition(x+width-rx/2,y)
                        p3t = coordinateTransformer.convertPosition(x+width,y+ry/2)
                        p4t = coordinateTransformer.convertPosition(x+width,y+height-ry/2)
                        p5t = coordinateTransformer.convertPosition(x+width-rx/2,y+height)
                        p6t = coordinateTransformer.convertPosition(x+rx/2,y+height)
                        p7t = coordinateTransformer.convertPosition(x,y+height-ry/2)
                        p8t = coordinateTransformer.convertPosition(x,y+ry/2)
                        p1t = coordinateTransformer.convertPosition(x+rx/2,y)
                        currentPath.appendSegment(StraightLineSegment(p1, p2))
                        currentPath.appendSegment(BezierSegment(p2, p3, p2t, p3t))
                        currentPath.appendSegment(StraightLineSegment(p3, p4))
                        currentPath.appendSegment(BezierSegment(p4, p5, p4t, p5t))
                        currentPath.appendSegment(StraightLineSegment(p5, p6))
                        currentPath.appendSegment(BezierSegment(p6, p7, p6t, p7t))
                        currentPath.appendSegment(StraightLineSegment(p7, p8))
                        currentPath.appendSegment(BezierSegment(p8, p1, p8t, p1t))
                        
                else:
                    # Straight edges
                    p1 = coordinateTransformer.convertPosition(x,y)
                    p2 = coordinateTransformer.convertPosition(x+width,y)
                    p3 = coordinateTransformer.convertPosition(x+width,y+height)
                    p4 = coordinateTransformer.convertPosition(x,y+height)
                    currentPath.appendSegment(StraightLineSegment(p1, p2))
                    currentPath.appendSegment(StraightLineSegment(p2, p3))
                    currentPath.appendSegment(StraightLineSegment(p3, p4))
                    currentPath.appendSegment(StraightLineSegment(p4, p1))
                returnPaths.append(currentPath)
            else:
                logError("Invalid rounding on rectangle") 
        else:
            logError("Cannot output invalid sized rectangle")
        return returnPaths

class BezierStraightLineConverter(object):
    """
    Converts BezierSegments into multiple StraightLineSegments.  Uses
    tolerance to determine how fine to break them up.
    """
    def __init__(self, tolerance):
        self.tolerance = tolerance
        self.toleranceSquared = tolerance**2
    def convertSegment(self, segment):
        if isinstance(segment, BezierSegment):
            logDebug("Examining bezier %s" % segment)
            return self.divideBezier( segment, 0.0, 1.0 )
        return None
    def divideBezier(self, segment, start, finish ):
        logDebug("Examining bezier between %.4f and %.4f" % (start, finish))
        returnSegments = []
        if self.isStraightBezier( segment, start, finish ):
            startPoint = interpolateBezier( segment, start )
            finishPoint = interpolateBezier( segment, finish )
            returnSegments.append(StraightLineSegment( startPoint, finishPoint ))
        else:
            # Simple binary subdivision using recursion
            middle = ( start + finish ) / 2.0
            firstHalf = self.divideBezier( segment, start, middle )
            secondHalf = self.divideBezier( segment, middle, finish )
            for newSegment in firstHalf:
                returnSegments.append(newSegment)
            for newSegment in secondHalf:
                returnSegments.append(newSegment)
        return returnSegments
    def isStraightBezier(self, segment, start, finish ):
        """Test to see if the portion of the bezier between start and finish is
        straight to within self.tolerance.  Note that we must test three points in
        order to avoid the case of a symmetric cubic curve."""
        startPoint = interpolateBezier( segment, start )
        finishPoint = interpolateBezier( segment, finish )
        straightLine = StraightLineSegment( startPoint, finishPoint )
        delta = (finish - start ) / 4
        current = start + delta
        while current < finish:
            currentPoint = interpolateBezier( segment, current )
            logDebug("Interpolated at %.4f as %s" % (current, currentPoint))
            if not isPointStraightLineDistanceLessThanTolerance( straightLine, currentPoint, self.toleranceSquared ):
                return False
            current += delta
        logDebug("Bezier is straight between %.4f and %.4f" % (start, finish) )
        return True
    
    
class ExportGcode(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-x", "--tab",
                        action="store", type="string",
                        dest="tab", default=None)
        self.OptionParser.add_option("--filenameBase",
                        action="store", type="string",
                        dest="filenameBase", default=option_filenameBase,
                        help="Base filename (layer names will be appended)")
        self.OptionParser.add_option("--machineTolerance",
                        action="store", type="float",
                        dest="machineTolerance", default=option_machineTolerance,
                        help="Machine tolerance/resolution (in drawing units)")
        self.OptionParser.add_option("--traverseZ",
                        action="store", type="float",
                        dest="traverseZ", default=option_traverseZ,
                        help="Height to traverse when not cutting (in drawing units)")
        self.OptionParser.add_option("--defaultCutZ",
                        action="store", type="float",
                        dest="defaultCutZ", default=option_defaultCutZ,
                        help="Default height to cut at (in drawing units)")
        self.OptionParser.add_option("--defaultFeedRate",
                        action="store", type="float",
                        dest="defaultFeedRate", default=option_defaultFeedRate,
                        help="Default tool feed rate (in drawing units per minute)")
        self.OptionParser.add_option("--returnToXYOrigin",
                        action="store", type="inkbool",
                        dest="returnToXYOrigin", default=option_returnToXYOrigin,
                        help="Return to X0.0 Y0.0 at end")
        self.OptionParser.add_option("--useCutVariables",
                        action="store", type="inkbool",
                        dest="useCutVariables", default=option_useCutVariables,
                        help="Include variables for cutZ, traverseZ and feedRate")
        self.OptionParser.add_option("--useScaleVariables",
                        action="store", type="inkbool",
                        dest="useScaleVariables", default=option_useScaleVariables,
                        help="Include scale variables in X Y coordinates")
        self.OptionParser.add_option("--reverseYAxis",
                        action="store", type="inkbool",
                        dest="reverseYAxis", default=option_reverseYAxis,
                        help="Reverse Y Axis (preserves convention for Inkscape/EMC)")
        self.OptionParser.add_option("--outputExtension",
                        action="store", type="string",
                        dest="outputExtension", default=option_outputExtension,
                        help="File extension for output files")

    def effect(self):
        logDebug("Args %s" % sys.argv)
        # Stage 1
        # -------
        # Set up all the helper objects.  This is a great place to extend the code
        # by implementing a new reader or converter.
        logStage("Creating helper readers and converters...")
        readers = {}
        readers[inkex.addNS("path","svg")] = PathElementReader()
        readers[inkex.addNS("line","svg")] = LineElementReader()
        readers[inkex.addNS("rect","svg")] = RectElementReader()
        readers[inkex.addNS("g","svg")] = NullElementReader()
        converters = []
        converters.append( BezierStraightLineConverter(self.options.machineTolerance) )
        
        # Stage 2
        # -------
        # Determine units for the document:
        logStage("Determining units...")
        unitAttr = self.document.xpath('//sodipodi:namedview/@inkscape:document-units', namespaces=inkex.NSS)
        units = "mm"
        if unitAttr:
            units = unitAttr[0]
        logMessage("Units: %s" % units )
        #Calculate conversion factor from drawing units to physical units
        try:
            self.scale = inkex.uuconv[units]
        except KeyError:
            logError("Invalid units - preserving scale")
            self.scale = 1.0
        logDebug("Using overall scale conversion of %.4f" % self.scale )

        # Stage 3
        # -------
        # See if an alternative Origin has been specified in the drawing (i.e.
        # a layer named 'Origin' that contains either a circle or two intersecting
        # lines.
        logStage("Trying to find alternative origin...")
        self.origin = None
        originLayers = self.document.xpath("svg:g[@inkscape:groupmode='layer' and @inkscape:label='Origin']", namespaces=inkex.NSS)
        logDebug("Found %s origin layers" % len(originLayers))
        if len(originLayers) == 1:
            logDebug("Attempting to find alternative origin")
            originLayer = originLayers[0]
            dummyCutList = CutList()
            self.convertElementAndChildren( readers, dummyCutList, originLayer )
            origin = None
            if len(dummyCutList.cutPaths) == 1:
                # Must be a single circle
                singlePath = dummyCutList.cutPaths[0]
                if len(singlePath.segments) == 1:
                    circle = singlePath.segments[0]
                    if isinstance(circle, ArcSegment):
                        origin = circle.center
            elif len(dummyCutList.cutPaths) == 2:
                # Must be two straight line segments
                path1 = dummyCutList.cutPaths[0]
                path2 = dummyCutList.cutPaths[1]
                if len(path1.segments) ==1 and len(path2.segments) == 1:
                    line1 = path1.segments[0]
                    line2 = path2.segments[0]
                    if isinstance(line1, StraightLineSegment) and isinstance(line2,StraightLineSegment):
                        origin = calculateLineSegmentIntersectionPoint(line1,line2)
            if origin is not None:
                logMessage("Found alternative origin")
                self.origin = origin
            else:
                logWarning("Origin must contain either one circle, or two intersecting straight lines")            
        else:
            if len(originLayers) > 1:
                logWarning("Found multiple origin layers - ignoring")
            
        cutLists = []
        # Stage 4
        # -------
        # Generate CutLists from the objects in the drawing.  Note that CutLists
        # may contain Beziers etc, but these will have to be converted prior to
        # exporting GCode     
        logStage("Generating cut lists...")
        if len(self.selected) > 0:
            logMessage("Converting only selected items to GCode")
            cutList = CutList()
            #Name attribute will be used as a hint for output filename
            cutList.attributes.setAttribute("Name","SelectedObjects")
            for id, selectedElement in self.selected.iteritems():
                self.convertElement(readers, cutList, selectedElement )
            if cutList.hasCutPaths():
                cutLists.append(cutList)
            else:
                logError("No valid paths were created from selected items")
        else:
            # Search for valid layers i.e. of the form 'ToolPath <name> <Fxx> <Tyy> <Zzz>'
            layers = []
            for node in self.document.xpath("svg:g[@inkscape:groupmode='layer']", namespaces=inkex.NSS):
                layerName = node.attrib.get(inkex.addNS("label","inkscape"))
                if layerName[:len(g_toolPath)]==g_toolPath:
                    layers.append(node)
            # If there are any ToolPath layers, then only output those
            if len(layers) > 0:
                toolre = re.compile(r"\s+([FTZ])(-?[0-9]*\.?[0-9]*)")
                for layer in layers:
                    layerName = layer.attrib.get(inkex.addNS("label","inkscape"))
                    cutList = CutList()
                    cutList.attributes.setAttribute("Name", layerName)
                    parameters = toolre.findall(layerName)
                    for parameter in parameters:
                        if parameter[0]=="F":
                            cutList.attributes.setAttribute("FeedRate",float(parameter[1]))
                        elif parameter[0]=="T":
                            cutList.attributes.setAttribute("ToolSelection",parameter[1])
                        elif parameter[0]=="Z":
                            cutList.attributes.setAttribute("CutZ",float(parameter[1]))
                    self.convertElementAndChildren( readers, cutList, layer )
                    if cutList.hasCutPaths():
                        cutLists.append(cutList)
                    else:
                        logWarning("No paths found for layer '%s'" % layerName)
            else:
                # No ToolPath layers found, so output all objects
                cutList = CutList()
                cutList.attributes.setAttribute("Name", "AllObjects")
                # add all the root elements (svg groups)
                for node in self.document.xpath("svg:g", namespaces=inkex.NSS):
                    self.convertElementAndChildren( readers, cutList, node )
                if cutList.hasCutPaths():
                    cutLists.append(cutList)
                else:
                    logWarning("No paths found in document")
                    
        # Stage 5
        # -------
        # We need to convert Beziers etc into straight lines or arcs.
        # Currently only one converter is implemented, but a great thing would
        # be to develop a converter that goes from Beziers to Arcs (TBD).   
        logStage("Converting to GCode primitives...")
        for converter in converters:
            for cutList in cutLists:
                for cutPath in cutList.cutPaths:
                    newSegments = []
                    for segment in cutPath.segments:
                        convertedSegments = converter.convertSegment( segment )
                        if convertedSegments is None:
                            # Nothing to do, preserve the original segment
                            newSegments.append(segment)
                        else:
                            # Must have converted it
                            logDebug("Replacing segment with %i" % len(convertedSegments))
                            assert len(convertedSegments) > 0
                            assert convertedSegments[0].source == segment.source
                            assert convertedSegments[-1].destination == segment.destination
                            for convertedSegment in convertedSegments:
                                newSegments.append( convertedSegment )
                    #Switch over the segment list
                    cutPath.segments = newSegments

        # Stage 6
        # -------
        # Optional step - we may want to include some form of optimisation here.
        # For example, a simulated annealing which would examine contiguous straight
        # line segments, and see if they can be aggregated into a longer segment
        # without deviating more than machineTolerance      
        logStage("Optimising toolpaths (not implemented)...")
              
        # Stage 7
        # -------
        # Each CutList should by now only contain straight lines and arcs,
        # which can be directly converted to GCode primitives      
        logStage("Writing GCode files...")
        for cutList in cutLists:
            self.writeCutList(cutList, units)

        logStage("Finished.")

    def convertElementAndChildren(self, readers, cutList, element):
        self.convertElement(readers, cutList, element )
        for child in element:
            self.convertElementAndChildren(readers, cutList, child )
            
    def appendTransforms(self, n, coordinateTransformer):
        """
        Given a node, works its way back up the tree looking for 'transform'
        nodes and appending them to the transformer
        """
        while (n is not None):
            t = None
            try:
                t = n.get('transform')
            except AttributeError:
                pass

            if t != None:
                transform = simpletransform.parseTransform(t)
                if transform != None:
                    coordinateTransformer.appendTransform(transform)
            n = n.getparent()
        
    def convertElement(self, readers, cutList, element):
        # Get the correct transform
        coordinateTransformer = CoordinateTransformer()
        if self.options.reverseYAxis:
            coordinateTransformer.setScale(self.scale, -self.scale)
        else:
            coordinateTransformer.setScale(self.scale, self.scale)
            
        if self.origin is not None:
            coordinateTransformer.origin = self.origin 
        self.appendTransforms(element, coordinateTransformer)
        if element.tag in readers:
            reader = readers[element.tag]
            logDebug("Found reader for '%s'" % element.tag)
            cutPaths = reader.convertElementToPaths(element, coordinateTransformer )
            assert cutPaths is not None
            for cutPath in cutPaths:
                cutList.appendCutPath(cutPath)
        else:
            logError("No valid reader for '%s'" % element.tag )

    def writeCutList( self, cutList, units):
        assert cutList.hasCutPaths()
        assert cutList.attributes.hasAttribute("Name")
        name = cutList.attributes.getAttribute("Name")
        mops = GcodeOperationList()
        mops.appendComment("Found %i paths for %s:" % (len(cutList.cutPaths), name) )
        if units == "mm":
            mops.setUnitsToMillimeters()
        elif units == "in":
            mops.setUnitsToInches()
        else:
            #leave units undefined
            logWarning("Invalid units '%s'" % units )
        # Work out the correct options, start from global settings
        traverseZ = self.options.traverseZ
        cutZ = self.options.defaultCutZ
        feedRate = self.options.defaultFeedRate
        toolSelection = None
        #override from cutlist
        traverseZ = cutList.attributes.getAttributeWithDefault("TraverseZ", traverseZ)
        cutZ = cutList.attributes.getAttributeWithDefault("CutZ", cutZ)
        feedRate = cutList.attributes.getAttributeWithDefault("FeedRate", feedRate)
        toolSelection = cutList.attributes.getAttributeWithDefault("ToolSelection", toolSelection)
        if self.options.useCutVariables:
            mops.defineVariable(1,traverseZ, "Height to traverse when not cutting")
            mops.defineVariable(2,cutZ, "Height to cut at")
            mops.defineVariable(3,feedRate, "Tool feed rate")
        if self.options.useScaleVariables:
            mops.defineVariable(4,1.0, "X axis scale factor")
            mops.defineVariable(5,1.0, "Y axis scale factor")
        if toolSelection is not None:
            mops.appendToolChange( toolSelection )
        allPathsValid = True
        for cutPath in cutList.cutPaths:
            description = self.describeCutPath(cutPath, "Path")
            for line in description:
                mops.appendComment( line )
                logDebug( line )
            if not self.millCutPath(cutList, cutPath, mops, traverseZ, cutZ, feedRate):
                allPathsValid = False
        if self.options.returnToXYOrigin:
            if self.options.useScaleVariables:
                mops.appendRapidMove(x="[0.0*#4]", y="[0.0*#5]")
            else:
                mops.appendRapidMove(x=0.0, y=0.0)
        mops.appendEnd()
        if allPathsValid:
            #Check if we are writing to a file or not
            logDebug("About to write output")
            if len(self.options.filenameBase) > 0:
                filename = self.options.filenameBase
                if not re.match(r'.*[\\/]$', filename):
                    filename += " "
                # Munge any bad characters out of the layer info
                filename += re.sub(r'[^a-zA-Z0-9]',"", name)
                self.writeMops(mops, filename + "." + self.options.outputExtension)
            else:
                # Being run as 'SaveAs', so just output to stdout
                for line in mops.lines:
                    print line
        else:
            logWarning("Skipping output for %s due to invalid paths" % name)
            
    def millCutPath(self, cutList, cutPath, mops, traverseZ, cutZ, feedRate):
        assert cutPath.hasSegments()

        if self.options.useCutVariables:
            mops.appendRapidMove(z="#1")
        else:
            mops.appendRapidMove(z=traverseZ)
        firstPoint = cutPath.segments[0].source
        if self.options.useScaleVariables:
            xformula = "[%.4f*#4]" % firstPoint.x
            yformula = "[%.4f*#5]" % firstPoint.y
            mops.appendRapidMove(x=xformula, y=yformula)
        else:
            mops.appendRapidMove(x=firstPoint.x, y=firstPoint.y)
        if self.options.useCutVariables:
            mops.appendRapidMove(z="#2",f="#3")
        else:
            mops.appendFeedMove(z=cutZ, f=feedRate) 
        #Now output all the line paths (all at the same cut depth)
        for segment in cutPath.segments:
            #we now have a line segment that goes from firstPoint to secondPoint
            secondPoint = segment.destination
            if isinstance(segment, StraightLineSegment):
                if self.options.useScaleVariables:
                    xformula = "[%.4f*#4]" % firstPoint.x
                    yformula = "[%.4f*#5]" % firstPoint.y
                    mops.appendFeedMove(x=xformula, y=yformula)
                else:
                    mops.appendFeedMove(x=secondPoint.x, y=secondPoint.y)
            elif isinstance(segment, ArcSegment):
                mops.appendFeedArc(segment.clockwise, x=secondPoint.x, y=secondPoint.y, i=segment.center.x-firstPoint.x,j=segment.center.y-firstPoint.y)
            else:
                logError("Unable to export GCode for %s" % type(segment))
                return False
            firstPoint = secondPoint
        #Retract out of the work piece
        if self.options.useCutVariables:
            mops.appendRapidMove(z="#1")
        else:
            mops.appendRapidMove(z=traverseZ)
        return True
    
    def describeCutPath(self, cutPath, name, showPoints=False):
        description = []
        bounds = calculateBounds( cutPath )
        description.append(name + " has size %.3f x %.3f" % (bounds.destination.x - bounds.source.x, bounds.destination.y -bounds.source.y) )
        description.append(name + " has bounding box %.3f-%.3f wide by %.3f-%.3f high" % (bounds.source.x, bounds.destination.x, bounds.source.y, bounds.destination.y) )
        if showPoints:
            for point in cutPath.points:
                description.append("Point " % point )
        return description
       
    def writeMops(self, mops, filename):
        logDebug("Writing to %s" % filename)
        outputfile = open(filename, 'wb')
        for line in mops.lines:
            outputfile.write(line)
            outputfile.write("\r\n")
        outputfile.close()
           
    def output(self):
        if len(self.options.filenameBase) > 0:
            # Run as an effect, so output, otherwise do nothing (run as output filter,
            # so will have already written to stdout)
            # BUG:  This line sometimes hangs (due to writing large string
            # to stdout.  Not sure how to fix.
            inkex.Effect.output(self)
         
effect = ExportGcode()
# If you are running from Eclipse, go to 'Run Configurations' and ensure that
# on the environment tab there is an environment variable named 'RUN_FROM_ECLIPSE'
# and it is set to 'True'
import os
try:
    isRunFromEclipse = os.environ["RUN_FROM_ECLIPSE"]
except KeyError:
    isRunFromEclipse = False
if isRunFromEclipse:
    #TODO: Need a unit testing framework.  Should run each test .svg file
    # and compare to 'known-good' output.  The only problem is that small changes
    # in algorithms may affect the resultant GCode quite a bit, which will result
    # in a lot of work verifying that the new output is in-fact correct.
    #effect.getoptions(["C:/Documents and Settings/General/workspace/InkscapeGcode/test/rectangle.svg"])
    #effect.getoptions(["C:/Documents and Settings/General/workspace/InkscapeGcode/test/transformrectangle.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/rect.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/rectrx.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/rectrxry.svg"])
    #effect.getoptions(["C:/Documents and Settings/General/workspace/InkscapeGcode/test/bezierstraight.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/jonata_cricket.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/ec_error1.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/rectorigin.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "Y:/Polaris/Property/Tutukaka/New Barn/Fitout/Temperature Pad/Temperature v1.3.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/text.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/paths.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/compoundpaths.svg"])
    #effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/compoundpaths2.svg"])
    effect.getoptions(["--filenameBase", "C:/Local/nc", "C:/Documents and Settings/General/workspace/InkscapeGcode/test/ec_error2.svg"])
    logDebug(effect.args)
    effect.parse()
    effect.getposinlayer()
    effect.getselected()
    effect.getdocids()
    effect.effect()
    effect.output()
else:
    effect.affect()

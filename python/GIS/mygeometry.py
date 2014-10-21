"""
Extensions to shapely geometry library
"""

from shapely.geometry import LineString, Point
import numpy as np

class MyLine(LineString):
    def __init__(self,xy):
        LineString.__init__(self,xy)

    def unitnormal(self,dist,dx=0.001,direction='left'):
        """
        Unit normal at dist along line
        """
        if dist+dx > 1.:
            P2 = self.interpolate(dist-dx,normalized=True)
            P1 = self.interpolate(dist,normalized=True)
        else:
            P1 = self.interpolate(dist,normalized=True)
            P2 = self.interpolate(dist+dx,normalized=True)
        dx = P2.x - P1.x
        dy = P2.y - P1.y
        mag = P1.distance(P2)
        if direction=='left':
            return -dy/mag,dx/mag
        elif direction=='right':
            return dy/mag,-dx/mag

    def perpendicular(self,dist,mag,**kwargs):
        """
        Returns the location of a point that is located 'mag' distance
        perpendicular to a line at point at 'dist'
        """
        nx,ny = self.unitnormal(dist,**kwargs)
        P1 = self.interpolate(dist,normalized=True)
        return Point([P1.x+nx*mag, P1.y+ny*mag])

    def perpline(self,dist,mag):
        """
        Returns a line that intersects the current line at perpendicularly
        """
        P1 = self.perpendicular(dist,mag,direction='right')
        P2 = self.perpendicular(dist,mag,direction='left')
        return MyLine([[P1.x,P1.y],[P2.x,P2.y]])

    def multiperpline(self,npoints,mag):
        """
        Returns a list of multiple equally spaced perpendicular lines
        """
        pts = np.linspace(0,1,npoints)
        return [self.perpline(pts[ii],mag) for ii in range(npoints)]


    def multipoint(self,npoints):
        """
        Returns a list of equally spaced points along a line
        """
        pts = np.linspace(0,1,npoints)
        return [self.interpolate(pts[ii],normalized=True) for ii in\
            range(npoints)]



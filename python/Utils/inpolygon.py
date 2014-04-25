"""
Wrapper function for finding points inside polygon

This provides backwards compatibility as 
matplotlib changed the method as of v1.3.x 
"""

import matplotlib

# Need to check the version of matplotlib as they got rid of the
# points_in_polygon tool

if matplotlib.__version__<='1.2.x':
    from matplotlib.nxutils import points_inside_poly 
    version = 'old'
else:
    from matplotlib.path import Path
    version = 'new'

def inpolygon(xy,xypoly):
    if version=='old':
        return points_inside_poly(xy,xypoly)
    else:
        Poly = Path(xypoly)
        return Poly.contains_points(xy)
 

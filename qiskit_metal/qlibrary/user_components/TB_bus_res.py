# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""File contains dictionary for NSquareSpiral and the make()."""
"""This file originated from rectangle spiral resonator component"""

from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np


class TBbusRes(QComponent):
    """A rectangle spiral resonator based on length input. The X dimension is
    modified by the code based on the total length inputed.

    Inherits `QComponent` class

    A rectangular spiral resonator. The width of the spiral is modified based on inputted values
    and given total length of the spiral.

        ::

            <--------X-------->
            __________________
            |   ___________   |
            |  |           |  |
            |  |           |  |
            |  |______________|
            |

    .. image::
        ResonatorCoilRect.png

    .. meta::
        Resonator Coil Rectangle

    Default Options:
        Convention: Values (unless noted) are strings with units included,
        (e.g., '30um')

        * n: '3' -- Number of turns of the spiral
        * length: '2000um' -- Total length of the spiral
        * line_width: '1um' -- The width of the line of the spiral
        * height: '40um' -- The height of the inner portion of the spiral
        * gap: '4um' -- The distance between each layer of the spiral
        * coupler_distance: '10um' -- The pin position from the grounded termination of the spiral
    """
    component_metadata = Dict(short_name='res')
    """Component metadata"""

    default_options = Dict(length='2000um',
                           line_width='20um',
                           width='400um',
                           gap='6um',
                           neck = '100um',
                           do_arc= True)
    """Default drawing options"""

    TOOLTIP = """A rectangle spiral resonator based on length input. The X dimension is
    modified by the code based on the total length inputed."""
    
    def make(self):
        """The make function implements the logic that creates the geoemtry
        (poly, path, etc.) from the qcomponent.options dictionary of
        parameters, and the adds them to the design, using
        qcomponent.add_qgeometry(...), adding in extra needed information, such
        as layer, subtract, etc."""
        p = self.p  # p for parsed parameters. Access to the parsed options.
        # Create the geometry

        meander_list = []
        pmeander_list = []

        def make_arc(x, y, radius, start_angle, end_angle):
            npts = abs(end_angle - start_angle)//5 +1
            angle_range = np.linspace(start_angle, end_angle, npts)
            arc_points = [(x + radius * np.cos(np.deg2rad(angle)),
                           y + radius * np.sin(np.deg2rad(angle))) for angle in angle_range]
            return arc_points
        arc_radius1 = p.line_width * 0.5
        arc_radius2 = arc_radius1 + p.line_width + p.gap

        # Calculate the number of turns based on the length and width
        turns = (p.length - p.width - 2 * p.neck) / 0.8 / p.width
        if turns < 1:
            self._error_message = f'Inputted values results in the number of turns being less than 1.'
            self.logger.warning(self._error_message)
            return
        
        # Define the fisrt, last, neck points of bus resonator
        x_first = -p.width/2
        y_first = 0
        x_last = p.width/2
        y_last = 0
        x_neck1 = -(p.gap + p.line_width)/2
        y_neck1 = 0
        x_neck2 = -(p.gap + p.line_width)/2
        y_neck2 = p.neck - (p.gap + p.line_width)/2
        x_neck3 = (p.gap + p.line_width)/2
        y_neck3 = p.neck + (p.gap + p.line_width)/2
        x_neck4 = (p.gap + p.line_width)/2
        y_neck4 = 0
        
        # add the first, neck 1 2 points to the meander list. other points will be added after the loop
        meander_list.append((x_first, y_first))
        pmeander_list.append((x_last, y_last))
        
        if p.do_arc:
            arc1 = make_arc(x_neck1 - arc_radius1*1, y_neck1 + arc_radius1*1, arc_radius1, 270, 360)
            arc2 = make_arc(x_neck2 - arc_radius1*1, y_neck2 - arc_radius1*1, arc_radius1, 0, 90)
            arc3 = make_arc(x_neck3 - arc_radius2*1, y_neck3 - arc_radius2*1, arc_radius2, 0, 90)
            arc4 = make_arc(x_neck4 + arc_radius1*1, y_neck4 + arc_radius1*1, arc_radius1, 270, 180)
            meander_list.extend(arc1)
            meander_list.extend(arc2)
            pmeander_list.extend(arc4)
            pmeander_list.extend(arc3)
        else:
            meander_list.append((x_neck1, y_neck1))
            meander_list.append((x_neck2, y_neck2))
            pmeander_list.append((x_neck4, y_neck4))
            pmeander_list.append((x_neck3, y_neck3))


        
        # Define points for the dual meander
        # i_point#p is the point on opposite side.
        for step in range(int(turns)):
            if step % 2 == 0:
                x_point1 = -p.width*0.4 - (p.gap + p.line_width)/2
                y_point1 = p.neck - (p.gap + p.line_width)/2 + step * 2*(p.line_width + p.gap)
                x_point2 = -p.width*0.4 - (p.gap + p.line_width)/2
                y_point2 = p.neck - (p.gap + p.line_width)/2 + (step+1.5) * 2*(p.line_width + p.gap)
                y_point1p = y_point1 + (p.gap + p.line_width)
                y_point2p = y_point2 - (p.gap + p.line_width)
                

            elif step % 2 == 1:
                x_point1 = p.width*0.4 - (p.gap + p.line_width)/2
                y_point1 = p.neck + (p.gap + p.line_width)/2 + (step) * 2*(p.line_width + p.gap)
                x_point2 = p.width*0.4 - (p.gap + p.line_width)/2
                y_point2 = p.neck + (p.gap + p.line_width)/2 + (step+0.5) * 2*(p.line_width + p.gap)
                y_point1p = y_point1 - (p.gap + p.line_width)
                y_point2p = y_point2 + (p.gap + p.line_width)

            x_point1p = x_point1 + (p.gap + p.line_width)
            x_point2p = x_point2 + (p.gap + p.line_width)

            if p.do_arc:
                # If do_arc is True, add arcs instead of straight lines
                # Use consistent radius for maintaining gap spacing
                if step % 2 == 0:
                    # For even steps, create arcs that maintain proper spacing
                    # Arc centers are adjusted to maintain the gap distance
                    arc1 = make_arc(x_point1 + arc_radius2*1, y_point1 + arc_radius2*1, arc_radius2, 270, 180)
                    arc2 = make_arc(x_point2 + arc_radius2*1, y_point2 - arc_radius2*1, arc_radius2, 180, 90)
                    arc1p = make_arc(x_point1p + arc_radius1*1, y_point1p + arc_radius1*1, arc_radius1, 270, 180)
                    arc2p = make_arc(x_point2p + arc_radius1*1, y_point2p - arc_radius1*1, arc_radius1, 180, 90)
                elif step % 2 == 1:
                    # For odd steps, create arcs that maintain proper spacing
                    arc1 = make_arc(x_point1 - arc_radius1*1, y_point1 + arc_radius1*1, arc_radius1, 270, 360)
                    arc2 = make_arc(x_point2 - arc_radius1*1, y_point2 - arc_radius1*1, arc_radius1, 0, 90)
                    arc1p = make_arc(x_point1p - arc_radius2*1, y_point1p + arc_radius2*1, arc_radius2, 270, 360)
                    arc2p = make_arc(x_point2p - arc_radius2*1, y_point2p - arc_radius2*1, arc_radius2, 0, 90)

                meander_list.extend(arc1)
                meander_list.extend(arc2)
                pmeander_list.extend(arc1p)
                pmeander_list.extend(arc2p)
            else:
                meander_list.append((x_point1, y_point1))
                meander_list.append((x_point2, y_point2))
                pmeander_list.append((x_point1p, y_point1p))
                pmeander_list.append((x_point2p, y_point2p))

        tail_length = (turns- int(turns))*0.8 * p.width
        if step % 2 == 0: # if the last step is even, add the tail to the right
            x_tail1 = x_point2 + tail_length
            y_tail1 = y_point2
            x_tail2 = x_point2 + tail_length
            y_tail2 = y_point2 - (p.line_width + p.gap)
        elif step % 2 == 1: # if the last step is odd, add the tail to the left
            x_tail1 = x_point2 - tail_length
            y_tail1 = y_point2
            x_tail2 = x_point2 - tail_length
            y_tail2 = y_point2 + (p.line_width + p.gap)

        if p.do_arc:
            if step % 2 == 0:
                arc_tail1 = make_arc(x_tail1 - arc_radius1*1, y_tail1 - arc_radius1*1, arc_radius1, 90, 0)
                arc_tail2 = make_arc(x_tail2 - arc_radius1*1, y_tail2 + arc_radius1*1, arc_radius1, 360, 270)
            elif step % 2 == 1:
                arc_tail1 = make_arc(x_tail1 + arc_radius1*1, y_tail1 + arc_radius1*1, arc_radius1, 270, 180)
                arc_tail2 = make_arc(x_tail2 + arc_radius1*1, y_tail2 - arc_radius1*1, arc_radius1, 180, 90)
            meander_list.extend(arc_tail1)
            meander_list.extend(arc_tail2)
        else:
            meander_list.append((x_tail1, y_tail1))
            meander_list.append((x_tail2, y_tail2))

        

        # reverse the pmeander_list and merge lists
        pmeander_list.reverse()
        meander_list.extend(pmeander_list)

        meander_list = draw.LineString(meander_list)
        meander_etch = draw.shapely.geometry.box(
            x_first, -(p.gap + p.line_width) / 2, 
            x_last, max(y_tail1, y_tail2) + (p.gap + p.line_width) / 2)
        
        c_items = [meander_list, meander_etch]
        c_items = draw.rotate(c_items, p.orientation, origin=(0, 0))
        c_items = draw.translate(c_items, p.pos_x, p.pos_y)
        [meander_list, meander_etch] = c_items
        ##############################################
        # add elements
        self.add_qgeometry('path', {'n_meander': meander_list},
                           width=p.line_width)
        self.add_qgeometry('poly', {'n_meander_etch': meander_etch},
                           subtract=True)
        

        # NEW PIN SPOT
        pin_list = meander_list.coords
        self.add_pin('busPin1',
                     points=np.array(pin_list[1::-1]),
                     width=p.line_width,
                     input_as_norm=True)
        self.add_pin('busPin2',
                     points=np.array(pin_list[-2:]),
                     width=p.line_width,
                     input_as_norm=True)


        

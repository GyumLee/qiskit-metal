# -*- coding: utf-8 -*-

# This code originate from the part of Qiskit.
# Modified by Heok

from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit
import numpy as np
class TransmonPocket_sqnl(BaseQubit):
    default_options = Dict(
        pad_gap='30um',
        inductor_width='20um',
        pad_width='455um',
        pad_height='90um',
        pocket_width='650um',
        pocket_height='650um',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        _default_connection_pads=Dict(
            connector_type='0',  # 0 = Claw type, 1 = T-shape
            pad_width='125um',
            pad_height='30um',
            pad_cpw_shift='0um',
            pad_cpw_extent='25um',
            claw_length='30um',
            claw_width='10um',
            claw_gap='6um',
            claw_cpw_length='40um',
            claw_cpw_width='10um',
            ground_spacing='5um',
            connector_location= '0'  # 0~5 means six locations with rotated CCW starting from right-top site
        ))
    """Default drawing options"""

    component_metadata = Dict(short_name='Pocket_claw',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')
    """Component metadata"""

    TOOLTIP = """The base `TransmonPocket_sqnl` class."""

    def make(self):
        """Define the way the options are turned into QGeometry.

        The make function implements the logic that creates the geoemtry
        (poly, path, etc.) from the qcomponent.options dictionary of
        parameters, and the adds them to the design, using
        qcomponent.add_qgeometry(...), adding in extra needed
        information, such as layer, subtract, etc.
        """
        self.make_pocket()
        self.make_connection_pads()

    def make_pocket(self):
        """Makes standard transmon in a pocket."""

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        # extract chip name
        chip = p.chip

        # since we will reuse these options, parse them once and define them as varaibles
        pad_width = p.pad_width
        pad_height = p.pad_height
        pad_gap = p.pad_gap

        # make the pads as rectangles (shapely polygons)
        pad = draw.rectangle(pad_width, pad_height)
        pad_top = draw.translate(pad, 0, +(pad_height + pad_gap) / 2.)
        pad_bot = draw.translate(pad, 0, -(pad_height + pad_gap) / 2.)

        rect_jj = draw.LineString([(0, -pad_gap / 2), (0, +pad_gap / 2)])
        # the draw.rectangle representing the josephson junction
        # rect_jj = draw.rectangle(p.inductor_width, pad_gap)

        rect_pk = draw.rectangle(p.pocket_width, p.pocket_height)

        # Rotate and translate all qgeometry as needed.
        # Done with utility functions in Metal 'draw_utility' for easy rotation/translation
        # NOTE: Should modify so rotate/translate accepts qgeometry, would allow for
        # smoother implementation.
        polys = [rect_jj, pad_top, pad_bot, rect_pk]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk] = polys

        # Use the geometry to create Metal qgeometry
        self.add_qgeometry('poly',
                           dict(pad_top=pad_top, pad_bot=pad_bot),
                           chip=chip)
        self.add_qgeometry('poly',
                           dict(rect_pk=rect_pk),
                           subtract=True,
                           chip=chip)
        # self.add_qgeometry('poly', dict(
        #     rect_jj=rect_jj), helper=True)
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=p.inductor_width,
                           chip=chip)

    def make_connection_pads(self):
        """Makes standard transmon in a pocket."""
        for name in self.options.connection_pads:
            self.make_connection_pad(name)

    def make_connection_pad(self, name: str):
        """Makes n individual connector.

        Args:
            name (str) : Name of the connector
        """
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.connection_pads[name]  # parser on connector options
        # extract chip name
        chip = p.chip
    
        if pc.connector_type == 0:  # Claw connector
            
            # define commonly used variables once
            pad_width = p.pad_width
            pad_height = p.pad_height
            pad_gap = p.pad_gap
            pocket_width = p.pocket_width
            pocket_height = p.pocket_height
    
            c_g = pc.claw_gap
            c_l = pc.claw_length
            c_w = pc.claw_width
            c_c_w = pc.claw_cpw_width
            c_c_l = pc.claw_cpw_length
            g_s = pc.ground_spacing
            con_loc = pc.connector_location
    
            claw_cpw = draw.box(-c_w, -c_c_w / 2, -c_c_l - c_w, c_c_w / 2)
            if con_loc == 0 or  con_loc == 2 or  con_loc == 3 or  con_loc == 5:
                t_claw_height = 2*c_g + 2 * c_w + pad_height  # temp value
    
                claw_base     = draw.box(-c_w, -(t_claw_height) / 2, c_l, t_claw_height / 2)
                claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l, t_claw_height / 2 - c_w)
                claw_base = claw_base.difference(claw_subtract)
    
                connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
                connector_etcher = draw.buffer(connector_arm, c_g)
            elif con_loc == 1 or  con_loc == 4:
                t_claw_height = 2*c_g + 2 * c_w + pad_width  # temp value
    
                claw_base     = draw.box(-c_w, -(t_claw_height) / 2, c_l, t_claw_height / 2)
                claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l, t_claw_height / 2 - c_w)
                claw_base = claw_base.difference(claw_subtract)
    
                connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
                connector_etcher = draw.buffer(connector_arm, g_s)
            else:
                print('Not appropriate position of connector')
            
            port_line = draw.LineString([(- c_w, 0),
                                     (-c_c_l - c_w, 0)])

            claw_rotate = 0
            if con_loc == 0:
                claw_rotate = -180
            elif con_loc == 1:
                claw_rotate = -90
            elif con_loc == 2:
                claw_rotate = 0
            elif con_loc == 3:
                claw_rotate = 0
            elif con_loc == 4:
                claw_rotate = 90
            elif con_loc == 5:
                claw_rotate = 180
    
            claw_translate_x = (con_loc==0 or con_loc==5)*(pad_width/2 + c_g) + (con_loc==2 or con_loc==3)*(-pad_width/2 - c_g)
            claw_translate_y = (con_loc==0 or con_loc==2)*(pad_height + pad_gap)/2 + (con_loc==5 or con_loc==3)*(-pad_height - pad_gap)/2 + (con_loc==1)*(pad_height + pad_gap/2 + c_g) - (con_loc==4)*(pad_height + pad_gap/2 + c_g)
            
            # Rotates and translates the connector polygons (and temporary port_line)
            polys = [connector_arm, connector_etcher, port_line]
            polys = draw.rotate(polys, claw_rotate, origin=(0, 0))
            polys = draw.translate(polys, claw_translate_x, claw_translate_y)
            #polys = draw.rotate(polys, p.orientation, origin=(0, 0))
            polys = draw.translate(polys, p.pos_x, p.pos_y)
            [connector_arm, connector_etcher, port_line] = polys
    
            # Generates qgeometry for the connector pads
            self.add_qgeometry('poly', {f'{name}_connector_arm': connector_arm},
                               chip=chip)
            self.add_qgeometry('poly',
                               {f'{name}_connector_etcher': connector_etcher},
                               subtract=True,
                               chip=chip)
    
            self.add_pin(name, points=port_line.coords, width=c_c_w,input_as_norm=True)
        elif pc.connector_type == 1:  # T-shape connector
            
            # define commonly used variables once
            pad_width = p.pad_width
            pad_height = p.pad_height
            pad_gap = p.pad_gap
            pocket_width = p.pocket_width
            pocket_height = p.pocket_height
    
            c_g = pc.claw_gap
            c_l = 0
            c_w = pc.claw_width
            c_c_w = pc.claw_cpw_width
            c_c_l = pc.claw_cpw_length
            g_s = pc.ground_spacing
            con_loc = pc.connector_location
    
            claw_cpw = draw.box(-c_w, -c_c_w / 2, -c_c_l - c_w, c_c_w / 2)
            if con_loc == 0 or  con_loc == 2 or  con_loc == 3 or  con_loc == 5:
                t_claw_height = pc.t_claw_height
    
                claw_base     = draw.box(-c_w, -(t_claw_height) / 2, c_l, t_claw_height / 2)
                claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l, t_claw_height / 2 - c_w)
                claw_base = claw_base.difference(claw_subtract)
    
                connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
                connector_etcher = draw.buffer(connector_arm, c_g)
            elif con_loc == 1 or  con_loc == 4:
                t_claw_height = pc.t_claw_height
    
                claw_base     = draw.box(-c_w, -(t_claw_height) / 2, c_l, t_claw_height / 2)
                claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l, t_claw_height / 2 - c_w)
                claw_base = claw_base.difference(claw_subtract)
    
                connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
                connector_etcher = draw.buffer(connector_arm, g_s)
            else:
                print('Not appropriate position of connector')
            
            port_line = draw.LineString([(- c_w, 0),
                                     (-c_c_l - c_w, 0)])

            claw_rotate = 0
            if con_loc == 0:
                claw_rotate = -180
            elif con_loc == 1:
                claw_rotate = -90
            elif con_loc == 2:
                claw_rotate = 0
            elif con_loc == 3:
                claw_rotate = 0
            elif con_loc == 4:
                claw_rotate = 90
            elif con_loc == 5:
                claw_rotate = 180
    
            claw_translate_x = (con_loc==0 or con_loc==5)*(pad_width/2 + c_g) + (con_loc==2 or con_loc==3)*(-pad_width/2 - c_g)
            claw_translate_y = (con_loc==0 or con_loc==2)*(pad_height + pad_gap)/2 + (con_loc==5 or con_loc==3)*(-pad_height - pad_gap)/2 + (con_loc==1)*(pad_height + pad_gap/2 + c_g) - (con_loc==4)*(pad_height + pad_gap/2 + c_g)
            
            # Rotates and translates the connector polygons (and temporary port_line)
            polys = [connector_arm, connector_etcher, port_line]
            polys = draw.rotate(polys, claw_rotate, origin=(0, 0))
            polys = draw.translate(polys, claw_translate_x, claw_translate_y)
            polys = draw.rotate(polys, p.orientation, origin=(0, 0))
            polys = draw.translate(polys, p.pos_x, p.pos_y)
            [connector_arm, connector_etcher, port_line] = polys
    
            # Generates qgeometry for the connector pads
            self.add_qgeometry('poly', {f'{name}_connector_arm': connector_arm},
                               chip=chip)
            self.add_qgeometry('poly',
                               {f'{name}_connector_etcher': connector_etcher},
                               subtract=True,
                               chip=chip)
    
            self.add_pin(name, points=port_line.coords, width=c_c_w,input_as_norm=True)

        else:
            print('Incorrect connector type')
            
        
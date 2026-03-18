# NEW OpenDrift module for bluebottle drift.
# Amandine  simplified 2023
# 12/4/24 Update equation drift angle 

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift
from opendrift.elements.passivetracer import PassiveTracer
from opendrift.readers import reader_global_landmask, reader_netCDF_CF_generic# , reader_netCDF_CF_unstructured # for corners
from opendrift.models.oceandrift import Lagrangian3DArray

# We add bluebottle properties to the existing element class PassiveTracer.
# AMANDINE careful that drag coef depends on angle of attack
PassiveTracer.variables = PassiveTracer.add_variables([
    ('Area_ratio', {'dtype': np.float32,
                           'units': 'm',
                     'description': 'Area above water / area under water in the long bluebottle axis',
                          'default': 2}),
    ('Drag_coef_ratio', {'dtype': np.float32,
                           'units': 'm',
                     'description': 'Drag coefficient above water / drag coefficient under water in the long bluebottle axis',
                          'default': 0.7}),
    ('Wind_perc', {'dtype': np.float32,
                           'units': 'm',
                     'description': 'Percentage of wind speed that contributes to the bluebottle movement',
                          'default': 0}),
    ('Wind_perc_calc', {'dtype': np.float32,
                           'units': 'm',
                     'description': 'Calculated Percentage of wind speed that contributes to the bluebottle movement',
                          'default': 0}),
    ('Angle_course_fromW', {'dtype': np.float32,
                           'units': 'deg',
                     'description': 'Course angle relative from wind. Positive anticlockwise. If not np.radians(40 *  np.exp(-0.4*Wind_speed)',
                          'default': 0}),
    ('Angle_course_fromW_calc', {'dtype': np.float32,
                           'units': 'deg',
                     'description': 'Calculated Course angle relative from wind. Positive anticlockwise. If not np.radians(40 * np.exp(-0.4*Wind_speed)',
                          'default': 0}),
    ('Orientation', {'dtype': np.float32,
                         'units': '1',
                     'description': '1 for right-handed (left-sailing). -1 for left-handed (right sailing)',
                         'default': 1}),
    ])

class BBDrift(OceanDrift):
    """Trajectory model based on the OceanDrift framework.
    """

    ElementType = PassiveTracer

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        #'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        #'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
#         'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
#         'latmin': {'fallback':-40},
        }

    # Default colors for plotting
#     status_colors = {'initial': 'green', 'active': 'blue',
#                      'outside domain': 'gray', 'stranded': 'red',
#                      'left-handed': 'yellow', 'right-handed': 'magenta'}
    
    # Configuration
    def __init__(self, *args, **kwargs):

        # Call parent constructor.
        super(BBDrift, self).__init__(*args, **kwargs)
        
#         self._add_config({
#             'seed:object_type': {'type': 'enum', 'enum': descriptions,
#                 'default': descriptions[0],
#                 'description': 'Leeway object category for this simulation',
#                 'level': self.CONFIG_LEVEL_ESSENTIAL},
#             'seed:jibe_probability': {'type': 'float',
#                 'default': 0.04, 'min': 0, 'max': 1,
#                 'description': 'Probability per hour for jibing (objects changing orientation)',
#                 'units': 'probability', 'level': self.CONFIG_LEVEL_BASIC},
#             })

    def update(self):
        """Update positions of bluebottles at each timestep."""

        
        # Calculate wind speed and direction
        Wind_speed = np.sqrt(self.environment.x_wind**2 + self.environment.y_wind**2)
        Wind_angle = np.arctan2(self.environment.y_wind,self.environment.x_wind)        

        # Densities of air and water.
        Rho_A = 1.225
        Rho_H = 1025        
        
        
        # Course Angle from wind. If not defined, use a function that depends on wind speed. 40deg for weak wind, 0 for strong winds 
        Angle_course_fromW_final = np.zeros(len(self.elements.Angle_course_fromW))
#         print("len angle", len(self.elements.Angle_course_fromW))
#         print('len wind', len(Wind_speed))
        for i in range(len(self.elements.Angle_course_fromW)):
            if self.elements.Angle_course_fromW[i] == 0: 
                Angle_course_fromW_final[i] = (50.8 * np.exp(-0.15*Wind_speed[i])-.5) # degrees
#                 print('wind speed', Wind_speed[i])
            else:     
                Angle_course_fromW_final[i] = self.elements.Angle_course_fromW[i] # degrees
        # Course depends on orientation of the bluebottle.(1 for right-handed (left-sailing). -1 for left-handed (right sailing)',)
        alpha = np.degrees(Wind_angle) + self.elements.Orientation * Angle_course_fromW_final
        alpha_rad = np.radians(alpha)

        
         # Shape parameter (Lambda). 
         # - either defined 
         # - either, if drift angle not defined, scale drag_coef_ratio to 1 if drift=0deg, to 0.5 if drift=40deg, then calculate from equation with area_ratio and the new drag_coef_ratio
        Wind_perc_final = 0 + np.zeros(len(self.elements.Wind_perc))
#         print("len perc", len(self.elements.Wind_perc))
        scale_ratio = 1+np.zeros(len(self.elements.Wind_perc))
        for i in range(len(self.elements.Wind_perc)):
            if self.elements.Wind_perc[i] == 0:
                scale_ratio[i] = 1 - (Angle_course_fromW_final[i]) /80
                Wind_perc_final[i] = np.sqrt(Rho_A/Rho_H * (self.elements.Area_ratio[i] * scale_ratio[i]) * self.elements.Drag_coef_ratio[i])
#                 print('perc', scale_ratio[i])
#                 print('perc', Wind_perc_final[i])
            else:     
                Wind_perc_final[i] = self.elements.Wind_perc[i]
#         print('Angle_course_fromW_final',  np.degrees(Angle_course_fromW_final))
#         print('Wind_perc_final', Wind_perc_final)
#         print('scale_ratio', scale_ratio)
        
        # Bluebottle final speed and course.
        V_bbx = Wind_perc_final*Wind_speed*np.cos(alpha_rad)
        V_bby = Wind_perc_final*Wind_speed*np.sin(alpha_rad)
#         print("vbb",V_bby)

        
        # Save
        self.elements.Wind_perc_calc  = Wind_perc_final
        self.elements.Angle_course_fromW_calc = Angle_course_fromW_final
        self.elements.current_drift_factor = 1.0 + np.zeros(len(self.elements.Wind_perc))
#         self.elements.mass = self.elements.mass - np.sqrt(self.environment.x_wind**2 + self.environment.y_wind**2)*.01

        # Update particle positions based on final speed and course.
        self.update_positions(V_bbx, V_bby)

         # Simply move particles with ambient current.
        self.advect_ocean_current()
       

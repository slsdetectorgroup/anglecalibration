
from angcal import AngleCalibration, EpicsMythenFileReader, PlotHelper
import numpy as np
import matplotlib.pyplot as plt

def select_base_peak(anglecalibration,
                      mythen_file_reader,
                      filelist : list[str],
                      module_index = 0):
    
    """ 
        select base peak by plotting module 0 redistributed to fixed
        angle-width bins for each frame in filelist with a detector angle between
        6° and 33°
    """
    
    plotter = PlotHelper(anglecalibration)

    fig, axis = plt.subplots()
    plotter.overwrite_plot = True

    detector_angle_range = (6.0, 33.0) 

    left_module_strip_angle = anglecalibration.diffraction_angle_from_DG_parameters(module_index, 0.0, 0, -0.5)
      
    print("strip angle: ", left_module_strip_angle)

    print("adjusted detector angle range for module {}: {} to {}".format(module_index,
        detector_angle_range[0], detector_angle_range[1]))
    
    for file in filelist:

        detector_angle = mythen_file_reader.read_detector_angle(file)

        # -5.0, 5.0 adjust for movement
        if (detector_angle + left_module_strip_angle - 5.0 >
                detector_angle_range[0] and
            detector_angle + left_module_strip_angle + 5.0 <
                detector_angle_range[1]):
            print("plotting file {} with detector angle {}".format(file, detector_angle))
              
            #plot module
            module_redistributed_to_fixed_angle_bins = anglecalibration.convert(
                    [file], module_index) # always plots diffraction pattern independant of detector location.

            plotter.plot_diffraction_pattern(
                module_redistributed_to_fixed_angle_bins.view(), module_index,
                detector_angle, axis=axis)
    
    

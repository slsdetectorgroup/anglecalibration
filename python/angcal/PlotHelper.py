import matplotlib.pyplot as plt 
import numpy as np 

class PlotHelper: 

    """ Helper class for plotting diffraction patterns. """
    
    def __init__(self, anglecalibration): 
        # maybe add freeze 
        self._anglecalibration = anglecalibration
        self._left_module_boundary_as_bin_index = lambda module_index, detector_angle : int((self._anglecalibration.diffraction_angle_from_DG_parameters( module_index, detector_angle, 0, -0.5) -  self._anglecalibration.angular_range[0]) // self._anglecalibration.histogram_bin_width)

        self._right_module_boundary_as_bin_index = lambda module_index, detector_angle : int((self._anglecalibration.diffraction_angle_from_DG_parameters( module_index, detector_angle, self._anglecalibration.MythenDetectorSpecifications.strips_per_module, +0.5) -  self._anglecalibration.angular_range[0]) // self._anglecalibration.histogram_bin_width)
        
        self._bin_to_diffraction_angle = lambda bin_index : bin_index * self._anglecalibration.histogram_bin_width + self._anglecalibration.angular_range[0]
    
        self._bin_to_diffraction_angle_base_peak_ROI_only = lambda bin_index : (bin_index * self._anglecalibration.histogram_bin_width - self._anglecalibration.base_peak_ROI_width + self._anglecalibration.base_peak_angle)

    @property 
    def overwrite_plot(self): 
        return self._overwrite_plot 
    
    @overwrite_plot.setter
    def overwrite_plot(self, value : bool):
        """if true, the next plot will overwrite the previous one"""
        self._overwrite_plot = value
        
    def plot_diffraction_pattern(self, photon_counts : np.ndarray, module_index=None, motor_position=None, axis=None): 
        """plot diffraction pattern, if module_index is given, only plot the part of the diffraction pattern that corresponds to the module, otherwise plot the whole diffraction pattern."""

        if(module_index is not None and motor_position is None):
            raise ValueError("if module_index is given, motor_position must also be given to determine the correct module region to plot.")

        plot_title = f"Diffraction Pattern for module {module_index}" if module_index is not None else "Diffraction Pattern"

        if axis is None:
            plt.figure(figsize=(8, 6))
            plt.xlabel("Diffraction Angle [degree]")
            plt.ylabel("Photon Counts")
            plt.title(plot_title)
        else:
            axis.set_xlabel("Diffraction Angle [degree]")
            axis.set_ylabel("Photon Counts")
            axis.set_title(plot_title)

        left_bin_boundary = 0
        right_bin_boundary = photon_counts.size
        if module_index is not None and (self._anglecalibration.num_fixed_angle_width_bins == photon_counts.size):
            left_bin_boundary = self._left_module_boundary_as_bin_index(
                module_index, motor_position)
            right_bin_boundary = self._right_module_boundary_as_bin_index(
                module_index, motor_position)
    
        bins = np.arange(left_bin_boundary, right_bin_boundary) 
        bins = np.vectorize(self._bin_to_diffraction_angle)(bins)

        if axis is None:
            plt.plot(bins, photon_counts[left_bin_boundary:right_bin_boundary])
            plt.show(block=False)
            input("Press Enter to continue...")
        else:
            if self._overwrite_plot:
                print("clearing axis for next plot")
                axis.clear()

            axis.plot(bins, photon_counts[left_bin_boundary:right_bin_boundary])
            if self._overwrite_plot:
                fig = axis.figure
                fig.canvas.draw()
                fig.show()  
                input("Press Enter to continue...")
            
    def plot_base_peak(self, photon_counts : np.ndarray, module_index=None, axis=None): 

        """plot base peak"""
        
        plot_title = f"Base Peak for module {module_index}" if module_index is not None else "Base Peak"

        if axis is None:
            plt.figure(figsize=(8, 6))
            plt.xlabel("Diffraction Angle [degree]")
            plt.ylabel("Photon Counts")
            plt.title(plot_title)
        else:
            axis.set_xlabel("Diffraction Angle [degree]")
            axis.set_ylabel("Photon Counts")
            axis.set_title(plot_title)

        left_bin_boundary = int((self._anglecalibration.base_peak_angle -
                                self._anglecalibration.base_peak_ROI_width -
                                self._anglecalibration.angular_range[0]) // self._anglecalibration.histogram_bin_width)
        right_bin_boundary = int((self._anglecalibration.base_peak_angle +
                                self._anglecalibration.base_peak_ROI_width -
                                self._anglecalibration.angular_range[0]) // self._anglecalibration.histogram_bin_width + 1)

        bins = np.arange(0, self._anglecalibration.base_peak_ROI_num_bins) 
        bins = np.vectorize(self._bin_to_diffraction_angle_base_peak_ROI_only)(bins)

        if axis is None:
            plt.plot(bins, photon_counts[left_bin_boundary:right_bin_boundary])
            plt.show(block=False)
            input("Press Enter to continue...")
        else:
            if self._overwrite_plot:
                axis.clear()
            axis.plot(bins, photon_counts[left_bin_boundary:right_bin_boundary])
            if self._overwrite_plot:
                plt.show(block=False)
                input("Press Enter to continue...")
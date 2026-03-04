from numpy import random 

from angcal import AngleCalibration, MythenDetectorSpecifications, FlatField, EpicsMythenFileReader

from angcal import DGParameters, BCParameters

anglecalibration = AngleCalibration(MythenDetectorSpecifications(), FlatField(MythenDetectorSpecifications()), EpicsMythenFileReader())

anglecalibration.read_initial_calibration_from_file("/home/mazzol_a/ANGCALDATA/angcal22_Jul25/angcal_Jul2025_P12_0p0105_original.off")

bcparameters = anglecalibration.BCparameters
bc_parameters_array = bcparameters.parameters()

#print("my base: ", bc_parameters_array.base)

bc_parameters_array[:,0] += random.choice([-1.0,1.0], size=bc_parameters_array.shape[0])*0.01 # angle_center_module_normal 

bc_parameters_array[:,1] += random.choice([-1.0,1.0], size=bc_parameters_array.shape[0])*0.05*bc_parameters_array[:,1] #modify module_center_sample_distances 

bc_parameters_array[:,2] += random.choice([-1.0,1.0], size=bc_parameters_array.shape[0])*0.01 # modify angle_center_beam

modified_bc_parameters =  bcparameters.parameters()# get modified parameters back as BCParameters object

dgparameters = DGParameters(bcparameters.num_modules()) 

bcparameters.convert_to_DGParameters(dgparameters) # update DG parameters accordingly

anglecalibration.write_DG_parameters_to_file("/home/mazzol_a/ANGCALDATA/angcal22_Jul25/angcal_Jul2025_P12_0p0105.off", dgparameters)



import numpy as np
import pandas as pd
from seidart.routines.definitions import *
from seidart.routines import prjrun, sourcefunction
from seidart.routines.definitions import Domain, Material, Model
from seidart.routines.arraybuild import Array
from glob2 import glob 
import os

# =============================================================================
class CommonOffset(Array):
    """
    Utilize the Array class to build a common offset image from a .prj file and 
    a source location file. If offsets are not provided then a receiver file 
    must be provided that coincides with the source file.

    The 
    """
    def __init__(
            self, source_file: str, 
            channel: str,
            prjfile: str, 
            receiver_file: str,
            receiver_indices: bool = False, 
            single_precision: bool = True,
            is_complex: bool = False,
            status_check: bool = False
        ):
        """
        Initialize the CommonOffset object. A receiver file must be provided 
        that coincides with the source file.

        :param source_file: The source file that contains the source locations
        :type source_file: str
        :param channel: The field direction to plot in either V[x,y,z] or E[x,y,z].  
        :type channel: str
        :param prjfile: The associated project file. 
        :type prjfile: str
        :param receiver_file: The receiver file that contains the receivers locations.
        :type receiver_file: str
        :param receiver_indices: If the receiver and source file are provided as indices then this needs to be flagged as True. Default to is in cartesian coordinates.
        :type receiver_indices: bool
        :param single_precision: 
        :type single_precision: bool
        :param status_check: Check to make sure the project file has all necessary parameters filled out. This will also recompute the tensor coefficients. 
        :type status_check: bool 

        """
        # super().__init__(source_file, *args, **kwargs)
        self.source_file = source_file
        self.prjfile = prjfile
        self.channel = channel
        self.receiver_file = receiver_file
        self.receiver_indices = receiver_indices
        self.single_precision = single_precision    
        self.is_complex = is_complex
        self.exaggeration = 0.5
        self.status_check = status_check
        self.output_basefile = None
        self.simulation_type = 'CO'
        self.build()

    def build(self) -> None:
        """
        Constructs domain and models based on project file and sets up source 
        and receiver configurations.
        """
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.prjfile,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        # Let's put a default gain value to not apply any gain to the 
        if 'E' in self.channel:
            self.gain = int(self.electromag.time_steps) 
        else:
            self.gain = int(self.seismic.time_steps)


        if self.channel in ['Vx','Vy','Vz']:
            self.is_seismic = True
            self.dt = self.seismic.dt 
            if self.status_check:
                prjrun.status_check(
                    self.seismic, self.material, self.domain, 
                    self.prjfile, append_to_prjfile = True
                )
        else:
            self.is_seismic = False
            self.dt = self.electromag.dt
            if self.status_check:
                prjrun.status_check(
                    self.electromag, self.material, self.domain, 
                    self.prjfile, append_to_prjfile = True
                )
        
         
        # Load the source and receiver locations
        self.source_receiver_xyz()

        # Copy the receiver_xyz to the variable receiver_xyz_all because we 
        # will modify the receiver_xyz variable
        self.receiver_xyz_all = self.receiver_xyz.copy()
        self.receiver_xyz = np.empty((0,3))

    def source_receiver_xyz(self):
        """
        Loads and sorts receiver and source locations from a file and adjusts them 
        according to the domain and CPML layer. 
        """
        self.source_xyz = pd.read_csv(self.source_file)
        self.receiver_xyz = pd.read_csv(self.receiver_file)
        
        # We need to make sure the receivers are ordered correctly and the 
        # absorbing boundary is corrected for
        # First check to see if the inputs are indices or
        cpml = int(self.domain.cpml)
        
        self.source_xyz = self.source_xyz.to_numpy() 
        self.receiver_xyz = self.receiver_xyz.to_numpy()
        if self.receiver_xyz.shape[1] == 1:
            self.receiver_xyz = self.receiver_xyz.T
        
        # We want to make sure the shape of the self.receiver_xyz array is in 
        # the correct shape. We won't be able to differentiate if the array is 
        # in correct shape if it is 3x3
        if self.receiver_xyz.shape[0] == 3 and np.prod(self.receiver_xyz.shape) > 9:
            self.receiver_xyz = self.receiver_xyz.T
        if self.receiver_xyz.shape[0] == 3 and np.prod(self.receiver_xyz.shape) == 6:
            self.receiver_xyz = self.receiver_xyz.T
        
        # If the receiver file contains values that are relative to the 
        # cartesian space of the domain, we want to change them to the indices
        # of the 
        if not self.receiver_indices:
            self.receiver_xyz = self.receiver_xyz / \
                np.array(
                    [
                        float(self.domain.dx), 
                        float(self.domain.dy), 
                        float(self.domain.dz) 
                    ]
                )
            self.receiver_xyz.astype(int)
        
        self.source_xyz = self.source_xyz
        self.receiver_xyz = self.receiver_xyz + cpml
    
    # -------------------------------------------------------------------------
    def co_run(self) -> None:
        """
        Run the electromagnetic or seismic model and extract the receiver 
        time series

        """
        self.timevec, self.fx, self.fy, self.fz, self.srcfn = sourcefunction(
            self.electromag, 1e7, 'gaus1'
        )
        n = len(self.source_xyz)
        self.co_image = np.zeros([len(self.timevec), n])
        self.co_image_complex = self.co_image.copy()

        # Loop through the source locations
        for i in range(n):
            self.receiver_xyz = self.receiver_xyz_all[i,:]
            self.source = self.source_xyz[i,:]

            print(f'Running shot gather for source location {self.source}')
            # Run the seismic or electromagnetic model
            if self.is_seismic:
                self.seismic.x = self.source[0]
                self.seismic.y = self.source[1]
                self.seismic.z = self.source[2]
                prjrun.runseismic(self.seismic, self.material, self.domain)
            else:
                self.electromag.x = self.source[0]
                self.electromag.y = self.source[1]
                self.electromag.z = self.source[2]
                prjrun.runelectromag(
                    self.electromag, 
                    self.material, 
                    self.domain, 
                    use_complex_equations = self.is_complex
                )
            
            # Extract the receiver time series
            self.domain.nx = self.domain.nx + 2 * int(self.domain.cpml)
            self.domain.nz = self.domain.nz + 2 * int(self.domain.cpml)
            if self.domain.dim == 2.5:
                self.domain.ny = self.domain.ny + 2 * int(self.domain.cpml)
            
            self.getrcx() 
            # Put a bandaid on the domain values 
            self.domain.nx = self.domain.nx - 2 * int(self.domain.cpml)
            self.domain.nz = self.domain.nz - 2 * int(self.domain.cpml)
            if self.domain.dim == 2.5:
                self.domain.ny = self.domain.ny - 2 * int(self.domain.cpml)
            
            # Remove all of the files 
            for file in glob('E*.0*.dat'):
                try:
                    os.remove(file) 
                except OSError as e:
                    print(f"Error removing file: {e.strerror}")
                    

            self.co_image[:,i] = self.timeseries
            # If the desired model is complex, extract the complex time series
            if self.is_complex:
                self.co_image_complex[:,i] = self.timeseries_complex
            
            # Overwrite the self.timeseries and self.timeseries_complex
            # variables wth the self.co_image variable
            self.timeseries = self.co_image
            self.timeseries_complex = self.co_image_complex
        
        if self.output_basefile:
            self.save(output_basefile = self.output_basefile)
        else:
            self.save() 


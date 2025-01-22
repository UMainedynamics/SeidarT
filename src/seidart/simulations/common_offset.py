import numpy as np
import pandas as pd
from seidart.routines.definitions import *
from seidart.routines.classes import Domain, Material, Model
from seidart.routines.arraybuild import Array
from glob2 import glob 
import subprocess
import os

from concurrent.futures import ThreadPoolExecutor

# =============================================================================
class CommonOffset(Array):
    """
    Utilize the Array class to build a common offset image from a .prj file and 
    a source location file. If offsets are not provided then a receiver file 
    must be provided that coincides with the source file.

    The 
    """
    def __init__(
            self, 
            source_file: str, 
            channel: str,
            project_file: str, 
            receiver_file: str,
            receiver_indices: bool = False, 
            single_precision: bool = True,
        ):
        """
        Initialize the CommonOffset object. A receiver file must be provided 
        that coincides with the source file.

        :param source_file: The source file that contains the source locations
        :type source_file: str
        :param channel: The field direction to plot in either V[x,y,z] or E[x,y,z].  
        :type channel: str
        :param project_file: The associated project file. 
        :type project_file: str
        :param receiver_file: The receiver file that contains the receivers locations.
        :type receiver_file: str
        :param receiver_indices: If the receiver and source file are provided as indices then this needs to be flagged as True. Default to is in cartesian coordinates.
        :type receiver_indices: bool
        :param single_precision: 
        :type single_precision: bool

        """
        # super().__init__(source_file, *args, **kwargs)
        self.source_file = source_file
        self.project_file = project_file
        self.channel = channel
        self.receiver_file = receiver_file
        self.receiver_indices = receiver_indices
        self.single_precision = single_precision    
        self.exaggeration = 0.5
        self.output_basefile = None
        self.simulation_type = 'CO'
        self.build()

    def build(self) -> None:
        """
        Constructs domain and models based on project file and sets up source 
        and receiver configurations.
        """
        self.domain, self.material, self.seismic, self.electromag = loadproject(
            self.project_file,
            Domain(), 
            Material(),
            Model(),
            Model()
        )
        # Let's put a default gain value to not apply any gain to the time series
        if 'E' in self.channel:
            self.gain = int(self.electromag.time_steps) 
        else:
            self.gain = int(self.seismic.time_steps)


        if self.channel in ['Vx','Vy','Vz']:
            self.is_seismic = True
            self.dt = self.seismic.dt 
            self.time_steps = self.seismic.time_steps
        else:
            self.is_seismic = False
            self.dt = self.electromag.dt
            self.time_steps = self.electromag.time_steps
        
        # Load the source and receiver locations
        self.source_receiver_xyz()

        # Copy the receiver_xyz to the variable receiver_xyz_all because we 
        # will modify the receiver_xyz variable
        self.receiver_xyz_all = self.receiver_xyz.copy()
        self.receiver_xyz = np.empty((0,3))
        
        self.source_xyz_all = self.source_xyz.copy() 
        self.source_xyz = np.empty((0,3))

    def source_receiver_xyz(self):
        """
        Loads and sorts receiver and source locations from a file and adjusts them 
        according to the domain and CPML layer. 
        """
        self.source_xyz = pd.read_csv(self.source_file)
        self.receiver_xyz = pd.read_csv(self.receiver_file)
        
        # We need to make sure the receivers are ordered correctly and the 
        # absorbing boundary is corrected for        
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
        
        # The source value is shifted in the Fortran file to accommodate the 
        # cpml, but the receivers are not so we must shift them now. 
        # self.source_xyz = self.source_xyz
        self.receiver_xyz = self.receiver_xyz + int(self.domain.cpml)
    
    
    # -------------------------------------------------------------------------
    def process_source(self, i):
        """
        Process a single source location.
        """
        self.receiver_xyz = self.receiver_xyz_all[i, :]
        self.source_xyz = self.source_xyz_all[i, :]
        
        json_parameterfile = '.'.join(self.project_file.split('.')[:-1])
        json_parameterfile = '.'.join([json_parameterfile,'-'.join(self.source_xyz.astype(str)), 'json'])
        
        print(f'Saving model inputs to {json_parameterfile}.')
        # Run the seismic or electromagnetic model
        if self.is_seismic:
            self.seismic.x = self.source_xyz[0]
            self.seismic.y = self.source_xyz[1]
            self.seismic.z = self.source_xyz[2]
            self.seismic.build(
                self.material, 
                self.domain, 
                recompute_tensors=False, 
                write_tensor = False
            )
            self.seismic.run(jsonfile = json_parameterfile)
        else:
            self.electromag.x = self.source_xyz[0]
            self.electromag.y = self.source_xyz[1]
            self.electromag.z = self.source_xyz[2]
            self.electromag.build(
                self.material, 
                self.domain, 
                recompute_tensors=False, 
                write_tensor = False
            )
            self.electromag.run(jsonfile = json_parameterfile)

        # Extract the receiver time series
        self.domain.nx = self.domain.nx + 2 * int(self.domain.cpml)
        self.domain.nz = self.domain.nz + 2 * int(self.domain.cpml)
        if self.domain.dim == 2.5:
            self.domain.ny = self.domain.ny + 2 * int(self.domain.cpml)
        
        self.receiver_xyz = self.receiver_xyz.reshape(1,3)
        self.getrcx()

        # Put a bandaid on the domain values
        self.domain.nx = self.domain.nx - 2 * int(self.domain.cpml)
        self.domain.nz = self.domain.nz - 2 * int(self.domain.cpml)
        if self.domain.dim == 2.5:
            self.domain.ny = self.domain.ny - 2 * int(self.domain.cpml)

        # Remove all of the files 
        source_ind = (self.source_xyz/np.array([self.domain.dx, self.domain.dy, self.domain.dz]) ).astype(int)
        if self.domain.dim == 2:
            pattern_input = '.'.join(source_ind.astype(str)[np.array([0,2])])
        else:
            pattern_input = '.'.join(source_ind.astype(str))
        
        if self.is_seismic:
            file_pattern = f'V[xyz].0*.{pattern_input}.dat'
        else:
            file_pattern = f'E[xyz].0*.{pattern_input}.dat'
        
        files_to_remove = glob(file_pattern)
        files_to_remove.append(json_parameterfile)
        
        if files_to_remove:
            print(f'Cleaning up files for {file_pattern}')
            subprocess.run(['rm', '-f'] + files_to_remove)
        
        # for file in glob(file_pattern):
        #     try:
        #         os.remove(file)
        #         os.remove(json_parameterfile)
        #     except OSError as e:
        #         print(f"Error removing file: {e.strerror}")

        return self.timeseries

    def co_run(self, parallel=False) -> None:
        """
        Run the electromagnetic or seismic model and extract the receiver time 
        series.

        This function loops through all source locations and runs the simulation
        for each source. It supports parallel execution to speed up processing.

        Parameters
        ----------
        parallel : bool, optional
            If True, runs the loop in parallel using a thread pool. If False,
            executes the loop sequentially. Default is True.

        Notes
        -----
        - In parallel mode, the function uses 
          `concurrent.futures.ThreadPoolExecutor` to process multiple source 
          locations concurrently.
        - Results from each simulation are stored in the `self.co_image` array,
          where each column corresponds to the time series for a specific 
          source.
        - File cleanup is performed after processing each source.

        See Also
        --------
        process_source : The function handling processing for each source 
        location.

        """
        n = len(self.source_xyz_all)
        self.co_image = np.zeros([self.time_steps, n])

        if parallel:
            # Parallel execution using ThreadPoolExecutor
            from concurrent.futures import ThreadPoolExecutor

            def process_source_wrapper(i):
                return self.process_source(i)

            with ThreadPoolExecutor() as executor:
                results = list(executor.map(process_source_wrapper, range(n)))

            # Update co_image using Numba for speed
            update_co_image(self.co_image, results)

        else:
            # Single-threaded execution
            for i in range(n):
                timeseries = self.process_source(i)
                self.co_image[:, i] = timeseries

        # Update self.timeseries and save results
        self.timeseries = self.co_image

        if self.output_basefile:
            self.save(output_basefile=self.output_basefile)
        else:
            self.save()


import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import mplstereonet 
from scipy.stats import norm, uniform, poisson, skewnorm
from scipy.spatial.transform import Rotation as R
import seaborn as sns

class Fabric:
    def __init__(self, params, plot = False):
        self.strikes = np.array([0]) 
        self.dips = np.array([0]) 
        self.density = None 
        self.xedges = None 
        self.yedges = None 
        self.contour_levels = None
        self.euler_angles = None 
        self.params = params
        self.npts = None
        self.plot = plot
        self.cmap_name = 'pink_r'
        self.cmap = None
        self.alpha = 0.3
        self.marker_size = 2
        self.figsize = (8,8)
        self.build()
    
    def build(self):
        self.custom_cmap()
        self.generate_strikes_dips()
        self.strike_dip_density()
        self.bunge_compute()
    
    def custom_cmap(self, n_colors = 100):
        palette = sns.color_palette(self.cmap_name, n_colors)
        self.cmap = sns.blend_palette(palette, as_cmap = True)
        self.cmap.set_under('white', 1.)
    
    # --------------------------------------------------------------------------
    def strike_dip_to_rotation_matrix(self, strike, dip):
        # Convert strike and dip to radians
        strike_rad = np.radians(strike)
        dip_rad = np.radians(dip)
        
        # Rotation matrix for strike (rotation around Z-axis)
        Rz = np.array([
            [np.cos(strike_rad), -np.sin(strike_rad), 0],
            [np.sin(strike_rad), np.cos(strike_rad), 0],
            [0, 0, 1]
        ])
        
        # Rotation matrix for dip (rotation around X-axis)
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(dip_rad), -np.sin(dip_rad)],
            [0, np.sin(dip_rad), np.cos(dip_rad)]
        ])
    
        # Combined rotation matrix
        rotation_matrix = Rz @ Rx
        return rotation_matrix
    
    # --------------------------------------------------------------------------
    def rotation_matrix_to_euler_angles(self, rotation_matrix):
        # Convert rotation matrix to Euler angles
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Input matrix must be 3x3")
        
        # Ensure the matrix is a valid rotation matrix
        if not np.allclose(
                np.dot(rotation_matrix, rotation_matrix.T), np.eye(3)
            ) or not np.isclose(np.linalg.det(rotation_matrix), 1):
            raise ValueError("Input matrix is not a valid rotation matrix")
        
        # Extract the Euler angles
        Phi = np.arccos(rotation_matrix[2, 2])
        phi1 = np.arctan2(rotation_matrix[0, 2], -rotation_matrix[1, 2])
        phi2 = np.arctan2(rotation_matrix[2, 0], rotation_matrix[2, 1])
        
        return phi1, Phi, phi2
    
    # ------------------------------------------------------------------------------
    # Function to generate strikes and dips based on a specified distribution
    def generate_strikes_dips(self):
        """
        Create a set of strikes and dips according to a distribution. 
        
        :param params: The parameters of the distribution.
        :type params: dict 
        """        
        n = len(self.params['npts'])
        for ind in range(n):
            if self.params['distribution'] == 'normal':
                self.strikes = np.hstack((
                    self.strikes,
                    skewnorm.rvs(
                        a=self.params['skew_strike'][ind], 
                        loc=self.params['strike'][ind], 
                        scale=self.params['strike_std'][ind], 
                        size=self.params['npts'][ind]
                    )
                ))
                self.dips = np.hstack((
                    self.dips, 
                    skewnorm.rvs(
                        a=self.params['skew_dip'][ind], 
                        loc=self.params['dip'][ind], 
                        scale=self.params['dip_std'][ind], 
                        size=self.params['npts'][ind]
                    )
                ))
            elif self.params['distribution'] == 'uniform':
                self.strikes = np.hstack((
                    self.strikes,
                    np.random.uniform(
                        self.params['strike_low'][ind], 
                        self.params['strike_high'][ind], 
                        self.params['npts'][ind]
                    )
                ))
                self.dips = np.hstack((
                    self.strikes,
                    np.random.uniform(
                        self.params['dip_low'][ind], 
                        self.params['dip_high'][ind], 
                        self.params['npts'][ind]
                    )
                ))
            elif self.params['distribution'] == 'poisson':
                self.strikes = np.hstack((
                    self.strikes,
                    poisson.rvs(
                        mu = self.params['lambda_strike'][ind], 
                        size = self.params['npts'][ind]
                    )
                ))
                self.dips = np.hstack((
                    self.strikes,
                    poisson.rvs(
                        mu = self.params['lambda_dip'][ind], 
                        size = self.params['npts'][ind]
                    )
                ))
        
        self.strikes = self.strikes[1:]
        self.dips = self.dips[1:]
    
    def strike_dip_density(self):
        self.xedges, self.yedges, self.density = mplstereonet.density_grid(
            self.strikes, self.dips, measurement = 'poles'
        )
        self.contour_levels = np.arange(
            0, np.ceil(self.density.max())
        )
    #     self.extent = [
    #         -np.pi/2, np.pi/2, -np.pi, np.pi
    #     ]

    #     self.extent = [
    #         self.xedges.min(), self.xedges.max(), 
    #         self.yedges.min(), self.yedges.max()
    #     ]
        
    # --------------------------------------------------------------------------
    def projection_plot(self, density_sigma = 3):
        """
        Create the stereonet plot 
        """
        # Plot the data using mplstereonet
        # self.fig = plt.figure(figsize = self.figsize) 
        # self.ax = self.fig.add_subplot(111, projection = 'stereonet')
        self.fig, self.ax = mplstereonet.subplots(figsize = self.figsize)
        self.cax = self.ax.density_contourf(
            self.strikes, self.dips, levels = self.contour_levels,
            measurement = 'poles', 
            cmap = self.cmap,
            sigma = density_sigma,
            gridsize = (100,100)
        )
        self.ax.pole(
            self.strikes, self.dips, 
            '.', c = '#050505', markersize = self.marker_size, 
            alpha = self.alpha 
        )
        # self.ax.set_azimuth_ticks(np.arange(0, 360, 10) )
        # self.ax.set_longitude_grid(10)
        self.ax.set_azimuth_ticklabels([])
        self.ax.set_xticks(np.radians(np.arange(-80, 90, 10)))
        self.ax.set_yticks(np.radians(np.arange(0,360,10)))
        self.ax.set_longitude_grid_ends(90)
        self.ax.set_rotation(0)
        self.ax._polar.set_position(self.ax.get_position() )
        # self.cax = self.ax.imshow(
        #     self.density#, extent = self.extent, origin = 'lower', cmap = self.cmap
        # )
        
        self.cbar = self.fig.colorbar(self.cax, shrink = 0.5, location = 'bottom' )
        self.ax.grid()
        plt.show()
    
    # --------------------------------------------------------------------------
    def bunge_compute(self):
        """
        
        """
        n = len(self.strikes)
        euler_zxz = np.zeros([n, 3])
        for ind in range(n):
            rotmat = self.strike_dip_to_rotation_matrix(self.strikes[ind], self.dips[ind])
            euler_zxz[ind,:] = self.rotation_matrix_to_euler_angles(rotmat)
        
        self.euler_angles = pd.DataFrame(euler_zxz, columns = ['z', 'x', 'z'])
        self.euler_angles.to_csv('euler_angles.csv')
        
        if self.plot:
            self.projection_plot()

# ==============================================================================
if __name__ == "__main__": 
    parser = argparse.ArgumentParser(
        description = """"""
    )
    parser.add_argument(
        '-t', '--distribution_type', nargs = 1, type = str, required = False, 
        default = ['uniform'], 
        help = """The name of the distribution. Options are: normal, uniform, 
        poisson, bimodal, multimodal. Default is normal."""
    )
    parser.add_argument(
        '-n', '--npts', nargs = '+', type = int, required = False, 
        default = [1000],
        help = """The number of points in the distribution. """
    )
    parser.add_argument(
        '-s', '--strike', nargs = '+', type = float, required = False, 
        default = [0], 
        help = """The strike mean(s) of the distribution (conditionally 
        optional). Bimodal and multimodal distributions require multiple 
        arguments. If the distribution is uniform, then this value becomes the 
        low strike limit to the distribution. Default is 0."""
    )
    parser.add_argument(
        '-S', '--strike2', nargs = 1, type = float, required = False,
        default = [360],
        help = """The high limit of the strike for the uniform distribution or 
        the standard deviation of the strike(s) for a normal distribution.
        Default is 360."""
    )
    parser.add_argument(
        '-d', '--dip', nargs = '+', type = float, required = False, 
        default = [0],
        help = """The dip means(s) of the distribution (conditionally optional). 
        Bimodal and multimodal distributions require multiple arguments. 
        Default is 0.""" 
    )
    parser.add_argument(
        '-D', '--dip2', nargs = 1, type = float, required = False,
        default = [90],
        help = """The high limit of the dip for the uniform distribution or the
        standard deviation of the dip(s) for a normal distribution. 
        Default is 90."""
    )
    parser.add_argument(
        '-k', '--strike_skewness_strike', nargs = '+', type = float, 
        required = False, default = [0], help = """The skewness of the 
        distribution of strike(s). Default is 0."""
    )
    parser.add_argument(
        '-K', '--dip_skewness', nargs = '+', type = float, required = False, 
        default = [0], 
        help = """The skewness of the distribution of dip(s). Default is 0."""
    ) 
    parser.add_argument(
        '-l', '--lambda_strike', nargs = '+', type = float, required = False, 
        default = [180],
        help = """The lambda parameter for the Poisson distribution of the 
        strike. Default is 180."""
    )
    parser.add_argument(
        '-L', '--lambda_dip', nargs = '+', type = float, required = False, 
        default = [45],
        help = """The lambda parameter for the Poisson distribution of the dip. 
        Default is 45."""
    )
    parser.add_argument(
        '-m', '--multimodal', required = False, action = 'store_true', 
        help = """Flag whether this is a mixed mode model. The number of modes 
        is determined by the length of the strike values."""
    )
    parser.add_argument('-P', '--plot', action = 'store_true', required = False,
        help = """Flag to plot the stereonet of strikes and dips."""
    )
    
    args = parser.parse_args()
    distribution = args.distribution_type[0]
    npts = args.npts
    strike = args.strike 
    dip = args.dip 
    strike2 = args.strike2
    dip2 = args.dip2
    skew_strike = args.skewness_strike 
    skew_dip = args.skewness_dip
    lambda_strike = args.lambda_strike 
    lambda_dip = args.lambda_dip 
    multimodal = args.multimodal 
    plot = args.plot 
    params = {
        'distribution':distribution,
        'multimodal': multimodal,
        'npts': npts
    }
    if distribution == 'normal':
        params['strike'] = strike 
        params['strike_std'] = strike2 
        params['dip'] = dip 
        params['dip_std'] = dip2
        params['skew_strike'] = skew_strike 
        params['skew_dip'] = skew_dip
    elif distribution == 'uniform':
        params['strike_low'] = strike 
        params['strike_high'] = strike2 
        params['dip_low'] = dip 
        params['dip_high'] = dip2 
    elif distribution == 'poisson':
        params['lambda_strike'] = lambda_strike 
        params['lambda_dip'] = lambda_dip
    else:
        raise ValueError("Unsupported distribution type")  
    
    Fabric(params, plot)
    
# Example usage
# distribution = 'bimodal'  # Change to 'normal', 'uniform', 'poisson', or 'bimodal'

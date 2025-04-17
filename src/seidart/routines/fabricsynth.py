import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap, BoundaryNorm
import mplstereonet 
from scipy.stats import norm, uniform, poisson, skewnorm
from scipy.spatial.transform import Rotation as R
import seaborn as sns
import copy

class Fabric:
    def __init__(self, params, output_filename = 'euler_angles.txt', plot = False):
        """
        Initialize the object. 
        
        :param params: A dictionary of the input parameters. 
        :type params: dict 
        :param output_filename: The name of the text file that will contain the 
            Euler angles from the given parameters. Default is 'euler_angles.txt'
        :type output_fielname: str
        :param plot: Generate a plot. If plotting is suppressed (False), it can 
            be plotted later using the class function 'projection_plot'
        
        """
        self.output_filename = output_filename
        self.trends = np.array([0]) # c-axes trends
        self.plunges = np.array([0]) # c-axes plunges
        self.a_trends = None # a-axes trends
        self.a_plunges = None # a-axes plunges
        self.orientations = np.array([0])
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
        '''
        Build the object. 
        '''
        self.custom_cmap()
        self.generate_trends_plunges()
        self.bunge_compute() 

    def custom_cmap(self, n_colors = 100):
        '''
        Compute a custom colormap for the density plot.
        
        :param n_colors: The number of colors to be used to compute the 
            colormap. Default is 100.
        :type n_colors: int
        '''
        palette = sns.color_palette(self.cmap_name, n_colors)
        self.cmap = sns.blend_palette(palette, as_cmap = True)
        self.cmap.set_under('white', 1.0)
    
    # --------------------------------------------------------------------------
    def trend_plunge_to_rotation_matrix(self, trend, plunge, orientation): 
        '''
        Compute the rotation matrix from a trend and plunge pair. 
        
        :param trend: The trend value in degrees 
        :type trend: float 
        :param plunge: The plunge value in degrees  
        :type plunge: float 
        :param orientation: The orientation/rotation angle around the axis of 
            the plunge.
        :type orientation: float
        
        :return: rotation_matrix
        :rtype: np.ndarray
        '''
        # Convert trend and plunge to radians
        trend_rad = np.radians(trend)
        plunge_rad = np.radians(plunge)
        orient_rad = np.radians(orientation)
        
        # Rotation matrix for trend (rotation around Z-axis)
        Rt = np.array([
            [np.cos(trend_rad), -np.sin(trend_rad), 0],
            [np.sin(trend_rad), np.cos(trend_rad), 0],
            [0, 0, 1]
        ])
        
        # Rotation matrix for plunge (rotation around X-axis)
        Rp = np.array([
            [1, 0, 0],
            [0, np.cos(plunge_rad), -np.sin(plunge_rad)],
            [0, np.sin(plunge_rad), np.cos(plunge_rad)]
        ])

        Ro = np.array([
            [np.cos(orient_rad), -np.sin(orient_rad), 0],
            [np.sin(orient_rad), np.cos(orient_rad), 0],
            [0, 0, 1]
        ])
        
        # Combined rotation matrix
        rotation_matrix = Rt @ Rp @ Ro
        return rotation_matrix
    
    # --------------------------------------------------------------------------
    def rotation_matrix_to_euler_angles(self, rotation_matrix):
        """
        Convert rotation matrix to Euler angles
        
        :param rotation_matrix: The 3-by-3 orthogonal rotation matrix
        :type rotation_matrix: np.ndarray 
        
        """
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
    # Function to generate trends and plunges based on a specified distribution
    def generate_trends_plunges(self):
        """
        Create a set of trends and plunges according to a distribution. 
        """        
        n = len(self.params['npts'])
        for ind in range(n):
            if self.params['distribution'] == 'normal':
                self.trends = np.hstack((
                    self.trends,
                    skewnorm.rvs(
                        a=self.params['skew_trend'][ind], 
                        loc=self.params['trend'][ind], 
                        scale=self.params['trend_std'][ind], 
                        size=self.params['npts'][ind]
                    )
                ))
                self.plunges = np.hstack((
                    self.plunges, 
                    skewnorm.rvs(
                        a=self.params['skew_plunge'][ind], 
                        loc=self.params['plunge'][ind], 
                        scale=self.params['plunge_std'][ind], 
                        size=self.params['npts'][ind]
                    )
                ))
                self.orientations = np.hstack((
                    self.orientations, 
                    skewnorm.rvs(
                        a=self.params['skew_orientation'][ind], 
                        loc=self.params['orientation'][ind], 
                        scale=self.params['orientation_std'][ind], 
                        size=self.params['npts'][ind]
                    )
                ))
            elif self.params['distribution'] == 'uniform':
                self.trends = np.hstack((
                    self.trends,
                    np.random.uniform(
                        self.params['trend_low'][ind], 
                        self.params['trend_high'][ind], 
                        self.params['npts'][ind]
                    )
                ))
                theta_low = np.radians(self.params['plunge_low'][ind])
                theta_high = np.radians(self.params['plunge_high'][ind])
                u = np.random.uniform(0, 1, self.params['npts'][ind])
                cos_theta = np.cos(theta_low) - u * (np.cos(theta_low) - np.cos(theta_high) )
                plunge_samples = np.degrees(np.arccos(cos_theta) )
                self.plunges = np.hstack((self.plunges, plunge_samples))
                self.orientations = np.hstack((
                    self.orientations,
                    np.random.uniform(
                        self.params['orientation_low'][ind], 
                        self.params['orientation_high'][ind], 
                        self.params['npts'][ind]
                    )
                ))
            elif self.params['distribution'] == 'poisson':
                self.trends = np.hstack((
                    self.trends,
                    poisson.rvs(
                        mu = self.params['lambda_trend'][ind], 
                        size = self.params['npts'][ind]
                    )
                ))
                self.plunges = np.hstack((
                    self.plunges,
                    poisson.rvs(
                        mu = self.params['lambda_plunge'][ind], 
                        size = self.params['npts'][ind]
                    )
                ))
                self.orientations = np.hstack((
                    self.orientations,
                    poisson.rvs(
                        mu = self.params['lambda_orientation'][ind], 
                        size = self.params['npts'][ind]
                    )
                ))
        
        self.trends = self.trends[1:]
        self.plunges = self.plunges[1:]
        self.orientations = self.orientations[1:]
    

    # --------------------------------------------------------------------------
    def projection_plot(
            self, 
            density_sigma = 3, 
            vmin = 1, 
            a_axes = False, 
            colorbar_location = 'bottom'
        ):
        """
        Create the stereonet plot. The axes and figure inputs are useful when 
        creating subplots since replacing empty axes objects with mplstereonet 
        axes objects is not straightforward. 
        
        :param density_sigma: The bandwidth parameter for the density plot. A 
            larger number uses a broader kernel and will smooth the data. 
            Default is 3.
        :type density_sigma: float 
         
        """
        if a_axes:
            trends = self.a_trends.copy() 
            plunges = self.a_plunges.copy()
        else:
            trends = self.trends.copy() 
            plunges = self.plunges.copy() 
        
        # Plot the data using mplstereonet
        self.fig, self.ax = mplstereonet.subplots(figsize = self.figsize)
        self.cax = self.ax.density_contourf(
            trends, plunges, levels = self.contour_levels,
            measurement = 'poles', 
            cmap = self.cmap,
            sigma = density_sigma,
            gridsize = (100,100),
            vmin = vmin
        )
        self.ax.pole(
            trends, plunges, 
            '.', c = '#050505', markersize = self.marker_size, 
            alpha = self.alpha 
        )
        
        # Customize ticks and grid as needed
        self.ax.set_azimuth_ticklabels([])
        self.ax.set_xticks(np.radians(np.arange(-80, 90, 10)))
        self.ax.set_yticks(np.radians(np.arange(0,360,10)))
        self.ax.set_longitude_grid_ends(90)
        self.ax.set_rotation(0)
        self.ax._polar.set_position(self.ax.get_position() )
        
        # Add the colorbar
        self.cbar = self.fig.colorbar(
            self.cax, shrink = 0.5, location = colorbar_location 
        )
        self.ax.grid()
        
        plt.show()
    
    # --------------------------------------------------------------------------
    def bunge_compute(self):
        """
        Compute the euler angles for the z-x-z transform (Bunge's notation). A 
        space delimited file will be produced from the output. 
        """
        
        n = len(self.trends)
        self.a_plunges = np.zeros([n])
        self.a_trends = np.zeros([n])
        a_local = np.array([1, 0, 0])
        euler_zxz = np.zeros([n, 3])
        
        for ind in range(n):
            rotmat = self.trend_plunge_to_rotation_matrix(self.trends[ind], self.plunges[ind], self.orientations[ind])
            euler_zxz[ind,:] = self.rotation_matrix_to_euler_angles(rotmat)
            a_global = rotmat @ a_local
            # Calculate trend and plunge from global coordinates
            x, y, z = a_global
            self.a_plunges[ind] = np.arcsin(z)
            self.a_trends[ind] = np.arctan2(y, x)
            
        # all a-axis plunges and trends need to be in degrees.
        self.a_plunges = np.degrees(self.a_plunges)
        self.a_trends = np.degrees(self.a_trends)
        
        self.euler_angles = pd.DataFrame(euler_zxz, columns = ['z', 'x', 'z'])
        # Write the space delimited text file. 
        self.euler_angles.to_csv(
            self.output_filename, index = False, header = True, sep = " "
        )
        
        if self.plot:
            self.projection_plot()

# --------------------------------------------------------------------------
def main():
    """
    Main function for console entry point. See help with 'fabricsynth -h'.
    """
    parser = argparse.ArgumentParser(
        description = """"""
    )
    parser.add_argument(
        '-t', '--distribution_type', nargs = 1, type = str, required = False, 
        default = ['uniform'], 
        help = """The name of the distribution. Options are: normal, uniform, 
        poisson. Bi/Multi-modal fabrics require multiple inputs, but the 
        distibution type is constant. Default is normal."""
    )
    parser.add_argument(
        '-n', '--npts', nargs = '+', type = int, required = False, 
        default = [1000],
        help = """The number of points in the distribution. """
    )
    parser.add_argument(
        '-s', '--trend', nargs = '+', type = float, required = False, 
        default = [0], 
        help = """The trend mean(s) of the distribution (conditionally 
        optional). Bimodal and multimodal distributions require multiple 
        arguments. If the distribution is uniform, then this value becomes the 
        low trend limit to the distribution. Default is 0."""
    )
    parser.add_argument(
        '-S', '--trend2', nargs = 1, type = float, required = False,
        default = [360],
        help = """The high limit of the trend for the uniform distribution or 
        the standard deviation of the trend(s) for a normal distribution.
        Default is 360."""
    )
    parser.add_argument(
        '-d', '--plunge', nargs = '+', type = float, required = False, 
        default = [0],
        help = """The plunge means(s) of the distribution (conditionally optional). 
        Bimodal and multimodal distributions require multiple arguments. 
        Default is 0.""" 
    )
    parser.add_argument(
        '-D', '--plunge2', nargs = 1, type = float, required = False,
        default = [90],
        help = """The high limit of the plunge for the uniform distribution or the
        standard deviation of the plunge(s) for a normal distribution. 
        Default is 90."""
    )
    parser.add_argument(
        '-k', '--trend_skewness_trend', nargs = '+', type = float, 
        required = False, default = [0], help = """The skewness of the 
        distribution of trend(s). Default is 0."""
    )
    parser.add_argument(
        '-K', '--plunge_skewness', nargs = '+', type = float, required = False, 
        default = [0], 
        help = """The skewness of the distribution of plunge(s). Default is 0."""
    ) 
    parser.add_argument(
        '-l', '--lambda_trend', nargs = '+', type = float, required = False, 
        default = [180],
        help = """The lambda parameter for the Poisson distribution of the 
        trend. Default is 180."""
    )
    parser.add_argument(
        '-L', '--lambda_plunge', nargs = '+', type = float, required = False, 
        default = [45],
        help = """The lambda parameter for the Poisson distribution of the plunge. 
        Default is 45."""
    )

    parser.add_argument('-P', '--plot', action = 'store_true', required = False,
        help = """Flag to plot the stereonet of trends and plunges."""
    )
    
    args = parser.parse_args()
    distribution = args.distribution_type[0]
    npts = args.npts
    trend = args.trend 
    plunge = args.plunge 
    trend2 = args.trend2
    plunge2 = args.plunge2
    skew_trend = args.skewness_trend 
    skew_plunge = args.skewness_plunge
    lambda_trend = args.lambda_trend 
    lambda_plunge = args.lambda_plunge 
    plot = args.plot 
    params = {
        'distribution':distribution,
        'npts': npts
    }
    if distribution == 'normal':
        params['trend'] = trend 
        params['trend_std'] = trend2 
        params['plunge'] = plunge 
        params['plunge_std'] = plunge2
        params['skew_trend'] = skew_trend 
        params['skew_plunge'] = skew_plunge
    elif distribution == 'uniform':
        params['trend_low'] = trend 
        params['trend_high'] = trend2 
        params['plunge_low'] = plunge 
        params['plunge_high'] = plunge2
        params['orientation_low'] = orientation 
        params['orientation_high'] = orientation2
    elif distribution == 'poisson':
        params['lambda_trend'] = lambda_trend 
        params['lambda_plunge'] = lambda_plunge
    else:
        raise ValueError("Unsupported distribution type")  
    
    Fabric(params, plot)

# ==============================================================================
if __name__ == "__main__": 
    main()
    
# Example usage
# distribution = 'bimodal'  # Change to 'normal', 'uniform', 'poisson', or 'bimodal'

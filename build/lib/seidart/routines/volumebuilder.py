import numpy as np 
from PIL import Image
from scipy.io import FortranFile
from scipy import ndimage as ndi
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D # noqa: F401
import struct 

n_slices = 200
files = np.array(
    [
        'slice1.png', 'slice2.png', 'slice3.png', 
        'slice2.png', 'slice3.png', 'slice4.png', 
        'slice3.png', 'slice4.png', 'slice5.png'
    ]
)
locs = np.array(
    [
        1, 20, 30, 
        50, 100, 110, 
        140, 170, 200
    ]
)


class VolumeBuilder:
    """
    Build a 3D labeled volume from a set of 2D images lying in a given plane
    (xy, yz, xz) and located at specific coordinates along the remaining axis.
    
    Conventions
    -----------
    - Internal 3D layout is (nx, ny, nz), corresponding to (x, y, z).
    - User provides 'locations' as 1-based slice indices along the orthogonal axis:
        plane='xy' -> locations are z-slices in [1..nz]
        plane='yz' -> locations are x-slices in [1..nx]
        plane='xz' -> locations are y-slices in [1..ny]
    - PNG images are assumed to be in (ny, nx, 3) pixel layout (row=y, col=x).
    """
    
    def __init__(self, 
                files, locations, 
                power = 1.0, plane='yz'
        ):
        """
        Parameters
        ----------
        files : list of str
            List of PNG file paths, in any order.
        locations : list or array of float
            1-based slice indices along the axis orthogonal to 'plane'.
            First must be 1.0, last must be n_slices.
        n_slices : int
            Number of slices along the orthogonal axis (including endpoints).
        plane : {'xy', 'yz', 'xz'}, optional
            Plane in which images lie. Default is 'yz':
              'xy' -> images define x,y; n_slices spans z  !! Height of image is the y-direction !!
              'yz' -> images define y,z; n_slices spans x
              'xz' -> images define x,z; n_slices spans y
        """
        self.files = list(files)
        self.locations = np.asarray(locations, dtype=int)
        self.n_slices = locations.max() 
        self.plane = plane.lower()
        self.weighted_average_power = float(power)
        
        if self.plane not in ('xy', 'yz', 'xz'):
            raise ValueError("plane must be one of 'xy', 'yz', 'xz'.")
        
        if len(self.files) != len(self.locations):
            raise ValueError("Files and locations must have the same length.")
        
        # Sort by user-provided locations
        order = np.argsort(self.locations)
        self.files = [self.files[i] for i in order]
        self.locations = self.locations[order]
        
        # Check endpoints in user convention: 1 .. n_slices
        if self.locations[0] != 1:
            raise ValueError("First location must be 1 (start of axis).")
        if self.locations[-1] != self.n_slices:
            raise ValueError("Last location must be n_slices (end of axis).")
        
        self.locations -= 1 # User input of locations starts with 1 and n_slices
        
        # Build integer label stack from images: rgb_int[height, width, Nimgs]
        self._image2int()
        
        # Interpolate along slices. Dims are now [height, weight, n_slices]
        self._weighted_average()
        
        # We need to reorder the arrays so that they are built in the correct 
        # directions
        if self.plane == 'xy':
            self.labels = np.transpose(self.labels, (2, 1, 0) )      # already (nx,ny,nz)
        elif self.plane == 'yz':
            pass
        elif self.plane == 'xz':
            self.labels = np.transpose(self.labels, (0, 2, 1) )
        else:
            raise ValueError("plane must be one of 'xy', 'yz', 'xz'.")
        
        f = FortranFile('geometry.dat', 'w')
        f.write_record( self.labels)
        f.close() 
    
    def _image2int(self):
        """
        Read all images and convert RGB combinations to integer labels.
        
        Builds:
        - self.rgb_int: shape (height, width, N_imgs) with integer codes
        - self.rgb:     array of unique RGB triplets (0-255)
        """
        # Insert each image at the correct internal slice index
        img = mpimg.imread(self.files[0])
        N = len(self.files)
        h, w, __ = img.shape 
        rgb_float = np.zeros((h, w, N))
        
        for ind in range(N):
            # read the image
            img = mpimg.imread(self.files[ind])
            # Convert RGB to a single value
            rgb_float[:,:,ind] = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])
        
        # Get the unique values of the image
        rgb_uni = np.unique(rgb_float)
        # We want the unique rgb values too
        rgb = np.zeros( [len(rgb_uni), 3] )
        rgb_int = np.zeros(rgb_float.shape)
        
        for ii in range(N):
            arr = rgb_float[:,:,ii]
            # reshape the image. We know it's three channels
            img_vect = np.zeros( [np.prod(arr.shape), 3] )
            img_vect[:,0] = np.reshape(img[:, :, 0], np.prod(np.shape(img[:, :, 0]) ) )
            img_vect[:,1] =	np.reshape(img[:, :, 1], np.prod(np.shape(img[:, :, 1]) ) )
            img_vect[:,2] =	np.reshape(img[:, :, 2], np.prod(np.shape(img[:, :, 2]) ) )
            
            for jj in range(0, len(rgb_uni) ):
                rgb_ind = np.reshape(arr == rgb_uni[jj], [np.prod(arr.shape)])
                rgb[jj,:] = (img_vect[rgb_ind,:])[0,:]
                rgb_int[ arr == rgb_uni[jj], ii ] = jj
        
        if np.max(rgb) <= 1.0:
            rgb = rgb * 255
            rgb = rgb.astype(int)
        
        rgb_int = rgb_int.astype(int)
        self.rgb_int = rgb_int 
        self.rgb = rgb
        
    def _weighted_average(self):
        nx, ny, __ = self.rgb_int.shape 
        labels = np.full((nx, ny, self.n_slices), -1, dtype=float)
        # labels = np.zeros((nx, ny, self.n_slices), dtype = int) - 1
        
        for ind in range(len(self.locations)-1 ):
            N = self.locations[ind+1] - self.locations[ind] 
            dist = np.linspace(0, N, N) / N 
            dist = dist**self.weighted_average_power
            image1 = self.rgb_int[:,:,ind]
            image2 = self.rgb_int[:,:,ind +1]
            for ii in range(self.locations[ind], self.locations[ind+1]):
                jj = ii - self.locations[ind]
                labels[:,:,ii] = self.interpolate_distance_transform(
                    image1, image2, dist[jj],
                )
        #         # labels[:,:,ii] = self.rgb_int[:,:,ind] * (1 - dist[jj]) + \
        #         #         self.rgb_int[:,:,ind+1] * dist[jj]
        
        labels[:,:,-1] = self.rgb_int[:,:,-1]
        self.labels = np.rint(labels).astype(int)
    
    def interpolate_distance_transform(self, image1, image2, t, spacing=(1.0, 1.0)):
        """
        Interpolate between two 2D labeled slices using signed distance transforms.
        
        Parameters
        ----------
        slice_A, slice_B : (H, W) int arrays
            Label images at positions z=0 and z=1.
        t : float in [0, 1]
            Interpolation parameter; 0 -> slice_A, 1 -> slice_B.
        spacing : tuple(float, float)
            Pixel spacing (dy, dx) for anisotropic grids, passed to EDT.
        
        Returns
        -------
        out : (H, W) int array
            Interpolated label image at position t.
        """
        if image1.shape != image2.shape:
            raise ValueError("image1 and image2 must have same shape.")
        
        labels = np.union1d(image1.ravel(), image2.ravel())
        labels = labels[labels != 0]   # optionally treat 0 as background
        
        H, W = image1.shape
        dt_A = np.zeros((len(labels), H, W), dtype=float)
        dt_B = np.zeros_like(dt_A)
        
        # build signed distance for each label in each slice
        for k, lab in enumerate(labels):
            mask_A = (image1 == lab)
            mask_B = (image2 == lab)
            
            # inside distances (positive)
            d_in_A = ndi.distance_transform_edt(mask_A, sampling=spacing)      # [web:45]
            d_in_B = ndi.distance_transform_edt(mask_B, sampling=spacing)      # [web:45]
            
            # outside distances (negative)
            d_out_A = ndi.distance_transform_edt(~mask_A, sampling=spacing)    # [web:45]
            d_out_B = ndi.distance_transform_edt(~mask_B, sampling=spacing)    # [web:45]
            
            dt_A[k] = d_in_A - d_out_A
            dt_B[k] = d_in_B - d_out_B
        
        # interpolate signed distance fields
        dt_t = (1.0 - t) * dt_A + t * dt_B
        
        # assign each pixel to the label with maximum signed distance
        max_idx = np.argmax(dt_t, axis=0)
        out = np.zeros((H, W), dtype=int)
        out[:, :] = labels[max_idx]
        
        return out
    
    def _label_from_rgb(self, rgb):
        """
        Map an (R,G,B) triplet to the corresponding integer label,
        using self.rgb (shape [n_unique, 3]) if available.
        """
        if not hasattr(self, 'rgb'):
            raise AttributeError("rgb not found; build rgb/rgb_int first.")
        rgb = tuple(int(c) for c in rgb)
        matches = np.all(self.rgb == np.array(rgb)[None, :], axis=1)
        idx = np.where(matches)[0]
        if idx.size == 0:
            raise ValueError(f"RGB {rgb} not found in unique color table.")
        return int(idx[0])
    
    def plot_scatter(self, label=None, rgb=None, dx=1.0, dy=1.0, dz=1.0, stride=1):
        """
        Plot a subset of voxels as a 3D scatter for a given label or RGB.
        
        Parameters
        ----------
        label : int or None
            Integer label to plot. If None, 'rgb' must be given.
        rgb : tuple or None
            (R,G,B) color to plot. If None, 'label' must be given.
        dx, dy, dz : float
            Physical spacing in x, y, z directions.
        stride : int
            Subsampling step along each axis to reduce plotted points.
        """
        if label is None and rgb is None:
            raise ValueError("Provide either 'label' or 'rgb'.")
        
        if label is None:
            label = self._label_from_rgb(rgb)
        
        vol = self.labels  # assumed shape (nx, ny, nz) or the Fortran shape you wrote
        nx, ny, nz = vol.shape
        
        xs = np.arange(0, nx, stride)
        ys = np.arange(0, ny, stride)
        zs = np.arange(0, nz, stride)
        
        sub = vol[np.ix_(xs, ys, zs)]
        mask = (sub == label)
        if not np.any(mask):
            print("No voxels with that label in the sampled grid.")
            return
        
        xi, yi, zi = np.where(mask)
        x = xs[xi] * dx
        y = ys[yi] * dy
        z = zs[zi] * dz
        
        # Color for scatter
        if hasattr(self, 'rgb'):
            rgb_val = self.rgb[label] / 255.0
            c = np.tile(rgb_val, (xi.size, 1))
        else:
            c = np.tile(np.array([[1.0, 0.0, 0.0]]), (xi.size, 1))
        
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=c, marker='s', s=5, linewidths=0)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        max_range = max((nx-1)*dx, (ny-1)*dy, (nz-1)*dz) / 2.0
        mid_x = (nx-1)*dx/2.0
        mid_y = (ny-1)*dy/2.0
        mid_z = (nz-1)*dz/2.0
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        plt.tight_layout()
        plt.show()
    
    def plot_voxels(self, label=None, rgb=None, dx=1.0, dy=1.0, dz=1.0, stride=1):
        """
        Plot a 3D voxel volume for a given label or RGB, downsampled.
        
        Parameters
        ----------
        label : int or None
            Integer label to plot. If None, 'rgb' must be given.
        rgb : tuple or None
            (R,G,B) color to plot. If None, 'label' must be given.
        dx, dy, dz : float
            Physical spacing in x, y, z directions.
        stride : int
            Subsampling step; >1 reduces number of voxels drawn.
        """
        if label is None and rgb is None:
            raise ValueError("Provide either 'label' or 'rgb'.")
        
        if label is None:
            label = self._label_from_rgb(rgb)
        
        vol = self.labels  # shape (nx, ny, nz)
        nx, ny, nz = vol.shape
        
        xs = np.arange(0, nx, stride)
        ys = np.arange(0, ny, stride)
        zs = np.arange(0, nz, stride)
        
        sub = vol[np.ix_(xs, ys, zs)]
        mask = (sub == label)
        if not np.any(mask):
            print("No voxels with that label in the sampled grid.")
            return
        
        nxp, nyp, nzp = sub.shape
        filled = mask
        
        # voxel colors: per-voxel RGBA
        if hasattr(self, 'rgb'):
            rgb_val = self.rgb[label] / 255.0
            facecolors = np.zeros((nxp, nyp, nzp, 4), dtype=float)
            facecolors[..., 0] = rgb_val[0]
            facecolors[..., 1] = rgb_val[1]
            facecolors[..., 2] = rgb_val[2]
            facecolors[..., 3] = 0.8  # alpha
        else:
            facecolors = None
        
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.voxels(filled, facecolors=facecolors, edgecolor='k', linewidth=0.1)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Note: dx,dy,dz scaling is visual, not applied to axes here;
        # for rough shape checking, this is usually OK. If you want
        # exact physical axes, you can rescale coordinates instead.
        ax.set_box_aspect((nxp*dx, nyp*dy, nzp*dz))
        
        plt.tight_layout()
        plt.show() 
    
    def plot_slice(self, label=None, rgb=None, ix=None, iy=None, iz=None,
               dx = 1.0, dy = 1.0, dz = 1.0, figure_size = (6,6), cmap='tab20'):
        """
        Plot orthogonal slices (x=const, y=const, z=const) through the 3D label volume.
        
        Parameters
        ----------
        label : int or None
            Integer label to highlight. If None and rgb is given, label is inferred
            from rgb. If both are None, all labels are shown.
        rgb : tuple or None
            (R,G,B) color to map to a label. Ignored if label is provided.
        ix, iy, iz : int or None
            Slice indices along x, y, z. If None, mid-plane is used.
            Indices are 0-based (Python convention).
        dx, dy, dz : float
            Physical spacing in x, y, z directions (used only for axis ticks).
        cmap : str
            Matplotlib colormap for imshow.
        """
        if not hasattr(self, 'labels'):
            raise AttributeError("self.labels not found; build/interpolate volume first.")
        
        vol = self.labels  # expected shape (nx, ny, nz) in your convention
        nx, ny, nz = vol.shape
        
        # Determine which label to show (optional highlighting)
        if label is None and rgb is not None:
            # Map RGB -> label using self.rgb
            if not hasattr(self, 'rgb'):
                raise AttributeError("rgb not found; run _image2int first.")
            rgb_arr = np.array(rgb, dtype=int)
            matches = np.all(self.rgb == rgb_arr[None, :], axis=1)
            idx = np.where(matches)[0]
            if idx.size == 0:
                raise ValueError(f"RGB {rgb} not found in unique color table.")
            label = int(idx[0])
        
        nz, ny, nx = self.labels.shape
        
        # Default slice indices: mid-planes
        if ix is not None:
            image_slice = self.labels[:,:,ix] 
            ratio = nz / ny
            ylab = 'Z'
            xlab = 'Y'
            extent=[0, ny*dy, 0, nz*dz]
            
        if iy is not None:
            image_slice = self.labels[:,iy,:]
            ratio = nz / nx
            ylab = 'Z'
            xlab = 'X'
            extent=[0, nx*dx, 0, nz*dz]
            
        if iz is not None:
            image_slice = self.labels[iz,:,:]
            ratio = ny / nx
            ylab = 'Y'
            xlab = 'X'
            extent=[0, nx*dx, 0, nz*dz]
            
        # Optionally convert to boolean mask for a single label
        if label is not None:
            image_slice = (image_sice == label).astype(int)
        
        fig, ax = plt.subplots( figsize=figure_size)
        im = ax.imshow(image_slice, cmap = cmap, extent = extent)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_aspect('auto')
        
        plt.tight_layout()
        plt.show()

vb = VolumeBuilder(files, locs)
vb.plot_slice(iy = 20)


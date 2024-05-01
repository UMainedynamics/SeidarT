#!/usr/bin/env python3
#
# This script will read an image and build the template project file template
# to be used in the seisarT program
#
# -----------------------------------------------------------------------------

import argparse
import numpy as np
import matplotlib.image as mpimg
from typing import Tuple

__all__ = ['prjbuild']

# ------------------------ Some Necessary Definitions -------------------------

def image2int(imfilename: str) -> Tuple[np.ndarray, np.ndarray]:
	"""
	Converts an image file to a 2D array of integer values representing unique
	RGB combinations and returns the unique RGB values.

	:param imfilename: The path to the image file.
	:type imfilename: str
	:return: A tuple containing the 2D array of integer values and the array of unique RGB values.
	:rtype: Tuple[np.ndarray, np.ndarray]
	"""
	# read the image
	img = mpimg.imread(imfilename)
	# Convert RGB to a single value
	rgb_int = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])
	# Get the unique values of the image
	rgb_uni = np.unique(rgb_int)
	# We want the unique rgb values too
	rgb = np.zeros( [len(rgb_uni), 3] )
	# reshape the image. We know it's three channels
	img_vect = np.zeros( [np.prod(rgb_int.shape), 3] )
	img_vect[:,0] = np.reshape(img[:, :, 0], np.prod(np.shape(img[:, :, 0]) ) )
	img_vect[:,1] =	np.reshape(img[:, :, 1], np.prod(np.shape(img[:, :, 1]) ) )
	img_vect[:,2] =	np.reshape(img[:, :, 2], np.prod(np.shape(img[:, :, 2]) ) )

	for ind in range(0, len(rgb_uni) ):
		rgb_ind = np.reshape(rgb_int == rgb_uni[ind], [np.prod(rgb_int.shape)])
		rgb[ind,:] = (img_vect[rgb_ind,:])[0,:]
		rgb_int[ rgb_int == rgb_uni[ind] ] = ind

	if np.max(rgb) <= 1.0:
		rgb = rgb * 255
		rgb = rgb.astype(int)

	return( rgb_int.astype(int), rgb)

def prjbuild(image_file: str, prjfile: str) -> None:
	"""
	Generates a project file (.prj) based on an input image file. This file
	contains domain parameters, material parameters, attenuation parameters,
	seismic parameters, and electromagnetic parameters derived from the image.

	:param image_file: The path to the input image file.
	:param prjfile: The path where the project file is to be saved.
	:type image_file: str
	:type prjfile: str
	"""
	# ----- Read the image file
	im, rgb = image2int(image_file)
	im = im.transpose()
	mat_id = np.unique(im)
	# Start writing the project file. To allow for headers we will start all
	# pertinant information after
	with open(prjfile, 'w') as prj:
		prj.write(header_comment)
		prj.write(new_line)
		prj.write('I,'+image_file+new_line)
		prj.write(new_line)
		
	# -------------------------------------------------------------------------
	# ------ Write domain parameters
	dim = 'D,dim,2'
	nx = 'D,nx,' + str(np.shape(im)[0])
	ny = 'D,ny,n/a'
	nz = 'D,nz,' + str(np.shape(im)[1])
	dx = 'D,dx,'
	dy = 'D,dy,n/a'
	dz = 'D,dz,'
	cpml = 'D,cpml,20'
	nmat = 'D,nmats,' + str(len( np.unique(im) ))
	tfile = 'D,tfile,'
	with open(prjfile, 'a') as prj:
		prj.write(dim+new_line)
		prj.write(nx+new_line)
		prj.write(ny+new_line)
		prj.write(nz+new_line)
		prj.write(dx+new_line)
		prj.write(dy+new_line)
		prj.write(dz+new_line)
		prj.write(cpml+new_line)
		prj.write(nmat+new_line)
		prj.write(tfile + new_line)
		prj.write(new_line)
	
	# -------------------------------------------------------------------------
	# ----- Write material parameters
	header = ("# number, id, R/G/B, Temperature, Density, Porosity, "
					"Water_Content, Anisotropic, ANG_File")
	
	with open(prjfile, 'a') as prj:
		i = 0
		prj.write(header + new_line )
		for x in mat_id:
			ln = ('M,' + str(x) + ',,' + str(rgb[x,0])  + '/' +
				str(rgb[x,1]) + '/' + str(rgb[x,2]) +
				',,,,,,,')
			prj.write( ln + new_line)
		prj.write(new_line)
	
	# -------------------------------------------------------------------------
	# ----- Write the attenuation parameters
	header = ("# number, Alpha1, Alpha2, Alpha3, fref")
	with open(prjfile, 'a') as prj:
		i = 0
		prj.write(header + new_line)
		for x in mat_id:
			ln = ('A,' + str(x) + ',,,,')
			prj.write(ln + new_line)
		prj.write(new_line)
	
	# -------------------------------------------------------------------------
	# ----- Write seismic parameters
	dt = 'dt,'
	steps = 'time_steps,'
	x = 'x,'
	y = 'y,'
	z = 'z,'
	f0 = 'f0,'
	theta = 'theta,0'
	phi = 'phi,0'
	source_file='source_file,'
	#
	comm = '# The source parameters for the seismic model'
	header = '# id, C11, C12, C13, C22, C23, C33, C44, C55, C66, rho'
	#
	with open(prjfile, 'a') as prj:
		i = 0
		prj.write(comm + new_line)
		prj.write('S,' + dt + new_line)
		prj.write('S,' + steps + new_line)
		prj.write('S,' + x + new_line)
		prj.write('S,' + y + new_line)
		prj.write('S,' + z + new_line)
		prj.write('S,' + f0 + new_line)
		prj.write('S,' + theta + new_line)
		prj.write('S,' + phi + new_line)
		prj.write(new_line)

		prj.write(header + new_line )
		for ind in mat_id:
			prj.write( 'C,' + str(ind) + ',,,,,,,,,,' + new_line)

		prj.write(new_line)
	
	# -------------------------------------------------------------------------
	# ----- Write EM Parameters 
	comm = '# The source parameters for the electromagnetic model'
	header = '# id, e11, e22, e33, s11, s22, s33'
	#
	with open(prjfile, 'a') as prj:
		i = 0
		prj.write(comm + new_line)
		prj.write('E,' + dt + new_line)
		prj.write('E,' + steps + new_line)
		prj.write('E,' + x + new_line)
		prj.write('E,' + y + new_line)
		prj.write('E,' + z + new_line)
		prj.write('E,' + f0 + new_line)
		prj.write('E,' + theta + new_line)
		prj.write('E,' + phi + new_line)

		prj.write(new_line)

		prj.write(header + new_line )
		for ind in mat_id:
			prj.write( 'P,' + str(ind) + ',,,,,,,,,,' + new_line)

		prj.write(new_line)

# -------------------------------- String Variables -------------------------------
new_line = '\n'
header_comment = """
# This is a project file template for the SeidarT software. In order to run the
# model for seismic, electromagnetic or both, the required inputs must be
#
# Domain Input Values:
#	dim 		- STR; either '2' or '2.5'; default is '2'
#	nx,ny,nz 	- INT; the dimensions of the image. If dim = 2.5, and ny is
#			  empty then default ny=1
#	dx,dy,dz	- REAL; the spatial step size for each dimension in meters. If
#			  dim = 2.5 and dy is empty then default dy=min(dx,dz)
#
# Material Input Values:
#	id 		- INT; the identifier given to each unique rgb value as it
#			  is read into the computer. It's recommended to use this
#			  script to make sure it's sorted correctly.
#	R/G/B 		- STR; the 0-255 values for each color code.
#	Temperature 	- REAL; temperature in Celsius.
#	Attenuation 	- REAL; (placeholder) will be attenuation length soon.
#	Density 	- REAL; density in kg/m^3
#	Porosity 	- REAL; percent porosity
#	Water_Content 	- REAL; percent of pores that contain water
#	Anisotropic 	- BOOL; whether the material is anisotropic (True) or
#			  isotropic (False).
#	ANG_File 	- STR; if Anisotrpic is True then the full path to the
#			  .ang file is supplied. The .ang file is a delimited text
#			  file that contains the 3-by-n array of euler rotation
#			  angles in radians.
#
#		or alternatively...
#	C11-C66 	- REAL; the stiffness coefficients with the appropriate id
#	E11-E33,S11-S33	- REAL; the permittivity and conductivity coefficients and
#			  'id' value corresponding to the coefficients along the diagonal
#			  of their respective tensors.
#
#
# Source Input Values:
#	dt 		- REAL; dx/(2*maxvelcity)
#	steps 		- INT; the total number of time steps
#	x,y,z 		- REAL; locations in meters, +x is to the right, +z is down, +y is into the screen
#	f0 		- REAL; center frequency for the guassian pulse function if
#			  'source_file' isn't supplied
#	theta 		- REAL; source orientation in the x-z plane,
#	phi 		- REAL; source orientation in the x-y plane for 2.5/3D only,
#	source_file	- STR; the pointer to the text file that contains the source
#			  timeseries as a steps-by-1 vector.
#
# 	**phi and theta are the rotation angles for spherical coordinates so
#		x = r sin(theta)cos(phi)
#		y = r sin(theta)sin(phi)
#		z = r cos(theta)
#
#	Theta is the angle from the z-axis (+ down relative to image), phi is the
#	angle from x-axis in the x-y plane (+ counterclockwise when viewed from above)
#
# Written by Steven Bernsen
# University of Maine
# -----------------------------------------------------------------------------

"""
def main(image_file, prjfile):
	prjbuild(image_file, prjfile)

# -------------------------- Command Line Arguments ---------------------------
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""The SeidarT software requires a
		.PNG image that is used to construct the model domain for seismic and
		electromagnetic wave propagation. Given the image file, a project file
		will be constructed which contains all the necessary parameters to be
		read in to the finite differences time domain modeling schemes.""" )
	
	parser.add_argument(
		'-i','--imagefile', 
		nargs=1, type=str, required = True,
		help='the full file path for the image', default=None
	)

	parser.add_argument(
		'-p', '--prjfile',
		nargs=1, type=str, required = False, default = 'jordan_downs.prj',
		help = """name of output file path with extension .prj and excluding
		the full path directory"""
	)

	# Get the arguments
	args = parser.parse_args()
	image_file = ''.join(args.imagefile)
	prjfile = ''.join(args.prjfile)
	main(image_file, prjfile)



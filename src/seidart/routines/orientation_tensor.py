#!/usr/bin/env python3

# Wrapper to generate the euler angles for the plunge and trend 
# then plot the results

import numpy as np
import argparse
import matplotlib.pyplot as plt
import mplstereonet

import seidart.fortran.orientsynth as ot



# -----------------------------------------------------------------------------
def	main(
		output_file: str, 
		npts: int, 
		plunge: float, trend: float, 
		angle_min: float, angle_max: float, 
		view_bool: bool
	) -> None:
	"""
	This function wraps the FORTRAN module orientsynth which synthesizes a 
	material fabric then saves outputs to a space delimited text file containing 
	the euler angles for z-x-z (bunge notation) that is correctly formatted for 
	input into the seidart program.
	
	:param output_file: Specify the output file that will be created to save the
		euler angles
	:type output_file: str
	:param npts: The number of crystal orientations within the substrate
	:type npts: int 
	:param plunge: the average plunge of the fabric in degrees
	:type plunge: float 
	:param trend: the average trend of the fabric in degrees
	:type trend: float 
	:param view_plot: Specify if a stereonet needs to be created.
	:type view_plot: bool
	"""
	euler_list= np.array([])
	orten = np.array([])

	euler_list, orten = ot.orientsynth(trend, plunge, amin, amax, npts)

	# Save the euler angles
	np.savetxt(output_file, euler_list, delimiter = " ")

	if view_plot:
		fig = plt.figure( figsize = (7, 5) )
		ax = fig.add_subplot(111,projection='stereonet')
		ax.pole(euler_list[:,0]*180/np.pi, euler_list[:,1]*180/np.pi, 'g^', markersize=6)
		ax.grid()

		# plt.rc('text', usetex=True)
		plt.figtext(0.01, 0.95, 'trend = ' + str(trend)   )
		plt.figtext(0.01, 0.9, 'plunge = ' + str(plunge) )
		plt.figtext(0.01, 0.85, 'angle min/max = ' + str(amin) + '/' + str(amax) )
		plt.figtext(0.01, 0.8, 'N = ' + str(npts) )

		# plt.figtext(0.01, 0.7, r'\underline{Tensor Coefficients}')
		plt.figtext(0.01, 0.7, 'Tensor Coefficients')
		plt.figtext(0.01, 0.65, '$a_{11} = $' + str( round( orten[0,0], 5) ) )
		plt.figtext(0.01, 0.6, '$a_{22} = $' + str( round( orten[1,1], 5) ) )
		plt.figtext(0.01, 0.55, '$a_{33} = $' + str( round( orten[2,2], 5) ) )

		plt.figtext(0.01, 0.5, '$a_{12} = $' +str( round( orten[0,1],5 ) ) )
		plt.figtext(0.01, 0.45, '$a_{13} = $' +str( round( orten[0,2], 5) ) )
		plt.figtext(0.01, 0.4, '$a_{23} = $' +str( round( orten[1,2], 5) ) )

		plt.show()




# -------------------------- Command Line Arguments ---------------------------
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=""" """ )

	parser.add_argument( '-o', '--outputfile', nargs=1, type=str,
							help='Specify the file to save the outputs', 
							default=None)

	parser.add_argument( '-n', '--npts', nargs = 1, type = int, required = False,
		help = """ Total number of grains in synthetic sample.""", default = [100])

	parser.add_argument( '-P', '--plunge', nargs = 1, type = float, 
		required = True, help = """ Plunge angle in degrees.""")

	parser.add_argument( '-t', '--trend', nargs =1, type = float, required = True,
		help = """ Trend angle in degrees.""")

	parser.add_argument( '-a', '--anglemin', nargs = 1, type = float, 
		required = True, help = """ Minimum angle deviation.""")

	parser.add_argument( '-A', '--anglemax', nargs = 1, type = float, 
		required = True, help = """Maximum angle deviation. """)

	parser.add_argument('-v', '--view_plot', action='store_true', 
		help = """Flag if you would like a plot to be shown""")
	
	# Get the arguments
	args = parser.parse_args()
	output_file = ''.join(args.outputfile)
	npts = args.npts[0]
	plunge=args.plunge[0]
	trend=args.trend[0]
	angle_min=args.anglemin[0]
	angle_max=args.anglemax[0]
	view_plot = args.view_plot
	
	main(output_file, npts, plunge, trend, angle_min, angle_max, view_plot)
	
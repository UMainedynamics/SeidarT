import numpy as np 
import pandas as pd 
from seidart.routines.definitions import * 
from seidart.routines import prjrun, sourcefunction 
from seidart.routines.arraybuild import Array
from seidart.simulations import CommonOffset
from glob2 import glob
import os

import numpy as np
from sklearn.cluster import DBSCAN
import numpy as np
from scipy.optimize import minimize

# =============================================================================
class MultiOffset(CommonOffset):
    def multioffset_run(self):
        '''
        Simulate a multi-offset radar survey. The outputs for each source location are saved in a m-by-n-by-p where m is the length of the time series, n is the number of source locations, and p is the number of receivers. 
        
        The receiver file and source file must be equal length. The number of survey points is going to be the number of source locations minus the number of receivers in the streamer. 
        '''
        pass

    def compute_nmo_correction(data, dt, offsets, velocities, n_estimators=10):
        """
        Use an ensemble method to compute the normal moveout (NMO) correction of seismic data.

        Parameters:
        data (np.array): The seismic data to correct. This should be a 2D array where each row is a trace and each column is a time sample.
        dt (float): The sampling interval in seconds.
        offsets (np.array): The offsets of the traces in meters.
        velocities (np.array): The NMO velocities in m/s.
        n_estimators (int): The number of base estimators in the ensemble.

        Returns:
        corrected_data (np.array): The NMO corrected data.
        """

        # Calculate the time samples
        t = np.arange(data.shape[1]) * dt

        # Initialize the corrected data
        corrected_data = np.zeros_like(data)

        # Loop over each trace
        for i in range(data.shape[0]):
            # Calculate the NMO correction
            t0 = np.sqrt(t**2 + (offsets[i]/velocities[i])**2)

            # Initialize the base estimator
            base_estimator = LinearRegression()

            # Initialize the ensemble method
            ensemble = BaggingRegressor(base_estimator=base_estimator, n_estimators=n_estimators)

            # Fit the ensemble to the data
            ensemble.fit(t.reshape(-1, 1), data[i])

            # Predict the corrected data
            corrected_data[i] = ensemble.predict(t0.reshape(-1, 1))

        return corrected_data


    def plot_nmo_performance_metrics(data, dt, offsets, velocities, ground_truth, n_estimators_values):
        """
        Plot the performance metrics of the compute_nmo_correction function.

        Parameters:
        data (np.array): The seismic data to correct. This should be a 2D array where each row is a trace and each column is a time sample.
        dt (float): The sampling interval in seconds.
        offsets (np.array): The offsets of the traces in meters.
        velocities (np.array): The NMO velocities in m/s.
        ground_truth (np.array): The ground truth data to compare the corrected data to.
        n_estimators_values (list of int): The values of n_estimators to calculate the performance metrics for.
        """

        # Initialize the RMSE values
        rmse_values = []

        # Loop over each value of n_estimators
        for n_estimators in n_estimators_values:
            # Compute the NMO correction
            corrected_data = compute_nmo_correction(data, dt, offsets, velocities, n_estimators)

            # Calculate the RMSE
            rmse = np.sqrt(mean_squared_error(ground_truth, corrected_data))

            # Append the RMSE to the RMSE values
            rmse_values.append(rmse)

        # Create a figure and an axis
        fig, ax = plt.subplots()

        # Plot the RMSE values
        ax.plot(n_estimators_values, rmse_values, label='RMSE')

        # Add labels and a legend
        ax.set_xlabel('Number of base estimators')
        ax.set_ylabel('RMSE')
        ax.legend()

        # Show the plot
        plt.show()

class LandStreamer(CommonOffset):
    def masw(self, data, sampling_rate):
        """
        Perform Multichannel Analysis of Surface Waves (MASW) on the given data.

        Parameters:
        data (np.array): The seismic data to analyze. This should be a 2D array where each row is a channel and each column is a time point.
        sampling_rate (float): The sampling rate of the data in Hz.

        Returns:
        dispersion_curve (np.array): The estimated dispersion curve. This is a 2D array where the first column is frequency and the second column is velocity.
        """

        # Number of channels and time points
        num_channels, num_time_points = data.shape

        # Calculate the Fourier transform of the data
        fft_data = np.fft.rfft(data, axis=1)

        # Calculate the frequencies corresponding to the Fourier transform
        frequencies = np.fft.rfftfreq(num_time_points, 1/sampling_rate)

        # Initialize the dispersion curve
        dispersion_curve = np.zeros((len(frequencies), 2))
        dispersion_curve[:, 0] = frequencies

        # Loop over frequencies
        for i, freq in enumerate(frequencies):
            # Calculate the cross-spectrum matrix
            cross_spectrum = np.outer(fft_data[:, i], fft_data[:, i].conj())

            # Calculate the eigenvalues and eigenvectors of the cross-spectrum matrix
            eigenvalues, eigenvectors = np.linalg.eig(cross_spectrum)

            # The velocity is the square root of the largest eigenvalue divided by the frequency
            velocity = np.sqrt(np.max(eigenvalues)) / freq

            # Store the velocity in the dispersion curve
            dispersion_curve[i, 1] = velocity

        return dispersion_curve


    def pick_dispersion_curves(self, dispersion_curve, eps=0.5, min_samples=5):
        """
        Use DBSCAN to pick the dispersion curve(s) from the MASW output.

        Parameters:
        dispersion_curve (np.array): The MASW output. This should be a 2D array where the first column is frequency and the second column is velocity.
        eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood.
        min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point.

        Returns:
        labels (np.array): The labels of the clusters. Noise points are given the label -1.
        """

        # Initialize the DBSCAN algorithm
        db = DBSCAN(eps=eps, min_samples=min_samples)

        # Fit the DBSCAN algorithm to the data
        db.fit(dispersion_curve)

        # Get the labels of the clusters
        labels = db.labels_

        return labels


    def tikhonov_inversion(self, dispersion_curve, alpha=0.1):
        """
        Use Tikhonov regularization to estimate the velocity model from the dispersion curve.

        Parameters:
        dispersion_curve (np.array): The dispersion curve. This should be a 2D array where the first column is frequency and the second column is velocity.
        alpha (float): The regularization parameter.

        Returns:
        velocity_model (np.array): The estimated velocity model.
        """

        # Define the forward model
        def forward_model(velocity_model):
            # In this simplified example, the forward model is a linear function of the velocity model
            return np.dot(velocity_model, dispersion_curve[:, 0])

        # Define the objective function
        def objective(velocity_model):
            # The objective function is the sum of the squared residuals plus the regularization term
            residuals = forward_model(velocity_model) - dispersion_curve[:, 1]
            regularization = alpha * np.sum(np.abs(np.gradient(velocity_model)))
            return np.sum(residuals**2) + regularization

        # Initialize the velocity model
        velocity_model = np.ones_like(dispersion_curve[:, 1])

        # Minimize the objective function
        result = minimize(objective, velocity_model)

        # Return the optimized velocity model
        return result.x

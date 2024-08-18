import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def kalman_filter_3d(csv_file, q, r, x0, p0, dt):
    """
    Parameters:
    csv_file (str): path to CSV file containing measurements
    q : process noise covariance
    r : measurement noise covariance
    x0: initial state estimate (3D)
    p0: initial state covariance (3x3)
    dt: time step

    Returns:
    x_upd : updated state estimate (3D)
    p_upd : updated state covariance (3x3)
    """
    # Load data from CSV file
    data = pd.read_csv("D:\Code\data.csv")
    
    

    # Initialize state estimate and covariance
    x_upd = x0
    p_upd = p0

    # Define system matrices
    F = np.array([[1, dt, 0.5 * dt**2],
                  [0, 1, dt],
                  [0, 0, 1]])
    H = np.eye(3)

    # Create arrays to store the updated state and covariance
    x_upd_array = np.zeros((len(data), 3))
    p_upd_array = np.zeros((len(data), 3, 3))

    # Process the data
    for i in range(len(data)):
        # Get the measurement
        z = np.array([data['ax'][i]*16384, data['ay'][i]*16384, data['az'][i]*16384])

        # Predict
        x_pred = F @ x_upd
        p_pred = F @ p_upd @ F.T + q * np.eye(3)

        # Update
        y = z - H @ x_pred
        s = H @ p_pred @ H.T + r * np.eye(3)
        k = p_pred @ H.T @ np.linalg.inv(s)
        x_upd = x_pred + k @ y
        p_upd = (np.eye(3) - k @ H) @ p_pred

        # Store the updated state and covariance
        x_upd_array[i] = x_upd
        p_upd_array[i] = p_upd

    # Plot the raw and filtered data
    plt.figure(figsize=(12, 6))
    plt.subplot(3, 1, 1)
    plt.plot(data['ax'], label='Raw X')
    plt.plot(x_upd_array[:, 0], label='Filtered X')
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(data['ay'], label='Raw Y')
    plt.plot(x_upd_array[:, 1], label='Filtered Y')
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(data['az'], label='Raw Z')
    plt.plot(x_upd_array[:, 2], label='Filtered Z')
    plt.legend()
    plt.tight_layout()
    plt.show()


    return x_upd_array, p_upd_array


csv_file = 'D:\Code\data.csv'
q = 1e-2  # process noise covariance
r = 1e-1  # measurement noise covariance
x0 = np.array([0, 0, 0])  # initial state estimate
p0 = np.eye(3)  # initial state covariance
dt = 0.0  # time step

x_upd, p_upd = kalman_filter_3d(csv_file, q, r, x0, p0, dt)
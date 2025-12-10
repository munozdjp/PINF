import numpy as np

def add_noise(data, noise_level=0.1):
    """
    This function adds Gaussian noise to the data. The noise is defined as a Gaussian distribution with mean 0 and
    standard deviation equal to noise_level times the standard deviation of the data.

    Parameters:
    data: numpy array, the original data
    noise_level: float, the level of noise to add, defined as a proportion of the data's standard deviation

    Returns:
    noisy_data: numpy array, the original data with added Gaussian noise
    """
    # Calculate the standard deviation of the data
    data_std = np.std(data)

    # Generate Gaussian noise
    noise = np.random.normal(loc=0, scale=noise_level * data_std, size=data.shape)

    # Add the noise to the data
    noisy_data = data + noise

    return noisy_data
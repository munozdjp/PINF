# Author: Juan Diaz Living System Lab, KAUST.
# License: BSD 3 clause
import scipy.io as sio
import os
# Standard scientific Python imports
import matplotlib.pyplot as plt
import numpy as np
# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, metrics
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import accuracy_score
import pickle
import time  # import time module


#Importing my Matrix Data
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier

# Define the classifiers to test
classifiers = {
    "MLP": MLPClassifier(random_state=1, max_iter=300),
    # "SVM": SVC(kernel='poly', gamma = 'auto'),  # probability=True to enable predict_proba method
    # "SVM": SVC(kernel='poly', gamma = 'auto', probability=True),  # probability=True to enable predict_proba method
    "Random Forest": RandomForestClassifier(max_depth=20, random_state=0),
    # "Logistic Regression": LogisticRegression(),
    "K-Nearest Neighbors": KNeighborsClassifier(n_neighbors=3)  # Hyperparam: number of neighbors
}

# Initialize a dictionary to store the accuracies of each classifier at each noise level
all_accuracies = {name: [] for name in classifiers.keys()}


# mat_contents=sio.loadmat('Data_And_labels.mat')
mat_contents=sio.loadmat('./DataSets/Data_And_labels_4sing_1Mixed_Biff.mat')
numpy_arrayData=mat_contents['matrix_with_data']

numpy_arrayLabels=mat_contents['allLabels']

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

# Define the dataset sizes to test
# dataset_sizes = [0.0005, 0.001, 0.002, 0.005, .01, .02]
dataset_sizes = [0.0005, 0.001, 0.002, 0.005, .01, .02]
# factor = 1/10
# dataset_sizes = [x * factor for x in dataset_sizes]
# dataset_sizes = [0.0005, 0.001, 0.002, 0.005, .01, .02]*1/100
# dataset_sizes = [0.0005, 0.001]
# dataset_sizes = [0.0005, 0.001, 0.002, 0.005, .01, .02]
# Define the noise levels to test
# noise_levels = np.linspace(0, 0.2, num=3)
noise_levels = np.linspace(0, 0, num=1)
# noise_levels = np.linspace(0, 0, num=1)
# Initialize a list to store the accuracy at each noise level
accuracies = []
# Loop over the dataset sizes



start_time = time.time()


for dataset_size in dataset_sizes:
    # Calculate the size of the reduced dataset
    reduced_size = int(dataset_size * numpy_arrayData.shape[0])

    # Randomly sample the data
    random_indices = np.random.choice(numpy_arrayData.shape[0], size=reduced_size, replace=False)
    reduced_data = numpy_arrayData[random_indices]
    reduced_labels = numpy_arrayLabels[random_indices]

    # Loop over the noise levels
    for noise_level in noise_levels:
        # Add noise to the data
        noisy_data = add_noise(reduced_data, noise_level=noise_level)

        # Split the noisy data into a training set and a test set
        X_train, X_test, y_train, y_test = train_test_split(noisy_data, reduced_labels, test_size=0.3, shuffle=True)

        # Loop over the classifiers
        for name, clf in classifiers.items():
            loop_start_time = time.time()

            # Train the classifier on the training set
            clf.fit(X_train, y_train.ravel())

            # Predict the labels of the test set
            predicted = clf.predict(X_test)

            # Calculate the accuracy of the predictions
            accuracy = metrics.accuracy_score(y_test, predicted)

            # Store the accuracy
            all_accuracies[name].append(accuracy)
            loop_end_time = time.time()
            print(f"Elapsed time for {name} at noise level {noise_level} and data size {dataset_size}: ",loop_end_time - loop_start_time, "seconds")
end_time = time.time()
print("Total elapsed time: ", end_time - start_time, "seconds")

# with open('accu_Downsampl_Compar_Noise.pkl', 'wb') as f:
#     pickle.dump(all_accuracies, f)
#for laoding accuracies
pickle_file_path = 'pickle_files/accuracies_Downsampling_Comparison_SmallDataset.pkl'

# Load the dictionary from the pickle file
with open(pickle_file_path, 'rb') as f:
    all_accuracies = pickle.load(f)
# Plot the relationship between noise level and accuracy for each classifier

for name, accuracies in all_accuracies.items():
    plt.plot(dataset_sizes, accuracies, marker='o', label=name)
plt.xlabel('Dataset size')
plt.ylabel('Accuracy')
plt.title('Effect of dataset size on classification accuracy')
plt.grid(True)
plt.legend()
# plt.savefig('Downsampling_Accuracy_2_plot.png', dpi=300)

save_directory = "figures_results_deep_learning"
# Ensure the directory exists
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# File path for saving the figure
save_path = os.path.join(save_directory, 'Downsampling_Accuracy_2_plot.png')

# Save the figure
plt.savefig(save_path)

# Show the plot
plt.show()
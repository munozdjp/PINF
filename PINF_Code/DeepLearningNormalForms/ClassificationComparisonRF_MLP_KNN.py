# Author: Juan Diaz Living System Lab, KAUST.

# License: BSD 3 clause
import scipy.io as sio
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
import sys  # Make sure to import sys
import os
#Importing my Matrix Data
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
sys.path.append('AdditionalpythonFunctions')
from descriptive_initials import get_descriptive_initials
from noise_level_addition import add_noise
import datetime
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


mat_contents=sio.loadmat('DataSets/Data_And_labels.mat')
numpy_arrayData=mat_contents['matrix_with_data']

numpy_arrayLabels=mat_contents['allLabels']
# Reducing the dataset size to 20%

reduced_size = int(0.02 * numpy_arrayData.shape[0])
# reduced_size = int(0.001 * numpy_arrayData.shape[0])

# Randomly sampling 20% of the data
random_indices = np.random.choice(numpy_arrayData.shape[0], size=reduced_size, replace=False)
reduced_data = numpy_arrayData[random_indices]
reduced_labels = numpy_arrayLabels[random_indices]



# Define the dataset sizes to test

# Define the noise levels to test
noise_levels = np.linspace(0, 0.5, num=6)
# noise_levels = np.linspace(0, 0.5, num=6)
# Initialize a list to store the accuracy at each noise level
accuracies = []

# Loop over the noise levels
# ...[previous code]...

# Loop over the noise levels
for noise_level in noise_levels:
    # Initialize a dictionary to store the accuracies of each classifier for each repetition at the current noise level
    accuracies_at_current_noise_level = {name: [] for name in classifiers.keys()}

    # Repeat the simulation 10 times for each noise level
    for _ in range(5):
        # Add noise to the data
        noisy_data = add_noise(reduced_data, noise_level=noise_level)

        # Split the noisy data into a training set and a test set
        X_train, X_test, y_train, y_test = train_test_split(noisy_data, reduced_labels, test_size=0.3, shuffle=True)

        # Loop over the classifiers
        for name, clf in classifiers.items():
            # Train the classifier on the training set
            clf.fit(X_train, y_train.ravel())

            # Predict the labels of the test set
            predicted = clf.predict(X_test)
            print(f"Classification using {name} finished at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            # Calculate the accuracy of the predictions
            accuracy = metrics.accuracy_score(y_test, predicted)

            # Store the accuracy
            accuracies_at_current_noise_level[name].append(accuracy)

    # Store the accuracies of each repetition at the current noise level
    for name in classifiers.keys():
        all_accuracies[name].append(accuracies_at_current_noise_level[name])

with open('accuraciesSlide0002.pickle', 'wb') as handle:
    pickle.dump(all_accuracies, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
# Plot the box plots for each classifier at each noise level
plt.figure(figsize=(10, 6))
for i, (name, accuracies) in enumerate(all_accuracies.items(), start=1):
    plt.subplot(1, len(classifiers), i)
    plt.boxplot(accuracies, labels=[f"{level:.1f}" for level in noise_levels])
    plt.title(name)
    plt.xlabel('Noise level')
    plt.ylabel('Accuracy')
    plt.grid(True)
plt.tight_layout()

# Directory where the figure should be saved
save_directory = "figures_results_deep_learning"
# Ensure the directory exists
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# File path for saving the figure
save_path = os.path.join(save_directory, 'AccuracyMLP_RF_KNN.png')

# Save the figure
plt.savefig(save_path)

# Show the plot
plt.show()


# n_classifiers = len(classifiers)
# n_noise_levels = len(noise_levels)
#
# # Define a small spacing to move boxplots of each classifier
# spacing = 0.15
#
# # Create figure
# plt.figure(figsize=(12, 6))
#
# # Adjust data format for boxplot and plot boxplots
# for i, (name, accuracies) in enumerate(all_accuracies.items()):
#     # For each classifier, we'll adjust the x positions of the boxplots
#     positions = [(j+1) + i*spacing for j in range(n_noise_levels)]
#     plt.boxplot(accuracies, positions=positions, widths=spacing)
#
# # Setting up the x-ticks and labels
# tick_positions = np.arange(1, n_noise_levels+1) + (spacing * (n_classifiers-1) / 2)  # Center the x-ticks for all classifiers
# tick_labels = [f"{level:.1f}" for level in noise_levels]
#
# # Extend the tick labels to incorporate classifier names
# extended_tick_labels = []
# for level in tick_labels:
#     for name in classifiers.keys():
#         extended_tick_labels.append(f"{name} ({level})")
#
# plt.xticks(np.arange(1, n_classifiers*n_noise_levels + 1), extended_tick_labels, rotation=45, ha='right')
#
# # Labels, title and grid
# plt.xlabel('Classifier (Noise level)')
# plt.ylabel('Accuracy')
# plt.grid(True, which='both', axis='y')
# plt.title('Classifier accuracy at different noise levels')
# plt.tight_layout()
# plt.show()


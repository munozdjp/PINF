# Author: Juan Diaz Living System Lab, KAUST.
# License: BSD 3 clause
import scipy.io as sio
import numpy as np
from sklearn import datasets, svm, metrics
from sklearn.model_selection import train_test_split
import os
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import accuracy_score
#Importing my Matrix Data
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier

import pickle
from matplotlib.backends.backend_pdf import PdfPages  # Add this line at the top of your script
import matplotlib.pyplot as plt
import sys
sys.path.append('./AdditionalpythonFunctions')
# Importing the function
from descriptive_initials import get_descriptive_initials
from noise_level_addition import add_noise
import datetime

plt.rcParams['figure.max_open_warning'] = 50


# Define the classifiers to test
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier

classifiers = {
    "Random Forest": RandomForestClassifier(max_depth=20, random_state=0)
}

noise_levels = np.linspace(0, 0.9, num=10)  # Change this line as needed.

def process_dataset(filename):
    mat_contents = sio.loadmat(filename)
    numpy_arrayLabels = mat_contents['allLabels']
    numpy_arrayData = mat_contents['matrix_with_data']
    # reduced_size = int(0.1 * numpy_arrayData.shape[0])
    reduced_size = int(0.1 * numpy_arrayData.shape[0])
    random_indices = np.random.choice(numpy_arrayData.shape[0], size=reduced_size, replace=False)
    reduced_data = numpy_arrayData[random_indices]
    reduced_labels = numpy_arrayLabels[random_indices]

    all_accuracies = {name: {'accuracy': {noise: [] for noise in noise_levels},
                             'f1_score': {noise: [] for noise in noise_levels}} for name in classifiers.keys()}

    for noise_level in noise_levels:
        for _ in range(5):  # 10 simulations for each noise level
            noisy_data = add_noise(reduced_data, noise_level=noise_level)
            X_train, X_test, y_train, y_test = train_test_split(noisy_data, reduced_labels, test_size=0.3, shuffle=True)
            for name, clf in classifiers.items():
                clf.fit(X_train, y_train.ravel())
                predicted = clf.predict(X_test)
                print(f"Classification using {name} finished at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                accuracy = metrics.accuracy_score(y_test, predicted)
                f1 = metrics.f1_score(y_test, predicted, average='micro')
                all_accuracies[name]['accuracy'][noise_level].append(accuracy)
                all_accuracies[name]['f1_score'][noise_level].append(f1)

    return all_accuracies


def plot_accuracies(all_accuracies, filename, pdf):
    descriptive_title = get_descriptive_initials(filename)
    fig, ax = plt.subplots(len(classifiers), 1, figsize=(10, len(classifiers) * 5))

    # If there's only one classifier, turn ax into a list for consistency
    if len(classifiers) == 1:
        ax = [ax]

    for idx, (name, scores) in enumerate(all_accuracies.items()):
        accuracies = [[acc for acc in scores['accuracy'][noise_level]] for noise_level in noise_levels]
        ax[idx].boxplot(accuracies, labels=[str(round(noise, 2)) for noise in noise_levels])
        ax[idx].set_title(f'{name} Accuracies at different noise levels')
        ax[idx].set_xlabel('Noise level')
        ax[idx].set_ylabel('Accuracy')

    plt.tight_layout()
    save_directory = "figures_results_deep_learning"
    # Ensure the directory exists
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    # Save the figure to the directory as a PNG
    png_save_path = os.path.join(save_directory, descriptive_title + '.png')
    plt.savefig(png_save_path)

    # Save the figure to the PDF passed to the function
    pdf.savefig(fig)

    plt.show()
    plt.close()

datasets_dir = './DataSets'

# Main code execution
with PdfPages('RForest_Noise_SH.pdf') as pdf:
    matrices = [
        'Data_And_labels_XSadd_final_XPitch_Final_XTrans_final_XHopf_Final_XHopfSaddBig_Final.mat'
    ]
    for matrix in matrices:
        matrix_path = os.path.join(datasets_dir, matrix)
        print(matrix_path)
        all_accuracies = process_dataset(matrix_path)
        plot_accuracies(all_accuracies, matrix, pdf)  # Pass the pdf

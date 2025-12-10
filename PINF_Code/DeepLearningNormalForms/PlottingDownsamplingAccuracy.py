import pickle
import matplotlib.pyplot as plt
pickle_file_path = './pickle_files/accuracies_Downsampling_Comparison_SmallDataset.pkl'
dataset_sizes = [0.0005, 0.001, 0.002, 0.005, .01, .02]
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
plt.show()
import pickle
import matplotlib.pyplot as plt

# Load accuracies
pickle_file_path = './pickle_files/accuracies_Downsampling_Comparison_SmallDataset.pkl'
with open(pickle_file_path, 'rb') as f:
    all_accuracies = pickle.load(f)

# Dataset sizes as fractions
dataset_sizes = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02]

# Convert to percentages
dataset_sizes_percent = [100 * s for s in dataset_sizes]

# Plot
for name, accuracies in all_accuracies.items():
    plt.plot(dataset_sizes_percent, accuracies, marker='o', label=name)

plt.xlabel('Dataset size (%)')
plt.ylabel('Accuracy')
plt.title('Effect of dataset size on classification accuracy')
plt.grid(True)
plt.legend()

plt.savefig('Downsampling_Accuracy_2_plot.pdf', dpi=300)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
noise_levels = np.linspace(0, 0.5, num=6)
pickle_dir = 'pickle_files'
pickle_file_path = os.path.join(pickle_dir, 'accuraciesSlide0002.pickle')

with open(pickle_file_path, 'rb') as handle:
    all_accuracies = pickle.load(handle)

#

# Define your noise levels
noise_levels = np.linspace(0, 0.5, num=6)

plt.figure(figsize=(10, 6))
for i, (name, accuracies) in enumerate(all_accuracies.items(), start=1):
    plt.subplot(1, len(all_accuracies), i)
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
save_path = os.path.join(save_directory, 'accuracy_vs_noise_level.png')

# Save the figure
plt.savefig(save_path)

# Show the plot
plt.show()

# # Define the relative path to the target folder
# relative_figure_path = './FiguresClassificationML/accuracy_vs_noise_level.png'
#
# # Save the figure using the relative path
# plt.savefig(relative_figure_path, dpi = 300)
#
# # Optionally, display the figure
# plt.show()
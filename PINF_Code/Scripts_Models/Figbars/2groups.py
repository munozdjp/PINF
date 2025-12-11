import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# 1. Locate script, data, and prepare output folder
script_path = os.path.abspath(__file__)
script_dir  = os.path.dirname(script_path)
script_name = os.path.splitext(os.path.basename(script_path))[0]

data_dir  = os.path.join(script_dir, '..', 'Data_from_scripts')
plots_dir = os.path.join(script_dir, f"{script_name}_plots")
os.makedirs(plots_dir, exist_ok=True)

# 2. Define your files and labels exactly as uploaded
files_and_labels = [
    ("Saddle_V4_violinMatrix.mat",        "Saddle"),
    ("Transcritical_V4_violinMatrix.mat", "Transcritical"),
    ("Pitchfork_V4_violinMatrix.mat",     "Pitchfork"),
    ("Hopf_V4_violinMatrix.mat",          "Hopf"),
    ("Lorentz_V4_violinMatrix.mat",       "Lorentz"),
    ("FitzHugh_V4_violinMatrix.mat",      "FitzHugh–Nagumo"),
    ("HH_V4_violinMatrix.mat",              "Hodgkin-Huxley"),
    ("vanderpol_V4_violinMatrix.mat",     "Van der Pol"),
    ("Volterra_V4_violinMatrix.mat",      "Lotka Volterra"),
]

# 3. Load all matrices into a dict: label → array (noise × repeats)
data_dict = {}
for fname, label in files_and_labels:
    mat = loadmat(os.path.join(data_dir, fname))
    data_dict[label] = mat['violinmatrix']

# 4. Common settings
noise_levels = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5])

# two groups as you requested
group1 = [ "Pitchfork", "Saddle", "Hopf", "Transcritical"]
group2 = ["Lorentz", "Van der Pol", "Hodgkin-Huxley", "FitzHugh–Nagumo", "Lotka Volterra"]


pastel9 = [
    '#8DD3C7',  # muted teal
    '#FFFFB3',  # pale butter
    '#BEBADA',  # lavender gray
    '#FB8072',  # dusty rose
    '#80B1D3',  # light sky blue
    '#FDB462',  # peach
    '#B3DE69',  # pastel lime
    '#FCCDE5',  # soft pink
    '#D9D9D9'   # light gray
]

# exact same colours you used before
labels = [
    "Pitchfork", "Saddle", "Hopf", "Transcritical",
    "Lorentz", "Van der Pol",
    "Hodgkin-Huxley", "FitzHugh–Nagumo",
    "Lotka Volterra"]

colors_map = { lab: pastel9[i] for i, lab in enumerate(labels) }


def plot_group_box(group, filename, title):
    fig, ax = plt.subplots(figsize=(14,7))
    ax.text(0.0, 1.06, 'a', transform=ax.transAxes,
            fontsize=20, va='top', ha='left', weight='bold')

    box_patches = []
    for i, label in enumerate(group):
        mat = data_dict[label]
        data = [mat[j,:] for j in range(len(noise_levels))]
        positions = [j*2 + i*0.2 for j in range(len(noise_levels))]
        box = ax.boxplot(data, positions=positions,
                         widths=0.15, patch_artist=True,
                        medianprops = dict(color='k', linewidth=1.5))
        for patch in box['boxes']:
            patch.set_facecolor(colors_map[label])
        box_patches.append(box['boxes'][0])

    ax.set_xticks([j*2 + 0.5 for j in range(len(noise_levels))])
    ax.set_xticklabels(noise_levels)
    ax.set_xlabel('Noise', fontsize=16)
    ax.set_ylabel('Error', fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(box_patches, group, loc='upper left', fontsize=16)

    plt.tight_layout()
    out_pdf = os.path.join(plots_dir, f"{script_name}_{filename}.pdf")
    fig.savefig(out_pdf, dpi=300, format="pdf", bbox_inches="tight")
    plt.close(fig)

# 5. Generate and save
plot_group_box(group1, "bifurcation_group", "Bifurcation Group — Error vs Noise")
plot_group_box(group2, "oscillator_group",  "Oscillator Group — Error vs Noise")


import numpy as np
import matplotlib.pyplot as plt

# Simulate ATAC-seq data
# 1000 regions, 80 bins
# Sparse signal: mostly zeros, some peaks
n_regions = 1000
n_bins = 80
matrix = np.zeros((n_regions, n_bins))

# Add some peaks
# Center is bin 40
# Signal intensity around 10-50
for i in range(n_regions):
    if np.random.random() > 0.5: # 50% of regions have signal
        height = np.random.uniform(10, 50)
        width = np.random.randint(2, 10)
        center = 40 + np.random.randint(-5, 5)
        for b in range(max(0, center-width), min(n_bins, center+width)):
            matrix[i, b] = height * np.exp(-0.5 * ((b - center) / (width/2))**2)

# Current logic in the script
all_values = matrix.flatten()
non_zero_values = np.array([v for v in all_values if v > 0])
percentile_90 = np.percentile(non_zero_values, 90)
vmax = percentile_90 * 1.2

print(f"90th percentile of non-zeros: {percentile_90:.4f}")
print(f"Calculated vmax (for heatmaps AND metaprofiles): {vmax:.4f}")

# Calculate mean profile
mean_profile = np.mean(matrix, axis=0)
max_mean = np.max(mean_profile)

print(f"Maximum value in mean profile: {max_mean:.4f}")
print(f"Ratio vmax / max_mean: {vmax/max_mean:.2f}")

if vmax > max_mean * 2:
    print("ISSUE CONFIRMED: vmax is much larger than the mean profile height.")
    print("Metaprofiles will look flat/compressed.")
else:
    print("Issue not reproduced with this simulation.")

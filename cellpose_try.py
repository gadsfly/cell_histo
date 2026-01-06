from skimage import io, exposure
import matplotlib.pyplot as plt
import numpy as np
from cellpose import models
from skimage import measure
import pandas as pd
from tifffile import imwrite


print("="*60)
print("Multi-Channel IF Analysis - Corrected Channel Mapping")
print("="*60)

# Load with CORRECT RGB channels
print("\n1. Loading images with correct channels...")
dapi_rgb = io.imread('PVN, MD, RE image61_R 91_Merged_ch00_SV.tif')
th_rgb = io.imread('PVN, MD, RE image61_R 91_Merged_ch01_SV.tif')
oxy_rgb = io.imread('PVN, MD, RE image61_R 91_Merged_ch02_SV.tif')
cfos_rgb = io.imread('PVN, MD, RE image61_R 91_Merged_ch03_SV.tif')

# Extract correct channels
dapi = dapi_rgb[..., 2]  # Blue channel
th = th_rgb[..., 1]      # Green channel
oxy = oxy_rgb[..., 0]    # Red channel
cfos = cfos_rgb[..., 0]  # Red channel (or use [..., 2], they're identical)

print(f"  DAPI: {dapi.shape}, range {dapi.min()}-{dapi.max()}, mean {dapi.mean():.1f}")
print(f"  TH:   {th.shape}, range {th.min()}-{th.max()}, mean {th.mean():.1f}")
print(f"  Oxy:  {oxy.shape}, range {oxy.min()}-{oxy.max()}, mean {oxy.mean():.1f}")
print(f"  cFos: {cfos.shape}, range {cfos.min()}-{cfos.max()}, mean {cfos.mean():.1f}")

# Quick preview
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for ax, img, title in zip(axes.flat, [dapi, th, oxy, cfos], 
                          ['DAPI', 'TH (Dopamine)', 'Oxytocin', 'cFos']):
    p2, p98 = np.percentile(img, (2, 98))
    ax.imshow(img, cmap='gray', vmin=p2, vmax=p98)
    ax.set_title(f'{title}\nmax={img.max()}, mean={img.mean():.1f}')
    ax.axis('off')
plt.tight_layout()
plt.savefig('preview_all_channels_correct.png', dpi=150)
print("\nSaved preview: preview_all_channels_correct.png")
plt.show()

# 2. Segment nuclei with Cellpose
print("\n2. Segmenting nuclei with Cellpose...")
print("   (This will take 2-5 minutes for large image...)")

model = models.CellposeModel(gpu=True, model_type='nuclei')
masks = model.eval(
    dapi,
    diameter=None,  # Auto-detect
    channels=[0, 0],
    flow_threshold=0.4,
    cellprob_threshold=0
)[0]

n_cells = masks.max()
print(f"   Found {n_cells} cells")

if n_cells == 0:
    print("\n❌ ERROR: No cells found! Check DAPI image quality.")
    import sys
    sys.exit(1)

# Save as 16-bit TIFF (supports up to 65535 unique cells)
if n_cells > 65535:
    print(f"Warning: Too many cells ({masks.max()}), saving as 32-bit")
    imwrite('masks_cellpose.tif', masks.astype(np.uint32))
else:
    imwrite('masks_cellpose.tif', masks.astype(np.uint16))

print(f"Saved masks: {masks.shape}, {masks.max()} cells")


# 3. Measure intensities
print("\n3. Measuring intensities in all channels...")
results = []

for region in measure.regionprops(masks):
    cell_id = region.label
    y, x = region.centroid
    mask = masks == cell_id
    
    results.append({
        'cell_id': cell_id,
        'x': x,
        'y': y,
        'area': region.area,
        'DAPI_mean': np.mean(dapi[mask]),
        'TH_mean': np.mean(th[mask]),
        'Oxy_mean': np.mean(oxy[mask]),
        'cFos_mean': np.mean(cfos[mask]),
        'TH_max': np.max(th[mask]),
        'Oxy_max': np.max(oxy[mask]),
        'cFos_max': np.max(cfos[mask]),
    })

df = pd.DataFrame(results)
print(f"   Measured {len(df)} cells")

# 4. Classify cells
print("\n4. Classifying cells...")
from skimage.filters import threshold_otsu

# Auto-threshold
th_thresh = threshold_otsu(df['TH_mean'].values)
oxy_thresh = threshold_otsu(df['Oxy_mean'].values)
cfos_thresh = threshold_otsu(df['cFos_mean'].values)

print(f"   Thresholds: TH={th_thresh:.1f}, Oxy={oxy_thresh:.1f}, cFos={cfos_thresh:.1f}")

df['TH_positive'] = df['TH_mean'] > th_thresh
df['Oxy_positive'] = df['Oxy_mean'] > oxy_thresh
df['cFos_positive'] = df['cFos_mean'] > cfos_thresh

# Colocalization
df['TH_cFos'] = df['TH_positive'] & df['cFos_positive']
df['Oxy_cFos'] = df['Oxy_positive'] & df['cFos_positive']
df['TH_Oxy'] = df['TH_positive'] & df['Oxy_positive']
df['triple'] = df['TH_positive'] & df['Oxy_positive'] & df['cFos_positive']

# 5. Results
print("\n" + "="*60)
print("RESULTS:")
print("="*60)
print(f"Total cells: {len(df)}")
print(f"\nSingle markers:")
print(f"  TH+ (dopamine):     {df['TH_positive'].sum():5d} ({df['TH_positive'].mean()*100:5.1f}%)")
print(f"  Oxy+ (oxytocin):    {df['Oxy_positive'].sum():5d} ({df['Oxy_positive'].mean()*100:5.1f}%)")
print(f"  cFos+ (active):     {df['cFos_positive'].sum():5d} ({df['cFos_positive'].mean()*100:5.1f}%)")

print(f"\nColocalization:")
if df['TH_positive'].sum() > 0:
    print(f"  TH+/cFos+ (active dopamine):    {df['TH_cFos'].sum():5d} ({df['TH_cFos'].sum()/df['TH_positive'].sum()*100:5.1f}% of TH+)")
else:
    print(f"  TH+/cFos+: N/A (no TH+ cells)")

if df['Oxy_positive'].sum() > 0:
    print(f"  Oxy+/cFos+ (active oxytocin):   {df['Oxy_cFos'].sum():5d} ({df['Oxy_cFos'].sum()/df['Oxy_positive'].sum()*100:5.1f}% of Oxy+)")
else:
    print(f"  Oxy+/cFos+: N/A (no Oxy+ cells)")

print(f"  TH+/Oxy+ (dual labeled):        {df['TH_Oxy'].sum():5d}")
print(f"  TH+/Oxy+/cFos+ (triple):        {df['triple'].sum():5d}")
print("="*60)

# 6. Save results
output_csv = 'pvn_cell_analysis.csv'
df.to_csv(output_csv, index=False)
print(f"\nSaved detailed results to: {output_csv}")

# 7. Visualize
print("\n5. Creating visualization...")
from skimage.transform import resize

# Downsample for display
scale = min(2000 / dapi.shape[0], 2000 / dapi.shape[1])
new_h = int(dapi.shape[0] * scale)
new_w = int(dapi.shape[1] * scale)

dapi_small = resize(dapi, (new_h, new_w), preserve_range=True, anti_aliasing=True)
masks_small = resize(masks.astype(float), (new_h, new_w), order=0, preserve_range=True)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: Segmentation
p2, p98 = np.percentile(dapi_small, (2, 98))
axes[0, 0].imshow(dapi_small, cmap='Blues', vmin=p2, vmax=p98)
axes[0, 0].set_title('DAPI', fontsize=14)
axes[0, 0].axis('off')

from matplotlib.colors import ListedColormap
np.random.seed(42)
n_cells = int(masks_small.max()) + 1
rand_cmap = ListedColormap(np.random.rand(n_cells, 3))
rand_cmap.colors[0] = [0, 0, 0]  # Background black
axes[0, 1].imshow(masks_small, cmap=rand_cmap)
axes[0, 1].set_title(f'Segmentation\n({n_cells} cells)', fontsize=14)
axes[0, 1].axis('off')

# Overlay
overlay = np.zeros((new_h, new_w, 3), dtype=float)
overlay[..., 0] = dapi_small / (dapi_small.max() + 1e-8) * 0.3

for _, row in df.iterrows():
    # Find corresponding location in downsampled image
    y_small = int(row['y'] * scale)
    x_small = int(row['x'] * scale)
    
    if 0 <= y_small < new_h and 0 <= x_small < new_w:
        cell_mask_small = masks_small == row['cell_id']
        if row['TH_positive']:
            overlay[cell_mask_small, 1] = 0.9  # Green
        if row['Oxy_positive']:
            overlay[cell_mask_small, 0] = 0.9  # Red
        if row['cFos_positive']:
            overlay[cell_mask_small, 2] = 0.9  # Blue

axes[0, 2].imshow(overlay)
axes[0, 2].set_title('Classification Overlay\nRed=Oxy, Green=TH, Blue=cFos', fontsize=14)
axes[0, 2].axis('off')

# Row 2: Intensity distributions
axes[1, 0].hist(df['TH_mean'], bins=50, alpha=0.7, color='green', edgecolor='black')
axes[1, 0].axvline(th_thresh, color='red', linestyle='--', linewidth=2)
axes[1, 0].set_xlabel('TH Intensity', fontsize=12)
axes[1, 0].set_ylabel('Number of Cells', fontsize=12)
axes[1, 0].set_title(f'TH Distribution\n{df["TH_positive"].sum()} positive', fontsize=12)
axes[1, 0].grid(alpha=0.3)

axes[1, 1].hist(df['Oxy_mean'], bins=50, alpha=0.7, color='red', edgecolor='black')
axes[1, 1].axvline(oxy_thresh, color='blue', linestyle='--', linewidth=2)
axes[1, 1].set_xlabel('Oxytocin Intensity', fontsize=12)
axes[1, 1].set_ylabel('Number of Cells', fontsize=12)
axes[1, 1].set_title(f'Oxytocin Distribution\n{df["Oxy_positive"].sum()} positive', fontsize=12)
axes[1, 1].grid(alpha=0.3)

axes[1, 2].hist(df['cFos_mean'], bins=50, alpha=0.7, color='purple', edgecolor='black')
axes[1, 2].axvline(cfos_thresh, color='orange', linestyle='--', linewidth=2)
axes[1, 2].set_xlabel('cFos Intensity', fontsize=12)
axes[1, 2].set_ylabel('Number of Cells', fontsize=12)
axes[1, 2].set_title(f'cFos Distribution\n{df["cFos_positive"].sum()} positive', fontsize=12)
axes[1, 2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('pvn_analysis_final.png', dpi=300, bbox_inches='tight')
print("Saved final figure to: pvn_analysis_final.png")
plt.show()

print("\n" + "="*60)
print("✓ ANALYSIS COMPLETE!")
print("="*60)
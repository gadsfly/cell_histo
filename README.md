# Multi-Channel Immunofluorescence Analysis Pipeline

Automated analysis for multi-channel IF images (DAPI, TH, Oxytocin, cFos) with ROI-based regional analysis.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/gadsfly/cell_histo.git
cd cell_histo
```

### 2. Create conda environment

```bash
conda env create -f cellpose.yml
conda activate cellpose
```

## Data Structure

Raw data location: `/data/big_rim/minji/Confocal`

Expected input files per section:
```
section_folder/
├── *_ch00_SV.tif    # DAPI (blue channel)
├── *_ch01_SV.tif    # TH (green channel)
├── *_ch02_SV.tif    # Oxytocin (red channel)
└── *_ch03_SV.tif    # cFos (red channel)
```

---

## Step 1: Run Cellpose Segmentation

Edit `cellpose_try.py` to set the correct paths for your section:

```python
dapi_rgb = io.imread('path/to/your_section_ch00_SV.tif')
th_rgb = io.imread('path/to/your_section_ch01_SV.tif')
oxy_rgb = io.imread('path/to/your_section_ch02_SV.tif')
cfos_rgb = io.imread('path/to/your_section_ch03_SV.tif')
```

Run:

```bash
python cellpose_try.py
```

**Outputs:**
- `preview_all_channels_correct.png` - Channel preview
- `masks_cellpose.tif` - Segmentation masks
- `pvn_cell_analysis.csv` - Cell measurements
- `pvn_analysis_final.png` - Initial analysis figure

---

## Step 2: Create ROI Masks

Create a region annotation image for your section:

1. Open your DAPI image in ImageJ/FIJI (convert to **grayscale** first)
2. Draw regions using distinct colors:
   - **Red** (RGB: 255, 0, 0) → PVN (Paraventricular Nucleus)
   - **Green** (RGB: 0, 255, 0) → MD (Mediodorsal)
   - **Blue** (RGB: 0, 0, 255) → RE (Reuniens)
3. Save as PNG (e.g., `ch00-1.png`)

> **Note:** Image dimensions should match your original DAPI image. Use pure colors without anti-aliasing.

---

## Step 3: Interactive Threshold Tuning

Choose one of two analysis approaches:

### Option A: Pixel Intensity Analysis

Analyzes raw pixel intensities within each ROI region. Faster, no additional segmentation needed.

Open `intensity_threshold.ipynb` and update the paths:

```python
DAPI_PATH = 'path/to/your_section_ch00_SV.tif'
TH_PATH = 'path/to/your_section_ch01_SV.tif'
OXY_PATH = 'path/to/your_section_ch02_SV.tif'
CFOS_PATH = 'path/to/your_section_ch03_SV.tif'
ROI_PATH = 'path/to/ch00-1.png'
```

Run all cells, then display the widget:

```python
display(widgets.HBox([controls, out]))
```

**Controls:**
| Control | Description |
|---------|-------------|
| **TH/Oxy/cFos sliders** | Adjust intensity thresholds |
| **Focus dropdown** | Select marker to highlight in spatial maps |
| **Zoom pad slider** | Adjust zoom level for detailed views |
| **💾 Save** | Export results |

**Saved Files:**
- `pixel_analysis_TIMESTAMP.json` - Thresholds, counts, coexpression stats, intensity statistics
- `pixel_analysis_TIMESTAMP.csv` - Flattened table for Excel/R
- `pixel_analysis_TIMESTAMP.png` - Summary figure

**Presentation Figures:**

After saving, run `generate_presentation_figures()` to create publication-ready figures:
- `fig_expression_by_region.png/svg`
- `fig_activation_rate.png/svg`
- `fig_mean_intensity.png/svg`
- `fig_coexpression_heatmap.png/svg`

---

### Option B: Cell-Based Analysis

Measures marker intensities per segmented cell. Better for cell counting and per-cell colocalization.

Open `cellpose_thres_sec.ipynb` and update paths:

```python
df = pd.read_csv('path/to/pvn_cell_analysis.csv')
masks = io.imread('path/to/masks_cellpose.tif')
dapi_rgb = io.imread('path/to/your_section_ch00_SV.tif')
th_rgb = io.imread('path/to/your_section_ch01_SV.tif')
oxy_rgb = io.imread('path/to/your_section_ch02_SV.tif')
cfos_rgb = io.imread('path/to/your_section_ch03_SV.tif')
roi_img = io.imread('path/to/ch00-1.png')
```

Run all cells, then display the widget:

```python
display(widgets.HBox([controls, out]))
```

**Controls:**
| Control | Description |
|---------|-------------|
| **TH/Oxy/cFos sliders** | Adjust intensity thresholds |
| **use max checkboxes** | Toggle mean vs max intensity |
| **Examples dropdown** | Select marker for example cells |
| **Region dropdown** | Filter example cells by region |
| **💾 Save Results** | Export results |

**Saved Files:**
- `thresholds_with_regions_TIMESTAMP.json` - Threshold settings and summary
- `pvn_analysis_with_regions_TIMESTAMP.csv` - Full cell data with classifications
- `pvn_analysis_with_regions_TIMESTAMP.png` - Summary figure

**Output CSV Columns:**
| Column | Description |
|--------|-------------|
| `cell_id` | Unique cell identifier |
| `x`, `y` | Cell centroid coordinates |
| `area` | Cell area in pixels |
| `DAPI_mean`, `TH_mean`, `Oxy_mean`, `cFos_mean` | Mean intensity per channel |
| `TH_max`, `Oxy_max`, `cFos_max` | Max intensity per channel |
| `region` | Assigned region (PVN, MD, RE, or outside) |
| `TH_positive`, `Oxy_positive`, `cFos_positive` | Boolean classification |
| `TH_cFos`, `Oxy_cFos` | Colocalization flags |

---

## File Descriptions

| File | Description |
|------|-------------|
| `cellpose_try.py` | Initial Cellpose segmentation script |
| `intensity_threshold.ipynb` | Pixel intensity-based threshold tuning |
| `cellpose_thres_sec.ipynb` | Cell-based threshold tuning |
| `cellpose.yml` | Conda environment specification |

---

## Citation

If using Cellpose:

> Stringer, C., Wang, T., Michaelos, M., & Pachitariu, M. (2021). Cellpose: a generalist algorithm for cellular segmentation. Nature Methods, 18(1), 100-106.

## License

MIT

## Contact

Mir Qi (mir.qi@duke.edu)
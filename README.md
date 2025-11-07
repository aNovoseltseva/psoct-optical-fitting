# PS-OCT Image Processing Pipeline

Automated processing pipeline for polarization-sensitive optical coherence tomography (PS-OCT) data with optical parameter fitting and birefringence quantification for brain tissue imaging.

## Overview

This pipeline processes PS-OCT data through three main stages:
1. **OCT Reconstruction** - Converts raw interferometric data to reconstructed volumes
2. **Calibration Fitting** - Extracts optical parameters (focal depth zf, Rayleigh range zr, scattering coefficients) from agar calibration tiles
3. **Tissue Fitting** - Applies calibration parameters to fit tissue optical properties (scattering μs, birefringence μb) and stitches tiles into full mosaic

## Requirements

### Software
- **MATLAB R2022a or later**
- **FIJI/ImageJ** - Required for BaSiC shading correction and stitching
  - Download from: https://fiji.sc/
  - Required plugin: BaSiC (for illumination correction)
- **xvfb-run** (Linux only) - For headless FIJI operation

### MATLAB Toolboxes
- Image Processing Toolbox (required)
- Signal Processing Toolbox (required for `findchangepts`)
- Statistics and Machine Learning Toolbox (optional, for advanced fitting)

### Custom Functions
All custom MATLAB functions are included in the `functions/` directory:
- Core optical fitting algorithms
- Image stitching and blending
- Data I/O utilities
- Surface detection algorithms

## Installation

1. **Install FIJI/ImageJ:**
   ```bash
   # Download from https://fiji.sc/
   # Or on Linux:
   wget https://downloads.imagej.net/fiji/latest/fiji-linux64.zip
   unzip fiji-linux64.zip
   ```

2. **Install BaSiC plugin in FIJI:**
   - Open FIJI
   - Go to Help > Update > Manage Update Sites
   - Check "BaSiC" in the list
   - Click "Close" and then "Apply Changes"
   - Restart FIJI

3. **Clone this repository:**
   ```bash
   git clone https://github.com/[your-username]/psoct-processing-pipeline.git
   cd psoct-processing-pipeline
   ```

4. **Configure paths:**

## Configuration

Edit `step0_info_d1.txt` with your sample parameters:

```
folder: Path to PS-OCT data directory
P2path: Path to two-photon reference data (if using for stitching)
nslice: Total number of tissue sections
stitch: 0 = use 2P coordinates, 1 = use PS-OCT coordinates
ds_factor: Downsampling factor (4 = 12×12 μm pixels, 16 = 48×48 μm pixels)
numX: Number of tiles in X direction
numY: Number of tiles in Y direction
Tile size: 1 = 1×1 mm tiles, 3 = 3×3 mm tiles
```

## Usage

### Step 1: OCT Reconstruction

Reconstruct raw interferometric data into co- and cross-polarized intensity volumes.

```bash
qsub step1_OCT_recon.sh
```

**Output:** `dist_corrected/co-[slice]-[tile]-*.dat` and `cross-[slice]-[tile]-*.dat` files

### Step 2: Agar Calibration Fitting

Extract focal depth (zf), Rayleigh range (zr), and scattering coefficients from agar tiles.

**Manual step:** In `step2_fitting_scatter_tile.m`, line 66, specify which tiles contain agar calibration:
```matlab
for iFile = [1:690]  % Update this range to match your agar tile numbers
```

Run fitting:
```bash
matlab -nodisplay -r "step2_fitting_scatter_tile; exit"
```

**Output:** `fitting/zf_zr_sim_fit.mat` containing calibration parameters

### Step 3: Tissue Fitting

Apply calibration parameters to fit tissue optical properties across all tiles.

Submit as array job (one job per slice):
```bash
qsub -t 1-[nslice] step4_fit.sh
```

Or run single slice:
```matlab
id = '1';  % slice number
step3_fitting_tissue
```

**Output:** 
- `fitting/vol[slice]/` - Individual tile fitting results
- Stitched parameter maps: `mus_[slice].mat`, `mub_[slice].mat`, `bfg_[slice].mat`

## Output Files

- **mus**: Scattering coefficient map (μs, mm⁻¹)
- **mub**: Backscattering coefficient map (μb, mm⁻¹)  
- **bfg**: Birefringence  map
- **zf_zr_sim_fit.mat**: System calibration parameters

## Notes

- Tile indices in file names: `co-[slice]-[tile]-[Z]-[X]-[Y].dat`
- Default fitting depth: 130 μm for agar, 80 μm for tissue (adjustable in function calls)
- Stitching uses bidirectional mosaic pattern with 15% overlap (1mm tiles) or 5% overlap (3mm tiles)



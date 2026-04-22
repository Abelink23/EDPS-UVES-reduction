# UVES Pipeline Data Management and Analysis

This repository contains Python tools for organizing, processing, and analyzing spectroscopic data from the VLT/UVES instrument, specifically designed to handle outputs from the ESO EDPS pipeline.

## Modules

### 1. `handle_data.py`
Focuses on file system organization and fixing FITS headers for compatibility.
- **`organize_EDPSdata`**: Matches pipeline outputs with final products using `ARCFILE` and `ESO PRO CATG` keywords.
- **`organize_fits_in_EDPSfolder_by_object`**: Scans directories and sorts FITS files into subfolders named after the `OBJECT` header keyword, while renaming files based on pipeline and archive identifiers.
- **`fake_master_response` / `fake_date_flat`**: Utility functions to manually correct or "spoof" header metadata (like dates or category tags) to ensure downstream pipeline compatibility.

### 2. `UVES.py`
A class-based approach to spectral reduction and visualization.
- **`ob` Class**: The core object for handling UVES spectra.
  - **Extraction**: Supports standard pipeline outputs, master reduced files, and ASCII formats.
  - **Corrections**:
    - Barycentric velocity correction.
    - Air-to-Vacuum conversion (Morton, Griesen, or IAU methods).
    - Automatic edge cutting based on instrument setup (Dichroic/Central Wavelength).
  - **Processing**:
    - Spike removal using Z-score or Kernel-based sigma clipping.
    - Flux normalization with iterative polynomial fitting and strong-line masking.
    - Negative flux cleaning and median "taring".
- **Visualization**:
  - `plot`: Standard spectrum plotting with support for flux, calibrated flux, or normalized flux.
  - `plot_orders`: Visualizes 2D resampled echelle orders.
  - `plt_diag_lines`: Multi-panel diagnostic plots for specific interstellar/stellar lines.
- **Master File Creation**:
  - `master_files`: Consolidates individual pipeline products (science, error, flux-calibrated) into a single multi-extension FITS "Master" file and a corresponding ASCII table.

## Requirements
- `numpy`
- `scipy`
- `matplotlib`
- `astropy`

## Usage Example
```python
from UVES import ob

# Load and process a spectrum
spec = ob("red_science_blue.fits", bar_corr=True, to_vac='Morton', cut_edges=True)

# Clean spikes and normalize
spec.clean_spikes(method='zscore')
spec.normalize()

# Plot the result
spec.plot(y='norm_flux')
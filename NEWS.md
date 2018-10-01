# News

## 1.0.2 - October 1, 2018
### Added
- Plots now search slot `@nmf.c1` as well.
- Additional `plotTree` plotting options (can plot +/- values, configure limits of color axis, turn off legend).
- Additional options in `edgesFromDM` for displaying diffusion map transitions in `plotDim`.
- Can fit onset/offset times against absolute value of curve in `impulseFit`.
- Added `treeForceStretchCoords` for additional refinement of force-directed layouts.

### Changed
- **Changed behavior**: Fixed bug where tip cells could end up being assigned to a segment other than the one they define.
- Fixed `treeForceDirectedLayout` to work without named segments.
- Fixed importing Seurat object when tSNE or PCA not calculated.
- Fixed bug in `floodPseudotimeCalc` when only one cell is used as the root. (Note: This is still not recommended usage, as in real data it is unlikely that a single cell can be definitely chosen as the starting point.)
- Reduced memory usage of `floodBuildTM` on large data
- Updated installation instructions with FAQ about udunits on Linux.
- Additional cleanup of documentation.

## 1.0.1 - May 20, 2018
### Added
- Added to repository: Supplementary Analysis

### Changed
- Fixed `plotTreeForce` with discrete labels.
- Fixed outdated references to deprecated `@group` slot.
# News

## 1.1.0 - July 25, 2019
### Added
- Added functions for using expression of NMF modules to identify cells that
are potentially doublets: `NMFDoubletsDefineModules`, `NMFDoubletsPlotModuleThresholds`, `NMFDoubletsDetermineCells`, `NMFDoubletsPlotModulesInCell`, and `NMFDoubletsPlotModuleCombos`. The idea is that many NMF modules encode cell type programs; if you identify modules that are not expressed in overlapping gradients (which might represent legitimate developmental transitions), and are strongly expressed in distinct sets of cells, then small numbers of cells that express both modules are likely to represent technical or biological doublets where two distinct cells (which each have their own cell type program expressed) have been detected as a single cell. Removing these can improve trajectory inference, as they often create 'short-circuits' between distinct portions of the developmental process.
- Added new functions for calculating smoothed expression and plotting expression curves: `geneSmoothFit`, `plotSmoothFit`, `cropSmoothFit`, `plotSmoothFitMultiCascade`. Provides alternative, less parametric options (LOESS smoothing and spline curves) to the previous impulse model
- Added `whichCells` to help identify cells that match particular criteria for use in plotting, subsetting, and differential expression testing.
- Added `plotDimDiscretized` and `plotTreeDiscretized` to allow plotting expression or other metadata in an on/off fashion (based on user-defined thresholds) to improve some visualizations.

### Changed
- Changed behavior of `pseudotimeWeightTransitionMatrix` to allow processing of more than 46,503 cells by processing the matrix in pieces.
- Fixed (another) bug in `treeForceDirectedLayout` that results from cells with duplicate random walk parameters
- Fixed an issue where `createURD` would fail because it would not find `method::new` for some reason.
- Structure of impulse fits has been modified to allow plotting with the original functions, but also the new `plotSmoothFit` and `plotSmoothFitMultiCascade` functions.
- Added additional input checking to `plotTree`, `plotTreeDual` to help identify bad inputs.
- Fixed bug in `plotTreeForce` when plotting TRUE/FALSE metadata
- Modified `treeForceDirectedLayout` to check for graph connectivity prior to calling some FDL routines that require connected graphs
- Fixed bug with segment names in `treeForceDirectedLayout`.

## 1.0.3 - May 28, 2019
### Added
- Added `plotTreeDual` function which can plot two labels as red/green on the tree dendrogram layout
### Changed
- Fixed bug in `processRandomWalks` when a cell occurs in in walks that is not in `@diff.data`
- Fixed bug in `plotTree` when using `color.limits` without explicitly setting `symmetric.scale`.
- Fixed bug in `plotTreeForce` when plotting discrete values.
- Fixed bugs in `treeForceDirectedLayout` that result from cells with duplicate random walk parameters or cells assigned NA by force-directed layout.

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
#' URD - Reconstruction of Branching Developmental Trajectories
#' 
#' URD reconstructs transcriptional trajectories underlying specification or differentiation 
#' processes in the form of a branching tree from single-cell RNAseq data. The URD
#' object class is documented under \code{\link{URDclass}}. There are many functions
#' in the package, but a starting guide to help you get started is:
#' 
#' @section Creating an URD object:
#' \itemize{
#'   \item \code{\link{createURD}} - Create an URD object
#'   \item \code{\link{seuratToURD}} - Import an object from Seurat
#' }
#' @section Basic data analysis:
#' \itemize{
#'   \item \code{\link{findVariableGenes}} - Find highly variable genes
#'   \item \code{\link{calcKNN}} - Calculate a k-nearest neighbor graph
#'   \item \code{\link{knnOutliers}} - Find outliers based on k-NN distance
#'   \item \code{\link{calcPCA}} - Principal Components Analysis
#'   \item \code{\link{calcTsne}} - tSNE project of data
#'   \item \code{\link{graphClustering}} - Graph-based clustering of cells
#' }
#' @section Calculate diffusion map:
#' \itemize{
#'   \item \code{\link{calcDM}} - Calculate a diffusion map
#' }
#' @section Calculate pseudotime:
#' \itemize{
#'   \item \code{\link{floodPseudotime}} - Perform 'flood' breadth-first graph search
#'   \item \code{\link{floodPseudotimeProcess}} - Convert 'floods' into cell pseudotime
#'   \item \code{\link{pseudotimePlotStabilityOverall}} - Verify that enough 'flood' simulations were performed
#' }
#' @section Find developmental trajectories:
#' \itemize{
#'   \item \code{\link{pseudotimeDetermineLogistic}} - Determine parameters of logistic for biasing transition probabilities
#'   \item \code{\link{pseudotimeWeightTransitionMatrix}} - Bias transition probabilities according to pseudotime
#'   \item \code{\link{simulateRandomWalksFromTips}} - Perform biased random walks starting at each 'tip'
#'   \item \code{\link{processRandomWalksFromTips}} - Convert biased random walks into visitation frequencies
#' }
#' @section Build developmental tree:
#' \itemize{
#'   \item \code{\link{loadTipCells}} - Load cells belonging to each tip into the URD object
#'   \item \code{\link{buildTree}} - Find branching tree structure from trajectory data
#'   \item \code{\link{treeForceDirectedLayout}} - Generate a force-directed layout of the tree
#' }
#' @section Plot with abandon:
#' \itemize{
#'   \item \code{\link{plotDim}} - Plot data on dimensionality reduction (tSNE, PCA, or diffusion map)
#'   \item \code{\link{plotDimDual}} - Plot multiple data on dimensionality reduction
#'   \item \code{\link{plotDimArray}} - Plot the same data on many pairs of dimensions
#'   \item \code{\link{plotDim3D}} - Plot dimensionality reduction in 3D
#'   \item \code{\link{plotDim3DStoreView}} - Store a plotDim3D orientation for making many similar plots
#'   \item \code{\link{plotTree}} - Plot data on dendrogram representation of tree
#'   \item \code{\link{plotTreeForce}} - Plot data on 3D force-directed layout representation of tree
#'   \item \code{\link{plotTreeForceStore3DView}} - Store a plotTreeForce orientation for making many similar plots
#' }
#' 
#' @docType package
#' @aliases URD-package
#' @name URD
NULL


# scClock
## Repository for B. Divergence scRNA seq

1. `util_funcs.R` contains the functions that are needed for the analysis.
1. `ReadScRNAseqData_BD.R`  read in the 0h 37 degree WT data and generate the Seurat object. Note that the data is down-sampled to include 800 points per cluster. This is done for memory management. Markers of each cluster are identified for later use.
1. `ReadScRNAseqData_TG.R`  read in scRNAseq from RH strain from ToxoAtlas (384 well). No need for down-sampling (~600 cells). This data comes with cell cycle phase info (DNA content assay). DEGs are performed based on cell cycle phase.
1. `ReadScRNAseqData_PB.R`  read in scRNAseq from Blood Stage of Plasmodilum Berghei from Malaria Atlas. Data is is down-sampled to include 800 points per cluster. This data comes with life-stage info. DEGs are performed based on Life Stage.
1. `DEAsyncTimeCourse.R` reads in the synchronized time-course data and performs a time-course DE analysis to include genes that show a significant change over time. Note that technical replicates are pooled, so there are 2 replicates per time-point. Outlier samples are excluded.
1. `fitPseudoTime.R` fits a principal curve to the first two PCA components of the sc data and diagonally projects the cells along the curve to order the cells and create a pseudo time. The pseudo time is scaled from 0 to 12 hours, binned every 20 min, and samples within each bin are considered as replicates. Genes are fitted against the pseudo time ina GAM model and genes that show no correlation with time are excluded.
1. `alignScWithSync.R` Fits smoothing splines to the gene curves and performs a cross-correlation analysis between the common genes in sc and sync data. The optimal shift is decided based on the distribution of the lag time between sc and sync common genes. The sc pseudo time is adjusted accordingly to re-map the start point.
1. `mixed_effect.R` Fits a mixed effect model to both sync and sc data to identify mean trends over time.
1. `clusterCurves.R` uses the mean-trend of Marker genes and performs a time-series clustering with dynamic time warping (dtw) to cluster the curves based on their similarity of shape, taking minor stretch and shifts into consideration. The curves are used to produce heatmaps.

### Glasso

To be filed
The Seurat objects for all three species are down-sampled to contain 800 cells/cluster and 2000 most variable genes. The code `genGlassoDat.R` generates the Seurat objects, uses network smoothing to reduce sparsity (params: alpha = 0.3, inter = 0.3) and generates the empirical contrivance matrices for Tg, Pb, and Bd. The contrivance matrices are used to generate Sparse covariance, precision, and partial correlation matrices (code `glasso_slurm_array.sh` on Chimera, parameter: rho = 1:999/1000). To see the effect of network smoothing, see exploratory code `estiamteCov.R`

Note: The slurm array has an upper limit of 1000 jobs, indexed by 0:1000. To scale lambda with step-size 0.0001 do:

```
#SBATCH --array=0-1000:1 
L_SCALE=$(echo "scale=4; x + ($SLURM_ARRAY_TASK_ID)/10000" | bc)
```
where x needs to change to `0, 0.1, 0.2`, etc. Keeping the upper limit of `x` to `0.5` seems sufficient. When `x = 0`, starting `x` from `10` maybe sufficient.

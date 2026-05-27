library(Seurat) # v4.3
library(harmony) # v1.2
library(qs)
`%nin%` = Negate(`%in%`)
`%!in%` = Negate(`%in%`)
set.seed(42)

# command line arguments passing
# args[1] - a list of file paths of Seurat objects (under the same time point) to be harmonized
# args[2] - the batch variable to correct for.
# args[3] - the output file name of the harmonized object to be saved
args <- commandArgs(trailingOnly=TRUE)
inSeurs <- as.vector(sapply(args[1], function(x){strsplit(x,",")[[1]]}))
group.by.vars <- args[2]
outFile <- args[3]

# load object into a list
myList <- list()
for(i in 1:length(inSeurs)){
  if (endsWith(inSeurs[i], ".RDS")) {
    myList[[i]] <- readRDS(inSeurs[i])
  } else {
    myList[[i]] <- qread(inSeurs[i])
  }
}

for (i in 1:length(myList)) {
  if (DefaultAssay(myList[[i]]) != "SCT") {
    stop(paste0("Object ", i, " needs to be preprocessed by SCTransform()."))
  }
  
  if (group.by.vars %nin% colnames(myList[[i]]@meta.data)) {
    stop(paste0("Object ", i, " is missing the required Harmony metadata column: ", group.by.vars, ". Please parse the information first."))
  }
  
  slot(myList[[i]]$SCT@SCTModel.list[[1]], 'median_umi') = median(myList[[i]]$SCT@SCTModel.list[[1]]@cell.attributes$umi)
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = myList)

# merge
merged_seurat <- merge(x = myList[[1]],
                       y = myList[2:length(myList)])
VariableFeatures(merged_seurat) <- features

# Run PCA()
merged_seurat <- merged_seurat %>%
  RunPCA()

# Run Harmony
harmonized_seurat <- RunHarmony(merged_seurat,
                                group.by.vars = group.by.vars,
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony",
                                project.dim = F)

harmonized_seurat <- harmonized_seurat %>%
  RunUMAP(reduction = "harmony", assay = "SCT", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims=1:30) %>%
  FindClusters(resolution=.8)

harmonized_seurat <- PrepSCTFindMarkers(object = harmonized_seurat, verbose = TRUE)
qsave(harmonized_seurat, outFile)



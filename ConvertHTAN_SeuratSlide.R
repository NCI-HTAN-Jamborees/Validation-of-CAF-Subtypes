library(Seurat) #v5.1.0
library(spacexr) #v2.1.1
ConvertSeuratToSpatial = function (spatialSeuratHTAN){
  spatialObject=CreateSeuratObject(counts=spatialSeuratHTAN@assays$RNA@counts,assay="Spatial")
  coord.df =data.frame(x= -Embeddings(spatialSeuratHTAN@reductions$spatial)[,2], y=Embeddings(spatialSeuratHTAN@reductions$spatial)[,1])
  spatialObject@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
  )
  spatialObject@meta.data= cbind(spatialObject@meta.data,spatialSeuratHTAN@meta.data)
  return(spatialObject)
}

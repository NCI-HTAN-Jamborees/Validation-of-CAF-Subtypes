library(Seurat) #v5.1.0
library(ggplot2) #v3.5.1
library(scCustomize) #v2.1.2


slidesList = list.files("../project-files/03.option3_CC/slide_seq/",pattern="RDS",full.names = T)
for (file in slidesList){
  sample = gsub("-","_",gsub("\\.RDS","",basename(file)))
  assign(sample,readRDS(file))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pdf("CAF_Subtype_Slides.pdf",width=12,height=6)
for (objName in ls(pattern="^HTAPP")){
  obj = eval(parse(text=objName))
  DimPlot(obj,group.by="cell_type",reduction="spatial",pt.size=0.01)+
    ggtitle(paste0(objName,":\nInitial Cell Types"))+
    coord_fixed(xlim = range(obj@reductions$spatial@cell.embeddings[,1]), ylim=range(obj@reductions$spatial@cell.embeddings[,2]))
  
  p1 = DimPlot(obj,group.by="celltype",reduction="spatial",pt.size=0.01)+
    ggtitle(paste0(objName,":\nCell Types: CAF Subtype"))+
    coord_fixed(xlim = range(obj@reductions$spatial@cell.embeddings[,1]), ylim=range(obj@reductions$spatial@cell.embeddings[,2]))
  
  fibroblast_highlight=list()
  for (cafSub in sort(unique(obj$celltype[which(obj$cell_type=="fibroblast")]))){
    fibroblast_highlight[[cafSub]] = colnames(obj)[which(obj$celltype==cafSub)]
  }
  
  p2=Cell_Highlight_Plot(seurat_object = obj, cells_highlight = fibroblast_highlight,reduction="spatial",highlight_color = gg_color_hue(length(fibroblast_highlight)))+
    ggtitle(paste0(objName,":\nCell Types: CAF Subtype"))+
    coord_fixed(xlim = range(obj@reductions$spatial@cell.embeddings[,1]), ylim=range(obj@reductions$spatial@cell.embeddings[,2]))
  
    
  subObj = subset(obj,cells=colnames(obj)[which(obj$cell_type=="fibroblast")])
  #subObj = FindNeighbors(subObj)
  #subObj = RunUMAP(subObj)
  p3 = DimPlot(subObj,reduction="spatial",group.by="celltype")+
    ggtitle(paste0(objName,":\nFibroblast CAF Subtypes"))+
    coord_fixed(xlim = range(obj@reductions$spatial@cell.embeddings[,1]), ylim=range(obj@reductions$spatial@cell.embeddings[,2]))
 print(p1|p2|p3)
  
}

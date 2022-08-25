
> folders=list.files('./',pattern='[12]$')
> folders
# [1] "G1"  "G2"  "L1"  "L2"  "NP1" "NP2" "PI1" "PI2"

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                               project = folder )
})
sce.big <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                       sceList[[6]],sceList[[7]],sceList[[8]]), 
                 add.cell.ids = folders, 
                 project = "mouse8")

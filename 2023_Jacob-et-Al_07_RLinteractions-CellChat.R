################################################################################
#####START SESSION##############################################################
################################################################################

setwd('C:/Users/User/Jacob-Et-Al/CellChat')
getwd()

library(Seurat)
library(dplyr)
library(scales)
library(RColorBrewer)
library(igraph)
library(parallel)
library(cowplot)
library(ggplot2)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)

data_path <- 'C:/Users/User/Jacob-Et-Al/CellChat'

################################################################################
#####RUN CELLCHAT on E12.5 EPI-FIB##############################################
################################################################################

###PART I: PREPARE THE REQUIRED DATA OBJECTS####################################

#####Load sc data###############################################################
all_cells_3_tp <- readRDS(file.path(data_path, 
                          'combined_filtered_normalized_seurat_object.rds'))

#####Read in metadata###########################################################
clustering <- read.csv(file.path(data_path, 'All_MetaData.csv'))

###Subset metadata info to filtered cells only##################################
all_cells_3_tp <- subset(all_cells_3_tp, cells = clustering$X)
rownames(clustering) <- clustering$X
all_cells_seurat <- CreateSeuratObject(counts=all_cells_3_tp@assays$RNA@counts,
                                       meta.data = clustering)

###Combine with second level clustering info####################################
sub_clustering <- read.csv('All-Cells_2nd-level-clustering.csv')
all_cells_seurat <- subset(all_cells_seurat, cells = sub_clustering$Cell_Barcode)

###Set cell identities as active ident##########################################
all_cells_seurat@active.ident <- factor(NULL)
cell_cluster_assignment <- sub_clustering$Inferred_Cell_Type
all_cells_seurat <- SetIdent(all_cells_seurat, 
                             value = as.character(cell_cluster_assignment))

###Create new seurat object prior to subsetting for re-use in later analyses####
all_cells_seurat2 <- all_cells_seurat

###Subset for E12.5#############################################################
all_cells_seurat <- subset(all_cells_seurat,
                           cells = 
                             which(all_cells_seurat@meta.data$embryonic_age == 'E12.5'))

###Subset for EPI and FIB clusters of interest##################################
all_cells_seurat <- subset(all_cells_seurat, 
                           cells = which(all_cells_seurat@active.ident == "EPI Periderm" | 
                                        all_cells_seurat@active.ident == "EPI Basal1" |
                                       all_cells_seurat@active.ident == "EPI BasalTagln" |
                                       all_cells_seurat@active.ident == "FIB Origin1" |
                                       all_cells_seurat@active.ident == "FIB Origin2" ))

###Create data and metadata object to make CellChat object######################
data.input = all_cells_seurat@assays$RNA@counts #normalized data matrix
meta = all_cells_seurat@meta.data #dataframe with meta data
cell.use = rownames(meta)
meta$new_clustering <- as.factor(all_cells_seurat@active.ident)
meta$labels <- meta$new_clustering
meta = data.frame(labels = meta$labels, row.names = colnames(data.input)) 
unique(meta$labels)

###Create CellChat Object#######################################################
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents)) #number of cells in each cell group

###Lig-Rec Database#############################################################
CellChatDB <- CellChatDB.mouse #use CellChatDB.mouse for mouse signalling
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

###Rename idents alphabetically#################################################
#Set desired cell order
cell.levels <- sort(levels(cellchat@meta$labels))

#Store cell labels in scRNA object metadata
all_cells_seurat$cellgroup <- Idents(all_cells_seurat)

#Add cell labels to cellchat metadata
identity = data.frame(group = all_cells_seurat$cellgroup, 
                      row.names = names(all_cells_seurat$cellgroup))
unique(identity$group)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell_type")
levels(cellchat@idents) #check idents are correct

#Reorder cells
cellchat <- setIdent(cellchat, ident.use = "cell_type", levels = cell.levels)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

###Preprocessing the expression data for cell-cell communication analysis#######

#Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

###PART II: INFERENCE OF CELL-CELL COMMUNICATION NETWORK########################

###Compute the communication probability and cellular communication network#####
cellchat <- computeCommunProb(cellchat)

###Filter out the cell-cell communication if there are only few cells###########
###in certain cell groups###
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat, slot.name = "netP")

###Infer the cell-cell communication at a signaling pathway level###############
cellchat <- computeCommunProbPathway(cellchat)
cellchat@net
cellchat@netP

###Calculate the aggregated cell-cell communication network#####################
cellchat <- aggregateNet(cellchat)

###Set colors###################################################################
colors_cell_groups <- c('#ADEBEB', '#29A3A3', '#1A75FF', '#666565', '#294D49')

###PART III: IDENTIFY GLOBAL COMMUNICAITON PATTERNS#############################

###Plot nr of interactions as circle############################################
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, color.use = colors_cell_groups, 
                 title.name = "Number of interactions")

###pathway circle plots SEMA3 and CXCL##########################################
for (pathway in c('SEMA3', 'CXCL')){
  pdf(paste0('EPI_FIB_EARLY_DIV_',pathway, '_circle_plot.pdf'))
  par(mfrow=c(1,1))
  print(netVisual_aggregate(cellchat, signaling = pathway, 
                            color.use = colors_cell_groups, layout = "circle"))
  dev.off()}

###Identify and visualize outgoing communication pattern of secreting cells#####
selectK(cellchat, pattern = "outgoing")

nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns)

#River plot
netAnalysis_river(cellchat, pattern = "outgoing", 
                  color.use = colors_cell_groups)

#Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing",  color.use = colors_cell_groups, 
                dot.size = c(1, 4), pathway.show = sort(cellchat@netP$pathways))

#####Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = nPatterns)

netAnalysis_river(cellchat, pattern = "incoming", 
                  color.use = colors_cell_groups)

netAnalysis_dot(cellchat,  pattern = "incoming", 
                color.use = colors_cell_groups, dot.size = c(1, 4),
                pathway.show = sort(cellchat@netP$pathways))





################################################################################
#####RUN CELLCHAT on E13.5 and 14.5  FIB INTER, FIB MUSCLE, MUSCLE##############
################################################################################

###PART I: PREPARE THE REQUIRED DATA OBJECTS####################################

###Use previously stored Seurat object and subset for  E13.5 and 14.5 ##########
all_cells_seurat <- subset(all_cells_seurat2, cells = 
                             which(all_cells_seurat2@meta.data$embryonic_age == 'E14.5'
                                   | all_cells_seurat2@meta.data$embryonic_age == 'E13.5'))

###Subset for FIB and MUSCLE clusters of interest###############################
all_cells_seurat <- subset(all_cells_seurat, 
                           cells = which(all_cells_seurat@active.ident == 'FIB Muscle1' | 
                                         all_cells_seurat@active.ident == 'FIB Muscle2'| 
                                         all_cells_seurat@active.ident == 'MUSCLE Early' |
                                         all_cells_seurat@active.ident == 'MUSCLE Mid' |
                                         all_cells_seurat@active.ident == 'FIB Inter1' |
                                         all_cells_seurat@active.ident == 'FIB Inter2' |
                                         all_cells_seurat@active.ident == 'FIB Inter3' |
                                         all_cells_seurat@active.ident == 'MUSCLE Late' ))
levels(all_cells_seurat)

###Create data and metadata object to make CellChat object######################
data.input = all_cells_seurat@assays$RNA@counts #normalized data matrix
meta = all_cells_seurat@meta.data #a dataframe with meta data
cell.use = rownames(meta)
meta$new_clustering <- as.factor(all_cells_seurat@active.ident)
meta$labels <- meta$new_clustering
meta = data.frame(labels = meta$labels, row.names = colnames(data.input)) 
unique(meta$labels)

###Create CellChat Object#######################################################
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))#number of cells in each cell group

###Lig-Rec Database#############################################################
cellchat@DB <- CellChatDB.use

###Rename idents alphabetically#################################################
#Set desired cell group order
cell.levels <- c("FIB Inter1","FIB Inter2","FIB Inter3", "FIB Muscle1", 
                 "FIB Muscle2", "MUSCLE Early", "MUSCLE Mid", "MUSCLE Late")

#Store cell labels in scRNA object metadata
all_cells_seurat$cellgroup <- Idents(all_cells_seurat)

#Add cell labels to cellchat metadata
identity = data.frame(group = all_cells_seurat$cellgroup, 
                      row.names = names(all_cells_seurat$cellgroup))
unique(identity$group)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell_type")
levels(cellchat@idents) #check idents are correct

#Reorder cells
cellchat <- setIdent(cellchat, ident.use = "cell_type", levels = cell.levels)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

###Preprocessing the expression data for cell-cell communication analysis#######

#Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

###PART II: INFERENCE OF CELL-CELL COMMUNICATION NETWORK########################

###Compute the communication probability and cellular communication network#####
cellchat <- computeCommunProb(cellchat)

###Filter out the cell-cell communication if there are only few cells###########
###in certain cell groups###
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat, slot.name = "netP")
levels(cellchat@idents)

###Infer the cell-cell communication at a signaling pathway level###############
cellchat <- computeCommunProbPathway(cellchat)

###Calculate the aggregated cell-cell communication network#####################
cellchat <- aggregateNet(cellchat)

###Set colors###################################################################
colors_cell_groups <- c('#32B6EA', '#1978E0', '#253EF2', '#02441C',
                        '#149914', '#29A3A3', '#47D1D1',  '#ADEBEB')

###Plot nr of interactions as circle############################################
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, color.use = colors_cell_groups,
                 title.name = "Number of interactions")

for (pathway in c('EDN', 'GRN')){
  pdf(paste0('MUSCLE_FIB_INTER_MUSCLE_',pathway, '_circle_plot.pdf'))
  par(mfrow=c(1,1))
  print(netVisual_aggregate(cellchat, signaling = pathway, 
                            color.use = colors_cell_groups, layout = "circle"))
  dev.off()
}

###PART III: IDENTIFY GLOBAL COMMUNICAITON PATTERNS#############################

###Identify and visualize outgoing communication pattern of secreting cells#####
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns)

#Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing", color.use = colors_cell_groups, 
                pathway.show= sort(c( 'SEMA3', 'ANGPTL', 'MIF',  'PSAP', 'TGFb', 
                                      'VISFATIN', 'CXCL', 'TWEAK','GRN', 'NRG',
                                      'NT', 'EDN','HGF','PTN', 'MK',  'WNT', 
                                      'PDGF', 'ncWNT', 'BMP', 'FGF', 'PERIOSTIN',
                                      'GALECTIN', 'GAS', 'NGF',  
                                      'EGF')), dot.size = c(1,4))

###Identify and visualize incoming communication pattern of target cells########
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

#Dot plot
netAnalysis_dot(cellchat, pattern = "incoming", pathway.show= 
                  sort(c( 'SEMA3', 'ANGPTL', 'MIF',  'PSAP', 'TGFb', 'VISFATIN',
                          'TWEAK','GRN', 'NRG','NT','CXCL', 'EDN', 'HGF','PTN',
                          'MK',  'WNT', 'PDGF', 'ncWNT', 'BMP', 'FGF', 'PERIOSTIN',
                          'GALECTIN', 'GAS', 'NGF',
                          'EGF')),
                color.use = colors_cell_groups, dot.size = c(1,4))





################################################################################
#####RUN CELLCHAT on E13.5 and 14.5  NEURO, VES, IMMU interactions##############
################################################################################

###PART I: PREPARE THE REQUIRED DATA OBJECTS####################################

###Use previously stored Seurat object and subset for  E13.5 and 14.5 ##########
all_cells_seurat <- subset(all_cells_seurat2, cells =
                             which(all_cells_seurat2@meta.data$embryonic_age == 'E13.5' |
                                  all_cells_seurat2@meta.data$embryonic_age == 'E14.5'))

###Subset for late populations of interest######################################
all_cells_seurat <- subset(all_cells_seurat, cells = which(all_cells_seurat@active.ident != "EPI BasalTagln" &
                                                             all_cells_seurat@active.ident != "EPI Basal1" &
                                                             all_cells_seurat@active.ident != "EPI Periderm" &
                                                             all_cells_seurat@active.ident != "FIB Origin1" &
                                                             all_cells_seurat@active.ident != "FIB Origin2" &
                                                             all_cells_seurat@active.ident != "FIB Origin3" &
                                                             all_cells_seurat@active.ident != "FIB Origin4" &
                                                             all_cells_seurat@active.ident != "FIB Origin5" &
                                                             all_cells_seurat@active.ident != "CHOND" &
                                                             all_cells_seurat@active.ident != "FIB Origin6" &
                                                             all_cells_seurat@active.ident != "FIB Deep1" &
                                                             all_cells_seurat@active.ident != "FIB Deep2" &
                                                             all_cells_seurat@active.ident != "FIB Deep3" 
))
levels(all_cells_seurat@active.ident)

###Create data and metadata object to make CellChat object######################
data.input = all_cells_seurat@assays$RNA@counts #normalized data matrix
meta = all_cells_seurat@meta.data #a dataframe with meta data
cell.use = rownames(meta)
meta$new_clustering <- as.factor(all_cells_seurat@active.ident)
meta$labels <- meta$new_clustering
meta = data.frame(labels = meta$labels, row.names = colnames(data.input)) 
unique(meta$labels)

###Create CellChat Object#######################################################
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents)) #number of cells in each cell group

###Lig-Rec Database#############################################################
cellchat@DB <- CellChatDB.use
###Rename idents alphabetically#################################################

#Set desired cell order
cell.levels <- c("EPI Basal2","EPI Basal3","EPI Basal4","EPI EarlyDiff","EPI LateDiff",
                 "EPI EarlyPlacode","EPI LatePlacode","FIB EarlyDC","FIB LateDC",
                 "FIB Upper1","FIB Upper2","FIB Upper3","FIB Upper4","FIB Lower",
                 "FIB Muscle1","FIB Muscle2","FIB Inter1","FIB Inter2","FIB Inter3",
                 "IMMU DendriticCells","IMMU Macrophages","IMMU MastCells","MUSCLE Early",
                 "MUSCLE Mid","MUSCLE Late","NC Melanocytes","NC SchwannCells",
                 "VESSEL BECs","VESSEL LECs","VESSEL MuralCells")

#Store cell labels in scRNA object metadata
all_cells_seurat$cellgroup <- Idents(all_cells_seurat)

#Add cell labels to cellchat metadata
identity = data.frame(group = all_cells_seurat$cellgroup, 
                      row.names = names(all_cells_seurat$cellgroup))
unique(identity$group)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell_type")
levels(cellchat@idents) #check idents are correct

#Reorder cells
cellchat <- setIdent(cellchat, ident.use = "cell_type", levels = cell.levels)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

###Preprocessing the expression data for cell-cell communication analysis#######

#Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

###PART II: INFERENCE OF CELL-CELL COMMUNICATION NETWORK########################

###Compute the communication probability and cellular communication network#####
cellchat <- computeCommunProb(cellchat)

###Filter out the cell-cell communication if there are only few cells###########
###in certain cell groups###
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat, slot.name = "netP")
cellchat@data.signaling
levels(cellchat@idents)

###Infer the cell-cell communication at a signaling pathway level###############
cellchat <- computeCommunProbPathway(cellchat)
cellchat@net
cellchat@netP

###Calculate the aggregated cell-cell communication network#####################
cellchat <- aggregateNet(cellchat)

###Set colors###################################################################
colors_cell_groups <- c('#248F24', '#D8D041', '#47D147', '#0047B3', '#B3D1Ff',
                        '#FF8c1A',  '#FF6666', '#A68CCE', '#37156B',
                        '#FABAA2',  '#F86B51', '#C12489','#CE73E0', '#FCC25A',
                        '#02441C','#149914',  '#32B6EA', '#1978E0','#253EF2',
                        '#B3D1Ff','#1A75FF', '#0047B3', '#29A3A3', '#47D1D1',
                        '#ADEBEB', '#FF6666', '#FF8c1A', '#DB70B8' , '#B82E8A',
                        '#6D3E91')

###PART III: IDENTIFY GLOBAL COMMUNICAITON PATTERNS#############################

###Identify and visualize outgoing communication pattern of secreting cells#####
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns)

###Dot plot (Vessel, NC, Immu together)#########################################
netAnalysis_dot(cellchat, pathway.show =
          c( 'GALECTIN', 'MK', 'PTN', 'ANGPTL','APELIN','APJ','CALCR','EGF','FGF',
            'GAS','HH','IGF','LIFR', 'ncWNT', 'NPR1','PDGF','PERIOSTIN','PROS',
            'PTH','SEMA3','TGFb','TNF','TWEAK','VEGF','VISFATIN','WNT','BMP','NGF',
            'ENHO','GRN','NRG','EDN','KIT','CCL','CHEMERIN','CSF','CX3C','FLT3',
            'IL2','IL4','MIF'), 
                pattern = "outgoing",  color.use = colors_cell_groups, 
          dot.size = c(1, 4))

#####Identify and visualize incoming communication pattern of target cells#####
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = nPatterns)
netAnalysis_dot(cellchat, 
              group.show = c('VESSEL BECs', 'VESSEL LECs', 'VESSEL MuralCells',
                             'NC Melanocytes', 'NC SchwannCells', 
                             'IMMU DendriticCells', 'IMMU Macrophages', 
                             'IMMU MastCells'), 
                color.use = c('#DB70B8', '#B82E8A', '#6D3E91', '#FF6666',
                              '#FF8c1A', '#B3D1Ff', '#1A75FF', '#0047B3'),
                pathway.show = c('GALECTIN', 'MK', 'PTN', 'ANGPTL','APELIN','APJ',
                                 'CALCR','EGF','FGF','GAS','HH','IGF','LIFR', 
                                 'ncWNT', 'NPR1','PDGF','PERIOSTIN','PROS','PTH',
                                 'SEMA3','TGFb','TNF','TWEAK','VEGF','VISFATIN',
                                 'WNT','BMP','NGF','ENHO','GRN','NRG','EDN','KIT',
                                 'CCL','CHEMERIN','CSF','CX3C','FLT3','IL2','IL4','MIF'), 
                pattern = 'incoming', 
                dot.size = c(1, 4))



sink(file = 'log_session_info_cellchat_plotting.txt')
sessionInfo()


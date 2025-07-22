#==============================================================================#
#            SCRIPT PARA ANÁLISE INTEGRATIVA DE DADOS DE scRNA-seq             #
#                                 com Seurat                                   #
#==============================================================================#

# --- 1. CARREGAR PACOTES ---
# Certifique-se de que todos os pacotes necessários estão instalados.
# install.packages(c("Seurat", "tidyverse", "patchwork", "harmony", "SeuratWrappers"))

library(Seurat)
library(tidyverse) # Para manipulação de dados (ex: %>% e read_tsv)
library(patchwork) # Para combinar plots
library(harmony)   # Para o método de integração Harmony
library(SeuratWrappers) # Para integração com o Seurat v5


#==============================================================================#
# --- 2. CARREGAMENTO E CONTROLE DE QUALIDADE (QC) AUTOMATIZADO POR AMOSTRA ---
#==============================================================================#

# Esta seção automatiza o carregamento de múltiplas amostras (formato 10x Genomics)
# e a aplicação de filtros de qualidade (QC) para cada uma delas.

# --- 2.1. Definir o diretório de trabalho ---
# Coloque aqui o caminho para a pasta que contém seus arquivos de amostra.
# Ex: setwd("/home/user/meu_projeto/dados_brutos")
setwd("/caminho/para/sua/pasta/de/dados")

# Criar subpastas para salvar os resultados do QC e os objetos RDS
dir.create("QC", showWarnings = FALSE)
dir.create("RDS", showWarnings = FALSE)


# --- 2.2. Função auxiliar para identificar outliers ---
# Esta função calcula outliers com base no Desvio Absoluto Mediano (MAD),
# que é mais robusto a valores extremos do que o desvio padrão.
mad_outlier <- function(seurat_object, metric, nmads) {
  M <- seurat_object@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  
  # Identifica células que estão 'nmads' desvios abaixo ou acima da mediana
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}

# --- 2.3. Função principal de QC ---
# Esta função itera sobre cada amostra, realiza o QC e salva os resultados.

# NOTA SOBRE GENES MITOCONDRIAIS:
# O Seurat busca por padrão o prefixo "mt-" (em minúsculo) para camundongo.
# Se seus dados usam "MT-" (humano) ou não têm prefixo, ajuste o `pattern` abaixo.
# Para camundongos sem prefixo, pode-se usar uma lista explícita de genes, como:
# pattern_mt <- "^Nd1|^Nd2|^Co1|^Co2|^Atp8|^Atp6|^Co3|^Nd3|^Nd4l|^Nd4|^Nd5|^Nd6|^Cytb"
# Para humano, o padrão é: pattern = "^MT-"

processar_amostras_qc <- function(samples_prefix, mt_pattern = "^mt-", perc_mt_threshold = 10) {
  
  for (sample_id in samples_prefix) {
    message(paste("Processando amostra:", sample_id))
    
    # Carregar dados do formato 10x (matrix, barcodes, features)
    counts <- ReadMtx(
      mtx = paste0(sample_id, "_matrix.mtx.gz"),
      cells = paste0(sample_id, "_barcodes.tsv.gz"),
      features = paste0(sample_id, "_features.tsv.gz")
    )
    
    # Criar o objeto Seurat
    # min.cells: Inclui features (genes) expressos em pelo menos 3 células.
    # min.features: Inclui células que expressam pelo menos 300 features.
    X <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 300)
    
    # Salvar o número de células antes do filtro
    cells_before <- ncol(X)
    
    # Calcular a porcentagem de genes mitocondriais
    X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = mt_pattern)
    
    # Filtros baseados em outliers
    # Remove células com contagens de RNA ou número de genes muito discrepantes (5 MADs da mediana)
    outlier_filter <- mad_outlier(X, 'nCount_RNA', 5) | mad_outlier(X, 'nFeature_RNA', 5)
    X <- subset(X, cells = colnames(X)[!outlier_filter])
    
    # Filtro de porcentagem mitocondrial
    # NOTA: Para single-nucleus (snRNA-seq), o limiar é bem mais baixo (ex: 1% ou 2%).
    # Altere o `perc_mt_threshold` se necessário.
    X <- subset(X, subset = percent.mt < perc_mt_threshold)
    
    cells_after <- ncol(X)
    
    # Salvar um resumo do QC
    qc_summary <- data.frame(
      Sample_ID = sample_id,
      Cells_Before_QC = cells_before,
      Cells_After_QC = cells_after
    )
    write.csv(qc_summary, file.path("QC", paste0(sample_id, "_QC_Summary.csv")), row.names = FALSE)
    
    # Salvar o objeto Seurat filtrado
    saveRDS(X, file.path("RDS", paste0(sample_id, "_AfterQC.rds")))
    
    message(paste("Amostra", sample_id, "finalizada. Células restantes:", cells_after))
  }
}

# --- 2.4. Executar o pipeline de QC ---
# Lista todos os arquivos da pasta que terminam com "_barcodes.tsv.gz"
# e extrai o prefixo de cada amostra.
file_list <- list.files(pattern = "_barcodes.tsv.gz$")
sample_prefixes <- sub("_barcodes.tsv.gz$", "", file_list)

# Executa a função de processamento
# Altere o padrão mitocondrial ou o limiar aqui, se necessário
processar_amostras_qc(sample_prefixes, mt_pattern = "^mt-", perc_mt_threshold = 10)


#==============================================================================#
# --- 3. VISUALIZAÇÃO DO QC (OPCIONAL, MAS RECOMENDADO) ---
#==============================================================================#

# Carregar os objetos RDS salvos após o QC para visualização
rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
seurat_list_qc <- lapply(rds_files, readRDS)
names(seurat_list_qc) <- sub("_AfterQC.rds", "", basename(rds_files))

# Função para criar Violin Plots padronizados
VlnPlot_Custom <- function(seurat_object, feature, ...) {
  VlnPlot(seurat_object, feature, pt.size = 0, ...) + # pt.size = 0 remove os pontos
    scale_y_continuous(limits = c(0, 20)) + # Padroniza o eixo Y
    theme(legend.position = "none", axis.title.x = element_blank()) +
    labs(title = seurat_object@project.name)
}

# Criar e combinar plots de `percent.mt` para todas as amostras
plot_list_mt <- lapply(seurat_list_qc, function(obj) VlnPlot_Custom(obj, "percent.mt"))
combined_plots <- wrap_plots(plotlist = plot_list_mt, ncol = 3) # Ajuste o ncol conforme o número de amostras

# Mostrar o plot combinado
print(combined_plots)
# ggsave("QC_ViolinPlots_PercentMT.png", combined_plots, width = 12, height = 8)


#==============================================================================#
# --- 4. INTEGRAÇÃO DE DADOS (SEURAT v5) ---
#==============================================================================#

# Esta é a abordagem moderna e recomendada no Seurat v5, usando `IntegrateLayers`.

# --- 4.1. Mergir objetos e preparar para integração ---
# Se os objetos não estiverem na memória, carregue-os
# rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
# seurat_list_qc <- lapply(rds_files, readRDS)

# Mergir todos os objetos da lista em um único objeto Seurat
merged_seurat <- merge(seurat_list_qc[[1]], y = seurat_list_qc[-1])

# Dividir o assay 'RNA' por amostra ('orig.ident'). Essencial para o `IntegrateLayers`.
merged_seurat[['RNA']] <- split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)

# --- 4.2. Pré-processamento padrão ANTES da integração ---
# Este pipeline é executado em cada "layer" (amostra) individualmente.
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)


# --- 4.3. Executar a Integração ---
# Escolha UM dos métodos abaixo. RPCA costuma ser mais rápido e robusto.

# Aumentar o limite de memória global para o processamento, se necessário
options(future.globals.maxSize = 8000 * 1024^2) # 8 GB

# MÉTODO A: RPCA (Reciprocal PCA) - Rápido e robusto
integrated_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

# MÉTODO B: CCA (Canonical Correlation Analysis)
# integrated_seurat <- IntegrateLayers(
#   object = merged_seurat,
#   method = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.cca",
#   verbose = TRUE
# )

# MÉTODO C: Harmony
# integrated_seurat <- IntegrateLayers(
#   object = merged_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = "pca",
#   new.reduction = "harmony",
#   verbose = TRUE
# )

# Re-juntar as layers do assay RNA após a integração para análises futuras
integrated_seurat[["RNA"]] <- JoinLayers(integrated_seurat[["RNA"]])

# Salvar o objeto integrado
# saveRDS(integrated_seurat, "integrated_seurat_rpca.rds")


#==============================================================================#
# --- 5. ANÁLISE DOWNSTREAM (CLUSTERIZAÇÃO E VISUALIZAÇÃO) ---
#==============================================================================#
# integrated_seurat <- readRDS("integrated_seurat_rpca.rds")

# --- 5.1. Análise de Componentes Principais (PCA) e Clusterização ---
# Use o "reduction" gerado pela sua integração (ex: "integrated.rpca")

# O ElbowPlot ajuda a determinar o número de dimensões (componentes principais) a usar
ElbowPlot(integrated_seurat, ndims = 30, reduction = "integrated.rpca")

# Defina o número de dimensões com base no "cotovelo" do gráfico
# Exemplo: dims_to_use <- 1:20
dims_to_use <- 1:20 
reduction_name <- "integrated.rpca" # Mude para "integrated.cca" ou "harmony" se usou outro método

# Encontrar vizinhos e clusters
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = reduction_name, dims = dims_to_use)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5) # Ajuste a resolução para mais/menos clusters

# --- 5.2. Visualização com UMAP ---
integrated_seurat <- RunUMAP(integrated_seurat, dims = dims_to_use, reduction = reduction_name)

# Plotar UMAP por cluster
p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")

# Plotar UMAP por amostra de origem (para verificar a qualidade da integração)
p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Amostra de Origem")

# Exibir plots lado a lado
p1 + p2


#==============================================================================#
# --- 6. IDENTIFICAÇÃO DE GENES MARCADORES ---
#==============================================================================#

# Mudar o Assay ativo para 'RNA' para encontrar marcadores no nível de expressão gênica
DefaultAssay(integrated_seurat) <- "RNA"

# Encontrar marcadores para cada cluster
# min.pct: detectado em pelo menos 25% das células do cluster
# logfc.threshold: log-fold change mínimo de 0.25
all_markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Salvar a lista de marcadores em um arquivo CSV
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

# Selecionar os 10 melhores marcadores de cada cluster com base no Log2 Fold Change
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# --- 6.2. Visualização com Heatmap ---
# NOTA: O ScaleData na etapa de integração foi executado apenas nos 'VariableFeatures'.
# Para um heatmap preciso com genes que podem não ser variáveis, é bom re-escalar os dados
# apenas com os genes de interesse.
integrated_seurat <- ScaleData(integrated_seurat, features = top10_markers$gene)

# Gerar o heatmap
DoHeatmap(
  subset(integrated_seurat, downsample = 300), # Downsample para performance
  features = top10_markers$gene,
  raster = TRUE # Usa rasterização para gerar um arquivo menor e mais rápido
)


#==============================================================================#
# --- 7. APÊNDICE: MÉTODO DE INTEGRAÇÃO LEGADO (SEURAT < v5) ---
#==============================================================================#

# Este método (FindIntegrationAnchors/IntegrateData) ainda funciona, mas o
# `IntegrateLayers` é o recomendado atualmente.

# # Carregar a lista de objetos QC'd
# rds_files <- list.files("RDS", pattern = "\\.rds$", full.names = TRUE)
# seurat_list_legacy <- lapply(rds_files, readRDS)
# 
# # Normalizar e encontrar features variáveis para cada objeto na lista
# for (i in 1:length(seurat_list_legacy)) {
#   seurat_list_legacy[[i]] <- NormalizeData(seurat_list_legacy[[i]], verbose = FALSE)
#   seurat_list_legacy[[i]] <- FindVariableFeatures(seurat_list_legacy[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# }
# 
# # Encontrar âncoras de integração
# integration_anchors <- FindIntegrationAnchors(object.list = seurat_list_legacy, dims = 1:30)
# 
# # Integrar os dados
# integrated_legacy <- IntegrateData(anchorset = integration_anchors, dims = 1:30)
# 
# # Definir o assay padrão para o integrado
# DefaultAssay(integrated_legacy) <- "integrated"
# 
# # Continuar com ScaleData, RunPCA, etc.
# integrated_legacy <- ScaleData(integrated_legacy, verbose = FALSE)
# integrated_legacy <- RunPCA(integrated_legacy, npcs = 30, verbose = FALSE)
# integrated_legacy <- RunUMAP(integrated_legacy, reduction = "pca", dims = 1:20)
# integrated_legacy <- FindNeighbors(integrated_legacy, reduction = "pca", dims = 1:20)
# integrated_legacy <- FindClusters(integrated_legacy, resolution = 0.5)
# 
# # Visualizar
# DimPlot(integrated_legacy, reduction = "umap")
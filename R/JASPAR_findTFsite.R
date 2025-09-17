# 安装所需的包
# install.packages(c("TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", "JASPAR2014"))
# BiocManager::install("Biostrings", force = TRUE)
# 
# BiocManager::install(c("TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", "JASPAR2014"))



# BiocManager::install("JASPAR2024", force = TRUE)

# 加载所需的包
library(biomaRt)
library(TFBSTools)
library(JASPAR2024)
library(BSgenome.Hsapiens.UCSC.hg38)

 # 安装所需的包
# install.packages(c("TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", "JASPAR2014"))
# BiocManager::install("Biostrings", force = TRUE)
# 
# BiocManager::install(c("TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", "JASPAR2014"))



# BiocManager::install("JASPAR2024", force = TRUE)

# 加载所需的包
library(biomaRt)
library(TFBSTools)
library(JASPAR2024)
library(BSgenome.Hsapiens.UCSC.hg38)



findTFsite <- function(genes = genes, species.code = species.code, maRtdataset = maRtdataset){
  
  library(biomaRt)
  library(TFBSTools)
  library(JASPAR2024)
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  # 指定Ensembl数据库和感兴趣的基因
  ensembl <- useMart("ensembl", dataset = maRtdataset)
  
  # genes = c("GLUD1") #输入要预测的基因，可以一次性输入多个基因，基因越多，运行时间越久。
  
  # 把基因名转换为ensembl_id
  gene_id<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                 filters = "external_gene_name", # 指定转化基因的格式
                 values = genes, # 基因列表
                 mart = ensembl)
  gene_id
  gene_id = gene_id$ensembl_gene_id
  
  # 获取基因组信息
  gene_info <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), 
                     filters = "ensembl_gene_id", 
                     values = gene_id, 
                     mart = ensembl)
  
  # 获取基因组数据
  genome <- BSgenome.Hsapiens.UCSC.hg38
  
  
  # 定义提取启动子序列的函数
  get_promoter_sequence <- function(chromosome, start, end, strand, upstream_distance = 2000) {
    if (strand == 1) {
      promoter_start <- max(start - upstream_distance, 1)
      promoter_end <- start - 1
    } else {
      promoter_start <- end + 1
      promoter_end <- min(end + upstream_distance, seqlengths(genome)[[chromosome]])
    }
    promoter_seq <- getSeq(genome, paste0("chr", chromosome), promoter_start, promoter_end)
    return(promoter_seq)
  }
  
  
  # 使用JASPAR数据库
  # db <- file.path(system.file("extdata", package="JASPAR2024"),
  #                 "JASPAR2024.sqlite")
  opts <- list()
  opts[["species"]] <- species.code
  # opts[["type"]] <- "ChIP-seq"
  opts[["all_versions"]] <- FALSE
  # pfm <- getMatrixSet(db, opts)
  
  library(RSQLite)
  JASPAR2024 <- JASPAR2024()
  JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
  
 
  
  pfm <- getMatrixSet(JASPARConnect, opts)
   head(pfm)
  # 可以查看getMatrixSet函数，查看这些opts()中各种参数的意义，根据不同目的可以使用不同的参数。
  
  
  # 打开一个csv文件进行写入
  file_conn <- file(paste0(genes,"_",species.code,"TFBS2024.csv"), "w")
  # 写入CSV文件头部
  writeLines("Gene,TF,Predicted_TFBS", file_conn)
  
  # 循环处理每个PWM对象
  for (name in names(pfm)) {
    # 提取基因的上游序列
    for (i in 1:dim(gene_info)[1]) {
      promoter_seq <- get_promoter_sequence(gene_info[i, "chromosome_name"], 
                                            gene_info[i, "start_position"], 
                                            gene_info[i, "end_position"], 
                                            gene_info[i, "strand"])
      
      # 预测转录因子结合位点
      predicted_TFBS <- matchPWM(pfm[[name]]@profileMatrix, promoter_seq)
      
      # 获取基因名
      gene_name <- getBM(attributes = "external_gene_name", 
                         filters = "ensembl_gene_id",
                         values = gene_info[i, "ensembl_gene_id"], 
                         mart = ensembl)
      
      # 如果预测结果不为空且基因名存在
      if (!is.null(predicted_TFBS) && nrow(gene_name) > 0) {
        gene_name <- gene_name$external_gene_name[1]
        
        # 将预测结果写入CSV文件
        writeLines(paste(gene_name, pfm[[name]]@name, predicted_TFBS, sep = ","), file_conn)
      }
    }
  }
  
  # 关闭文件连接
  close(file_conn)
}



#9606 hsapiens_gene_ensembl  human
findTFsite(genes = "CXCL9", species.code = 9606, maRtdataset = "hsapiens_gene_ensembl")


#10090 mmusculus_gene_ensembl mouse
findTFsite(genes = "CXCL9", species.code = 10090, maRtdataset = "mmusculus_gene_ensembl") 











setwd("D:\\CAUCLOUD\\Peacock\\02.test\\26.RNA")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
###读入数据
database <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
###设置分组
Feature <- factor(c("control","control","control","control","KD","KD","KD","KD"))
###分组与性状联系
coldata <- data.frame(row.names = colnames(database), Feature)
#countData <- countData[, rownames(colData)]
#database<-as.matrix(database)
###DEseq2读入数据
dds <- DESeqDataSetFromMatrix(countData = database , colData = coldata, design = ~Feature)
dds <- DESeq(dds)
res <- results(dds)
resordered <- res[order(res$padj),]
summary(res)
###写入基因表达结果
write.csv(as.data.frame(resordered), file="results.csv")
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)
###写入上调和下调基因
write.csv(as.data.frame(res1_up),file="upgene.csv")
write.csv(as.data.frame(res1_down),file="downgene.csv")
#3.3 绘制火山图
genes<- res1
# 根据上调、下调、不变为基因添加颜色信息
genes$color <- ifelse(genes$padj<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")
###要标注得基因
geneList0 <- c('MSTRG.10952|EDNRB2')
geneList <- subset(genes, rownames(genes) %in% geneList0)
p<- ggplot(# 数据、映射、颜色
  genes, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #突出表示差异基因
  geom_point(data=geneList,aes(x = log2FoldChange, y = -log10(padj)),colour="yellow",size=3.5)+
  #辅助线
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="Blue/White,log2(Fold change)",y="-log10 (p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))
p
###要标注得基因重命名
#geneList$genename <- c('EDNRB2')
geneList$genename <- c('CYP2J19')
library(ggrepel)
p + geom_label_repel(data = geneList, 
                     aes(x = log2FoldChange, y = -log10(padj), label = genename),
                     size = 4,color="black",
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",   #连线的颜色
                     segment.size = 0.4,  #连线的粗细
)
p


#3.1展示某个基因的表达量
#可以直接用DESeq2的plotCounts
#plotCounts(dds, gene = 'MSTRG.10952|EDNRB2', intgroup="Feature") # 指定某个基因
plotCounts(dds, gene = 'CYP2J19|CYP2J19', intgroup="Feature") 
  
#可用ggplot对图片样式进行修改，并用ggrepel进行标注
d <- data.frame(t(subset(database,rownames(database)=="MSTRG.10952|EDNRB2")))

d$Feature <- coldata$Feature

 ggplot(d, aes(x = factor(d$Feature) , y = d$MSTRG.10952.EDNRB2, color = "Feature"))+
   geom_point(position=position_jitter(w=0.2,h=0))+
   geom_text_repel(aes(label=rownames(d)))+
   theme_bw()+
   ggtitle("EDNRB2")+
  theme(plot.title=element_text(hjust=0.5))

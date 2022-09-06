library(foreach)
library(pheatmap)
library(metacell)
library(tgconfig)
library(tgstat)
#
set_param("mc_plot_device",'pdf', "metacell")
set_param("scm_spike_regexp","^ERCC-","metacell")
set_param("scm_mc_mark_k_per_clust",100,"metacell") #default: 5
set_param("scm_mc_mark_min_gene_cov",0.3,"metacell") # default: 0.25
set_param("scm_mc_mark_min_gene_fold",2,"metacell") # default: 1.5

set_param("mcell_mc2d_K",30,"metacell") # default: 20
set_param("mcell_mc2d_T_edge",0.02,"metacell") # default: 0.05
set_param("mcell_mc2d_max_confu_deg",4,"metacell") # default: 5
set_param("mcell_mc2d_edge_asym",FALSE,"metacell") # default: TRUE
set_param("mcell_mc2d_proj_blur",0.02,"metacell") # default: 0.02
grDevices::pdf.options(useDingbats = FALSE)




#Patient 1
setwd("./Project/LUAD")
if(!dir.exists("LUADdb")) dir.create("LUADdb/")
scdb_init("LUADdb/", force_reinit=T)
if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

mcell_import_multi_mars(mat_nm="LUAD",base_dir="./umi.tab/" ,dataset_table_fn="./umi.tab/Lung_pat_info_20210703.txt",force=TRUE,patch_cell_name=T)
mcell_plot_umis_per_cell("LUAD", min_umis_cutoff = 400)

nms = c(rownames(mat@mat), rownames(mat@ignore_gmat)) 
nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
pre_nr_term <- c("^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("MBALAT1", "XIST", "XIST_intron")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
pre_bad_genes
##
genes_RBC <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
cells_RBC <- names(which((apply(mat@mat[intersect(genes_RBC,rownames(mat@mat)),],2,sum))>=1))
genes_Mt<-foreach(i='^MT-', .combine = c) %do% grep(i, nms, v=T)
cells_mt= names(which((Matrix::colSums(mat@mat[genes_Mt,])/Matrix::colSums(mat@mat))>0.2))

mcell_mat_ignore_genes(new_mat_id='LUAD', mat_id='LUAD', pre_bad_genes, reverse=F)
ignore<-as.character(read.table("./ignore_cell.txt",sep='\t',header=T)$x)
cell=c(names(which(Matrix::colSums(mat@mat)>8000)),names(which(Matrix::colSums(mat@mat)<400)))
unique(c(cells_RBC,cell,cells_mt,ignore))->rmcells
mcell_mat_ignore_cells(new_mat_id='LUAD', mat_id='LUAD', ig_cells = rmcells, reverse = F)


genes_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','STMN1','TOP2A','MKI67','MX1','RSAD2')
tab_fn = "./lateral_gmods.txt"
mcell_mat_rpt_cor_anchors(mat_id='LUAD', gene_anchors = genes_anchors, cor_thresh = 0.1,
                          gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table('./lateral_gmods.txt', header=T)
foc_genes = apply(gcor_mat[, intersect(colnames(gcor_mat),genes_anchors)], 1, which.max)
#
mat = scdb_mat("LUAD")
mcell_add_gene_stat(gstat_id="LUAD", mat_id="LUAD", force=T)
mcell_gset_filter_varmean(gset_id="LUAD_feats", gstat_id="LUAD", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "LUAD_feats", gstat_id="LUAD", T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id="LUAD", gset_id="LUAD_feats")
gset <- scdb_gset("LUAD_feats")
pst_genes <- names(gset@gene_set)
pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                 "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                 "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
pst_ex_genes <- c()
pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
pst_add_genes <- c()
final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
final_genes

gset = gset_new_gset(sets = gset@gene_set[final_genes], desc = "final genes")
scdb_add_gset("LUAD_feats", gset)


mcell_add_cgraph_from_mat_bknn(mat_id="LUAD",
                gset_id = "LUAD_feats",
                graph_id="LUAD_graph",
                K=252,
                dsamp=T)###一般K为总细胞数的平方根的数目，若细胞数很少，K=20-40
mcell_add_cgraph_from_mat_bknn(mat_id="LUAD",
                gset_id = "LUAD_feats",
                graph_id="LUAD_graph",
                K=270,
                dsamp=T)


mcell_coclust_from_graph_resamp(
                coc_id="LUAD_coc500",
                graph_id="LUAD_graph",
                min_mc_size=70,
                p_resamp=0.75, n_resamp=500)




mcell_mc_from_coclust_balanced(
                coc_id="LUAD_coc500",
                mat_id= "LUAD",
                mc_id= "LUAD_mc",
                K=70, min_mc_size=70, alpha=2) 

mcell_mc_from_coclust_balanced(
                coc_id="LUAD_coc500",
                mat_id= "LUAD",
                mc_id= "LUAD_mc",
                K=90, min_mc_size=100, alpha=2) 




#
mc = scdb_mc("LUAD_mc")
mc@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc@mc_fp))
scdb_add_mc("LUAD_mc",mc)
mcell_gset_from_mc_markers(gset_id="LUAD_markers", mc_id="LUAD_mc")
mcell_mc_plot_marks(mc_id="LUAD_mc", gset_id="LUAD_markers", mat_id="LUAD")
##Projecting metacells and cells in 2D
mcell_mc2d_force_knn(mc2d_id="LUAD_2dproj",mc_id="LUAD_mc", graph_id="LUAD_graph")
#tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
#tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="LUAD_2dproj")#

mc_hc <-mcell_mc_hclust_confu(mc_id="LUAD_mc",graph_id="LUAD_graph")
mc_sup <- mcell_mc_hierarchy(mc_id="LUAD_mc",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="LUAD_mc",
                        graph_id="LUAD_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2,
                         show_mc_ids=TRUE)
#
lfp <- round(log2(mc@mc_fp),2)
write.table(lfp,file='./lfp.txt',sep='\t',quote=F)
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="TCL1A")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="CD79A")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="CD8A")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="CD27")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="FCRL4")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="TNFRSF14")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="RGS1")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="IGHG1")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="IGHG1")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="CD3D")
mcell_mc2d_plot_gene(mc2d_id="LUAD_2dproj",gene="CXCL13")

lfp <- log2(mc@mc_fp)
plt = function(gene1, gene2, lfp, colors)
{
    plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
    text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))

}
plt(gene1 = 'CD3D', gene2 = 'FOXP3', lfp = lfp, colors = mc@colors)



mc@annots=c(
  rep('GNLY+ cNK',17),rep('GNLY- cNK',6),rep('GNLY+ cNK',31),'CD8+ CTL','GNLY+ cNK',rep("CD8+ CTL",25),'resident-like NK',rep("CD8+ CTL",2),rep('resident-like NK',4),
  rep("B cell",20),'mDC','cDC1',rep('cDC2',5),rep("IL1B+ Macrophage",2),rep('APOC1+ Macrophage',4),rep("CD16+ Monocyte",3),rep("CD14+ Monocyte",5),'ICOS+ ILC2','KIT+ ILC3','Mast Cell','Basophil',
  'pDC','Plasma B',rep("CXCL13-CD8+ Dysfunctional Trm",7),'CD8+ISG15+ T',rep("CXCL13+CD8+ Dysfunctional Trm",2),rep("CD8+ Tem",3),"CXCL13-CD8+ Dysfunctional Trm",rep("CD8+ Tem",2),
  "CXCL13+CD8+ Dysfunctional Trm",rep("CXCL13-CD8+ Dysfunctional Trm",2),"CD8+ CTL",'CD8+ Tem','Naive-like','CD8+ MAIT','CD8+ MAIT',rep('CD4+ Tem',9),'CD4+ Treg','Naive-like',rep('CD4+ Treg',4),rep("Naive-like",7),
  rep("CD4+ Tem",12),'CD8+ MAIT',"CD4+ Tem",'CD8+ MAIT','Naive-like',"CD4+ Tem","CD4+ Tem","CD4+ Tem",'CD8+ MAIT','CD4+ Tem','CD4+ Treg','CD4+ Tfh','Cycling T','CD4+ Tem'
)
mc@annots<-factor(mc@annots)
mc@colors=mc@annots


library(chameleon)
levels(mc@colors)=distinct_colors(16,minimal_saturation = 12,minimal_lightness = 14,maximal_lightness = 100)$name
levels(mc@colors)=distinct_colors(28,minimal_saturation = 5,minimal_lightness = 15,maximal_lightness = 90)$name

mc@colors<-as.character(mc@colors)
mc@annots<-as.character(mc@annots)

scdb_add_mc("mc_annot",mc)

scdb_mc2d("LUAD_2dproj")->mc2d
pdf("./mc2d.pdf",useDingbats=F,width=8,height=7)
cols = mc@colors
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(mc2d@sc_x, mc2d@sc_y, pch = 19, col=cols[mc@mc[names(mc2d@sc_x)]],bg = cols[mc@mc[names(mc2d@sc_x)]],cex = 0.2, lwd = 0.2)
legend("topright",title="celltype",inset=c(-0.25,0),legend =unique(mc@annots),col=unique(cols),pch=19,cex=0.5,ncol=1) 
dev.off()



mcell_mc2d_plot_by_factor(
  mc2d_id='LUAD_2dproj',
  mat_id='LUAD',
  meta_field='stage',
  meta_data_vals = NULL,
  single_plot = F,
  filter_values = NULL,
  filter_name = NULL,
  ncols = NULL,
  neto_points = F,
  colors = mc@colors
)


###
source("/home/youwh/mouse_CRC/li_et_al_cell_2018_melanoma_scrna-master/util_funcs.r")

mel_plot_e_gc_barplots = function(mc_id, name, genes=NULL, ncolumns=2, panel_height=50, panel_width=300, ord_first_by_color=T, n_ideal_umi=1000) 
{
	mc = scdb_mc(mc_id)
	col2group = c(1:length(mc@mc))
	names(col2group) = as.character(mc@colors)

	if (is.null(genes)) {
		marks_gset = scdb_gset(mc_id)
		genes = names(marks_gset@gene_set)
	}
	e_gc = mc@e_gc[genes, ] * n_ideal_umi
	
	if (ord_first_by_color) {
		e_gc = e_gc[, ord_by_id[[mc_id]]]
	}
	
	.plot_start(scfigs_fn(mc_id, sprintf("mc_geom_mean_%s", name)), w=ncolumns * panel_width, h=panel_height * ceiling(length(genes) / ncolumns))
	layout(matrix(1:(length(genes) + length(genes) %% 2), ncol=ncolumns))
	par(mar=c(0.5, 12, 0.5, 1))
	
	for (g in genes) {
		barplot(e_gc[g, ], border=NA, col=mc@colors[as.numeric(colnames(e_gc))], xaxt='n', yaxt='n', space=0)
		yaxp = par("yaxp")
		axis(2, yaxp=c(yaxp[1], yaxp[2], 1), las=2, cex=1)
		mtext(g, 2, line=1.5, cex=1.2, las=2)
	}
	dev.off()
	
}
genes=c("CD3D",'CD4','CD8A','ITGAE','ZNF683','CXCL13','CXCR5','GZMK','TCF7','LEF1','CCR7','IL7R','CCR6','CD44','GZMB','PRF1','PDCD1','HAVCR2','CTLA4','TIGIT')
genes=c('CD3D','CD4','CD8A','TRAC','NKG7','PRF1','GZMB','GZMK','IFNG','CTLA4','FOXP3','TNFRSF4','TCF7','CCR7','IL7R','GNLY','IL32','CCL5','HAVCR2','LAG3','PDCD1','ISG15','CXCR6','KIT','ICOS','CD79A','IGHG1','IGHD','CD27','CD83','BCL6','CD1C','LAMP3','CLEC9A','CD14','FCGR3A','CD68','CD163','CXCR5','CPA3')
ord_by_id<-list()
ord_by_id[["mc_annot"]]=as.numeric(order(mc@annots))
mel_plot_e_gc_barplots("mc_annot",'Tgene',genes=genes,ncolumns=3)



ord_by_id<-list()
ord_by_id[["mc_annot"]]=as.numeric(order(mc@annots))
mel_plot_e_gc_barplots("mc_annot",'Tgene',genes=genes,ncolumns=2)


###等高线图
data=data.frame(metacell_1=mc2d@sc_x,metacell_2=mc2d@sc_y,stage=aa$stage)
ggplot(data[which(data$stage=='adjacent'),])+
  geom_point(aes(x=metacell_1, y=metacell_2),color=NA)+
  geom_density_2d(aes(x=metacell_1, y=metacell_2))+
  theme_bw()+
  #删除网格线
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(min(data$metacell_1) - 0.1*diff(range(data$metacell_1)),
                              max(data$metacell_1) + 0.1*diff(range(data$metacell_1))))+
  scale_y_continuous(limits=c(min(data$metacell_2) - 0.1*diff(range(data$metacell_2)),
                              max(data$metacell_2) + 0.1*diff(range(data$metacell_2))))




######
mat@cell_metadata[names(mc@mc),c('stage','PatientID')]->aa
aa$type=LUAD@meta.data[rownames(aa),'type_global']
aa$id=paste0(aa$stage,"_",aa$PatientID)
aa$stage<-NULL
aa$PatientID<-NULL
data.frame(table(aa))->aa
for(i in 1:nrow(aa)){
  aa$percentage[i]=aa$Freq[i]/sum(aa$Freq[which(aa$id==aa$id[i])])
}
library(stringr)
aa$stage=substr(aa$id,1,str_locate(as.character(aa$id),"_")[,1]-1)
aa$stage=factor(aa$stage,levels=c("adjacent","AIS","MIA","INV"))
p<-list()
color<-c('B cell'="#673997",'Basophil'="#B35A20",'CD4+ T'="#048D5E",'CD8+ T'="#ADDA81",'Cycling T'="#FD7906",'DC'="#FF9998",'ILC'="#1775B6",'Macrophage'="#EA2087",'Mast cell'="#F7F393",'Monocyte'="#FFBB6D",'NK'="#9ACDE6",'Plasma B'="#ED201E")
for(i in unique(aa$type)){
my_comparisons=list(c('MIA','INV'),c("adjacent",'MIA'),c('adjacent','INV'))
p[[i]]<-ggplot(aa[which(aa$type==i),],aes(x=stage,y=percentage,color=type))+
geom_boxplot(lwd=1)+
geom_point(aes(color=type),size=1,position='jitter')+
scale_colour_manual(values=color[i])+
scale_fill_manual(values=color[i])+
theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title =element_blank(),legend.position = 'none')+
        stat_compare_means(comparisons = my_comparisons) + #添加成对p值
        stat_compare_means()+
ylab("percentage of CD45+ immune cells")+ ggtitle(i)
}
pdf("./percentage_patient.pdf",useDingbats=F,width=15,height=6)
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]], nrow = 2)
dev.off()


color=c('Naive B'='#938A81','Memory B'='#FAD09F','Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED','Plasma B'='#FF913F')


pdf("percentage_pid_bCELL.pdf",width=8,height=5)
ggplot(aa,aes(x=id,y=percentage,fill=type))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=color)+
NoLegend()
dev.off()

###
p<-list()
for(i in unique(aa$type)){
p[[i]]<-ggplot(aa[which(aa$type==i),],aes(x=stage,y=percentage))+
geom_point(aes(fill=stage),color='black',pch=21,size=3,position='jitter')+
scale_fill_manual(values=c("#C5E4CB","#00AB1E","#FF9D9F","#FF0000"))+
stat_boxplot(geom = "errorbar",position = position_dodge(0.5))+
stat_compare_means(comparisons=list(c('adjacent','AIS'),c("adjacent",'INV'),c("adjacent",'MIA'),c('MIA','INV')))

}
pdf("./percentage_type1.pdf",width=20/6*8,height=3)
cowplot::plot_grid(p[[1]],p[[3]],p[[4]],p[[6]],p[[8]],p[[10]],p[[11]],p[[12]], nrow = 1)
dev.off()

pdf("./percentage_type1.pdf",width=20/6*8,height=3)
cowplot::plot_grid(p[[1]],p[[3]],p[[4]],p[[6]],p[[8]],p[[10]],p[[11]],p[[12]], nrow = 1)
dev.off()


####
pdf("percentage_pid.pdf",width=8,height=5)
ggplot(aa,aes(x=id,y=percentage,fill=type))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=color)+
NoLegend()
dev.off()
###
pdf("percentage_stage.pdf",width=1.9,height=1.7)
ggplot(bb,aes(x=stage,y=percentage,fill=type))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=color)+
NoLegend()
dev.off()



####相关性
porpation=list()

porpation[['MIA']]= matrix(nrow=22,ncol=31)
colnames(porpation[['MIA']])=as.character(unique(aa$type))
rownames(porpation[['MIA']])=as.character(unique(aa$id[which(aa$stage=='MIA')]))
for(i in 1:ncol(porpation[['MIA']])){
  porpation[['MIA']][,i]=aa[which(aa$stage=='MIA'& aa$type==colnames(porpation[['MIA']])[i]),]$percentage
}


porpation[['INV']]= matrix(nrow=33,ncol=32)
colnames(porpation[['INV']])=as.character(unique(aa$type))
rownames(porpation[['INV']])=as.character(unique(aa$id[which(aa$stage=='INV')]))
for(i in 1:ncol(porpation[['INV']])){
  porpation[['INV']][,i]=aa[which(aa$stage=='INV'& aa$type==colnames(porpation[['INV']])[i]),]$percentage
}



porpation[['AIS']]= matrix(nrow=4,ncol=32)
colnames(porpation[['AIS']])=as.character(unique(aa$type))
rownames(porpation[['AIS']])=as.character(unique(aa$id[which(aa$stage=='AIS')]))
for(i in 1:ncol(porpation[['AIS']])){
  porpation[['AIS']][,i]=aa[which(aa$stage=='AIS'& aa$type==colnames(porpation[['AIS']])[i]),]$percentage
}

porpation[['adjacent']]= matrix(nrow=34,ncol=32)
colnames(porpation[['adjacent']])=as.character(unique(aa$type))
rownames(porpation[['adjacent']])=as.character(unique(aa$id[which(aa$stage=='adjacent')]))
for(i in 1:ncol(porpation[['adjacent']])){
  porpation[['adjacent']][,i]=aa[which(aa$stage=='adjacent'& aa$type==colnames(porpation[['adjacent']])[i]),]$percentage
}

porpation[['total']]= matrix(nrow=93,ncol=32)
colnames(porpation[['total']])=as.character(unique(aa$type))
rownames(porpation[['total']])=as.character(unique(aa$id))
for(i in 1:ncol(porpation[['total']])){
  porpation[['total']][,i]=aa[which(aa$type==colnames(porpation[['total']])[i]),]$percentage
}


porpation[['tumor']]= matrix(nrow=33+4+22,ncol=32)
colnames(porpation[['tumor']])=as.character(unique(aa$type))
rownames(porpation[['tumor']])=as.character(unique(aa$id[which(aa$stage %in% c('AIS','MIA','INV'))]))
for(i in 1:ncol(porpation[['tumor']])){
  porpation[['tumor']][,i]=aa[which(aa$stage %in% c('AIS','MIA','INV') & aa$type==colnames(porpation[['tumor']])[i]),]$percentage
}


col = colorRampPalette(c("blue", "white", "red"))(100)

pdf("percentage_corrlation.pdf")
heatmap(cor(porpation[['adjacent']]),col=col,symm = TRUE)
heatmap(cor(porpation[['AIS']]),col=col,symm = TRUE)
heatmap(cor(porpation[['MIA']]),col=col,symm = TRUE)
heatmap(cor(porpation[['INV']]),col=col,symm = TRUE)
heatmap(cor(porpation[['total']]),col=col,symm = TRUE)
dev.off()
heatmap(cor(porpation[['tumor']]),col=col,symm = TRUE)

######
mat@cell_metadata[names(mc@mc),c('stage','PatientID')]->bb
bb$type=factor(mc@mc)
levels(bb$type)=mc@annots
bb$PatientID<-NULL
data.frame(table(bb))->bb
for(i in 1:nrow(bb)){
  bb$percentage[i]=bb$Freq[i]/sum(bb$Freq[which(bb$stage==bb$stage[i])])
}
bb[order(bb$stage),]->bb
bb$order=c(rep(1:16,4))
bb$stage=factor(bb$stage,levels=c("adjacent","AIS","MIA","INV"))
library(ggalluvial)
library(reshape)
pdf("./percentage_sankey.pdf")
ggplot(bb,aes(x = stage, stratum =type, 
        alluvium = order,y = percentage,
        fill = type, label = type)) +
  scale_x_discrete(expand = c(.1, .1)) +
  scale_fill_manual(values=unique(mc@colors))+
  geom_flow() +
  geom_stratum(alpha = .5) +
  #geom_text(stat = "stratum", size = 2.5) +
  theme(panel.border=element_rect(color='black', fill=NA), 
       panel.grid.major =element_blank(), 
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text=element_text(colour="black")) 
dev.off()

#B cell
#提取mc 
mc = scdb_mc('mc_annot')
mc@color_key=data.frame(group=mc@annots,color=mc@colors)
ucolkey =  unique(mc@color_key[, c('group', 'color')])
col2grp = ucolkey$group
names(col2grp) = ucolkey$color
cg = split(names(mc@mc), col2grp[mc@colors[mc@mc]])
nms = c(cg[['B cell']])
mc_map = 1:length(unique(mc@mc[nms]))
names(mc_map) = names(table(mc@mc[nms]))
dst_mc = mc_map[as.character(mc@mc[nms])]
names(dst_mc) = nms
new_id = paste("submc", "Bcell", sep="_")
mcell_mat_ignore_cells(new_id, "LUAD", nms, reverse=T)
mcell_add_gene_stat(new_id, new_id)
mcell_new_mc(mc_id = new_id,mc = dst_mc,outliers = character(0), #setdiff(c(mc@outliers, names(mc@mc)), nms),
									 scmat = scdb_mat(new_id))
mc_Bcell=scdb_mc("submc_Bcell")
mat_Bcell=scdb_mat("submc_Bcell")
mc_compute_fp(mc_Bcell ,mat_Bcell@mat , norm_by_mc_size  =  T , min_total_umi  =  10 )->lfp
colnames(lfp)=names(mc_map)

lat_genes=names(scdb_gset("LUAD_feats")@gene_set)
B_genes = intersect(names(which(apply(abs(lfp),1,max)>log2(3))),lat_genes)
mat_b_ds = scm_downsamp(mat_Bcell@mat, 1000)
library(tgstat)
B_cor_c = tgs_cor(t(as.matrix(mat_b_ds[B_genes, ])), spearman=T)
diag(B_cor_c) = NA
B_cor_mc=cor(t(lfp[B_genes,]))
diag(B_cor_mc) = NA
blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
pdf(file='./figs/submc_B_cor.pdf', width=max(700, 300 + length(B_genes) * 12)/72, height=max(700, 300 + length(B_genes) * 12)/72)
pheatmap(pmin(pmax(B_cor_mc, -0.8), 0.8), clustering_method="ward.D2", cutree_rows=15, cutree_cols=15, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols)
dev.off()

colnames(mc_Bcell@e_gc)=names(mc_map)
e_gc = mc_Bcell@e_gc[c(B_genes), ] * 1000
pdf(file='./figs/submc_B_exp.pdf', width=max(700, 300 + length(B_genes) * 12)/(72*2), height=max(700, 300 + length(B_genes) * 12)/(72*1.1))
pheatmap(e_gc,clustering_method='ward.D2',border_color=NA,col=blwtrd_cols,scale='row',cellwidth=10, cellheight=10, fontsize=12,breaks=unique(c(seq(-4,4, length=101))))
dev.off()
###画图
library(randomcoloR)

mc_Bcell@annots<-factor(names(mc_map))
levels(mc_Bcell@annots)=mc@annots[119:141]
mc_Bcell@colors=mc_Bcell@annots
levels(mc_Bcell@colors)=distinctColorPalette(k = 4)
mc_Bcell@colors<-as.character(mc_Bcell@colors)
mc_Bcell@annots<-as.character(mc_Bcell@annots)

pdf("./figs/mc2d_B.pdf",useDingbats=F,width=8,height=7)
cols = mc_Bcell@colors
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(mc2d@sc_x[names(mc_Bcell@mc)], mc2d@sc_y[names(mc_Bcell@mc)], pch = 19, col=cols[mc_Bcell@mc],bg = cols[mc_Bcell@mc],cex = 0.7, lwd = 0.5)
legend("topright",title="celltype",inset=c(-0.25,0),legend =unique(mc_Bcell@annots),col=unique(cols),pch=19,cex=0.5,ncol=1) 
dev.off()
##
mat@cell_metadata[c(names(mc_Bcell@mc),names(mc@mc[which(mc@mc=='171')])),c('stage','PatientID')]->bb
bb$type=factor(mc@mc[rownames(bb)])
levels(bb$type)=c(mc@annots[119:141],'plasma B')
bb$id=paste0(bb$stage,"_",bb$PatientID)
bb$stage<-NULL
bb$PatientID<-NULL
data.frame(table(bb))->bb
for(i in 1:nrow(bb)){
  bb$percentage[i]=bb$Freq[i]/sum(bb$Freq[which(bb$id==bb$id[i])])
}
library(stringr)
bb$stage=substr(bb$id,1,str_locate(as.character(bb$id),"_")[,1]-1)
bb$stage=factor(bb$stage,levels=c("adjacent","AIS","MIA","INV"))

pdf("./percentage_patient_B.pdf",useDingbats=F,width=5,height=5)
for(i in unique(bb$type)){
my_comparisons=list(c('adjacent','AIS'),c("AIS",'MIA'),c('MIA','INV'))
p<-ggplot(bb[which(bb$type==i),],aes(x=stage,y=percentage,color=stage))+
geom_boxplot(lwd=1)+
geom_point(aes(color=stage),size=2,position='jitter')+
scale_color_manual(values=c("#367DC0","#72B626","#F19598","#D11E1A"))+
theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title =element_blank(),legend.position = 'none')+
        stat_compare_means(comparisons = my_comparisons) + #添加成对p值
        stat_compare_means()+
ylab("percentage of total B cells")+ ggtitle(i)
print(p)
}
dev.off()


####T cell
mc = scdb_mc('mc_annot')
mc@color_key=data.frame(group=mc@annots,color=mc@colors)
ucolkey =  unique(mc@color_key[, c('group', 'color')])
col2grp = ucolkey$group
names(col2grp) = ucolkey$color
cg = split(names(mc@mc), col2grp[mc@colors[mc@mc]])
nms =c(cg[['T cell']])

mc_map = 1:length(unique(mc@mc[nms]))
names(mc_map) = names(table(mc@mc[nms]))
dst_mc = mc_map[as.character(mc@mc[nms])]
names(dst_mc) = nms
new_id = paste("submc", "T", sep="_")
mcell_mat_ignore_cells(new_id, "LUAD", nms, reverse=T)
mcell_add_gene_stat(new_id, new_id)

mcell_new_mc(mc_id = new_id,mc = dst_mc,outliers = character(0), #setdiff(c(mc@outliers, names(mc@mc)), nms),
									 scmat = scdb_mat(new_id))
mc_T=scdb_mc("submc_T")
mat_T=scdb_mat("submc_T")
mc_compute_fp(mc_T ,mat_T@mat , norm_by_mc_size  =  T , min_total_umi  =  10 )->lfp
colnames(lfp)=names(mc_map)

lat_genes=names(scdb_gset("LUAD_feats")@gene_set)
T_genes = intersect(names(which(apply(abs(lfp),1,max)>log2(5))),lat_genes)
mat_b_ds = scm_downsamp(mat_T@mat, 500)
library(tgstat)
T_cor_c = tgs_cor(t(as.matrix(mat_b_ds[T_genes, ])), spearman=T)
diag(T_cor_c) = NA
T_cor_mc=cor(t(lfp[T_genes,]))
diag(T_cor_mc) = NA
blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
pdf(file='./figs/submc_T_cor.pdf', width=max(700, 300 + length(T_genes) * 12)/72, height=max(700, 300 + length(T_genes) * 12)/72)
pheatmap(pmin(pmax(T_cor_mc, -0.7), 0.7), clustering_method="ward.D2", cutree_rows=15, cutree_cols=15, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols)
dev.off()


colnames(mc_T@e_gc)=names(mc_map)
e_gc = mc_T@e_gc[T_genes, ] * 1000
pdf(file='./figs/submc_T_exp.pdf', width=max(700, 300 + length(T_genes) * 12)/72, height=max(700, 300 + length(T_genes) * 12)/72)
pheatmap(e_gc,clustering_method='ward.D2',border_color=NA,col=blwtrd_cols,scale='row',cellwidth=10, cellheight=10, cutree_rows=15, cutree_cols=15,fontsize=12,breaks=unique(c(seq(-6,6, length=101))))
dev.off()

##
library(randomcoloR)
mc_T@annots=names(mc_map)
mc_T@annots=mc@annots[as.numeric(as.character(mc_T@annots))]
mc_T@colors=factor(mc_T@annots)
levels(mc_T@colors)=distinctColorPalette(k = 8)
mc_T@colors<-as.character(mc_T@colors)
mc_T@annots<-as.character(mc_T@annots)

pdf("./figs/mc2d_T.pdf",useDingbats=F,width=8,height=7)
cols = mc_T@colors
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(mc2d@sc_x[names(mc_T@mc)], mc2d@sc_y[names(mc_T@mc)], pch = 19, col=cols[mc_T@mc],bg = cols[mc_T@mc],cex = 0.3, lwd = 0.5)
legend("topright",title="celltype",inset=c(-0.25,0),legend =unique(mc_T@annots),col=unique(cols),pch=19,cex=0.5,ncol=1) 
dev.off()



#####着色图
cell_stats=mat@cell_metadata[names(mc2d@sc_x),]
fat <- data.frame(mc2d@sc_x, mc2d@sc_y,
                  cell_stats$PatientID,
                  cell_stats$stage,
                  cell_stats$amp_batch_id,
                  cell_stats$seq_batch_id,
                  cell_stats$location,
                  mc@mc[names(mc2d@sc_x)],
                  mc@mc[names(mc2d@sc_x)])  
colnames(fat) <- c("x", "y", "Donor","Stage", "amp_batch_id","seq_batch_id", "location","Type","color")
fat$Type <- factor(fat$Type)
levels(fat$Type)=mc@annots
fat$Type<-as.character(fat$Type)
fat$color <- factor(fat$color)
levels(fat$color)=mc@colors
fat$color<-as.character(fat$color)
dim(fat)
fat$CellID=names(mc2d@sc_x)


fat.scatter2 <- fat[which(fat$Type=='T cell'), ]
dim(fat.scatter2)

fat.scatter3 <- fat[which(fat$Type=='B cell'), ]

theme_publa <- function(base_size = 8, base_family = "sans", legend_position = "right",
                     title_size=8){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
  theme(legend.position = legend_position, legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000"), axis.ticks = element_line(colour = "#000000"), 
        legend.key = element_blank(), 
        axis.text = element_text(size = base_size, face="plain"), plot.title=element_text(face="plain", size = title_size),
        axis.title = element_text(face="plain", size = base_size), legend.text=element_text(size = base_size),
        legend.title=element_text(face="plain", size = base_size), strip.text=element_text(face="plain", size = base_size)
        )
}

plotgeneexp<-function(genes,mat_id,mc_id,mc2d_id,fat,spec=NULL,figwidth = 200/4, figheight = figwidth, figtextsize = 8, dotsize = 0.5, stroke = 0.05){
     mat=scdb_mat(mat_id)
     mc=scdb_mc(mc_id)
     mc2d=scdb_mc2d(mc2d_id)
     for(gene in genes){
    # or fat.scatter2 for only annotated celltypes
    spot <- data.frame(mc2d@sc_x[intersect(colnames(mat@mat),unique(fat$CellID))], 
                       mc2d@sc_y[intersect(colnames(mat@mat),unique(fat$CellID))],
                       mat@mat[gene,intersect(colnames(mat@mat),unique(fat$CellID))]
    )
    colnames(spot) <- c("x", "y", "value")
    spot <- spot[order(as.numeric(factor(spot$value))),]
    require(ggplot2)
    p1<-ggplot2::ggplot(spot,aes(x=x,y=y, colour=log2(value+1),fill=log2(value+1)))+
        geom_point(size=dotsize, shape=21, stroke = stroke*(figwidth/200))+
        scale_colour_gradientn(name=NULL, breaks = c(0, floor(max(log2(na.omit(spot$value+1))))),
                             colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)) +
        scale_fill_gradientn(name=NULL, breaks = c(0, floor(max(log2(na.omit(spot$value+1))))),
                           colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                             c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                               "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)) +
        theme_publa(base_size = figtextsize) +
      theme(legend.position = "bottom", legend.background = element_blank(), 
            legend.text = element_text(size=figtextsize), 
            legend.title = element_text(size=figtextsize), plot.title = element_text(size=figtextsize),
            axis.line.x = element_blank(), axis.line.y = element_blank(),
            axis.ticks = element_blank(), axis.text = element_blank(),
            legend.key.size = unit(1.5,"mm")
      ) + guides(colour = guide_colorbar(ticks = FALSE))
      p1
    
    ggsave(filename = paste0(getwd(),"/figs/Spotgenes_",mat_id, "_", spec, "_" , gene,".pdf"),
           width=figwidth, height=figheight, units="mm", dpi=150, useDingbats=FALSE)

     }  
}
plotgeneexp(c('TCF7','CCR7','LEF1','CCR6','GZMK'),'LUAD','mc_annot','LUAD_2dproj',fat=fat.scatter2,figwidth = 40, figheight=50, figtextsize = 8, dotsize = 0.2, stroke=0.05, spec = "Tcell")

plotgeneexp(genes,'LUAD','mc_annot','LUAD_2dproj',fat=fat,figwidth = 40, figheight=50, figtextsize = 8, dotsize = 0.4, stroke=0.05, spec = "cell")
mc@annots[152:181]=c(rep("Memory B cell",13),rep("GC B cell"),'Memory B cell','Memory B cell',rep("Naive B cell",5),rep("Activated B cell",9))


plotgeneexp(c('TCL1A','CD27','IGHD','BCL6','MZB1','IGHA1','IGHG1','IGHG2','IGHG4','IGHM','CD83','FCRL4'),'LUAD','mc_annot','LUAD_2dproj',fat=fat,figwidth = 35, figheight=50, figtextsize = 8, dotsize = 0.5, stroke=0.05, spec = "Bcell")

####
mctoseurat("LUAD",'mc_annot_test','LUAD_2dproj')->LUAD
subset(LUAD,type %in% c('CD8+ memory T',"CD4+ memory T","Cytotoxic T","Treg","Tfh","Dysfunctional T","Memory B","plasma B","Naive B","Translational T","MAIT","Activated B","GC B"))->TB
mypvals <- read.delim("./pvalues.txt", check.names = FALSE)
mymeans <- read.delim("./means.txt", check.names = FALSE)
mypvals[!duplicated(mypvals$interacting_pair),]->mypvals
mymeans[!duplicated(mymeans$interacting_pair),]->mymeans


rownames(mypvals)<-mypvals$interacting_pair
pval.f<-mypvals[,12:ncol(mypvals)]
RR=c("Activated B|Cytotoxic T","Activated B|Dysfunctional T","Activated B|MAIT","Activated B|Naive/Memory T","Activated B|Tfh","Activated B|Translational T","Activated B|Treg",
"Naive B|Cytotoxic T","Naive B|Dysfunctional T","Naive B|MAIT","Naive B|Naive/Memory T","Naive B|Tfh","Naive B|Translational T","Naive B|Treg",
"plasma B|Cytotoxic T","plasma B|Dysfunctional T","plasma B|MAIT","plasma B|Naive/Memory T","plasma B|Tfh","plasma B|Translational T","plasma B|Treg",
"GC B|Cytotoxic T","GC B|Dysfunctional T","GC B|MAIT","GC B|Naive/Memory T","GC B|Tfh","GC B|Translational T","GC B|Treg",
"Memory B|Cytotoxic T","Memory B|Dysfunctional T","Memory B|MAIT","Memory B|Naive/Memory T","Memory B|Tfh","Memory B|Translational T","Memory B|Treg")

PR=c("Cytotoxic T|Activated B","Cytotoxic T|Memory B","Cytotoxic T|Naive B","Cytotoxic T|GC B",
"Dysfunctional T|Activated B","Dysfunctional T|Memory B","Dysfunctional T|Naive B","Dysfunctional T|GC B",
"Translational T|Activated B","Translational T|Memory B","Translational T|Naive B","Translational T|GC B",
"Naive/Memory T|Activated B","Naive/Memory T|Memory B","Naive/Memory T|Naive B","Naive/Memory T|GC B",
"Treg|Activated B","Treg|Memory B","Treg|Naive B","Treg|GC B",
"MAIT|Activated B","MAIT|Memory B","MAIT|Naive B","MAIT|GC B",
"Tfh|Activated B","Tfh|Memory B","Tfh|Naive B","Tfh|GC B")

pval.f[,PR]->pval.f
index<-c()
for(i in 1:nrow(pval.f)){
	if(length(which(pval.f[i,]<0.05))>0){
		index[i]=i
	}else{
		index[i]=NA
	}
}
as.numeric(as.character(na.omit(index)))->index
pval.f[index,]->pval.f

rownames(mymeans)=mymeans$interacting_pair
means.f<-mymeans[,12:ncol(mymeans)]
means.f<-means.f[rownames(pval.f),PR]

df_names = expand.grid(rownames(pval.f), colnames(pval.f))
pval.f1 = unlist(pval.f)
pval.f1[pval.f1==0] = 0.00009
plot.data = cbind(df_names,pval.f1)
pr = unlist(as.data.frame(means.f))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')


pdf("./interaction_TB.pdf",useDingbats=F,width=22,height=16)
ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=scale(mean))) +
   scale_color_distiller(palette = "RdBu")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(size=7,angle = 90,colour="black", hjust = 1),
        axis.text.y = element_text(size=6, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
dev.off()


############
mc_id='mc_annot'
mat_id='LUAD'
mc = scdb_mc(mc_id)
mat = scdb_mat(mat_id)
dat=data.frame(patient=mat@cell_metadata[names(mc@mc), 'PatientID'],mc=mc@mc)
dat$type=factor(dat$mc)
levels(dat$type)=mc@annots
dat$type=as.character(dat$type)
table(dat[,c('patient','mc')])->pat_grp

pat_grp_n=pat_grp / rowSums(pat_grp)
pat_grp_n[,dat$mc[order(dat$type)]]->pat_grp_n


merge <- RunVelocity(merge, deltaT = 1, kCells = 25, fit.quantile = 0.02, 
        spliced.average = 0.2, unspliced.average = 0.05, ncores = 10)



ident.colors <- (scales::hue_pal())(n = length(x = levels(merge$cluster)))
names(x = ident.colors) <- levels(merge$cluster)
cell.colors <- ident.colors[merge$cluster]
names(x = cell.colors) <- colnames(x = merge)

emb = Embeddings(merge, reduction = "umap")
vel = Tool(merge, slot = "RunVelocity")
show.velocity.on.embedding.cor(emb = emb, vel = vel, n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)


#####

library(survival)
library(survminer)
t(score)->score
cbind(eLUAD$OS,score)->survival
for(j in colnames(score)){
  val=quantile(as.numeric(survival[,j]),1/3)
  for(i in rownames(survival)){
    if(survival[i,j]<val){
      survival[i,j]=0
    }else{
      survival[i,j]=1
    }
  }
}
covariates <-colnames(OS)[35:61]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(PFI.time/30, PFI)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = OS)})
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value <- signif(x$waldtest["pvalue"], digits = 2)
                         wald.test <- signif(x$waldtest["test"], digits = 2)
                         beta <- signif(x$coef[1], digits = 2);
                         HR <- signif(x$coef[2], digits = 2);
                         HR.confint.lower <- signif(x$conf.int[ ,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[ ,"upper .95"], 2)
                         HR <- paste0(HR, " (",
                                     HR.confint.lower, "-", HR.confint.upper, ") ")
                         res <- c(beta, HR, wald.test, p.value,HR.confint.lower,HR.confint.upper)
                         names(res) <- c("beta", "HR(95% CI for HR)", "wald.test",
                                         "p.value","HR.confint.lower","HR.confint.upper")
                         return(res)
                       })
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.name = F))
as.data.frame(res)->res
##
HR=gsub('[\\(\\)]','',res$'HR(95% CI for HR)')
HR=gsub("-",' ',HR)
HR=as.data.frame(do.call(cbind,strsplit(HR,' ')),stringsAsFactors=FALSE)
names(HR)=rownames(res)

pdf("./Project/LUAD/figs/forest.pdf")
par(mar=c(5,6,4,13))
plot(as.numeric(HR[1,]),1:dim(HR)[2],
pch=15,cex=2,col='blue',bty='n',yaxt='n',ylab=NA,xlab='HR',
xlim=range(as.numeric(unlist(HR))))
abline(v=1,col='grey',lwd=2,lty=2)

for(i in 1:ncol(HR)){
  x=as.numeric(HR[2:3,i])
  lines(x,c(i,i),col='blue')
  text(0.2,i,rownames(res)[i],xpd=T,adj=c(0,0))
  text(2.1,i,as.numeric(res[i,4]),xpd=T,adj=c(0,0))
  text(2.7,i,as.character(res[i,2]),xpd=T,adj=c(0,0))
}
text(2.1,ncol(HR)+0.5,'pvalue',xpd=T,adj=c(0,0))
text(2.7,ncol(HR)+0.5,'HR(CI)',xpd=T,adj=c(0,0))


dev.off()

##
 ggplot(group, aes(x=Pathological.stage, y=cNK,fill=Pathological.stage)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  chameleon::scale_fill_chameleon()+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  ylab("Value")+xlab("")+ #设置x轴和y轴的标题
  stat_compare_means(comparisons=list(c('AIS','stage I'),c('MIA','stage I')))

c("#673997","#B35A20","#048D5E","#ADDA81","#FD7906","#FF9998","#1775B6","#EA2087","#F7F393","#FFBB6D","#9ACDE6","#ED201E")

####
mel_mc2d_plot_genes_of_interest = function(mc_id, mc2d_id=mc_id, min_max_lfp=1)
{
	mc = scdb_mc(mc_id)
	lfp = log2(mc@mc_fp)
	genes_of_interest=c('CD3D','CD4','CD8A','FOXP3','LAG3','CTLA4','GZMB','IL7R','PRF1','GNLY','NKG7','CD79A','MZB1','C1QA','CD14','FCGR3A','LAMP3','CLEC10A','CLEC9A')

	for (gene in genes_of_interest) {
		mcell_mc2d_plot_gene(mc2d_id, gene, show_legend=T, neto_points=T)
	}
	
}
mel_mc2d_plot_genes_of_interest('mc_annot_global','LUAD_2dproj')
###
library(Biobase)
library(BisqueRNA)
bulk.eset <- Biobase::ExpressionSet(assayData = NC$exp)
sample.ids <- colnames(LUAD)
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=LUAD$PatientID,
                       cellType=LUAD$type)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=as_matrix(LUAD@assays$RNA@counts),
                                  phenoData=sc.pdata)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
ref.based.estimates <- res$bulk.props
 data.frame(prop=as.numeric(ref.based.estimates['CD27+ Memory B',]),group=NC$group)->ACB
 levels(ACB$group.Pathological.stage)=c('AIS','INV','INV','IIIA','MIA')



 ###TRUST4
 /home/youwh/bioapp/TRUST4/run-trust4 -b CHG020333Aligned.sortedByCoord.out.bam \
-f /home/youwh/bioapp/TRUST4/hg38_bcrtcr.fa \
--ref /home/youwh/ref/IMGT+C.fa


cat id.txt | while read line
do
samtools fastq $line/*.bam -@ 10  -1 _1.fq -2 _2.fq
done



ggplot(aa,aes(x=amp_batch_id,y=nCount_RNA,color=stage))+
geom_jitter(size=0.3)+
scale_color_manual(values=c('adjacent'='#C5E4CB','AIS'='#00AB1E','MIA'='#FF9D9F','INV'='#FF0000'))+
scale_y_continuous(limits=c(2,4),breaks=seq(2,4,0.5))+
theme_bw()+
theme(panel.grid=element_blank())+
theme(axis.text.x=element_blank())+
theme(axis.ticks.x=element_blank())+
ylab('Number of UMIs')+
xlab("")

p1<-ggplot(number,aes(x=Var1,y=percentage,color=stage))+
geom_jitter(size=0.4)+
scale_color_manual(values=c('adjacent'='#C5E4CB','AIS'='#00AB1E','MIA'='#FF9D9F','INV'='#FF0000'))+
theme_bw()+
theme(panel.grid=element_blank())+
theme(axis.text.x=element_blank())+
theme(axis.ticks.x=element_blank())+
ylab('')+
xlab("")

##DotPlot
library(tidyverse)
library(cowplot)
library(ggdendro)
library(ggtree)
markers=c('CD3D','CD8A','CD4','LEF1','TCF7','IL7R','GZMK','FOXP3','GZMB','PRF1','SLC4A10','KLRB1','CXCL13','ISG15','ZNF683','ITGAE','NKG7','GNLY','CX3CR1','ICOS','KIT','CD79A','MZB1','BCL6','CD27','CD83','CD14','FCGR3A','CD68','IL1B','APOC1','CLEC9A','CD1C','LAMP3','LILRA4','TPSAB1','CPA3')

LUAD <- NormalizeData(LUAD, normalization.method = "LogNormalize", scale.factor = 10000)
LUAD <- ScaleData(LUAD, features = rownames(LUAD))


ggplot(num,aes(x=num))+
geom_bar(stat='bin',bins=20)+
theme_bw()+
theme(panel.grid=element_blank())+
xlim(30, 60)
DimPlot(LUAD,reduction='metacell',group.by='type_global',split.by='stage',ncol=2)+scale_color_manual(values=c('B cell'="#673997",'Basophil'="#B35A20",'CD4+ T'="#048D5E",'CD8+ T'="#ADDA81",'Cycling T'="#FD7906",'DC'="#FF9998",'ILC'="#1775B6",'Macrophage'="#EA2087",'Mast cell'="#F7F393",'Monocyte'="#FFBB6D",'NK'="#9ACDE6",'Plasma B'="#ED201E"))
##

table(paste0(LUAD$PatientID,'_',LUAD$stage),LUAD$type)->tab
tab <- tab/rowSums(tab)

tab <- log10(tab+1e-3)
dists1 <- matrix(NA,nrow(tab),nrow(tab),dimnames=list(rownames(tab),rownames(tab)))

for(row in rownames(dists1)){
  for(col in colnames(dists1)){
    dists1[row,col] <- sqrt(sum((tab[row,]-tab[col,])^2))
  }
}

diag(dists1) <- NA

dists_normal <- dists1[grepl("adjacent",rownames(dists1)),grepl("adjacent",colnames(dists1))]
dists_normal <- dists_normal[upper.tri(dists_normal,diag=F)]
dists_AIS<- dists1[grepl("AIS",rownames(dists1)),grepl("AIS",colnames(dists1))]
dists_AIS <- dists_AIS[upper.tri(dists_AIS,diag=F)]
dists_MIA<- dists1[grepl("MIA",rownames(dists1)),grepl("MIA",colnames(dists1))]
dists_MIA <- dists_MIA[upper.tri(dists_MIA,diag=F)]
dists_INV<- dists1[grepl("INV",rownames(dists1)),grepl("INV",colnames(dists1))]
dists_INV <- dists_INV[upper.tri(dists_INV,diag=F)]
##INV vs adjacet
dists_cross1 <- dists1[grepl("INV",rownames(dists1)),grepl("adjacent",colnames(dists1))]
dists_cross2 <- dists1[grepl("MIA",rownames(dists1)),grepl("adjacent",colnames(dists1))]
dists_cross3 <- dists1[grepl("INV",rownames(dists1)),grepl("MIA",colnames(dists1))]


INV <- unlist(lapply(strsplit(rownames(dists_cross1),"_"),function(x){x[1]}))
normal <- unlist(lapply(strsplit(colnames(dists_cross1),"_"),function(x){x[1]}))
MIA <- unlist(lapply(strsplit(colnames(dists_cross3),"_"),function(x){x[1]}))


for(norm_iter in normal){
  dists_cross1[which(INV==norm_iter),which(normal==norm_iter)] <- NA
}
for(norm_iter in normal){
  dists_cross2[which(MIA==norm_iter),which(normal==norm_iter)] <- NA
}

dists_cross1 <- array(dists_cross1)
dists_cross2 <- array(dists_cross2)
dists_cross3 <- array(dists_cross3)


dists2 <- c(dists_normal,dists_AIS,dists_MIA,dists_INV,dists_cross1,dists_cross2,dists_cross3)
labs <- rep(c("normal","AIS",'MIA','INV',"cross1",'cross2','cross3'),times=c(length(dists_normal),length(dists_AIS),length(dists_MIA),length(dists_INV),length(dists_cross1),length(dists_cross2),length(dists_cross3)))
data.frame(dist=dists2,lab=labs)->dists3

dists3$lab=factor(dists3$lab,levels=c('normal',"AIS",'MIA','INV','cross2','cross3',"cross1"))
pdf("Euclideandistance_stage.pdf",useDingbats=F,width=9,height=7)
ggplot(dists3,aes(x=lab,y=dist))+
stat_boxplot(aes(x=lab,y=dist,fill=lab),geom='errorbar',width=0.15,position=position_dodge(0.8))+
geom_boxplot(aes(fill=lab),width=0.6,outlier.color = "white",position=position_dodge(0.8))+
chameleon::scale_fill_chameleon()+
geom_jitter(size=0.1)+
theme_bw()+
theme(panel.grid=element_blank())+
#stat_compare_means(comparisons = my_comparisons) + 
ylab('Euclidean distances')+
coord_flip()
dev.off()

my_comparisons=list(c('normal','cross1'),c("cross2",'cross1'),c("cross2",'normal'))
##
mc@n_bc->aa
aa[setdiff(rownames(aa),'AB862'),]->aa
cor(t(aa),method='pearson')->aa.cor
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)[25:1]
pheatmap(aa.cor,col=coul)


##
color<-c('B cell'="#673997",'Basophil'="#B35A20",'CD4+ T'="#048D5E",'CD8+ T'="#ADDA81",'Cycling T'="#FD7906",'DC'="#FF9998",'ILC'="#1775B6",'Macrophage'="#EA2087",'Mast cell'="#F7F393",'Monocyte'="#FFBB6D",'NK'="#9ACDE6",'Plasma B'="#ED201E")

for(i in 1:nrow(AIS_num)){
  AIS_num$percent[i]=AIS_num$Freq[i]/sum(AIS_num$Freq[which(AIS_num$Var1==AIS_num$Var1[i])])
}
AIS_num$Var2=factor(AIS_num$Var2,levels=c("CD4+ T","CD8+ T","Cycling T","Basophil","DC",'NK','ILC','Monocyte',"Macrophage","Mast cell","Plasma B","B cell"))
AIS_num$Var1=factor(AIS_num$Var1,levels=levels(AIS_num$Var1)[order(AIS_num[which(AIS_num$Var2=='B cell'),'percent'])])
pdf("./youwh_work/Project_eLUAD/LUAD/figs/percentage_AIS_patient.pdf")
ggplot(AIS_num,aes(x=Var1,y=percent,fill=Var2))+
geom_col(position="stack",width=0.6)+
scale_fill_manual(values=color)+
theme_bw()+
theme(panel.grid=element_blank())
dev.off()

data.frame(INV_num)->INV_num
for(i in 1:nrow(INV_num)){
  INV_num$percent[i]=INV_num$Freq[i]/sum(INV_num$Freq[which(INV_num$Var1==INV_num$Var1[i])])
}
INV_num$Var2=factor(INV_num$Var2,levels=c("CD4+ T","CD8+ T","Cycling T","Basophil","DC",'NK','ILC','Monocyte',"Macrophage","Mast cell","Plasma B","B cell"))
INV_num$Var1=factor(INV_num$Var1,levels=levels(INV_num$Var1)[order(INV_num[which(INV_num$Var2=='B cell'),'percent'])])
pdf("./youwh_work/Project_eLUAD/LUAD/figs/percentage_INV_patient.pdf")
ggplot(INV_num,aes(x=Var1,y=percent,fill=Var2))+
geom_col(position="stack",width=0.6)+
scale_fill_manual(values=color)+
theme_bw()+
theme(panel.grid=element_blank())
dev.off()

total_num=data.frame(table(paste0(LUAD$stage,'_',LUAD$PatientID),LUAD$type_global))
for(i in 1:nrow(total_num)){
  total_num$percent[i]=total_num$Freq[i]/sum(total_num$Freq[which(total_num$Var1==total_num$Var1[i])])
}
total_num$Var2=factor(total_num$Var2,levels=c("CD4+ T","CD8+ T","Cycling T","Basophil","DC",'NK','ILC','Monocyte',"Macrophage","Mast cell","Plasma B","B cell"))
total_num$Var1=factor(total_num$Var1,levels=c(paste0('adjacent_',levels(Normal_num$Var1)),paste0('AIS_',levels(AIS_num$Var1)),paste0('MIA_',levels(MIA_num$Var1)),paste0('INV_',levels(INV_num$Var1))))
pdf("./youwh_work/Project_eLUAD/LUAD/figs/percentage_total_patient.pdf",width=10)
ggplot(total_num,aes(x=Var1,y=percent,fill=Var2))+
geom_col(position="stack",width=0.6)+
scale_fill_manual(values=color)+
theme_bw()+
theme(panel.grid=element_blank())
dev.off()

too-many-cells make-tree \
    --matrix-path ./input/tumor_matrix.csv \
    --labels-file ./input/tumor_label.csv \
    --draw-collection "PieRing" \
    --output ./tumor_out/ 

too-many-cells make-tree --prior ./tumor_out/ --labels-file ./input/tumor_label.csv \
--draw-colors "[\"#673997\", \"#B35A20\", \"#048D5E\",\"#ADDA81\",\"#FD7906\",\"#FF9998\",\"#1775B6\",\"#EA2087\",\"#F7F393\",\"#FFBB6D\",\"#9ACDE6\",\"#ED201E\",\"#2FA62A\",\"#C8AFD4\"]" --smart-cutoff 2 -M 1 \
--draw-collection "PieChart" \
--dendrogram-output "tree_labeled_alternate.svg" \
--output ./tumor_out_pruned 
####
for(i in 1:nrow(aa)){
  aa$per[i]=aa$Freq[i]/sum(aa$Freq[which(aa$Var1==aa$Var1[i])])
}
library(ggalluvial)
library(reshape)
aa$order=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
aa$Var2=factor(aa$Var2,levels=c('B cell','Plasma B','CD4+ T','CD8+ T','Cycling T','DC','Macrophage','Mast cell','Monocyte','NK'))
pdf("")


ggplot(aa,aes(x = Var1, stratum =Var2, 
        alluvium = order,y = per,
        fill = Var2, label = Var2)) +
  scale_x_discrete(expand = c(.1, .1)) +
  scale_fill_manual(values=c(c('B cell'="#673997",'CD4+ T'="#048D5E",'CD8+ T'="#ADDA81",'Cycling T'="#FD7906",'DC'="#FF9998",'Macrophage'="#EA2087",'Mast cell'="#F7F393",'Monocyte'="#FFBB6D",'NK'="#9ACDE6",'Plasma B'="#ED201E")))+
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 2.5) +
  theme(panel.border=element_rect(color='black', fill=NA), 
       panel.grid.major =element_blank(), 
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text=element_text(colour="black")) 


##TCGA
ggplot(aa,aes(x=CD83,y=CD27))+
#geom_point()+
stat_density_2d(aes(fill = ..level..), geom="polygon")+
scale_fill_distiller(palette='YlOrRd',direction=1)+
theme_bw()+
theme(panel.grid=element_blank())

ggplot(bb,aes(x=RGS1,y=CD27))+
#geom_point()+
stat_density_2d(aes(fill = ..level..), geom="polygon")+
scale_fill_distiller(palette='YlOrRd',direction=1)+
theme_bw()+
theme(panel.grid=element_blank())

##
res_scDC_noClust1 <- scDC::scDC_noClustering(as.character(Bcell$type), as.character(Bcell$stage), calCI = TRUE, 
                                     calCI_method = c("BCa"),
                                     nboot = 1000,ncores=10)
res1=res_scDC_noClust1$result
res1$mean=apply(res_scDC_noClust1$thetastar,1,mean)

res1[,c('cellTypes','subject','mean')]->res1
res1<-data.frame(type=as.character(res1$cellTypes[1:33]),N_mean=res1$mean[1:33],T_mean=res1$mean[34:66])



per2<-data.frame(compare=c(rep('INVvsNormal',5),rep('MIAvsNormal',5),rep('INVvsMIA',5)))
per2$celltype=rep(res1$cellTypes[1:5],3)
per2$fc=1
per2$fc[1:5]=res1$mean[which(res1$subject=='INV')]/res1$mean[which(res1$subject=='adjacent')]
per2$fc[6:10]=res1$mean[which(res1$subject=='MIA')]/res1$mean[which(res1$subject=='adjacent')]
per2$fc[11:15]=res1$mean[which(res1$subject=='INV')]/res1$mean[which(res1$subject=='MIA')]
per2$logfc=log2(per2$fc)
##
for(i in 1:nrow(per2)){
    if(per2$logfc[i]>0.25){
        per2$reg[i]='UP'
    }else if(per2$logfc[i]<(-0.25)){
         per2$reg[i]='DOWN'
    }else{
        per2$reg[i]='NAN'
    }
}
pdf("./Project/LUAD/figs/scdc_compare_tissue_stage.pdf",useDingbats=F,width=3.5,height=5)
ggplot(per2,aes(x=compare,y=celltype))+
geom_point(aes(size=abs(logfc),color=reg))+
scale_color_manual(values=c('UP'='#ED2A28','DOWN'='#2C9685','NAN'='#B5B6B8'))+
theme_bw()+
theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL,use.overlap=F)
res1=res$bulk.props[c('CD27+ Memory B','CD83+ Activated B','TCL1A+ Naive B','GC B'),]
data.frame(t(res1))->res1
data.frame(prop=c(res1$CD27..Memory.B,res1$CD83..Activated.B,res1$TCL1A..Naive.B,res1$GC.B),type=c(rep('CD27 Memory B',278),rep('CD83 Memory B',278),rep('TCL1A..Naive.B',278),rep('GC B',278)))->type
library(plyr)
mu <- ddply(prop, "group", summarise, grp.mean=mean(prop))
ggplot(prop, aes(x=prop, color=group1)) +
  geom_density()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
   theme_classic()

for(i in 1:nrow(prop1)){
  if(prop1$prop[i]<quantile(prop1$prop,0.5)){
    prop1$group[i]='low'
}else{
  prop1$group[i]='high'
}}
fit<-survfit(Surv(PFS.time/30,PFS)~group,data=prop1)
pdf('./')
ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'),legend='none')
####
library(splatter)
library(Seurat)
library(msigdbr)
library(singleseqgset)
library(heatmap3)
h.human <- msigdbr(species="Homo sapiens",category="C7",subcategory='IMMUNESIGDB')



p1<-ggplot(aa,aes(x=mem,y=plasma,color=type))+geom_point(size=4)+
scale_color_manual(values=c('Naive B'='#938A81','Memory B'='#FAD09F','Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED','Plasma B'='#FF913F'))

p2<-gg.gap(plot = p1,
           segments = c(1.03, 1.09),
           tick_width = 0.01,
           rel_heights = c(0.25, 0, 0.1),
           ylim = c(0.97, 1.2)
           )
##
ggplot(up, aes(x=Term, y=log10P)) +
  geom_bar(stat='identity', aes(fill=type), width=.7) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  xlab("") +
  ylab("-log10(pvalue)")+
  coord_flip()

#
gene=c('CD8A','CD4','CXCL13','GZMK','GZMB','PDCD1','CTLA4','HAVCR2','LAG3','TIGIT','ITGAE','ICOS','TOX','TOX2','BHLHE40')
DotPlot(Tfr,features=gene[15:1],group.by='type')+
scale_color_distiller(palette='Reds',direction=1)+
coord_flip()
#
sc_info$type=TILC$type
p1<-ggplot(sc_info,aes(x=sc_x,y=sc_y,color=type))+
geom_point(size=0.2,alpha=0.6)+NoLegend()+
scale_color_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1','KIT+ ILC3'='#452C2C','ICOS+ ILC2'='#B4A8BD'))
#
mc_info$type=mc@annots[as.numeric(rownames(mc_info))]
p2<-ggplot(mc_info,aes(x=mc_x,y=mc_y))+
geom_point(aes(fill=type),color='black',pch=21,size=2)+NoLegend()+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1','KIT+ ILC3'='#452C2C','ICOS+ ILC2'='#B4A8BD'))+
scale_y_continuous(limits=c(0,1200), breaks=seq(0,1200,300))
##
pdf("./figs/TILC.pdf",useDingbats=F)
p<-ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,color=type),size=0.2,alpha=0.4)+NoLegend()+
scale_color_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1','KIT+ ILC3'='#452C2C','ICOS+ ILC2'='#B4A8BD'))+
geom_point(data=mc_info,aes(x=mc_x,y=mc_y,fill=type),color='black',pch=21,size=2)+NoLegend()+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1','KIT+ ILC3'='#452C2C','ICOS+ ILC2'='#B4A8BD'))
dev.off()

res_scDC_noClust1 <- scDC::scDC_noClustering(as.character(TILC$type), as.character(TILC$stage), calCI = TRUE, 
                                     calCI_method = c("BCa"),
                                     nboot = 1000,ncores=20)
res1=res_scDC_noClust1$result
res1$mean=apply(res_scDC_noClust1$thetastar,1,mean)

res1[,c('cellTypes','subject','mean')]->res1
res1<-data.frame(row.names=as.character(res1$cellTypes[1:16]),Normal=res1$mean[1:16], AIS=res1$mean[17:32],INV=res1$mean[33:48],MIA=res1$mean[49:64])


pheatmap(res1,scale='row',border_color=NA,cellwidth=15,cellheight=10,fontsize=10,color=colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu")[11:0])(100),cluster_cols=F)
####TCGA eLUAD
for(i in 1:nrow(group)){
  if(group$exp[i]<quantile(group$exp,0.5)){
    group$group[i]='low'
  }else{
    group$group[i]='high'
  }
}
library(limma)
TCGA$exp->exp
some_cluster<-group
some_cluster$group<-as.character(some_cluster$group)
some_cluster$group[which(some_cluster$group!= "high")]<-"control"
some_cluster$group[which(some_cluster$group== "high")]<-"case"
grouP<-as.factor(some_cluster$group)
desigN <- model.matrix(~ grouP + 0)
rownames(desigN)<-colnames(exp)
comparE <- makeContrasts(grouPcase-grouPcontrol,levels=desigN)
fiT <- lmFit(exp, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
Diff<-topTable(fiT3,p.value=1,num=Inf)

Diff$sig=''
Diff[rownames(Diff)[which(abs(Diff$logFC)>log2(1.5) & Diff$adj.P.Val<0.05)],]$sig=rownames(Diff)[which(abs(Diff$logFC)>log2(1.5) & Diff$adj.P.Val<0.05)]
pdf("./youwh_work/Project_eLUAD/LUAD/figs/eLUAD_TCGA_CD4CXCL13+vs-.pdf",useDingbats=F)
ggplot(Diff, aes(logFC, -log10(adj.P.Val), label = sig)) +
  geom_point(color = ifelse(Diff$sig == "", "grey", "red")) +
  geom_text_repel(data = Diff[Diff$sig != "",], col="blue",max.overlaps=35) +
  geom_vline(xintercept = log2(1.5) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(1.5), linetype="dashed", color="grey77") + ylab('(-log10 q value)') + xlab('log2 f')
dev.off()



p<-ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,color=type),size=0.5)+NoLegend()+
scale_color_manual(values=c('TCL1A+ Naive B'='#938A81','CD27+ Memory B'='#FAD09F','CD83+ Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED'))+
geom_point(data=mc_info,aes(x=mc_x,y=mc_y,fill=type),color='black',pch=21,size=3)+NoLegend()+
scale_fill_manual(values=c('TCL1A+ Naive B'='#938A81','CD27+ Memory B'='#FAD09F','CD83+ Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED'))
###
RR=c('CD4+ Tfh|CD83+ Activated B','CD4+ Tfh|GC B','CD4+ Tfh|FCRL4+ B','CD4+ Tfh|Plasma B','CD4+ Tfh|TCL1A+ Naive B','CD4+ Tfh|CD27+ Memory B',
'CD4+ Treg|CD83+ Activated B','CD4+ Treg|GC B','CD4+ Treg|FCRL4+ B','CD4+ Treg|Plasma B','CD4+ Treg|TCL1A+ Naive B','CD4+ Treg|CD27+ Memory B',
'CD8+ CTL|CD83+ Activated B','CD8+ CTL|GC B','CD8+ CTL|FCRL4+ B','CD8+ CTL|Plasma B','CD8+ CTL|TCL1A+ Naive B','CD8+ CTL|CD27+ Memory B',
'CXCL13+CD8+ Dysfunctional Trm|CD83+ Activated B','CXCL13+CD8+ Dysfunctional Trm|GC B','CXCL13+CD8+ Dysfunctional Trm|FCRL4+ B','CXCL13+CD8+ Dysfunctional Trm|Plasma B','CXCL13+CD8+ Dysfunctional Trm|TCL1A+ Naive B','CXCL13+CD8+ Dysfunctional Trm|CD27+ Memory B',
'CD8+ Tem|CD83+ Activated B','CD8+ Tem|GC B','CD8+ Tem|FCRL4+ B','CD8+ Tem|Plasma B','CD8+ Tem|TCL1A+ Naive B','CD8+ Tem|CD27+ Memory B'
)
##
pheatmap(cor.m,border_color='NA',breaks=unique(c(seq(-0.6,0.6, length=100))),color=colorRampPalette(c("#61C1F0", "white", "#EC3C32"))(100))


res.cut <- surv_cutpoint(OS, time = "OS.time", event = "OS",
   variables = c("score"))
  res.cat <- surv_categorize(res.cut)
  ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'),legend='none')
##T &ILC gene list
Tcell=subset(LUAD,type_global %in% c('CD4+ T','CD8+ T'))
T_ligand_receptor=c('CXCL13','IL21','XCL1','IL6ST','CD40LG','BTLA','ICAM1','TNFSF8','TNFSF4','TNFSF13B','CTLA4','LTB','ICOS','ITGB1','TNFRSF1B','CXCR3','CXCR6','CCR4')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- ScaleData(Tcell, features = rownames(Tcell))
AverageExpression(Tcell,features=T_ligand_receptor,group.by='type')->aa
pheatmap::pheatmap(aa$RNA[,c('CD4+ Tfh','CXCL13+CD8+ Dysfunctional Trm','CD8+ CTL','CXCL13-CD8+ Dysfunctional Trm','CD8+ Tem','CD4+ Treg','Naive-like')],scale='row',cluster_rows=F,cluster_cols=F,breaks=unique(c(seq(-3,3, length=100))),border=NA,color=colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(100),cellwidth=15,cellheight=13)
###
Bcell=subset(LUAD,type_global %in% c('B cell','Plasma B'))
Bcell <- NormalizeData(Bcell, normalization.method = "LogNormalize", scale.factor = 10000)
Bcell <- ScaleData(Bcell, features = rownames(Bcell))
B_receptor_ligand=c('CXCR5','IL21R','XCR1','EBI3','CD40','TNFRSF14','AREG','TNFRSF8','TNFRSF4','TNFRSF13B','CD86','LTBR','ICOSLG','VCAM1','LTA','CXCL9','CXCL16','CCL17')
AverageExpression(Bcell,features=B_receptor_ligand,group.by='type')->bb
pheatmap::pheatmap(bb$RNA,scale='row',cluster_rows=F,cluster_cols=F,breaks=unique(c(seq(-3,3, length=100))),border=NA,color=colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(100),cellwidth=15,cellheight=13)
###

pdf("./figs/eTCGA_TLS_CXCL13CD4.pdf")
 ggplot(score,aes(x=Tfh,y=TLS))+
 geom_point(aes(color=color,alpha=0.7),size=3)+
 scale_color_manual(values='#90D3F2')+
 geom_smooth(method = 'lm', se = F, color = 'black')+NoLegend()
dev.off()


pdf("./figs/eTCGA_CD8CXCL13_CXCL13CD4.pdf")
 ggplot(score,aes(x=Tfh,y=CD8))+
 geom_point(aes(color=color,alpha=0.7),size=3)+
 scale_color_manual(values='#90D3F2')+
 geom_smooth(method = 'lm', se = F, color = 'black')+NoLegend()
dev.off()
##Myeloid
subset(LUAD,type_global %in% c('DC','Macrophage','Monocyte'))->Myeloid
 mc=scdb_mc("mc_annot")
 mc2d=scdb_mc2d("LUAD_2dproj")

  sc_info=data.frame(sc_x=mc2d@sc_x[colnames(LUAD)],sc_y=mc2d@sc_y[colnames(LUAD)],row.names=names(mc2d@sc_y[colnames(LUAD)]))
sc_info$type=Myeloid$type

 sc_info=data.frame(sc_x=mc2d@sc_x[colnames(Myeloid)],sc_y=mc2d@sc_y[colnames(Myeloid)],row.names=names(mc2d@sc_y[colnames(Myeloid)]))
sc_info$type=Myeloid$type

ggplot(data=sc_info,aes(x=sc_x,y=sc_y,color=type))+
geom_point(size=3)+
scale_color_manual(values=c("APOC1+ Macrophage"='#C2FF99',"cDC2"='#00FECF','cDC1'='#575329','mDC'='#B05B6F','pDC'='#8CD0FF','IL1B+ Macrophage'='#00A6AA','CD14+ Monocyte'='#788D66','CD16+ Monocyte'='#885578'))+
NoLegend()


ggplot(data=sc_info[which(sc_info$stage=='MIA'),],aes(x=sc_x,y=sc_y,color=type))+
geom_point(size=3)+
scale_color_manual(values=c("APOC1+ Macrophage"='#C2FF99',"cDC2"='#00FECF','cDC1'='#575329','mDC'='#B05B6F','pDC'='#8CD0FF','IL1B+ Macrophage'='#00A6AA','CD14+ Monocyte'='#788D66','CD16+ Monocyte'='#885578'))+
NoLegend()


minnumbcells <- min(sapply(c(as.character(unique(sc_info$stage))), function(x) length(rownames(sc_info[sc_info$stage==x,]))))
set.seed(22)
plotcells <-  as.vector(sapply(c(as.character(unique(sc_info$stage))), function(x) sample(rownames(sc_info[sc_info$stage==x,]), size=minnumbcells, replace=F, prob=NULL )))


plotcells <-  as.vector(sapply(c('tLung','tL/B'), function(x) sample(rownames(sc_info[sc_info$origin==x,]), size=minnumbcells, replace=F, prob=NULL )))

pdf('./figs/density_2d.pdf',useDingbats=F)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='adjacent'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='adjacent'),],aes(x=sc_x,y=sc_y),color="black")


ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='AIS'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='AIS'),],aes(x=sc_x,y=sc_y),color="black")

ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='MIA'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='MIA'),],aes(x=sc_x,y=sc_y),color="black")

ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='INV'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$stage=='INV'),],aes(x=sc_x,y=sc_y),color="black")

dev.off()




geom_point(data=fat,aes(x=x,y=y,color=type),size=1,color='lightgrey')+
geom_point(data=fat[which(fat$tissue=='IA_N'),],aes(x=x,y=y,color=type),size=1)+scale_color_manual(values=c('cDC1'="#5A0007",'cDC2'="#809693",'mDC'="#6A3A4C",'pDC'="#1B4400",'Kupffer cell'="#A30059",'Macrophage'="#FFDBE5",'Monocyte'="#0000A6"))+
geom_density_2d(data=fat[which(fat$tissue=='IA_N'),],aes(x=x,y=y),colour="black")+Seurat::NoLegend()




####

res_scDC_noClust1 <- scDC::scDC_noClustering(as.character(Myeloid$type), as.character(Myeloid$stage), calCI = TRUE, 
                                     calCI_method = c("BCa"),
                                     nboot = 1000,ncores=10)
res1=res_scDC_noClust1$result
res1$mean=apply(res_scDC_noClust1$thetastar,1,mean)

res1[,c('cellTypes','subject','mean')]->res1
res1<-data.frame(type=as.character(res1$cellTypes[1:8]),normal=res1$mean[1:8],AIS=res1$mean[9:16],INV=res1$mean[17:24],MIA=res1$mean[25:32])

per2<-data.frame(compare=c(rep('INVvsNormal',5),rep('MIAvsNormal',5),rep('INVvsMIA',5)))
per2$celltype=rep(res1$cellTypes[1:5],3)
per2$fc=1
per2$fc[1:5]=res1$mean[which(res1$subject=='INV')]/res1$mean[which(res1$subject=='adjacent')]
per2$fc[6:10]=res1$mean[which(res1$subject=='MIA')]/res1$mean[which(res1$subject=='adjacent')]
per2$fc[11:15]=res1$mean[which(res1$subject=='INV')]/res1$mean[which(res1$subject=='MIA')]
per2$logfc=log2(per2$fc)
##
for(i in 1:nrow(per2)){
    if(per2$logfc[i]>0.25){
        per2$reg[i]='UP'
    }else if(per2$logfc[i]<(-0.25)){
         per2$reg[i]='DOWN'
    }else{
        per2$reg[i]='NAN'
    }
}
pdf("./Project/LUAD/figs/scdc_compare_tissue_stage.pdf",useDingbats=F,width=3.5,height=5)
ggplot(per2,aes(x=compare,y=celltype))+
geom_point(aes(size=abs(logfc),color=reg))+
scale_color_manual(values=c('UP'='#ED2A28','DOWN'='#2C9685','NAN'='#B5B6B8'))+
theme_bw()+
theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


##
gene_tilc<-c('CD4','CD8A','CXCR6','ITGAE','ZNF683','CXCL13','LAG3','FOXP3','CTLA4','TCF7','CCR7','LEF1','GZMB','IFNG','PRF1','GZMK','CD44','EOMES','SLC4A10','KLRB1','ISG15','MX1','IFIT3','MKI67','TOP2A','FCGR3A','CX3CR1','GNLY','XCL1','IL7R','ICOS','KIT','TNFRSF18')

c('CD4+ Treg','CD4+ Tfh','CD4+ Tem','Naive-like','Cycling T','CD8+ISG15+ T','CD8+ CTL','CD8+ Tem','CD8+ MAIT','CXCL13+CD8+ Dysfunctional Trm','CXCL13-CD8+ Dysfunctional Trm','GNLY+ cNK','GNLY- cNK','resident-like NK','ICOS+ ILC2','KIT+ ILC3')


ggplot(bb,aes(x=type,y=CXCL13,fill=type))+geom_boxplot(alpha=0.9)+scale_fill_manual(values=c('CXCL13+CD8+ Dysfunctional Trm'='#D157A0','CD4+ Tfh'='#6F0062'))+
NoLegend()

ggplot(bb,aes(x=type,y=GAPDH,fill=type))+geom_boxplot(alpha=0.9)+scale_fill_manual(values=c('CXCL13+CD8+ Dysfunctional Trm'='#D157A0','CD4+ Tfh'='#6F0062'))+
NoLegend()
####
Myeloid_gene<-c('APOC1','APOE','MARCO','CD14','FCGR3A','IL1B','C1QA','C1QB','S100A9','S100A8','XCR1','CLEC9A','CD1A','CD1C','CLEC10A','CD207','CCL17','HLA-DQA1','HLA-DQB1','CCR7','LAMP3','FSCN1','CCL12','CXCL9','CXCL10','GZMB','JCHAIN')

Myeloid_gene<-c('C1QA','C1QB','S100A9','S100A8','FCN1','CD14','FCGR3A','APOC1','APOE','IL1B','TREM2','MARCO','CCR7','LAMP3','XCR1','CLEC9A','CLEC10A','CD207','CD1C','CD1A','CXCL10','CXCL9','CCL17','GZMB','FSCN1','JCHAIN','HLA-DQA1','HLA-DQB1')
###
adata <- sc$AnnData(
  X   = Matrix::t((GetAssayData(LUAD,slot='counts'))), #scVI requires raw counts
  obs = LUAD[[]],
  var = GetAssay(LUAD)[[]]
)
scvi$model$SCVI$setup_anndata(adata,batch_key="orig.ident")
model = scvi$model$SCVI(adata)
model$train()
latent = model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(Bcell)
Bcell[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(Bcell))

Bcell <- FindNeighbors(Bcell, dims = 1:10, reduction = "scvi")
Bcell <- FindClusters(Bcell, resolution =0.5)

Bcell <- RunUMAP(Bcell, dims = 1:10, reduction = "scvi", n.components = 2)

###
mc_info=data.frame(mc_x=mc2d@mc_x,mc_y=mc2d@mc_y)
sc_info=data.frame(sc_x=mc2d@sc_x[colnames(LUAD)],sc_y=mc2d@sc_y[colnames(LUAD)],row.names=names(mc2d@sc_x[colnames(LUAD)]))
mc_info$type=mc@annots
sc_info$type_global=LUAD$type_global
p<-ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type_global),color='black',pch=21,size=1.8)+
scale_fill_manual(values=c('B cell'="#673997",'Basophil'="#B35A20",'CD4+ T'="#048D5E",'CD8+ T'="#ADDA81",'Cycling T'="#FD7906",'DC'="#FF9998",'ILC'="#1775B6",'Macrophage'="#EA2087",'Mast cell'="#F7F393",'Monocyte'="#FFBB6D",'NK'="#9ACDE6",'Plasma B'="#ED201E"))+NoLegend()



adata <- sc$AnnData(
  X   = Matrix::t((GetAssayData(Bcell,slot='counts'))), #scVI requires raw counts
  obs = Bcell[[]],
  var = GetAssay(Bcell)[[]]
)


data.frame(celltype=c(rep('B cell',2),rep("CD4+ T",2),rep("CD8+ T",2),rep("NK",2)),DOWN=c(67,73,25,19,54,6,4,3),UP=c(69,96,29,18,74,16,32,4),type=rep(c('IAC','MIA'),4))->count_deg
count_deg$DOWN=(-count_deg$DOWN)

count<-data.frame(celltype=rep(count_deg$celltype,2),stage=rep(count_deg$type,2),count=c(count_deg$UP,count_deg$DOWN))
ggplot(count,aes(x=celltype,y=count,fill=stage))+
geom_bar(stat = 'identity', position ='dodge')


ggplot()+
geom_point(data=sc_info,aes(x=metacell_1,y=metacell_2,fill=type),color='black',pch=21,size=1.8)+
scale_fill_manual(values=c('CD27+ Memory B'='#FAD09F','CD83+ Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED','Plasma B'='#FF913F','TCL1A+ Naive B'='#938A81'))+
NoLegend()

sc_info$PID=LUAD@meta.data[rownames(sc_info),'PatientID']
sc_info$stage=LUAD@meta.data[rownames(sc_info),'stage']
data.frame(table(sc_info$type,sc_info$stage))->Bcell_count
for(i in 1:nrow(Bcell_count)){
  Bcell_count$percentage[i]=Bcell_count$Freq[i]/sum(Bcell_count$Freq[which(Bcell_count$Var2==Bcell_count$Var2[i])])
}
pdf("percentage_stage.pdf",width=1.9,height=1.7)
ggplot(Bcell_count,aes(x=Var2,y=percentage,fill=Var1))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c('CD27+ Memory B'='#FAD09F','CD83+ Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED','Plasma B'='#FF913F','TCL1A+ Naive B'='#938A81'))+
NoLegend()
dev.off()

pdf("./Tcell_metacell.pdf",useDingbats=F)
ggplot()+
geom_point(data=sc_info,aes(x=metacell_1,y=metacell_2,fill=type),color='black',pch=21,size=1.8)+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1'))+
NoLegend()
dev.off()

data.frame(table(Tcell$type,Tcell$stage))->T_count
for(i in 1:nrow(T_count)){
  T_count$percentage[i]=T_count$Freq[i]/sum(T_count$Freq[which(T_count$Var2==T_count$Var2[i])])
}
T_count$Var2<-factor(T_count$Var2,levels=c("adjacent",'AIS','MIA','INV'))
pdf("percentage_stage_tcell.pdf",width=1.9,height=1.7)
ggplot(T_count,aes(x=Var2,y=percentage,fill=Var1))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1'))+NoLegend()
dev.off()


cluster_markers<-list()
deg<-deg[which(deg$p_val_adj<0.05 & deg$avg_log2FC>0.3),]
deg[substr(deg$gene,1,3) !='RPL',]->deg
deg[substr(deg$gene,1,3) !='RPS',]->deg

deg[substr(deg$gene,1,2) !='AP',]->deg
deg[substr(deg$gene,1,4) !='SNOR',]->deg

for(i in unique(deg$cluster)){
  cluster_i <- deg %>% filter(cluster == i)
  cluster_markers[[paste("C",i,sep='')]] <- as.numeric(bitr(cluster_i$gene, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID)
}
res<-compareCluster(geneCluster = cluster_markers, fun = "enrichGO",ont='BP',OrgDb='org.Hs.eg.db')
res@compareClusterResult->cres
data.frame(cres %>% group_by(Cluster) %>% top_n(n = 5, wt =-pvalue))->top5
cres[which(cres$Description %in% top5$Description),]->top5
 top5[,c('Cluster','Description','pvalue','qvalue','Count')]->top5
 top5$log10P=-log10(top5$pvalue)

library(tidyr)
matrix(nrow=39,ncol=10)->top5_matrix
rownames(top5_matrix)=unique(top5$Description)
colnames(top5_matrix)=unique(deg$cluster)
for(i in rownames(top5_matrix)){
  for(j in unique(deg$cluster)){
    if(nrow(top5[which(top5$Description==i &top5$Cluster==paste0('C',j)),])==1){
      top5_matrix[i,j]=top5[which(top5$Description==i &top5$Cluster==paste0('C',j)),]$log10P
    }else{
    top5_matrix[i,j]=0
    }
  }
}

bk <- c(seq(0,10,by=0.01))
pdf("heatmap_go_tcell.pdf",useDingbats=F,width=20,height=20)
pheatmap(top5_matrix,cellheight=15,breaks=bk,cellwidth=15,border_color='lightgrey',cluster_rows=F,cluster_cols=F,color =colorRampPalette(c('white','#FDEEEF','#FAB8BA','#F69EA1','#EF7377','#EC4F56','#E52733','#E2001D'))(1000))
dev.off()
###
Tcell$id=paste0(Tcell$stage,"_",Tcell$PatientID)
data.frame(table(Tcell$type,Tcell$id))->tt
for(i in 1:nrow(tt)){
  tt$percentage[i]=tt$Freq[i]/sum(tt$Freq[which(tt$Var2==tt$Var2[i])])
}
library(stringr)
tt$stage=substr(tt$Var2,1,str_locate(as.character(tt$Var2),"_")[,1]-1)
tt$stage=factor(tt$stage,levels=c("adjacent","AIS","MIA","INV"))


p<-list()
for(i in unique(tt$type)){
p[[i]]<-ggplot(tt[which(tt$type==i),],aes(x=stage,y=percentage))+
geom_point(aes(fill=stage),color='black',pch=21,size=3,position='jitter')+
scale_fill_manual(values=c("#C5E4CB","#00AB1E","#FF9D9F","#FF0000"))+
stat_boxplot(geom = "errorbar",position = position_dodge(0.5))+
stat_compare_means(comparisons=list(c('adjacent','AIS'),c("adjacent",'INV'),c("adjacent",'MIA'),c('MIA','INV')))

}
pdf("./percentage_type_tcell.pdf",width=20/6*9,height=3)
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]], nrow = 1)
dev.off()
pdf("./percentage_type_bcell1.pdf",width=20/6*6,height=3)
cowplot::plot_grid(p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[22]], nrow = 1)
dev.off()



tt$id=factor(tt$id,levels=c(levels(tt$id)[1:38],levels(tt$id)[72:93],levels(tt$id)[39:71]))
pdf("percentage_pid_tcell.pdf",width=8,height=5)
ggplot(tt,aes(x=id,y=percentage,fill=type))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1'))+NoLegend()
dev.off()

pdf("./TNK_metacell.pdf",useDingbats=F)
ggplot()+
geom_point(data=sc_info,aes(x=metacell_1,y=metacell_2,fill=type),color='black',pch=21,size=1.8)+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1','Cycling T'='#456648'))+
theme_classic()+
NoLegend()
dev.off()


###T&NK markergene
data.frame(table(NKT$type,NKT$stage))->NKT_count

for(i in 1:nrow(NKT_count)){
  NKT_count$percentage[i]=NKT_count$Freq[i]/sum(NKT_count$Freq[which(NKT_count$Var2==NKT_count$Var2[i])])
}
NKT_count$Var2<-factor(NKT_count$Var2,levels=c("adjacent",'AIS','MIA','INV'))
pdf("percentage_stage_t_NKcell.pdf",width=1.9,height=1.7)
ggplot(NKT_count,aes(x=Var2,y=percentage,fill=Var1))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1'))+
theme_classic()+
NoLegend()
dev.off()


NKT$id=paste0(NKT$stage,"_",NKT$PatientID)

data.frame(table(NKT$type,NKT$id))->NKT_count1
for(i in 1:nrow(NKT_count1)){
  NKT_count1$percentage[i]=NKT_count1$Freq[i]/sum(NKT_count1$Freq[which(NKT_count1$Var2==NKT_count1$Var2[i])])
}
library(stringr)
NKT_count1$stage=substr(NKT_count1$Var2,1,str_locate(as.character(NKT_count1$Var2),"_")[,1]-1)
NKT_count1$stage=factor(NKT_count1$stage,levels=c("adjacent","AIS","MIA","INV"))
NKT_count1$id=factor(NKT_count1$id,levels=c(levels(NKT_count1$id)[1:38],levels(NKT_count1$id)[72:93],levels(NKT_count1$id)[39:71]))
pdf("percentage_stage_t_NKcell_patient.pdf",width=10,height=5)
ggplot(NKT_count1,aes(x=id,y=percentage,fill=type))+geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c("CXCL13-CD8+ Dysfunctional Trm"='#BEC459', "Naive-like"='#A3C8C9','CD8+ CTL'='#EEC3FF',"resident-like NK"='#3B9700','CD4+ Treg'='#0CBD66',
'CD4+ Tem'='#00489C','GNLY+ cNK'='#886F4C','GNLY- cNK'='#34362D','CD4+ Tfh'='#6F0062','CXCL13+CD8+ Dysfunctional Trm'='#D157A0','Cycling T'='#456648','CD8+ Tem'='#B77B68',
'CD8+ MAIT'='#456D75','CD8+ISG15+ T'='#7A87A1'))+NoLegend()+theme_classic()
dev.off()


####


#### cell cycle score
mat=scdb_mat("LUAD")
mc=scdb_mc("mc_annot")
mc2d=scdb_mc2d("LUAD_2dproj")

height = tgconfig::get_param("mcell_mc2d_gene_height", "metacell")
width = tgconfig::get_param("mcell_mc2d_gene_width", "metacell")
mc_cex = tgconfig::get_param("mcell_mc2d_gene_mc_cex", "metacell")
sc_cex = tgconfig::get_param("mcell_mc2d_gene_cell_cex", "metacell")

cells = colnames(NKT)
cc_genes =  cc.genes$g2m.genes
diff_genes = setdiff(cc_genes, mat@genes)

cc_genes = setdiff(cc_genes, diff_genes)

tot_cc_c = Matrix::colSums(mat@mat[cc_genes, cells])
f_cc_c =  tot_cc_c / Matrix::colSums(mat@mat[, cells])
f_cc_cutoff=0.01
tot_cc = tapply(f_cc_c > f_cc_cutoff, mc@mc[cells], sum)
tara = tapply(f_cc_c, mc@mc[cells], length) - tot_cc
df = data.frame(row.names=names(tot_cc), tot_cc=tot_cc, tara=tara, f_cc=tot_cc/(tot_cc+tara))
df = df[order(df$tot_cc / (df$tot_cc + df$tara)), ]
#plot_mode == 'circle_size'

plot(mc2d@sc_x[cells], mc2d@sc_y[cells], pch=19, cex=sc_cex, col='lightgrey', xlab="", ylab="")
mcs = as.numeric(rownames(df))
points(mc2d@mc_x[mcs], mc2d@mc_y[mcs], pch=21, bg=mc@colors[mcs], cex=mc_cex * log10(1 + 100*df$f_cc)/2)
fs = seq(0.05, length=1 + ceiling(max(df$f_cc)/0.15), by=0.15)
legend("topright", legend=paste0(fs * 100, "%"), pch=21, pt.bg='darkgrey', col='black', pt.cex=mc_cex * log10(1 + 100*fs)/2, bty='n', y.intersp=2, inset=0.04)

##
ord_by_id<-list()
ord_by_id[["mc_annot"]]=as.numeric(order(mc@annots))

df_cc = data.frame(row.names=names(f_cc_c), f_cc=f_cc_c, tot_cc=tot_cc_c, stringsAsFactors = F)
df_cc$group=NKT@meta.data[rownames(df_cc),'type']
f_cc_by_grp = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$group, mean)
tot_cc_by_grp = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$group, length)
f_cc_by_grp = f_cc_by_grp
tot_cc_by_grp = tot_cc_by_grp[names(f_cc_by_grp)]
pdf("barplot_profiler.pdf",useDingbats=F)
par(mar=c(8,4,4,1))
barplot(f_cc_by_grp, col=mc@colors[names(f_cc_by_grp)], las=2, ylab='% prolif cells')
mtext(tot_cc_by_grp, 3, at=seq(0.6, by=1.2, length=length(tot_cc_by_grp)), las=2, line=0.5)
dev.off()


###

 sc_info=data.frame(sc_x=mc2d@sc_x,sc_y=mc2d@sc_y,row.names=names(mc2d@sc_y))

ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=1.8)+
scale_fill_manual(values=c("CD8-TRM"='#BEC459', "CD4-Naive"='#A3C8C9','CD4-Treg'='#0CBD66',
'CD4-CXCL13'='#6F0062','CD8-CXCL13'='#D157A0','CD4-Th17like'='#456648','CD8-TEM'='#B77B68','CD4-Tem'='#00489C'))+
theme_classic()+
NoLegend()


ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=clon,color=clon),pch=21,size=1)+
scale_fill_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000))+
 scale_color_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000))+      
theme_classic()
##
library(immunarch)
igdata <- repLoad('./vdj_t_info.csv', .mode = "paired")
igdata$bc_cluster=sc_info[vdj_t$barcode,'type']
igdata$bc_cluster=as.character(igdata$bc_cluster)
names(igdata$bc_cluster)=vdj_t$barcode
CD8_barcode=rownames(sc_info)[which(sc_info$type %in% c('CD8-CXCL13','CD8-TEM','CD8-TRM'))]
CD8_barcode=intersect(vdj_t$barcode,CD8_barcode)
vdj_t[which(vdj_t$barcode %in% CD8_barcode),]->CD8_df
write.csv(CD8_df,'./CD8_vdj.csv',quote=F)
CD8_data=repLoad('./CD8_vdj.csv', .mode = "paired")
CD8_data$bc_cluster=sc_info[CD8_barcode,'type']
CD8_data$bc_cluster=sc_info[CD8_barcode,'type']
CD8_data$bc_cluster=as.character(CD8_data$bc_cluster)
names(CD8_data$bc_cluster)=CD8_barcode

scdata_CD8<- select_clusters(CD8_data, CD8_data$bc_cluster, "Cluster")

repOverlap(scdata_CD8$data) %>% vis()
scdata<- select_clusters(igdata, igdata$bc_cluster, "Cluster")
repOverlap(scdata$data) %>% vis()


CD4_barcode=rownames(sc_info)[which(sc_info$type %in% c('CD4-Naive','CD4-CXCL13','CD4-Tem','CD4-Treg','CD4-Th17like'))]
CD4_barcode=intersect(vdj_t$barcode,CD4_barcode)
vdj_t[which(vdj_t$barcode %in% CD4_barcode),]->CD4_df
write.csv(CD4_df,'./CD4_vdj.csv',quote=F)
CD4_data=repLoad('./CD4_vdj.csv', .mode = "paired")
CD4_data$bc_cluster=sc_info[CD4_barcode,'type']
CD4_data$bc_cluster=sc_info[CD4_barcode,'type']
CD4_data$bc_cluster=as.character(CD4_data$bc_cluster)
names(CD4_data$bc_cluster)=CD4_barcode

scdata_CD4<- select_clusters(CD4_data, CD4_data$bc_cluster, "Cluster")
repOverlap(scdata_CD4$data) %>% vis()
pdf("percentage_10X.pdf",width=4,height=5)
ggplot(aa,aes(x=Var1,y=percentage,fill=Var2))+
geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c("CD8-TRM"='#BEC459', "CD4-Naive"='#A3C8C9','CD4-Treg'='#0CBD66',
'CD4-CXCL13'='#6F0062','CD8-CXCL13'='#D157A0','CD4-Th17like'='#456648','CD8-TEM'='#B77B68','CD4-Tem'='#00489C'))+
theme_classic()
dev.off()

mean_CD4=data.frame(row.names=c("CD4-CXCL13",'CD4-Treg','CD4-Naive','CD4-Tem','CD4-Th17like'),luo=c(10.14,1.77,1.84,2.51,1.64),zhang=c(1.5,1.78,1.01,1,1.45))
mean_CD4a=data.frame(meansize=c(mean_CD4$luo,mean_CD4$zhang),pid=c(rep('luo',5),rep('zhang',5)),type=rep(rownames(mean_CD4),2))

ggplot(mean_CD4a,aes(x=pid,y=meansize,fill=type))+
geom_bar(stat='identity',position = 'stack')+
scale_fill_manual(values=c( "CD4-Naive"='#A3C8C9','CD4-Treg'='#0CBD66','CD4-CXCL13'='#6F0062','CD4-Th17like'='#456648','CD4-Tem'='#00489C'))+
theme_classic()

mean_CD8=data.frame(row.names=c("CD8-CXCL13",'CD8-TRM','CD8-TEM'),luo=c(64.71,15.76,10.71),zhang=c(0,3.61,2.82))
mean_CD8a=data.frame(meansize=c(mean_CD8$luo,mean_CD8$zhang),pid=c(rep('luo',3),rep('zhang',3)),type=rep(rownames(mean_CD8),2))
ggplot(mean_CD8a,aes(x=pid,y=meansize,fill=type))+
geom_bar(stat='identity',position = 'stack')+
scale_fill_manual(values=c("CD8-TRM"='#BEC459','CD8-TEM'='#B77B68','CD8-CXCL13'='#D157A0'))+
theme_classic()

sc.pl.dotplot(adata,['PTPRC','C1QA','CD14','CD79A','NCR3','CD3D','CD4','CD8A','NCAM1','CD1C','CLEC9A'],groupby='louvain')

###
vdj_t[which(vdj_t$sample=='P1'),]->P1
P2=vdj_t[which(vdj_t$sample=='P2'),]
P1$sample<-NULL
P2$sample<-NULL
list(P1,P2)->vdj_t1
combined <- combineTCR(vdj_t1, 
           samples = c("P1", "P2"), 
                ID = c("P1",'P2'), cells="T-AB",removeMulti=TRUE)

vizGenes(combined, 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
aa$stage=substr(aa$Var1,1,str_locate(as.character(aa$Var1),"_")[,1]-1)
aa$stage=factor(aa$stage,levels=c("adjacent","AIS","MIA","INV"))

for(i in 1:nrow(aa)){
  aa$percentage[i]=aa$Freq[i]/sum(aa$Freq[which(aa$id==aa$id[i])])
}


p<-list()
for(i in unique(aa$type)){
p[[i]]<-ggplot(aa[which(aa$type==i),],aes(x=stage,y=percentage))+
geom_point(aes(fill=stage),color='black',pch=21,size=3,position='jitter')+
scale_fill_manual(values=c("#C5E4CB","#00AB1E","#FF9D9F","#FF0000"))+
stat_boxplot(geom = "errorbar",position = position_dodge(0.5))+
stat_compare_means(comparisons=list(c('adjacent','AIS'),c("adjacent",'INV'),c("adjacent",'MIA'),c('MIA','INV')))+
theme_classic()+
ylab(i)

}
pdf("./percentage_type_tnk.pdf",width=20/6*13,height=3)
cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[11]],p[[12]],p[[13]],p[[14]], nrow = 1)
dev.off()




##
library(AUCell)
library(GSEABase)
geneSets <- GeneSet(Inhibit, setName="geneSet1")
cells_AUC <- AUCell_run(NKT@assays$RNA@counts, geneSets)
NKT$Inhibit=as.numeric(cells_AUC@assays@data$AUC)
NKreceptor=c('KLRD1','FGFBP2','FCGR3A','S1RP5','KLRC1','KLRC3','KLRB1','KLRC2')
Inhibit=c("HAVCR2",'PDCD1','CTLA4','TIGIT','LAG3')
Cyto=c("PRF1",'GZMB','GZMA','GZMH','NKG7','GNLY')
TRM=c("CA10",'ITGA1','ITGAE','IL2','IL10','CXCR6','CXCL13','KCNK5','RGS1','CRTAM','DUSP6','IL23R','PDCD1')
FeaturePlot(NKT,reduction='metacell',features='inhibit')+
scale_colour_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                                 c('white','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410','#E52C22','#DD1D23','#C20030','#930039','#8C003A','#6F003D','#56033F'))(1000)) 


data.frame(pid=NKT$PatientID,stage=NKT$stage,Cyto=NKT$Cyto)->Cyto
Cyto$id=paste0(Cyto$stage,'_',Cyto$pid)
Cyto$type=NKT$type_global
CD8=Cyto[which(Cyto$type=='CD8+ T'),]
NK=Cyto[which(Cyto$type=='NK'),]

aggregate(CD8$Cyto, by=list(type=CD8$id),mean)->CD8_mean
colnames(CD8_mean)
colnames(CD8_mean)=c('id','score')
CD8_mean$stage=substr(CD8_mean$id,1,str_locate(as.character(CD8_mean$id),"_")[,1]-1)
CD8_mean$stage=factor(CD8_mean$stage,levels=c('adjacent','AIS','MIA','INV'))
CD8_mean$median=c("")
CD8_mean[which(CD8_mean$stage=='adjacent'),]$median=median(CD8_mean[which(CD8_mean$stage=='adjacent'),'score'])
CD8_mean[which(CD8_mean$stage=='AIS'),]$median=median(CD8_mean[which(CD8_mean$stage=='AIS'),'score'])
CD8_mean[which(CD8_mean$stage=='MIA'),]$median=median(CD8_mean[which(CD8_mean$stage=='MIA'),'score'])
CD8_mean[which(CD8_mean$stage=='INV'),]$median=median(CD8_mean[which(CD8_mean$stage=='INV'),'score'])
CD8_sce<-Rmisc::summarySE(CD8_mean,measurevar='score',groupvars='stage')
CD8_sce$type='CD8'
aggregate(NK$Cyto, by=list(type=NK$id),mean)->NK_mean
colnames(NK_mean)
colnames(NK_mean)=c('id','score')
NK_mean$stage=substr(NK_mean$id,1,str_locate(as.character(NK_mean$id),"_")[,1]-1)
NK_mean$stage=factor(NK_mean$stage,levels=c('adjacent','AIS','MIA','INV'))
NK_mean$median=c("")
NK_mean[which(NK_mean$stage=='adjacent'),]$median=median(NK_mean[which(NK_mean$stage=='adjacent'),'score'])
NK_mean[which(NK_mean$stage=='AIS'),]$median=median(NK_mean[which(NK_mean$stage=='AIS'),'score'])
NK_mean[which(NK_mean$stage=='MIA'),]$median=median(NK_mean[which(NK_mean$stage=='MIA'),'score'])
NK_mean[which(NK_mean$stage=='INV'),]$median=median(NK_mean[which(NK_mean$stage=='INV'),'score'])
NK_sce<-Rmisc::summarySE(NK_mean,measurevar='score',groupvars='stage')
NK_sce$type='NK'
sce=rbind(CD8_sce,NK_sce)


ggplot(data=sce,aes(x=stage,y=score,color='type'))+
geom_errorbar(aes(ymin=score-se,ymax=score+se),width=0.1,position=position_dodge(0.1))+
scale_color_manual(values=c("CD8"='#ADDA81','NK'="#9ACDE6"))+
geom_line()+
geom_point()



###
data.frame(pid=NKT$PatientID,stage=NKT$stage,Inhibit=NKT$Inhibit)->Inhibit
Inhibit$id=paste0(Inhibit$stage,'_',Inhibit$pid)
Inhibit$type=NKT$type_global
CD8=Inhibit[which(Inhibit$type=='CD8+ T'),]
NK=Inhibit[which(Inhibit$type=='NK'),]
CD4=Inhibit[which(Inhibit$type=='CD4+ T'),]

aggregate(CD8$Inhibit, by=list(type=CD8$id),mean)->CD8_mean
colnames(CD8_mean)
colnames(CD8_mean)=c('id','score')
CD8_mean$stage=substr(CD8_mean$id,1,str_locate(as.character(CD8_mean$id),"_")[,1]-1)
CD8_mean$stage=factor(CD8_mean$stage,levels=c('adjacent','AIS','MIA','INV'))
CD8_mean$median=c("")
CD8_mean[which(CD8_mean$stage=='adjacent'),]$median=median(CD8_mean[which(CD8_mean$stage=='adjacent'),'score'])
CD8_mean[which(CD8_mean$stage=='AIS'),]$median=median(CD8_mean[which(CD8_mean$stage=='AIS'),'score'])
CD8_mean[which(CD8_mean$stage=='MIA'),]$median=median(CD8_mean[which(CD8_mean$stage=='MIA'),'score'])
CD8_mean[which(CD8_mean$stage=='INV'),]$median=median(CD8_mean[which(CD8_mean$stage=='INV'),'score'])
CD8_sce<-Rmisc::summarySE(CD8_mean,measurevar='score',groupvars='stage')
CD8_sce$type='CD8'
aggregate(NK$Inhibit, by=list(type=NK$id),mean)->NK_mean
colnames(NK_mean)
colnames(NK_mean)=c('id','score')
NK_mean$stage=substr(NK_mean$id,1,str_locate(as.character(NK_mean$id),"_")[,1]-1)
NK_mean$stage=factor(NK_mean$stage,levels=c('adjacent','AIS','MIA','INV'))
NK_mean$median=c("")
NK_mean[which(NK_mean$stage=='adjacent'),]$median=median(NK_mean[which(NK_mean$stage=='adjacent'),'score'])
NK_mean[which(NK_mean$stage=='AIS'),]$median=median(NK_mean[which(NK_mean$stage=='AIS'),'score'])
NK_mean[which(NK_mean$stage=='MIA'),]$median=median(NK_mean[which(NK_mean$stage=='MIA'),'score'])
NK_mean[which(NK_mean$stage=='INV'),]$median=median(NK_mean[which(NK_mean$stage=='INV'),'score'])
NK_sce<-Rmisc::summarySE(NK_mean,measurevar='score',groupvars='stage')
NK_sce$type='NK'


aggregate(CD4$Inhibit, by=list(type=CD4$id),mean)->CD4_mean
colnames(CD4_mean)
colnames(CD4_mean)=c('id','score')
CD4_mean$stage=substr(CD4_mean$id,1,str_locate(as.character(CD4_mean$id),"_")[,1]-1)
CD4_mean$stage=factor(CD4_mean$stage,levels=c('adjacent','AIS','MIA','INV'))
CD4_mean$median=c("")
CD4_mean[which(CD4_mean$stage=='adjacent'),]$median=median(CD4_mean[which(CD4_mean$stage=='adjacent'),'score'])
CD4_mean[which(CD4_mean$stage=='AIS'),]$median=median(CD4_mean[which(CD4_mean$stage=='AIS'),'score'])
CD4_mean[which(CD4_mean$stage=='MIA'),]$median=median(CD4_mean[which(CD4_mean$stage=='MIA'),'score'])
CD4_mean[which(CD4_mean$stage=='INV'),]$median=median(CD4_mean[which(CD4_mean$stage=='INV'),'score'])
CD4_sce<-Rmisc::summarySE(CD4_mean,measurevar='score',groupvars='stage')
CD4_sce$type='CD4'

sce=rbind(CD8_sce,CD4_sce,NK_sce)

ggplot(data=sce,aes(x=stage,y=score,color=type))+
geom_errorbar(aes(ymin=score-se,ymax=score+se),width=0.1,position=position_dodge(0.1))+
scale_color_manual(values=c("CD8"='#ADDA81','NK'="#9ACDE6",'CD4'='#048D5E'))+
geom_point()
##
data.frame(celltype=c(rep('DC',2),rep("Macrophage",2),rep("Monocyte",2),rep("ILC",2)),DOWN=c(41,87,38,123,38,24,24,36),UP=c(199,215,1270,1118,37,9,173,123),type=rep(c('IAC','MIA'),4))->count_deg
count_deg$DOWN=(-count_deg$DOWN)

count<-data.frame(celltype=rep(count_deg$celltype,2),stage=rep(count_deg$type,2),count=c(count_deg$UP,count_deg$DOWN))
ggplot(count,aes(x=celltype,y=count,fill=stage))+
geom_bar(stat = 'identity', position ='dodge')

###
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=3)+
scale_fill_manual(values=c('SPP1+ B'='#DB5142','Naive B'='#938A81','CD83+ Memory B'='#FF8A9A','Plasma'='#FF913F'))+
theme_classic()+
NoLegend()

scale_color_manual(values=c('Naive B'='#938A81','Memory B'='#FAD09F','Activated B'='#FF8A9A','FCRL4+ B'='#DB5142','GC B'='#0086ED','Plasma B'='#FF913F'))
####
iTalk_res
LRPlot(iTalk_res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:20],link.arr.width=iTalk_res$cell_to_mean_exprs[1:20])
##
bitr(deg, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = 'org.Hs.eg.db')
edo <- enrichGO(id,OrgDb = 'org.Hs.eg.db',ont= "BP",readable=TRUE,pvalueCutoff  = 0.05)
edo <- pairwise_termsim(edo)
emapplot(edo, cex_category=1.5)
##
ggplot(data=tmp,aes(x=order))+
geom_smooth(aes(y=log(log2(CXCL13+1))+1),method='loess',se=F,color='#CB5052')+
geom_smooth(aes(y=log(log2(GZMB+1))+1),method='loess',se=F,color='#31729C')+
geom_smooth(aes(y=log(log2(LAG3+1))+1),method='loess',se=F,color='#EF9541')+
geom_smooth(aes(y=log(log2(ENTPD1+1))+1),method='loess',se=F,color='#8F6CAD')+
geom_smooth(aes(y=log(log2(HAVCR2+1))+1),method='loess',se=F,color='#805B50')+
geom_smooth(aes(y=log(log2(TOX+1))+1),method='loess',se=F,color='#57A457')+
theme_classic()

CreateSeuratObject(dat$'Gene Expression',min.features=200,min.cells=3)->exp
exp <- NormalizeData(exp)
# Find and scale variable features
exp  <- FindVariableFeatures(exp , selection.method = "mean.var.plot")
exp  <- ScaleData(exp , features = VariableFeatures(exp))
# Add HTO data as a new assay independent from RNA
hto=dat$Custom[c(paste0('HTO_',1:10)),colnames(exp)]
exp [["HTO"]] <- CreateAssayObject(counts = hto)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
exp <- NormalizeData(exp , assay = "HTO", normalization.method = "CLR")
exp <- HTODemux(exp, assay = "HTO", positive.quantile = 0.99)
exp<-subset(exp,HTO_classification.global=='Singlet')
subset(exp,HTO_classification %in% c(paste0('HTO-',1:5)))->Tumor
subset(exp,HTO_classification %in% c(paste0('HTO-',6:10)))->Normal

###
data.frame(CD79A=as.numeric(dat1['CD79A',]),CD27=as.numeric(dat1['CD27',]),group=group)->aa
pdf("./CD79A_CD27.pdf",useDingbats=F)
ggplot(aa,aes(x=CD79A,y=CD27))+
stat_density_2d(aes(fill = ..level..), geom="polygon")+
scale_fill_distiller(palette='YlOrRd',direction=1)+
geom_point(aes(color=group),alpha=0.75)+
scale_color_manual(values=c('T'='#2AFFE4','N'='#595757'))+
xlim(0,15)+
ylim(0,15)+
theme_bw()+
theme(panel.grid=element_blank())
dev.off()

data.frame(CD27=as.numeric(dat1['CD27',]),CD83=as.numeric(dat1['CD83',]),group=group)->bb
pdf("./CD83_CD27.pdf",useDingbats=F)
ggplot(bb,aes(x=CD27,y=CD83))+
stat_density_2d(aes(fill = ..level..), geom="polygon")+
scale_fill_distiller(palette='YlOrRd',direction=1)+
geom_point(aes(color=group),alpha=0.75)+
scale_color_manual(values=c('T'='#2AFFE4','N'='#595757'))+
xlim(0,15)+
ylim(0,15)+
theme_bw()+
theme(panel.grid=element_blank())
dev.off()

paste ./*/*.count | awk '{printf $1 "\t";for(i=2;i<=NF;i+=2) printf $i"\t";printf $i}'





fit<-survfit(Surv(OS.time,OS)~score,data=res.cat)
pdf('./')
ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'),legend='none')


data.frame(CD83=c(as.numeric(tumor['CD83',]),as.numeric(nlung['CD83',])),CD27=c(as.numeric(tumor['CD27',]),as.numeric(nlung['CD27',])),group=c(rep('T',278),rep('N',288)))->bb

res.cut <- surv_cutpoint(OS, time = "OS.time", event = "OS",
   variables = c("score"))
  res.cat <- surv_categorize(res.cut)
  ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'),legend='none')


all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(all, features = rownames(all))
all <- RunPCA(all, features = VariableFeatures(object = all))
all<-harmony::RunHarmony(all,'PID')
all <- FindNeighbors(all, dims = 1:50,reduction='harmony')
all <- FindClusters(all, resolution =1)
all <- RunUMAP(all, dims = 1:50,reduction='harmony')
pdf("./all_umap.pdf",useDingbats=F)
DimPlot(all,reduction="umap",label=T,raster=T)
dev.off()


ggplot(aa,aes(x=group,y=CXCR5,fill=group))+
 geom_boxplot()+
 scale_fill_manual(values=c('high'='#e883b4','low'='#f6cbde'))+
 theme_classic()+
 NoLegend()


 pdf('./figures/density_2d.pdf',useDingbats=F)
ggplot()+
geom_point(data=sc_info,aes(x=UMAP_1,y=UMAP_2),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$origin=='tLung'),],aes(x=UMAP_1,y=UMAP_2,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$origin=='tLung'),],aes(x=UMAP_1,y=UMAP_2),color="black")+
theme_classic()


ggplot()+
geom_point(data=sc_info,aes(x=UMAP_1,y=UMAP_2),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$origin=='tL/B'),],aes(x=UMAP_1,y=UMAP_2,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$origin=='tL/B'),],aes(x=UMAP_1,y=UMAP_2),color="black")+
theme_classic()

dev.off()




ggplot()+
geom_smooth(data=pre,aes(x=order,y=CXCL13),method='loess',se=F,color='#CB5052')+
geom_smooth(data=IAC,aes(x=order,y=CXCL13),method='loess',se=F,color='#CB5052',linetype = "dashed")+
geom_smooth(data=pre,aes(x=order,y=GZMB),method='loess',se=F,color='#31729C')+
geom_smooth(data=IAC,aes(x=order,y=GZMB),method='loess',se=F,color='#31729C',linetype = "dashed")+
geom_smooth(data=pre,aes(x=order,y=PRF1),method='loess',se=F,color='#8F6CAD')+
geom_smooth(data=IAC,aes(x=order,y=PRF1),method='loess',se=F,color='#8F6CAD',linetype = "dashed")+
geom_smooth(data=pre,aes(x=order,y=HAVCR2),method='loess',se=F,color='#EF9541')+
geom_smooth(data=IAC,aes(x=order,y=HAVCR2),method='loess',se=F,color='#EF9541',linetype = "dashed")+
geom_smooth(data=pre,aes(x=order,y=TOX),method='loess',se=F,color='#57A457')+
geom_smooth(data=IAC,aes(x=order,y=TOX),method='loess',se=F,color='#57A457',linetype = "dashed")+

theme_classic()



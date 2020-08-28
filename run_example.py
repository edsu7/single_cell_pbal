## Import functions
from pdclust_expanded import *

## Initialize variables
out_dir="/out_dir/"
index_file=sys.argv[1]#sys.argv[1]#out_dir+"/metadata.csv"
project_name="PX0740"#sys.argv[2]
ref="hg38_no_alt"#sys.argv[3]
ref_dir='/ref/'
print("output directory:"+out_dir)
print("output index file:"+index_file)
difference_type='man_dist_scaled'
core_count=16
###
subprocess.run(["mkdir","-p",out_dir+"/results/"])
subprocess.run(["mkdir","-p",out_dir+"/tmp/"])
runPipeline(index_file,out_dir,ref,ref_dir,project_name,jobs=4,threads=core_count)

### Intialize chromosomes to do operations
chr_list=["chr"+str(x) for x in list(range(1,23))]+["chrX","chrY"]

### Specify bam files. Looks for associated bam file features i.e. json,cnv,cpg meth
targets=[x for x in glob.iglob(out_dir+"/mapping/**/*.bam", recursive=True)]
target_df=findFiles(targets,out_dir,project_name)

### Pull stats ; Set custom annotation
stats,CNV_array=pullStatistics(target_df,chr_list)

annotation=pd.read_csv(out_dir+"/annotations.tsv",sep='\t').set_index("sample").rename(columns={'anno':'annotation'})
stats=stats.merge(annotation.loc[:,['delta_ct','annotation']],left_index=True,right_index=True)\
.replace("single-cell","sc")\
.replace("negative","neg")\
.replace("positive","pos")

###Initialize figure counts
figure_count=65

### Items of interest
stuff_to_plot=[
["total_reads","mapped","mapped_minus_dup"],
["cpg_count"],
["cpg_count_cov3"],
["average_meth","average_meth_cov3"],
["mapped%","dup_rate%","mapped_minus_dup%"],
["T7_conversion"],
["lambda_conversion"],
["general_conversion"],
["CNV"]
]

list_shared_axes=[
False,
False,
False,
True,
True,
True,
True,
True,
True
]

### Plot Items of interest
for x,shared_axes in zip(stuff_to_plot,list_shared_axes):
        fig=plotBoxplot(x,stats,['annotation'],project_name,shared_axes)
        plot_figure(fig,out_dir,"fig"+chr(figure_count))
        figure_count+=1

### Ready annotation for downstream plotting in heatmaps
all_annotations_category_colored=ready_annotations(stats,['annotation','average_meth'])
QC_samples=stats.query("delta_ct<=12 and CNV<=350 and annotation=='sc'").index.values.tolist()
qc_annotations_category_colored=ready_annotations(stats.loc[QC_samples,:],['annotation','average_meth'])

### Cluster by All CNV and plot
fig,cnv_clusters=plot_dendrogram_CNV(CNV_array,stats,['annotation','average_meth'],all_annotations_category_colored,3)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1
### Cluster by SC CNV only
fig,cnv_clusters=plot_dendrogram_CNV(CNV_array.loc[:,["Chromosome","Start"]+QC_samples],
                                     stats,
                                     ['annotation','average_meth'],
                                     qc_annotations_category_colored,3)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1
stats['cnv_clusters']=cnv_clusters

### Cluster by pdclust and plot
pairwise = pool_pairwise_combination(target_df.cpg.values.tolist(),core_count,chr_list,difference_type,out_dir,project_name)

pairwise_array = pairwise\
.assign(sample_1 = lambda row : row['sample_1'].str.split("/").str[-1].str.split(".").str[0])\
.assign(sample_2 = lambda row : row['sample_2'].str.split("/").str[-1].str.split(".").str[0])\
.drop_duplicates(keep='first')[['sample_1','sample_2',difference_type]]\
.pivot(index='sample_1',columns='sample_2',values=difference_type)\
.loc[QC_samples,QC_samples]


pairwise_array.to_pickle("/out_dir/results/pairwise_array.pkl")
fig,pdclust_clusters=plotPairwise_heatmap(pairwise_array,stats,
                                          ['annotation','average_meth',"cnv_clusters"],
                                          qc_annotations_category_colored,3,difference_type)

plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1
stats['pdclust_clusters']=pdclust_clusters


### Plot CpG pairwise via MDS/PCA with annotations
fig=plot_scatter_MDS(pairwise_array,stats,'pdclust_clusters',qc_annotations_category_colored)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1

fig=plot_scatter_MDS(pairwise_array,stats,'average_meth',qc_annotations_category_colored)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1

fig=plot_scatter_pca(pairwise_array,stats,'pdclust_clusters',qc_annotations_category_colored)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1

fig=plot_scatter_pca(pairwise_array,stats,'average_meth',qc_annotations_category_colored)
plot_figure(fig,out_dir,"fig"+chr(figure_count));figure_count+=1
## Merge and find DMRs by specifying groups

CpGa=target_df.loc[stats.query("pdclust_clusters=='A'").index.values.tolist(),'cpg'].values.tolist()
CpGb=target_df.loc[stats.query("pdclust_clusters=='B'").index.values.tolist(),'cpg'].values.tolist()

smoothed_python_df=merge_cpgs([CpGa,CpGb],core_count)

min_cpg_cov=3
min_cpg_in_window=3
fdr_cutoff=0.01
cpg_window=200
dm_CpGs,filtered_dmrs,fig_diff,fig_dist,fig_dmr=find_DMRs(smoothed_python_df,min_cpg_cov,min_cpg_in_window,fdr_cutoff,cpg_window)
plot_figure(fig_diff,out_dir,"fig"+chr(figure_count));figure_count+=1
plot_figure(fig_dist,out_dir,"fig"+chr(figure_count));figure_count+=1
plot_figure(fig_dmr,out_dir,"fig"+chr(figure_count));figure_count+=1

subprocess.run(["magick",out_dir+"/results/"+"*.png",out_dir+"/results/"+"out.pdf"])
stats.to_csv(out_dir+"/results/stats.csv")


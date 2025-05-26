#coding: utf-8

import subprocess as sb
import os
import sys

assert len(sys.argv) in [1,2,3]

path="data/"
batch=sys.argv[-1] if (len(sys.argv) == 2) else "" #"P7 (inas P5)"
batchnames={"P7 (inas P5)":"P7-inasP5",
        "results_user_run409":"P7ivP5",
        "results_user_run412":"P12inasP10",
        "results_user_run413":"P12ivP10",
        "results_user_run415":"P22ivP20",
        "results_user_run416":"P22inasP20"}
batchname=batchnames[batch] if (len(batch) > 0) else ""

if (not os.path.exists("results "+batch+"/")):
	sb.call("mkdir -p results "+batch+"/", shell=True)

try:
	import pandas
	import xlsx2csv
	import matplotlib_venn
	import matplotlib.pyplot
except:
	sb.check_output("python3 -m pip install pandas==1.1.3 xlsx2csv==0.7.7 matplotlib-venn==0.11.6 matplotlib==2.1.1 seaborn==0.11.0 scipy==1.5.4 scikit-learn==0.24.0 numpy==1.19.5", shell=True)
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import pandas as pd

##parameters
all_batches = list(batchnames.keys()) if (len(batch) == 0) else [batch]

thres_fc = 1e-2
alpha_p = 0.05
if (len(sys.argv)<3):
	list_genes=["CD", "DESeq2"][1]
else:
	list_genes=sys.argv[2]
	assert list_genes in ["CD", "DESeq2"]

## Plot Venn diagrams
if (len(batch) > 0):
    conditions = ["IL1_Dose1", "IL1_Dose2", "IL1_Dose3"]
    compare_condition = ["PBS"]
    control = "IL1"
    runname = "run404" if (batch == "P7 (inas P5)") else batch.split("_")[-1]
    path_ = path+batch+"/"
    for condition in conditions+compare_condition:
        fname = "diffanaresultsannotation_deseq2_PREMSTEM_A2020_"+runname+"-diffana_"+batchname+"_"+condition+"_vs_"+batchname+"_"+control
        if (not os.path.exists(path_+fname+".csv")):
            assert os.path.exists(path_+fname+".xlsx")
            sb.call("xlsx2csv \""+path_+fname+".xlsx\" > \""+path_+fname+".csv\"", shell=True)

    if (list_genes == "DESeq2"):
        ## Compute DEGs from DeSeq2 DE analysis
        from utils import BH_adjust
        degs_di = {}
        for condition in conditions+compare_condition:
            ## get data
            fname = "diffanaresultsannotation_deseq2_PREMSTEM_A2020_"+runname+"-diffana_"+batchname+"_"+condition+"_vs_"+batchname+"_"+control
            df=pd.read_csv(path_+fname+".csv")
            df.index = df["Gene name"] #df["Id"]
            fc_col = 'log2foldchange '+batchname+'_'+condition+' vs '+batchname+'_'+control
            cols = ['baseMean', fc_col, 'Wald test p-value', 'BH adjusted p-values']
            df=df[cols].dropna().groupby(level=0).mean()
            pvalues = df[['Wald test p-value']]
            fcs = df[fc_col]
            #adj_pvalues = df[['BH adjusted p-values']]
            cpvalues, significant = BH_adjust(pvalues.values.flatten(), alpha=alpha_p)
            corr_pvalues = pd.DataFrame(cpvalues[significant], index=pvalues.index[significant], columns=["corr_pvalues"]).sort_values(by="corr_pvalues", ascending=True)
            degs_di.setdefault(condition, [g for g in list(corr_pvalues.index) if (abs(fcs.loc[g]) > thres_fc)])

    elif (list_genes == "CD"):
        degs_di = {}
        for condition in conditions+compare_condition:
            fname = "results "+batch+"/"+batch+condition+"-"+control+".csv"
            df=pd.read_csv(fname)
            ## already selected by p-value < 0.05
            assert alpha_p == 0.05
            degs_di.setdefault(condition, list(df[df.columns[0]]))

    else:
        raise ValueError("No source for DEG lists provided.")

    ## Venn diagrams
    figure, axes = plt.subplots(1, len(conditions), figsize=(15, 5))
    for sc, condition in enumerate(conditions):
        degs0, degs1 = [set(degs_di[c]) for c in [condition]+compare_condition]
        intersec = degs0.intersection(degs1)
        diff0, diff1 = degs0.difference(degs1), degs1.difference(degs0)
        venn2(subsets = (len(diff0), len(diff1), len(intersec)), set_labels=[""]*2, ax=axes[sc])
        axes[sc].legend(labels=[condition]+compare_condition)
    axes[1].set_title("Venn diagram of DE genes (compared to "+control+", batch "+batch+")")
    plt.savefig("results "+batch+"/venn_diagram_DEG_"+"_".join(conditions)+"_VS_"+compare_condition[0]+".png", bbox_inches="tight")

## Otherwise build a confusion matrix
else:
    conditions = ["IL1_Dose1", "IL1_Dose2", "IL1_Dose3", "PBS"]
    control = "IL1"
    modes = ["inas", "iv"]
    times = ["P7-P5", "P12-P10", "P22-P20"]
    signed_degs = {}
    genes = set()
    for batch in batchnames:
        runname = "run404" if (batch == "P7 (inas P5)") else batch.split("_")[-1]
        path_ = path+batch+"/"
        for condition in conditions:
            mode = "inas" if ("inas" in batchnames[batch]) else "iv"
            time = "P7-P5" if ("P7" in batchnames[batch]) else ("P12-P10" if ("P12" in batchnames[batch]) else "P22-P20")
            analysis = ",".join([condition, mode, time])
            if (list_genes == "DESeq2"):
                fname = "diffanaresultsannotation_deseq2_PREMSTEM_A2020_"+runname+"-diffana_"+batchnames[batch]+"_"+condition+"_vs_"+batchnames[batch]+"_"+control
                df=pd.read_csv(path_+fname+".csv")
                df.index = df["Gene name"] #df["Id"]
                fc_col = 'log2foldchange '+batchnames[batch]+'_'+condition+' vs '+batchnames[batch]+'_'+control
                cols = ['baseMean', fc_col, 'Wald test p-value', 'BH adjusted p-values']
                df=df[cols].dropna().groupby(level=0).mean()
                pvalues = df[['Wald test p-value']]
                fcs = df[fc_col]
                from utils import BH_adjust
                cpvalues, significant = BH_adjust(pvalues.values.flatten(), alpha=alpha_p)
                corr_pvalues = pd.DataFrame(cpvalues[significant], index=pvalues.index[significant], columns=["corr_pvalues"]).sort_values(by="corr_pvalues", ascending=True)
                up_degs = set([g for g in list(corr_pvalues.index) if (fcs.loc[g] > thres_fc)])
                dn_degs = set([g for g in list(corr_pvalues.index) if (fcs.loc[g] < -thres_fc)])
                #fcs = [[g,float(fcs.loc[g])] for g in list(corr_pvalues.index) if (abs(fcs.loc[g]) > thres_fc)]
            elif (list_genes == "CD"):
                fname = "results "+batch+"/"+batch+condition+"-"+control+".csv"
                df=pd.read_csv(fname)
                df.index = df[df.columns[0]]
                ## already selected by p-value < 0.05
                assert alpha_p == 0.05
                up_degs = set([g for g in list(df.index) if (df.loc[g][batch+condition+"||"+control] > thres_fc)])
                dn_degs = set([g for g in list(df.index) if (df.loc[g][batch+condition+"||"+control] < -thres_fc)])
                #fcs = [[g, float(df.loc[g][condition+"||"+control])] for g in list(df.index) if (abs(df.loc[g][condition+"||"+control]) > thres_fc)]
            else:
                raise ValueError("No source for DEG lists provided.")
            total = up_degs.union(dn_degs)
            genes = genes.union(total)
            signed_degs.setdefault(analysis, {'up': up_degs, "dn": dn_degs, "total": total})#, "fcs": fcs})
    analyses = list(signed_degs.keys())
    signed_degs_df = pd.DataFrame([["|".join(signed_degs[analysis][idx]) for idx in ["up", "dn", "total"]] for analysis in analyses], columns=["up", "dn", "total"], index=analyses)
    signed_degs_df.to_csv("signed_degs_df"+("" if (list_genes == "DESeq2") else "_"+list_genes)+".csv")
    init_df = lambda _ : pd.DataFrame([], index=analyses, columns=analyses)
    confusion_df, annot_df, annot_df_up, annot_df_dn, perc_df, perc_df_up, perc_df_dn = [init_df(0) for i in range(7)]
    import numpy as np
    for ia, a in enumerate(analyses):
        for b in analyses[ia:]:
            nup, ndn, ntotal = [len(list(signed_degs[a][idx].intersection(signed_degs[b][idx]))) for idx in ["up", "dn", "total"]]
            nbkup, nbkdn, nbktotal = [len(list(signed_degs[a][idx].union(signed_degs[b][idx]))) for idx in ["up", "dn", "total"]]
            confusion_df.loc[b][a] = "|".join(["#up="+str(nup)+"/"+str(nbkup), "#dn="+str(ndn)+"/"+str(nbkdn), "#tot="+str(ntotal)+"/"+str(nbktotal)])
            annot_df.loc[b][a] = str(ntotal)+"/"+str(nbktotal)
            annot_df_up.loc[b][a] = str(nup)+"/"+str(nbkup)
            annot_df_dn.loc[b][a] = str(ndn)+"/"+str(nbkdn)
            perc_df.loc[b][a] = 100*ntotal/float(nbktotal) if (nbktotal > 0) else np.nan
            perc_df_up.loc[b][a] = 100*nup/float(nbkup) if (nbkup > 0) else np.nan
            perc_df_dn.loc[b][a] = 100*ndn/float(nbkdn) if (nbkdn > 0) else np.nan
            perc_df.loc[a][b] = perc_df.loc[b][a]
            perc_df_up.loc[a][b] = perc_df_up.loc[b][a]
            perc_df_dn.loc[a][b] = perc_df_dn.loc[b][a]
            annot_df.loc[a][b] = annot_df.loc[b][a]
            annot_df_up.loc[a][b] = annot_df_up.loc[b][a]
            annot_df_dn.loc[a][b] = annot_df_dn.loc[b][a]

    confusion_df.to_csv("confusion_df"+("" if (list_genes == "DESeq2") else "_"+list_genes)+".csv")

    print(confusion_df.head())

    genes = list(genes)
    confusion_gene = np.matrix(np.zeros((len(genes), len(analyses))))
    for ia, a in enumerate(analyses):
        #fcs, fc_genes = [x[1] for x in signed_degs[a]["fcs"]], [x[0] for x in signed_degs[a]["fcs"]]
        #confusion_gene[:, ia] = np.array([(0. if (not g in fc_genes) else fcs[fc_genes.index(g)]) for g in genes]).reshape((len(genes), 1))
        confusion_gene[:, ia] = np.array([int(g in list(signed_degs[a]["total"])) for g in genes]).reshape((len(genes), 1))
    confusion_gene_df = pd.DataFrame(confusion_gene, index=genes, columns=analyses)
    genes = list(sorted(genes, key=lambda g : np.sum(confusion_gene_df.loc[g].values), reverse=True))
    confusion_gene_df = confusion_gene_df.loc[genes]

    confusion_gene_df.to_csv("confusion_gene_df"+("" if (list_genes == "DESeq2") else "_"+list_genes)+".csv")

    print(confusion_gene_df.head())

    import seaborn as sns
    import numpy as np
    from matplotlib.patches import Patch

    ## draw a heatmap
    def draw_heatmap(df, annot=None, type_="", order_row=True, order_col=True, plot_colorbar=True, row_colors=None, col_colors=None, lut_row=None, lut_col=None):
        fig, ax = plt.subplots(1, 1, figsize=(20, 20))
        mask = np.zeros_like(df.values)
        mask[np.triu_indices_from(mask, k=1)] = True
        method, metric = "average", "euclidean"#"cosine" if (list_genes == "CD") else "euclidean"
        with sns.axes_style("white"):
            if (str(annot) != "None"):
                #cns = sns.heatmap(df.fillna(0), annot=annot, ax=ax, fmt="", mask=mask, square=True, row_colors=row_colors, col_colors=col_colors)
                cns = sns.clustermap(df.fillna(0), annot=annot, fmt="", method=method, metric=metric, row_colors=row_colors, col_colors=col_colors, z_score=None, standard_scale=None, row_cluster=order_row, col_cluster=order_col, dendrogram_ratio=(.1, .2), figsize=(20,20))
            else:
                cns = sns.clustermap(df.fillna(0), method="average", metric="cosine", z_score=None, row_colors=row_colors, col_colors=col_colors, standard_scale=None, row_cluster=order_row, col_cluster=order_col, dendrogram_ratio=(.1, .2), figsize=(20,20))
            if (not plot_colorbar):
                cns.cax.set_visible(False)
            if (str(row_colors) != "None" or str(col_colors) != "None"):
                handles = []
                lut_all = {}
                for ix, x in enumerate([row_colors, col_colors]):
                    lut = [lut_row, lut_col][ix]
                    if (str(x) != "None" and str(lut) != "None"):
                        handles += [Patch(facecolor=lut[name]) for name in lut]
                        lut_all.update(lut)
                plt.legend(handles, lut_all, bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')#title
        ax.set_title("Confusion matrix of lists of "+(type_+" " if (len(type_) > 0) else "")+"DEGs for all doses, batches (compared to "+control+")")
        plt.savefig("results/heatmap_confusion_matrix"+("_"+type_ if (len(type_) > 0) else "")+("" if (list_genes == "DESeq2") else "_"+list_genes)+".png", bbox_inches="tight")

    modes = [a.split(",")[1] for a in analyses]
    times = [a.split(",")[-1] for a in analyses]
    lut_mode = dict(zip(list(set(modes)), sns.color_palette("Paired")))
    lut_time = dict(zip(list(set(times)), sns.color_palette("Set2")))
    row_colors = list(map(lambda x : lut_mode[x], modes))
    col_colors = list(map(lambda x : lut_time[x], times))
    draw_heatmap(perc_df, annot_df, row_colors=row_colors, col_colors=col_colors, lut_row=lut_mode, lut_col=lut_time)
    draw_heatmap(perc_df_up, annot_df_up, type_="up-regulated", row_colors=row_colors, col_colors=col_colors, lut_row=lut_mode, lut_col=lut_time)
    draw_heatmap(perc_df_dn, annot_df_dn, type_="down-regulated", row_colors=row_colors, col_colors=col_colors, lut_row=lut_mode, lut_col=lut_time)
    ## Histogram for distribution of genes in terms of presence in DEG lists
    batches = [str(x) for x in range(1, len(confusion_gene_df.columns)+1)]
    pd.DataFrame(confusion_gene_df.sum(axis=1), columns=["#groupes où DE"]).to_csv("histogramme_file_"+list_genes+".csv")
    gene_presence = confusion_gene_df.sum(axis=1).values.T.tolist()
    plt.figure(figsize=(10,8))
    kwargs = dict(alpha=0.85, density=False, histtype='bar', ec='black')
    plt.hist(gene_presence, bins=len(batches), **kwargs)
    #plt.xticks(range(1,len(batches)), batches, rotation=80)
    locs, labels = plt.xticks()
    for sb, b in enumerate(locs.tolist()):
        ngenes = len(np.argwhere(np.array(gene_presence)==sb+1))
        plt.text(sb+1+0.1*int(len(str(ngenes)) < 3 and sb < 14), ngenes+5, str(ngenes))
    plt.xlabel("#batches in which gene belongs to DEG list")
    plt.ylabel("#genes")
    plt.savefig("histogram_confusion_gene_df"+("" if (list_genes == "DESeq2") else "_"+list_genes)+".png")
    plt.close()
    exit()
    from scipy.spatial.distance import pdist, squareform
    N=200 #readable
    most_de_genes = confusion_gene_df.sum(axis=1).sort_values(ascending=False).index
    confusion_gene_df_ = confusion_gene_df.loc[most_de_genes[:N]].T
    dist_gene = np.round(100*(1-squareform(pdist(confusion_gene_df_, metric="euclidean"))/float(len(analyses))), 1)
    dist_gene_df = pd.DataFrame(dist_gene, index=confusion_gene_df_.index, columns=confusion_gene_df_.index)
    draw_heatmap(dist_gene_df, annot=dist_gene_df, type_="genes (distance)", row_colors=row_colors, col_colors=col_colors, lut_row=lut_mode, lut_col=lut_time)
    ## Trim most DE genes
    N=200 #readable
    confusion_gene_df_ = confusion_gene_df.loc[most_de_genes[:N]]
    draw_heatmap(confusion_gene_df_, annot=np.asarray(confusion_gene_df_, dtype=int), type_="genes", order_row=False, plot_colorbar=False)

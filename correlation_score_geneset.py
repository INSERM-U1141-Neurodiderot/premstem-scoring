#coding: utf-8

## Compute correlation between cosine score and intersection of genes
from scipy.stats import pearsonr, spearmanr
from sklearn import linear_model
from sklearn.metrics import r2_score
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

batchnames={"P7 (inas P5)":"P7-inasP5",
        "results_user_run409":"P7ivP5",
        "results_user_run412":"P12inasP10",
        "results_user_run413":"P12ivP10",
        "results_user_run415":"P22ivP20",
        "results_user_run416":"P22inasP20"}

test = ["pearsonr", "spearmanr"][0]

alpha_p=0.05 #OK for CD genes
thres_fc=0#1e-2

print("Filter genes fdr<=%.2f thres_fc>=%.2f" % (alpha_p,thres_fc))

cosine_scores, intersection_genes = {}, {}
conditions, control = ["IL1_Dose1", "IL1_Dose2", "IL1_Dose3"], "PBS"
for batch in batchnames:
    scores = pd.read_csv("results "+batch+"/dose_ranking.csv")
    scores.index = scores[scores.columns[0]]
    scores = scores[scores.columns[1]]
    signatures = [None]*len(conditions+[control])
    for ic, c in enumerate(conditions+[control]):
        df = pd.read_csv("results "+batch+"/"+batch+""+c+"-IL1.csv")
        df.index = df[df.columns[0]]
        df = df.drop(columns=df.columns[:1])
        df = df.loc[np.abs(df[df.columns[0]]) > thres_fc]
        signatures[ic] = df
    for ic, c in enumerate(conditions):
        gene_list = list(set(list(signatures[ic].dropna().index)).intersection(set(list(signatures[-1].dropna().index))))
        size = len(gene_list)
        nm = batch+":"+c+"_IL1:"+control+"_IL1"
        intersection_genes.setdefault(nm, size)
        cosine_scores.setdefault(nm, float(scores.loc[batch+c+"||IL1"]))
scores, sizes = [cosine_scores[k] for k in cosine_scores], [intersection_genes[k] for k in cosine_scores]

doses = [k.split(":")[1].split("_")[1] for k in cosine_scores]
ages = [batchnames[k.split(':')[0]][:3] for k in cosine_scores]
ages = ["P7" if ("P7" in a) else a for a in ages]
modes = ["iv" if ("iv" in batchnames[k.split(":")[0]]) else "inas" for k in cosine_scores]
excel_file = pd.DataFrame([scores,sizes,doses,modes,ages], index=["score","#genes in common","dose","mode","age"], columns=[k for k in cosine_scores])
excel_file.to_csv("scores_vs_geneset_fdr=%.2f_fc=%.2f.csv" % (alpha_p,thres_fc))

stat, p = eval(test)(scores, sizes)
linreg = linear_model.LinearRegression().fit(np.ravel(sizes).reshape(-1, 1), np.ravel(scores).reshape(-1, 1))
pred_scores = linreg.predict(np.ravel(sizes).reshape(-1,1))

matplotlib.rcParams.update({'font.size': 22})
ffsize=25
plt.figure(figsize=(15,10))
import seaborn as sns
edges = {"iv": "none", "inas":""}
colors = {"Dose1": "c", "Dose2": "b", "Dose3": "k"}#, "r", "g", "y", "m", "c"]
markers = {'P7': "D","P12": "*","P22": "."}
orders = {5: "1", 10: "2", 20: "3"}
for ik, k in enumerate(list(sorted(list(batchnames.keys()), reverse=False, key=lambda k : orders[int("-".join(batchnames[k].split("-")).split('P')[-1])]+("iv" if ("iv" in batchnames[k]) else "inas")))):
    ks = [ky for ky in cosine_scores if (k in ky)]
    sizes_ = [sz for sz in [intersection_genes[ky] for ky in ks]]
    scores_ = [sc for sc in [cosine_scores[ky] for ky in ks]]
    doses_ = ["".join(":".join(ky.split(":")[1:]).split("_IL1")[0].split("_")[-1].split(":")[0]) for ky in ks]
    plt.plot(sizes_[0], scores_[0], colors[doses_[0]]+(markers["P7"] if ("P7" in batchnames[k]) else (markers["P12"] if ("P12" in batchnames[k]) else markers["P22"])), markerfacecolor=(edges["iv"] if ("iv" in batchnames[k]) else colors[doses_[0]]), 
            markersize=40 if ("P7" not in batchnames[k]) else 25, label="".join(batchnames[k].split("-")))#k.split("_")[-1])
    for iky,ky in enumerate(ks[1:]):
        iky += 1
        plt.plot(sizes_[iky], scores_[iky], colors[doses_[iky]]+(markers["P7"] if ("P7" in batchnames[k]) else (markers["P12"] if ("P12" in batchnames[k]) else markers["P22"])), markerfacecolor=(edges["iv"] if ("iv" in batchnames[k]) else colors[doses_[iky]]), 
                markersize=40 if ("P7" not in batchnames[k]) else 25)#, label="".join(batchnames[k].split("-")))#k.split("_")[-1])
    for iky,ky in enumerate(ks):
        dose = doses_[iky]
        if ((batchnames[k]=="P22ivP20" and dose=="Dose2") or (batchnames[k]=="P7-inasP5" and dose=="Dose2") or (batchnames[k]=="P12inasP10" and dose=="Dose3")):
            plt.text(intersection_genes[ky]-1600, cosine_scores[ky]-0.01, s=dose, fontsize=int(ffsize*3/4))
        else:
            plt.text(intersection_genes[ky]-500*int(thres_fc==0)+1000, cosine_scores[ky]-0.01, s=dose, fontsize=int(ffsize*3/4))
plt.plot(sizes, pred_scores, "r-", linewidth=2, label=r"$R^2="+str(round(r2_score(scores,pred_scores), 2))+"$")
plt.xlabel("Number of common DE genes with the control signature",fontsize=ffsize)
plt.ylabel("Score",fontsize=ffsize)
plt.xticks(fontsize=ffsize)
plt.yticks(fontsize=ffsize)
plt.legend(loc="upper right")
#plt.title("Correlation cosine score ~ gene intersection\n("+test+" = "+str(round(stat, 2))+", p="+str(round(p, 3))+", thres_fc="+str(thres_fc)+")")
plt.savefig("correlation_score_intersections"+("" if (thres_fc==0) else "thres_fc="+str(thres_fc))+".png", bbox_inches="tight")

print("r = %0.2f, p-value = %0.2f, N=%d" % (stat, p, len(sizes)))
print("R^2 = %0.2f" % r2_score(scores,pred_scores))
ids = [k.split(":")[0].split("_")[-1]+"|"+"".join(":".join(k.split(":")[1:]).split("_IL1")).split("_")[-1] for k in cosine_scores]
scs, szs = [cosine_scores[k] for k in cosine_scores], [intersection_genes[k] for k in cosine_scores]
res = pd.DataFrame([scs, szs], index=["Score", "Sizes"], columns=ids)
print(res.T.sort_values("Score", ascending=False))

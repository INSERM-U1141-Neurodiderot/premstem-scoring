#coding: utf-8

## Needs normalised counts matrix (if performed on a single batch) or raw counts matrices for each batch (if performed on several batches)

import sys
assert len(sys.argv) in [1,2]

### modifiable parameters for signature computation
alpha=0.05
thr_var=0.8
thr_grp_var=float("inf")
sig_only=True
nperm=100
gene_id="Gene name"
trim_beforehand=False
path="data/"
batch="" if (len(sys.argv) == 1) else sys.argv[-1] #"P7 (inas P5)"
batchnames={"P7 (inas P5)":"P7-inasP5", 
        "results_user_run409":"P7ivP5", 
        "results_user_run412":"P12inasP10", 
        "results_user_run413":"P12ivP10", 
        "results_user_run415":"P22ivP20", 
        "results_user_run416":"P22inasP20"}
batchname=(batchnames[batch]+"_") if (len(batch) > 0) else ""
all_batches=list(batchnames.keys())
assert batch in all_batches or len(batch) == 0

import os
import subprocess as sb

from utils import load_dataset, batch_normalize, compute_CD_signature, cosine_score

if (len(batch) > 0):
    assert os.path.exists(path+batch)

results_folder="results"+(" "+batch if (len(batch) > 0) else "")+"/"
if (not os.path.exists(results_folder)):
	sb.call("mkdir \""+results_folder+"\"", shell=True)

try:
	from geode import chdir
except:
	sb.check_output("python3 -m pip install git+https://github.com/Maayanlab/geode.git pandas==1.1.3 xlsx2csv==0.7.7", shell=True)
import pandas as pd

#get data
if (len(batch) > 0):
   fname = load_dataset(batch, "normalised", verbose=False)
   df=pd.read_csv(fname+".csv")
   df.index = df[gene_id]
   df=df[list(filter(lambda x : batchname in x, df.columns))].groupby(level=0).mean()
else:
    ## Collect raw expression data for each batch, build matrix + metadata file, quantile-normalize and run COMBAT
    fname = "combat_adjusted_data"
    if (not os.path.exists(fname+".csv")):
        samples_nested, genes = [None]*len(all_batches), []
        for ib, batch_ in enumerate(all_batches):
            fname_ = load_dataset(batch_, "rawPooled", verbose=False)
            df_=pd.read_csv(fname_+".csv")
            batchname_=batchnames[batch_]
            samples_ = list(filter(lambda x : batchname_ in x, df_.columns))
            samples_nested[ib] = samples_
            genes = list(set(genes+list(df_[gene_id])))
        samples = [si for s in samples_nested for si in s]
        raw_df = pd.DataFrame([], index=genes, columns=samples)
        batches = [spp for sp in [[batch_]*len(samples_nested[ib]) for ib, batch_ in enumerate(all_batches)] for spp in sp]
        sample_names = ["".join(sample_.split("_")[0].split("-")) for sample_ in samples]
        modes = ["iv" if ("iv" in sample_name_) else "inas" for sample_name_ in sample_names]
        times = ["-".join(sample_name_.split(modes[isample])) for isample, sample_name_ in enumerate(sample_names)]
        doses = ["Control" if ("Dose" not in sample_name_) else sample_name_.split("_")[-2] for sample_name_ in sample_names]
        metadata = pd.DataFrame([batches, times, modes, doses], columns=samples, index=["Batch", "Time", "Mode", "Dose"])
        for ib, batch_ in enumerate(all_batches):
            fname_ = load_dataset(batch_, "rawPooled")
            raw_df_ = pd.read_csv(fname_+".csv")
            raw_df_.index = raw_df_[gene_id]
            raw_df_ = raw_df_.loc[genes]
            raw_df_=raw_df_[samples_nested[ib]].groupby(level=0).mean()
            raw_df[samples_nested[ib]] = raw_df_.values
        df = batch_normalize(raw_df, metadata, batch_var="Batch")
    else:
        df=pd.read_csv(fname+".csv").dropna()
        df.index=df[df.columns[0]]
        df=df.drop(columns=df.columns[:1])

print(df.head())

if (trim_beforehand):
	## trim out low variance genes across all samples
	gene_var = df.var(axis=1).values.flatten()
	print("#genes = "+str(len(gene_var)))
	gene_var = gene_var[gene_var > thr_var]
	print("#genes with expression variance across samples > "+str(thr_var)+" = "+str(len(gene_var)))
	df = df.loc[df.var(axis=1) > thr_var]

## IL1 (diseased samples), PBS (healthy samples), IL1_Dose* (diseased samples treated with dose *)
df.loc["annotation"] = ["_".join(x.split("_")[1:-1]) if (len(x.split("_")) > 3) else x.split("_")[1] for x in df.columns]

from time import time
start_time = time()

## Signature building:
#- Signature associated with healthy phenotype is the one built using CD[PBS||IL1]
#- Signature associated with dose*-treated phenotype is the one built using CD[IL1_Dose*||IL1]
if (not os.path.exists(results_folder+"CD_signatures.csv")):
    from functools import reduce
    sig_di = {}
    all_batches = [batch] if (len(batch) > 0) else all_batches
    for batch_ in all_batches:
        df_ = df[list(filter(lambda x : batchnames[batch_] in x, df.columns))]
        print(batch_)
        for sign in ["PBS||IL1", "IL1_Dose1||IL1", "IL1_Dose2||IL1", "IL1_Dose3||IL1"]:
            subdf = df_.iloc[:,(df_.loc["annotation"].values==sign.split("||")[0])|(df_.loc["annotation"].values==sign.split("||")[-1])]
            annot = [int(x == sign.split("||")[0])+1 for x in subdf.loc["annotation"]]
            sig = compute_CD_signature(subdf.loc[list(filter(lambda x : x != "annotation", subdf.index))], samples=annot, name=sign, alpha=alpha, thr_var=-float("inf") if (trim_beforehand) else thr_var, thr_grp_var=float("inf"), sig_only=True, nperm=100)
            sig.columns = [batch_+sig.columns[0]]
            sig.to_csv(results_folder+batch_+("-".join(sign.split("||")))+".csv")
            sig_di.setdefault(batch_+sign, sig)
    sig_df = reduce(lambda x,y: x.join(y, how="outer"), list(sig_di.values()))
    sig_df.to_csv(results_folder+"CD_signatures.csv")
else:
    sig_df = pd.read_csv(results_folder+"CD_signatures.csv")
    sig_df.index = sig_df[sig_df.columns[0]]
    sig_df = sig_df.drop(columns=sig_df.columns[:1])

print(sig_df.head())

## Ranking
#- For each dose *, compute the cosine similarity score between CD[IL1_Dose*||IL1] and CD[PBS||IL1]
#- The higher the cosine score is, the most similar to the healthy phenotype is the treated phenotype
#- Using these cosine scores, we can then rank doses by increasing similarity to the healthy phenotype
scores_di = {}
all_batches = [batch] if (len(batch) > 0) else all_batches
for batch_ in all_batches:
    for sign in list(filter(lambda x : "PBS||IL1" not in x and batch_ in x, sig_df.columns)):#sig_di.keys())):
        #score = cosine_score(sig_di[sign], sig_di[batch_+"PBS||IL1"], scale=False)
        #score = cosine_score(sig_df[[sign]].fillna(0), sig_df[[batch_+"PBS||IL1"]].fillna(0), scale=False)
        score = cosine_score(sig_df[[sign]].dropna(), sig_df[[batch_+"PBS||IL1"]].dropna(), scale=False)
        scores_di.setdefault(sign, score)

scores = pd.DataFrame(scores_di, index=["rank score"]).T.sort_values("rank score", ascending=False)
print(scores)
scores.to_csv(results_folder+"dose_ranking.csv")

end_time = time()
print("Elapsed time = "+str(round(end_time-start_time, 2))+" sec.")

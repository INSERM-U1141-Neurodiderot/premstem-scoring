# coding:utf-8

import numpy as np
import pandas as pd
import subprocess as sb
import os

def load_dataset(batch, type_, path="data/", verbose=True):
    runname = "run404" if (batch == "P7 (inas P5)") else batch.split("_")[-1]
    fname = path+batch+"/diffanaresultsannotation_deseq2_PREMSTEM_A2020_"+runname+"-normalisation_"+type_+"CountMatrix"
    #res = sb.check_output("[ -d \""+fname+".csv\" ] && echo 1 || echo 0", shell=True).decode("utf-8").split("\n")[0]
    if (not os.path.exists(fname+".csv")):
        assert os.path.exists(fname+".xlsx")
        sb.call("xlsx2csv \""+fname+".xlsx\" > \""+fname+".csv\"", shell=True)
    if (verbose):
        print("\""+fname+".csv\"")
        sb.call("head -n5 \""+fname+".csv\"", shell=True)
    return fname

#' @param x NumPy array of float of size n
#' @param P NumPy array of float of size n
#' @return Python float: cosine score value
def cosine_score(x, P, scale=False, verbose=True):
    '''Computes cosine score S : x, P -> (x.P)/(||x||.||P||)'''
    df = x.join(P, how="inner")
    if (verbose):
        print("Size of intersection |x & P| = "+str(len(df.index))+" (|x|="+str(len(x.index))+", |P|="+str(len(P.index))+")")
    if (scale):
        scale_func = lambda col : (np.array(col)-np.mean(col))/np.std(col)
        df = df.apply(lambda col : scale_func(col.values.flatten().tolist()))
    x, P = df[[df.columns[0]]], df[[df.columns[1]]]
    cos = float(np.dot(np.transpose(x), P)/(np.linalg.norm(x, 2)*np.linalg.norm(P, 2)))
    return cos

## https://www.statsmodels.org/stable/_modules/statsmodels/stats/multitest.html#fdrcorrection
## adjusted for FDR using classic Benjamini & Hochberg (1995) method
## from fdrcorrection function (method="indep") in package statsmodels
def BH_adjust(pvalues, alpha=0.05):
	_ecdf = lambda x : np.arange(1,len(x)+1)/float(len(x))
	sorted_ind = np.argsort(pvalues)
	sorted_pvalues = pvalues[sorted_ind]
	ecdffactor = _ecdf(sorted_pvalues)
	reject = sorted_pvalues <= ecdffactor*alpha
	if reject.any():
		rejectmax = max(np.nonzero(reject)[0])
		reject[:rejectmax] = True
	pvals_corrected_raw = sorted_pvalues / ecdffactor
	pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
	del pvals_corrected_raw
	pvals_corrected[pvals_corrected>1] = 1
	pvals_corrected_ = np.empty_like(pvals_corrected)
	pvals_corrected_[sorted_ind] = pvals_corrected
	del pvals_corrected
	reject_ = np.empty_like(reject)
	reject_[sorted_ind] = reject
	adj_pvalues = pvals_corrected_
	significance = reject_
	return adj_pvalues, significance

## Implementation of CD in http://www.maayanlab.net/CD/
#' @param df Pandas DataFrame of expression profiles
#' @param samples Python list of integers (1: control or 2: treated)
#' @param thres threshold number of significant genes
#' @param nperm number of iterations for p-value computation
#' @return sig signature CD[treated||control] restricted to significantly DE genes
def compute_CD_signature(df, samples=[], thres=None, nperm=100, name="signature", thr_var=0., thr_grp_var=float("inf"), sig_only=True, alpha=0.05, binarize=False, verbose=True):
    from geode import chdir
    df = df.dropna()
    ## remove genes with low expression variance across samples and high variance among same-group samples
    df = df.loc[(np.var(df.values, axis=1)>thr_var)&(np.var(df.iloc[:,np.array(samples)==1].values, axis=1)<thr_grp_var)&(np.var(df.iloc[:,np.array(samples)==2].values, axis=1)<thr_grp_var)]
    ## returns chdir_res: a list of tuples (signed magnitude CD, gene name, p-value) sorted by decreasing magnitude
    ## automatically computes significant genes (calculate_sig=True computes p-values, sig_only=True returns only significant genes)
    chdir_res = chdir(np.asarray(df.values, dtype=np.float64), samples, df.index, calculate_sig=True, nnull=nperm, sig_only=sig_only, gamma=1.)
    if (not sig_only):
    	pvalues = np.array(list(map(lambda x : x[2], chdir_res)))
    	_, significant = BH_adjust(pvalues, alpha=alpha)
    else:
    	_, significant = [], [True]*len(chdir_res)
    magnitudes = np.array(list(map(lambda x : x[0], chdir_res)))[significant]
    if (binarize):
        magnitudes = np.array([int(m > 0)-int(m < 0) for m in magnitudes.tolist()])
    genes = np.array(list(map(lambda x : x[1], chdir_res)))[significant]
    sig = pd.DataFrame(magnitudes, index=genes, columns=[name])
    if (verbose):
        print(sig.head())
    return sig

## https://github.com/Maarten-vd-Sande/qnorm/blob/master/qnorm/quantile_normalize.py
## adapted for Python 2.7
def quantile_normalize(df):
	n_rows, n_cols = np.shape(df.values)
	qnorm = np.zeros((n_rows, n_cols), dtype=np.float64)
	sorted_val = np.zeros((n_rows, n_cols), dtype=np.float64)
	sorted_idx = np.zeros((n_rows, n_cols), dtype=np.uint32)
	sorted_rowmeans = np.zeros(n_rows, dtype=np.float64)
	for col_i in range(n_cols):
		argsort = np.argsort(df.values[:, col_i])
		sorted_idx[:, col_i] = argsort
		sorted_val[:, col_i] = np.array([df.values[i, col_i] for i in argsort])
	for row in range(n_rows):
		sorted_rowmeans[row] = np.mean(sorted_val[row, :])
	# we quantile normalize separately per column
	for col_i in range(n_cols):
		i = 0
		# we fill out a column not from lowest index to highest index,
		# but we fill out qnorm from lowest value to highest value
		while i < n_rows:
			n = 0
			val = 0.0
			# since there might be duplicate numbers in a column, we search for
			# all the indices that have these duplicate numbers. Then we take
			# the mean of their rowmeans.
			while ((i + n < n_rows) and sorted_val[i, col_i] == sorted_val[i + n, col_i]):
				val += sorted_rowmeans[i + n]
				n += 1
			# fill out qnorm with our new value
			if n > 0:
				val /= n
				for j in range(n):
					idx = sorted_idx[i + j, col_i]
					qnorm[idx, col_i] = val
			i += n
	return pd.DataFrame(qnorm, columns=df.columns, index=df.index)

def batch_normalize(df, metadata, batch_var="Batch", path="", verbose=True):
    df = quantile_normalize(df)
    ## 1. Remove Genes with low standard deviation
    thres = 0.8
    keep_ids = np.sqrt(np.var(df.values, axis=1)) > thres 
    df = df.loc[keep_ids]
    ## 2. Make sure none of your batches contain genes with only zeroes
    for batch in list(set(metadata.loc[batch_var])):
    	ids = metadata.loc[batch_var] == batch
    	sum_expr = np.sum(df.values[:,ids], axis=1)
    	df = df.loc[sum_expr > 0]
    ## 3. Run ComBat
    if (not os.path.exists("combat.py")):
    	sb.call("wget -O combat.py https://raw.githubusercontent.com/brentp/combat.py/master/combat.py", shell=True)
    from combat import combat
    batch = metadata.loc[batch_var]
    import time
    t = time.time()
    data_fit = combat(df, batch, None, None)
    if (verbose):
        print("%.2f seconds\n" % (time.time() - t))
    data_fit = pd.DataFrame(data_fit.values.tolist(), index=list(data_fit.index), columns=list(data_fit.columns))
    if (verbose):
        print(data_fit.iloc[:5, :5])
    data_fit.to_csv(path+"combat_adjusted_data"+("_"+batch_var if (batch_var != "Batch") else "")+".csv")
    if (verbose):
        print("")
    sb.call("rm -rf *.pyc", shell=True)
    return data_fit

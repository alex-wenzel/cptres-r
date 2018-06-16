"""
This is a simple driver script for RGES based
on the pyRGES.ipynb notebook
"""

import json
from multiprocessing import Pool
import sklearn
from scipy.stats import binom_test
import sys
sys.path.insert(0, '')

from RGES.RGES.DiffEx import DiffEx
from RGES.RGES.L1KGCT import L1KGCTX
from RGES.RGES.Score import score

PHENOTYPE_PATH = "RGES/deep_untreated_markers.tsv"
DRUG_PROFILE_PATH = "LINCS_SIGS/GSE70138_2017-03-06_landmarks_ranked_n118050x972.gctx"

DE = DiffEx(PHENOTYPE_PATH)
DE.data.rename(columns={"avg_logFC":"log2FoldChange"}, inplace=True)

LINCS = L1KGCTX(DRUG_PROFILE_PATH)

OUTPATH = "RGES/deep_untreated_scores.json"
PERMS_PATH = "RGES/deep_untreated_permutations.json"

def mt_score_CHILD(signame):
    """Returns the RGES score for signame based on DE and LINCS"""
    return ((signame, score(DE, LINCS, signame)))

def mt_score(procs):
    """Returns a dictionary of {drug_profile: score}"""
    p = Pool(processes=procs)
    res = p.map(mt_score_CHILD, list(LINCS.data))
    p.close()
    p.join()
    return {r[0]: r[1] for r in res}

def shuffle_sigs():
    LINCS.data = LINCS.data.apply(sklearn.utils.shuffle, axis=0)
    LINCS.data.index = map(str, sorted(map(int, list(LINCS.data.index))))

def stop_permuting(true_score, perm_vals):
    if true_score < 0:
        n_more_extreme = len([s for s in perm_vals if s < true_score])
    else:
        n_more_extreme = len([s for s in perm_vals if s > true_score])
    return binom_test(n_more_extreme, len(perm_vals)+1, 0.05, alternative='greater')*2 < 0.05

PERMUTATIONS = 100000
PROCESSES = 16

try:
    true_scores = json.loads(open(OUTPATH).read())
except FileNotFoundError:
    print("Calculating true scores...")
    true_scores = mt_score(PROCESSES)
    open(OUTPATH, 'w').write(json.dumps(true_scores))

try:
    perms_d = json.loads(open(PERMS_PATH).read())
    to_del = []
    for signame in list(LINCS.data):
        if stop_permuting(true_scores[signame], perms_d[signame]):
            to_del.append(signame)
    LINCS.data.drop(to_del, axis=1, inplace=True)
except FileNotFoundError:
    perms_d = {signame: [] for signame in list(LINCS.data)}

for i in range(1, PERMUTATIONS+1):
    if len(list(LINCS.data)) == 0:
        break
    print("Starting permutation "+str(i)+" with "+str(len(list(LINCS.data)))+" profiles remaining...")
    print("\tShuffling drug signature rankings...")
    shuffle_sigs()
    print("\tCalculating RGES for all profiles...")
    p_res = mt_score(PROCESSES)
    to_del = []
    for signame in list(LINCS.data):
        perms_d[signame].append(p_res[signame])
        if stop_permuting(true_scores[signame], perms_d[signame]):
            to_del.append(signame)
    LINCS.data.drop(to_del, axis=1, inplace=True)
    open(PERMS_PATH, 'w').write(json.dumps(perms_d))
open(PERMS_PATH, 'w').write(json.dumps(perms_d))

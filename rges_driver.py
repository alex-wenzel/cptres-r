"""
This is a simple driver script for RGES based
on the pyRGES.ipynb notebook
"""

import json
from multiprocessing import Pool
import sklearn
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
print(LINCS.data)  #Debug
sys.exit(1)  #Debug

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

PERMUTATIONS = 100
PROCESSES = 16

print("Calculating true scores...")
true_scores = mt_score(PROCESSES)
open(OUTPATH, 'w').write(json.dumps(true_scores))

perms_d = {signame: [] for signame in list(LINCS.data)}

for i in range(PERMUTATIONS):
    print("Starting permutation "+str(i))
    print("\tShuffling drug signature rankings...")
    shuffle_sigs()
    print("\tCalculating RGES for all profiles...")
    p_res = mt_score(PROCESSES)
    for signame in p_res.keys():
        perms_d[signame].append(p_res[signame])
    open(PERMS_PATH, 'w').write(json.dumps(perms_d))
open(PERMS_PATH, 'w').write(json.dumps(perms_d))

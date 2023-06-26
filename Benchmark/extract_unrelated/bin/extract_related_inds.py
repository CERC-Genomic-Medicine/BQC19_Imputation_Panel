import pandas as pd

### Should be modified:)) 
### Test example
related_inds = pd.read_csv("path to/plink2.kin0", sep = "\t")
IID1 = list(related_inds['#IID1'])
IID2 = list(related_inds['IID2'])
IID = IID1 + IID2
IID = list(set(IID))
dict = {'' : IID}
res_df = pd.DataFrame(dict)
res_df.to_csv("./related_inds_based_on_make_rel.txt", index = False, header = False)
related_inds = pd.read_csv("./related_inds_based_on_make_rel.txt", header = None)
all_inds = pd.read_csv("./samples.txt", header = None)
all_inds_list = list(all_inds[0])
related_inds_list = list(related_inds[0])
unrelated = [ x for x in all_inds_list if (x not in related_inds_list)]
unrelated_dict = {'':unrelated}
unrelated_df = pd.DataFrame(unrelated_dict)
unrelated.to_csv('./unrelated_individuals_based_on_new_rel_model.txt', sep = "\t", header = None, index = None)
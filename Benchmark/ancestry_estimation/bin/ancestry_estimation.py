import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn import metrics

argparser = argparse.ArgumentParser(description = 'Ancestry Estimation')
argparser.add_argument('-r', '--reference_PCs', required = True, help = 'Reference PCs')
argparser.add_argument('-s', '--study_PCs', required = True, help = 'Study PCs')
argparser.add_argument('-rl', '--reference_labels', required = True, help = 'Reference ancestry labels')
argparser.add_argument('-np', '--number_PCs', required = True, help = 'number of PCs')

ref_PC = pd.read_csv(args.reference_PCs, sep='\t')
study_PC = pd.read_csv(args.study_PCs, sep='\t')
ref_ancestry = pd.read_csv(args.ref_ancestry_labels)
    
ref = ref_PC.merge(ref_ancestry, left_on='indivID', right_on='ID', how='right')
feature_list = ["PC"+str(i) for i in range(1,args.number_PCs)]
labels = np.array(ref['genetic_region'])
features = ref[feature_list]
features = np.array(features)
    

train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.20, random_state = 42)

#random forest classifier
random_forest_clf = RandomForestClassifier(n_estimators=1000, random_state=42)
random_forest_clf.fit(train_features, train_labels)
test_pred=random_forest_clf.predict(test_features)

study_PC['random_forest_predicted_ancestry'] = random_forest_clf.predict(study_PC[feature_list].values)
rf_CVSs_mean = np.mean(cross_val_score(random_forest_clf, train_features, train_labels, cv=5))
rf_accuracy = metrics.accuracy_score(test_labels, test_pred)
# confusion matrix
cm = metrics.confusion_matrix(test_labels, test_pred, labels=random_forest_clf.classes_)
disp = metrics.ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=random_forest_clf.classes_)
disp.plot()
plt.savefig('./confusion_matrix_random_forest.png', bbox_inches='tight')
plt.clf()

#k-nearest neighbours classifier
knn_clf = KNeighborsClassifier(n_neighbors=10)
knn_clf.fit(train_features, train_labels)
test_pred=knn_clf.predict(test_features)

study_PC['k_nearest_neighbours_predicted_ancestry'] = knn_clf.predict(study_PC[feature_list].values)
study_PC[['indivID','k_nearest_neighbours_predicted_ancestry','random_forest_predicted_ancestry']].to_csv('./predicted_ancestry.txt', index=False, sep='\t')
knn_CVSs_mean = np.mean(cross_val_score(knn_clf, train_features, train_labels, cv=5))
knn_accuracy = metrics.accuracy_score(test_labels, test_pred)

# confusion matrix
cm = metrics.confusion_matrix(test_labels, test_pred, labels=knn_clf.classes_)
disp = metrics.ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=knn_clf.classes_)
disp.plot()
plt.savefig('./confusion_matrix_rknn.png', bbox_inches='tight')
plt.clf()



###### 3D PCA plot

fig = plt.figure(figsize = (15,8))
ax = fig.add_subplot(111, projection='3d')
ref_PCA['color'] = ref_PCA.genetic_region.map({'EUR':'orange', 'EAS':'hotpink', 'CSA':'green', 'AMR':'blue', 'AFR':'yellow', 'MID': 'red', 'OCE':'brown'})

ax.scatter(ref_PC.PC1, ref_PCA.PC2, ref_PCA.PC3, s=20,  c=ref_PCA.color, alpha = 0.3)
ax.scatter(study_PC.PC1, study_PC.PC2, study_PC.PC3, s=20,  c='black', alpha = 0.8, marker='*')

matplotlib.rcParams.update({'font.size': 10})
ax.set_xlabel('PC1', fontsize=13, fontweight='bold', fontfamily='serif')
ax.set_ylabel('PC2', fontsize=13, fontweight='bold', fontfamily='serif')
ax.set_zlabel('PC3', fontsize=13, fontweight='bold', fontfamily='serif')
e = mpatches.Patch(color='orange', label="EUR")
ea = mpatches.Patch(color='hotpink', label="EAS")
cs = mpatches.Patch(color='green', label="CSA" )
am = mpatches.Patch(color='blue', label="AMR" )
af = mpatches.Patch(color='yellow', label="AFR" )
mi = mpatches.Patch(color='red', label="MID" )
oc = mpatches.Patch(color='brown', label="OCE" )


plt.legend(handles=[e, ea, cs, am, af, mi, oc], loc ='lower left', fontsize = 7)
fig.savefig('./PCA3D.png', bbox_inches='tight')
plt.clf()

###### PC1 vc PC2 plot

fig = plt.figure(figsize = (8,8), dpi=100)
plt.scatter(ref_PCA.PC1, ref_PCA.PC2,  s=10,  c=ref_PCA.color, alpha = 0.4)
plt.scatter(study_PC.PC1, study_PC.PC2,  s=10,  c='black', marker='*', alpha = 0.5)


e = mpatches.Patch(color='orange', label="EUR")
ea = mpatches.Patch(color='hotpink', label="EAS")
cs = mpatches.Patch(color='green', label="CSA" )
am = mpatches.Patch(color='blue', label="AMR" )
af = mpatches.Patch(color='yellow', label="AFR" )
mi = mpatches.Patch(color='red', label="MID" )
oc = mpatches.Patch(color='brown', label="OCE" )

plt.grid(ls='--', lw=2, alpha=0.1)
plt.xticks(fontsize=12, fontweight='bold', fontfamily='serif')
plt.yticks(fontsize=12, fontweight='bold', fontfamily='serif')

plt.ylabel('PC1', fontsize=20, fontweight='bold', fontfamily='serif')
plt.xlabel('PC2', fontsize=20, fontweight='bold', fontfamily='serif')
plt.legend(handles=[e, ea, cs, am, af, mi, oc], loc ='lower left', fontsize = 7)
fig.savefig('./PC1_PC2.png', bbox_inches='tight')
plt.clf()


##### PC1 vs PC3

fig = plt.figure(figsize = (8,8), dpi=100)
plt.scatter(ref_PCA.PC1, ref_PCA.PC3,  s=10,  c=ref_PCA.color, alpha = 0.4)
plt.scatter(Bstudy_PC.PC1, study_PC.PC3,  s=10,  c='black', marker='*', alpha = 0.5)


e = mpatches.Patch(color='orange', label="EUR")
ea = mpatches.Patch(color='hotpink', label="EAS")
cs = mpatches.Patch(color='green', label="CSA" )
am = mpatches.Patch(color='blue', label="AMR" )
af = mpatches.Patch(color='yellow', label="AFR" )
mi = mpatches.Patch(color='red', label="MID" )
oc = mpatches.Patch(color='brown', label="OCE" )

plt.grid(ls='--', lw=2, alpha=0.1)
plt.xticks(fontsize=12, fontweight='bold', fontfamily='serif')
plt.yticks(fontsize=12, fontweight='bold', fontfamily='serif')

plt.ylabel('PC1', fontsize=20, fontweight='bold', fontfamily='serif')
plt.xlabel('PC3', fontsize=20, fontweight='bold', fontfamily='serif')
plt.legend(handles=[e, ea, cs, am, af, mi, oc], loc ='lower left', fontsize = 7)
fig.savefig('./PC1_PC3.png', bbox_inches='tight')
plt.clf()

##### PC2 vs PC3

fig = plt.figure(figsize = (8,8), dpi=100)
plt.scatter(ref_PCA.PC2, ref_PCA.PC3,  s=10,  c=ref_PCA.color, alpha = 0.4)
plt.scatter(Bstudy_PC.PC2, study_PC.PC3,  s=10,  c='black', marker='*', alpha = 0.5)


e = mpatches.Patch(color='orange', label="EUR")
ea = mpatches.Patch(color='hotpink', label="EAS")
cs = mpatches.Patch(color='green', label="CSA" )
am = mpatches.Patch(color='blue', label="AMR" )
af = mpatches.Patch(color='yellow', label="AFR" )
mi = mpatches.Patch(color='red', label="MID" )
oc = mpatches.Patch(color='brown', label="OCE" )

plt.grid(ls='--', lw=2, alpha=0.1)
plt.xticks(fontsize=12, fontweight='bold', fontfamily='serif')
plt.yticks(fontsize=12, fontweight='bold', fontfamily='serif')

plt.ylabel('PC2', fontsize=20, fontweight='bold', fontfamily='serif')
plt.xlabel('PC3', fontsize=20, fontweight='bold', fontfamily='serif')
plt.legend(handles=[e, ea, cs, am, af, mi, oc], loc ='lower left', fontsize = 7)
fig.savefig('./PC2_PC3.png', bbox_inches='tight')
plt.clf()

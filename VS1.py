from openbabel import pybel
from pandas import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, auc
#from matplotlib import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
#%matplotlib inline

scores = []
for mol in pybel.readfile('sdf','gnina_scored_vinardo.sdf.gz'):
    scores.append({'title': mol.title, 
                   'CNNscore': float(mol.data['CNNscore']), 
                   'CNNaffinity': float(mol.data['CNNaffinity']),
                   'Vinardo': float(mol.data['minimizedAffinity'])})
scores = pd.DataFrame(scores)  
scores['label'] = scores.title.str.contains('active')

scores

plt.figure(dpi=150)
plt.plot([0,1],[0,1],'k--',alpha=0.5,linewidth=1)
fpr,tpr,_ = roc_curve(scores.label,-scores.Vinardo)
plt.plot(fpr,tpr,label="Vinardo (AUC = %.2f)"%auc(fpr,tpr))
fpr,tpr,_ = roc_curve(scores.label,scores.CNNaffinity)
plt.plot(fpr,tpr,label="CNNaffinity (AUC = %.2f)"%auc(fpr,tpr))
fpr,tpr,_ = roc_curve(scores.label, scores.CNNaffinity.rank() + (-scores.Vinardo).rank())
plt.plot(fpr,tpr,label="Consensus (AUC = %.2f)"%auc(fpr,tpr))
plt.legend(loc='lower right')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.gca().set_aspect('equal')
plt.savefig('AUC.png')

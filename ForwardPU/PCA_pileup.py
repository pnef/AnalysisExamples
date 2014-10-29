#!/usr/bin/env ipython
""" 
Principle Component Analysis of calorimeter-based forward pileup jet feature space. 

"""
import sys, os

sys.path.append(os.path.abspath("../base/"))
import matplotlib.pyplot as plt
import analysisbase as ana
from analysisbase import wait
from rootpy.io import root_open
import array

import ROOT as R
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc
import numpy as np

# Main ------------------------------------------------------------------------
if __name__ == "__main__":
    tfile = "/Users/pascal/Data/ATLAS/ForwardPileup/20141010.14.24_ForwardPileupJets.Zmumu_PowhegPythia8_MC12_COMMON.forward.root"

    myfile = root_open(tfile)
    tree = myfile.EventTree
    
    analysis = ana.AnalysisBase()

    Jpt           = array.array('f',[0])
    Jeta          = array.array('f',[0])
    Jm            = array.array('f',[0])
    isHardScatter = array.array('i',[0])
    isPileup      = array.array('i',[0])
    Jwidth        = array.array('f',[0])
    delRsqr       = array.array('f',[0])
    delRSkewness  = array.array('f',[0])
    delR_34       = array.array('f',[0])
    delR_01       = array.array('f',[0])


    # set branch address
    tree.SetBranchAddress("Jpt",          Jpt);
    tree.SetBranchAddress("Jeta",         Jeta);
    tree.SetBranchAddress("Jm",           Jm);
    tree.SetBranchAddress("isHardScatter",isHardScatter);
    tree.SetBranchAddress("isPileup",     isPileup);
    tree.SetBranchAddress("Jwidth",       Jwidth);
    tree.SetBranchAddress("delRsqr",      delRsqr);
    tree.SetBranchAddress("delRSkewness", delRSkewness);
    tree.SetBranchAddress("delR_34",      delR_34);
    tree.SetBranchAddress("delR_01",      delR_01);


    # X will be the data matrix: column = different events, row = different variable
    Xsig      = []
    Xbkg      = []

    for entry in range(tree.GetEntries()):
        if entry > 10000:
            break
        tree.GetEntry(entry)
        if(np.abs(Jeta[0])<2.4 or (Jpt[0]<20 and Jpt[0]>40) or Jwidth[0]==0):
            continue
        if(isPileup[0]):
            Xbkg.append([Jwidth[0], delRsqr[0], delRSkewness[0], delR_34[0], delR_01[0]])
        if(isHardScatter[0]):
            Xsig.append([Jwidth[0], delRsqr[0], delRSkewness[0], delR_34[0], delR_01[0]]) 


    length = min(len(Xbkg), len(Xsig))

    X         = Xsig[:length]+Xbkg[:length]
    SignalBit = [0 for i in range(length)] + [1 for i in range(length)]


    # scale variables to mean zero and unit variance
    X_scaled = preprocessing.scale(X)

    pca = PCA(n_components=3)
    X_r = pca.fit(X_scaled).transform(X_scaled)
    print(pca.explained_variance_ratio_)
    

    width_sig =[ii[0] for ii in X[:length]]
    width_bkg =[ii[0] for ii in X[length:]]
    PC1_sig   =[ii[0] for ii in X_r[:length]]
    PC1_bkg   =[ii[0] for ii in X_r[length:]]
    PC2_sig   =[ii[1] for ii in X_r[:length]]
    PC2_bkg   =[ii[1] for ii in X_r[length:]]
    PC3_sig   =[ii[2] for ii in X_r[:length]]
    PC3_bkg   =[ii[2] for ii in X_r[length:]]

    plt.figure()
    plt.scatter(PC1_sig, width_sig, c='b', marker='o', alpha=0.5, label='sig')
    plt.scatter(PC1_bkg, width_bkg, c='r', marker='x', alpha=0.5, label='bkg')
    plt.xlabel('PC1')
    plt.ylabel('width')
    plt.legend()
    plt.show()

    # Construct projective likelihod
    analysis.gethistofromarray("widthsig", width_sig , bins=100, range=(-0.1,0.5))
    analysis.gethistofromarray("widthbkg", width_bkg , bins=100, range=(-0.1,0.5))
    analysis.gethistofromarray("PC1sig",   PC1_sig   , bins=100, range=(-8,6))
    analysis.gethistofromarray("PC1bkg",   PC1_bkg   , bins=100, range=(-8,6))
    analysis.gethistofromarray("PC2sig",   PC2_sig   , bins=100, range=(-3,6))
    analysis.gethistofromarray("PC2bkg",   PC2_bkg   , bins=100, range=(-3,6))
    analysis.gethistofromarray("PC3sig",   PC3_sig   , bins=100, range=(-5,5))
    analysis.gethistofromarray("PC3bkg",   PC3_bkg   , bins=100, range=(-5,5))

    plt.figure()
    plt.bar(analysis.H["PC1sig"].bins[:-1], analysis.H["PC1sig"].hist, width=(analysis.H["PC1sig"].bins[0]-analysis.H["PC1sig"].bins[1]), color='b', alpha=0.5, label="signal")     # gets all but last items in list
    plt.bar(analysis.H["PC1bkg"].bins[:-1], analysis.H["PC1bkg"].hist, width=(analysis.H["PC1bkg"].bins[0]-analysis.H["PC1bkg"].bins[1]), color='r', alpha=0.5, label="background") # gets all but last items in list
    plt.xlabel('PC1')
    plt.ylabel('Entries')
    plt.legend()
    plt.show()
    
    analysis.makelikelihood("widthsig", "widthbkg", "Widthlikelihood")
    analysis.makelikelihood("PC1sig"  , "PC1bkg"  , "PC1likelihood")
    analysis.makelikelihood("PC2sig"  , "PC2bkg"  , "PC2likelihood")
    analysis.makelikelihood("PC3sig"  , "PC3bkg"  , "PC3likelihood")
#    plt.figure()
#    plt.bar(analysis.H["PC1likelihood"].bins[:-1], analysis.H["PC1likelihood"].hist, width=(analysis.H["PC1likelihood"].bins[0]-analysis.H["PC1likelihood"].bins[1]), color='g', alpha=0.5) # gets all but last items in list
#    plt.show()

    LL_sig       = analysis.evaluate_likelihood(zip(PC1_sig,PC2_sig,PC3_sig), ("PC1likelihood","PC2likelihood","PC3likelihood") )
    LL_bkg       = analysis.evaluate_likelihood(zip(PC1_bkg,PC2_bkg,PC3_bkg), ("PC1likelihood","PC2likelihood","PC3likelihood") )
    LL_PC1_sig   = analysis.evaluate_likelihood(zip(PC1_sig),                 ("PC1likelihood",) )
    LL_PC1_bkg   = analysis.evaluate_likelihood(zip(PC1_bkg),                 ("PC1likelihood",) )
    LL_PC2_sig   = analysis.evaluate_likelihood(zip(PC2_sig),                 ("PC2likelihood",) )
    LL_PC2_bkg   = analysis.evaluate_likelihood(zip(PC2_bkg),                 ("PC2likelihood",) )
    LL_PC3_sig   = analysis.evaluate_likelihood(zip(PC3_sig),                 ("PC3likelihood",) )
    LL_PC3_bkg   = analysis.evaluate_likelihood(zip(PC3_bkg),                 ("PC3likelihood",) )
    LL_width_sig = analysis.evaluate_likelihood(zip(width_sig),               ("Widthlikelihood",) )
    LL_width_bkg = analysis.evaluate_likelihood(zip(width_bkg),               ("Widthlikelihood",) )
    plt.figure()
    plt.hist(LL_sig, bins=100, range=(0,1), label='sig', color='b', alpha=0.5,normed=True)
    plt.hist(LL_bkg, bins=100, range=(0,1), label='bkg', color='r', alpha=0.5,normed=True)
    plt.xlabel('Likelihood')
    plt.ylabel('Norm Entries')
    plt.show()



    # make ROC curve for first PC
    fpr_L,     tpr_L,     _ = roc_curve(np.append(np.zeros(len(LL_sig))       ,np.ones(len(LL_bkg))),       np.append(LL_sig, LL_bkg))
    fpr_PC1,   tpr_PC1,   _ = roc_curve(np.append(np.zeros(len(LL_PC1_sig))   ,np.ones(len(LL_PC1_bkg))),   np.append(LL_PC1_sig,LL_PC1_bkg))
    fpr_PC2,   tpr_PC2,   _ = roc_curve(np.append(np.zeros(len(LL_PC2_sig))   ,np.ones(len(LL_PC2_bkg))),   np.append(LL_PC2_sig,LL_PC2_bkg))
    fpr_PC3,   tpr_PC3,   _ = roc_curve(np.append(np.zeros(len(LL_PC3_sig))   ,np.ones(len(LL_PC3_bkg))),   np.append(LL_PC3_sig,LL_PC3_bkg))
    fpr_width, tpr_width, _ = roc_curve(np.append(np.zeros(len(LL_width_sig)) ,np.ones(len(LL_width_bkg))), np.append(LL_width_sig,LL_width_bkg))

    plt.plot(fpr_L,      tpr_L,     label='likelihood(PC1,2,3)', color='black',   linestyle='-', )
    plt.plot(fpr_PC1,    tpr_PC1,   label='PC1',                 color='b',       linestyle='-', )
    plt.plot(fpr_PC2,    tpr_PC2,   label='PC2',                 color='g',       linestyle='-', )
    plt.plot(fpr_PC3,    tpr_PC3,   label='PC3',                 color='r',       linestyle='-', )
    plt.plot(fpr_width,  tpr_width, label='width',               color='magenta', linestyle='--', )
    plt.xlabel('Background Efficiency')
    plt.ylabel('Signal Efficiency')
    plt.legend(loc='upper left',frameon=False)
    plt.title('ROC curve')
    plt.show()
    

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math as mt
import numpy as np
import pandas as pd
import ROOT
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# In[3]:


ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/Delphes/libDelphes.so")


# In[4]:


ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/Delphes/libDelphes")
#ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/ExRootAnalysis/libExRootAnalysis.so")

ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/ExRootAnalysis/ExRootAnalysis/ExRootTreeReader.h"')
ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/DelphesClasses.h"')
#ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/DelphesModule.h"')


# # Gluons Jets image construction

# In[5]:


data_path = '/Users/bcoder/Bonesini_qq_gg/roots/gluons1.root'


# In[6]:


f = ROOT.TFile.Open(data_path)
chain = ROOT.TChain("Delphes")
chain.Add(data_path)


# In[8]:


treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()


# In[9]:


numberOfEntries


# In[10]:


branchJet = treeReader.UseBranch("Jet")
branchFatJet = treeReader.UseBranch("FatJet")
branchMuon = treeReader.UseBranch("Muon")
branchElectron = treeReader.UseBranch("Electron")
branchMET = treeReader.UseBranch("MissingET")
branchPhoton = treeReader.UseBranch("Photon")
branchParticle = treeReader.UseBranch("Particle");
branchEFlowTrack = treeReader.UseBranch("EFlowTrack");
branchEFlowTower = treeReader.UseBranch("EFlowTower");
branchEFlowMuon = treeReader.UseBranch("EFlowMuon");
branchEFlowPhoton = treeReader.UseBranch("EFlowPhoton");
branchEFlowNeutralHadron = treeReader.UseBranch("EFlowNeutralHadron");
branchTower = treeReader.UseBranch("Tower")


# # Parton jet match

# In[11]:


from tqdm import tqdm


# In[12]:


def Delta_R(parton, jet):
    return mt.sqrt((jet.Eta-parton.Eta)**2+(jet.Phi-parton.Phi)**2)


# In[13]:


def high_pt_jet_match(partons, jets):
    
    #selecting two highest PT jets
    jets = sorted(jets, key= lambda x: x.PT)
    jets = jets[0:2]
    
    p_idx = list(np.arange(0, len(partons), 1))
    j_idx = list(np.arange(0, len(jets), 1))
    parton_jet = []
    
    #checking for DeltaR compatibility
    for i in range(len(p_idx)):
        parton_jet_R = []
        for p in p_idx:
            parton = partons[p]
            for j in j_idx:
                jet = jets[j]
                parton_jet_R.append([p, j, Delta_R(parton,jet)])
                
        parton_jet_R = sorted(parton_jet_R, key= lambda x:x[2])
     
        if (parton_jet_R[0][2] < .7 ):
            parton_jet.append([partons[parton_jet_R[0][0]], jets[parton_jet_R[0][1]]])
            if len(p_idx) != 1:
                p_idx.pop(parton_jet_R[0][0])
                j_idx.pop(parton_jet_R[0][1])
            
    if len(parton_jet) == len(partons):
        return parton_jet
        
    else:
        return False


# In[14]:


images_quark = []
images_qcd = []

for entry in tqdm(range(0, numberOfEntries)):
    
    
    treeReader.ReadEntry(entry)
    number_of_jets = branchJet.GetEntries()
    #print("Number of jets: ", number_of_jets)
    
    #building partons
    partons = []
    for p in branchParticle:
        if (p.Status == 23) and (abs(p.PID) == 21):
            partons.append(p)
    
    if branchJet.GetEntries() >= 2:
        
        jets = [i for i in branchJet]
        
        #checking if match max DeltaR with highest Pt jets
        association = high_pt_jet_match(partons, jets)
        
        #if no complete match between highest PT and partons
        #move on
        if not association:
            continue
            
        associated_jets = [row[1] for row in association]
                    
        
        #building jets images
        for j in range(branchJet.GetEntries()):

            jet_im = []
            particles = branchJet.At(j).Constituents
            jet_Eta = branchJet.At(j).Eta
            jet_Phi = branchJet.At(j).Phi
            #print("Num of particles: ", particles.GetEntries())

            for p in range(particles.GetEntries()):

                if(type(particles.At(p)) == ROOT.GenParticle ): 
                    four_mom_p = particles.At(p).P4()
                    pt = four_mom_p.Pt()
                    jet_im.append([jet_Eta-particles.At(p).Eta, jet_Phi-particles.At(p).Phi, pt ])

                if(type(particles.At(p)) == ROOT.Track ): 
                    four_mom_p = particles.At(p).P4()
                    pt = four_mom_p.Pt()
                    jet_im.append([jet_Eta-particles.At(p).Eta, jet_Phi-particles.At(p).Phi, pt ])

                if(type(particles.At(p)) == ROOT.Tower ): 
                    four_mom_p = particles.At(p).P4()
                    pt = four_mom_p.Pt()
                    jet_im.append([jet_Eta-particles.At(p).Eta, jet_Phi-particles.At(p).Phi, pt ])
            
            
            if branchJet.At(j) in associated_jets:
                images_qcd.append(jet_im)
            else: 
                images_quark.append(jet_im)


# In[15]:


print(len(images_quark))
print(len(images_qcd))


# In[16]:


etas = np.linspace(-0.8, 0.8, 100)
phis = np.linspace(-0.8, 0.8, 100)

im_to_plt = []
count = 0
for jet_im in images_qcd:
    if count%100 == 0:
        print("It: ", count)
        
    im = np.zeros((100,100))
    for i in range(100-1):
        for j in range(100-1):
            eta_inf = etas[i]
            eta_sup = etas[i+1]
            phi_inf = phis[j]
            phi_sup = phis[j+1]
            for el in jet_im:
                if (el[0] > eta_inf) & (el[0] < eta_sup) & (el[1] > phi_inf) & (el[1] < phi_sup):
                    im[i,j] += el[2]
    im_to_plt.append(im)
    count += 1
    
im_to_plt = np.array(im_to_plt)
im_to_plt.shape      
        


# In[18]:


plt.figure(figsize=(15,15))
SUM_Image = np.sum(im_to_plt, axis = 0)
plt.imshow(SUM_Image/float(im_to_plt.shape[0]), origin='lower',norm=LogNorm(vmin=0.01))
plt.colorbar()
plt.xlabel("Δη", fontsize=20)
plt.ylabel("ΔΦ", fontsize=20)
plt.title("Pt Density Gluon jets", fontsize=20)
plt.savefig("../graphs/Pt_density_gluons_jets.pdf")


# In[19]:


plt.figure(figsize=(15,15))
plt.imshow(im_to_plt[7], origin='lower',norm=LogNorm(vmin=0.01))
plt.colorbar()
plt.xlabel("Δη", fontsize=20)
plt.ylabel("ΔΦ", fontsize=20)
plt.title("Single Gluon Image", fontsize=20)
plt.savefig("../graphs/Pt_density_single_gluons.pdf")


# In[20]:


np.save('./Gluons_images_NEW.npy', im_to_plt)


# In[ ]:





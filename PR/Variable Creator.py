#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math as mt
import numpy as np
import pandas as pd
import ROOT
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# # Creating jets FeedForward vectors

# In[3]:


ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/Delphes/libDelphes.so")


# In[4]:


ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/Delphes/libDelphes")
#ROOT.gSystem.Load("/Users/bcoder/MG5_aMC_v2_6_6/ExRootAnalysis/libExRootAnalysis.so")

ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/ExRootAnalysis/ExRootAnalysis/ExRootTreeReader.h"')
ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/DelphesClasses.h"')
#ROOT.gInterpreter.Declare('#include "/Users/bcoder/MG5_aMC_v2_6_6/Delphes/classes/DelphesModule.h"')


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


g_array = []

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
        
        for j in associated_jets:
            j_att = []
            j_att.append(j.PT) #T mom
            j_att.append(j.Eta) #Eta
            j_att.append(j.Phi) #Phi
            j_att.append(branchJet.GetEntries()) #number of jets
            j_att.append(j.Mass) #invariant mass
            j_att.append(j.Mass/j.PT) #rateo
            j_att.append(j.NCharged) #charged multiplicity
            j_att.append(j.NNeutrals) #neutral multiplicity
            j_att.append(j.EhadOverEem) #Rateo of energy in ECAL and HCAL
            j_att.append(j.MeanSqDeltaR)
            j_att.append(j.PTD)
            j_att.append(j.Tau[0]) #starting with N-Subjettiness
            if j.Tau[0] > 0:
                j_att.append(j.Tau[1]/j.Tau[0])
            else:
                j_att.append(0)
                
            j_att.append(j.Tau[1])
            if j.Tau[1] > 0:
                j_att.append(j.Tau[2]/j.Tau[1])
            else:
                j_att.append(0)
                
            j_att.append(j.Tau[2])
            if j.Tau[2] > 0:
                j_att.append(j.Tau[3]/j.Tau[2])
            else:
                j_att.append(0)
            
            j_att.append(j.Tau[3])
            if j.Tau[3] > 0:
                j_att.append(j.Tau[4]/j.Tau[3])
            else:
                j_att.append(0)
            
            
            g_array.append(j_att)


# In[17]:


g_array = np.array(g_array)
print(g_array.shape)
np.save("./g_var_arr.npy", g_array)


# In[18]:


data_path = '/Users/bcoder/Bonesini_qq_gg/roots/quarks1.root'

f = ROOT.TFile.Open(data_path)
chain = ROOT.TChain("Delphes")
chain.Add(data_path)


# In[19]:


treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()


# In[20]:


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


# In[21]:


q_array = []

for entry in tqdm(range(0, numberOfEntries)):
    
    
    treeReader.ReadEntry(entry)
    number_of_jets = branchJet.GetEntries()
    
    #building partons
    partons = []
    for p in branchParticle:
        if (p.Status == 23) and (abs(p.PID) <= 8):
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
        
        for j in associated_jets:
            j_att = []
            j_att.append(j.PT) #T mom
            j_att.append(j.Eta) #Eta
            j_att.append(j.Phi) #Phi
            j_att.append(branchJet.GetEntries()) #number of jets
            j_att.append(j.Mass) #invariant mass
            j_att.append(j.Mass/j.PT) #rateo
            j_att.append(j.NCharged) #charged multiplicity
            j_att.append(j.NNeutrals) #neutral multiplicity
            j_att.append(j.EhadOverEem) #Rateo of energy in ECAL and HCAL
            j_att.append(j.MeanSqDeltaR)
            j_att.append(j.PTD)
            j_att.append(j.Tau[0]) #starting with N-Subjettiness
            if j.Tau[0] > 0:
                j_att.append(j.Tau[1]/j.Tau[0])
            else:
                j_att.append(0)
                
            j_att.append(j.Tau[1])
            if j.Tau[1] > 0:
                j_att.append(j.Tau[2]/j.Tau[1])
            else:
                j_att.append(0)
                
            j_att.append(j.Tau[2])
            if j.Tau[2] > 0:
                j_att.append(j.Tau[3]/j.Tau[2])
            else:
                j_att.append(0)
            
            j_att.append(j.Tau[3])
            if j.Tau[3] > 0:
                j_att.append(j.Tau[4]/j.Tau[3])
            else:
                j_att.append(0)
            
            q_array.append(j_att)
    


# In[22]:


q_array_2 = np.array(q_array)


# In[23]:


q_array_2 = q_array[:g_array.shape[0]]


# In[24]:


q_array_2 = np.array(q_array_2)
q_array_2.shape


# In[28]:


np.save("./q_var_arr.npy", np.array(q_array))


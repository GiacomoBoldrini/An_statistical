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


# # Quark Jets image construction

# In[5]:


data_path = '/Users/bcoder/Bonesini_qq_gg/roots/quarks1.root'


# In[6]:


f = ROOT.TFile.Open(data_path)


# In[7]:


chain = ROOT.TChain("Delphes")
chain.Add(data_path)


# In[9]:


treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()


# In[10]:


numberOfEntries


# In[11]:


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

# In[12]:


from tqdm import tqdm


# In[13]:


def Delta_R(parton, jet):
    return mt.sqrt((jet.Eta-parton.Eta)**2+(jet.Phi-parton.Phi)**2)


# In[14]:


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
         

    


# In[16]:


hard_p = []
hard_p_id = []
count_qcd_jets = 0
event_counter = 0
for entry in range(0, 5000):
    
    
    treeReader.ReadEntry(entry)
    
    #building partons
    partons = []
    for p in branchParticle:
        if (p.Status == 23) and (abs(p.PID) <= 8):
            partons.append(p)
        
    #building jets for only events with more than 2 jets

    if branchJet.GetEntries() >= 2:
        
        count_qcd_jets+= branchJet.GetEntries() - len(partons)
        jets = [i for i in branchJet]
        
        high_pt_p = sorted(partons, key= lambda x: x.PT)
        association = high_pt_jet_match(high_pt_p, jets)
    
        if len(association) < 2:
            association = parton_jet_match(partons, jets)
        else:
            #print("Passed")
            event_counter+=1
        


# In[35]:


event_counter


# In[17]:


images_quark = []
images_qcd = []

for entry in tqdm(range(0, numberOfEntries)):
    
    
    treeReader.ReadEntry(entry)
    number_of_jets = branchJet.GetEntries()
    #print("Number of jets: ", number_of_jets)
    
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
                images_quark.append(jet_im)
            else: 
                images_qcd.append(jet_im)


# In[18]:


print(len(images_quark))
print(len(images_qcd))


# # Plotting Total Quark Events and building images as numpy

# In[19]:


etas = np.linspace(-0.8, 0.8, 100)
phis = np.linspace(-0.8, 0.8, 100)

im_to_plt_q = []
for jet_im in images_quark:
        
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
    im_to_plt_q.append(im)
    
im_to_plt_q = np.array(im_to_plt_q)
im_to_plt_q.shape      
np.save('./Quark_images.npy', im_to_plt_q)   


# In[20]:


plt.figure(figsize=(15,15))
SUM_Image = np.sum(im_to_plt_q, axis = 0)
plt.imshow(SUM_Image/float(im_to_plt_q.shape[0]), origin='lower',norm=LogNorm(vmin=0.01))
plt.colorbar()
plt.xlabel("Δη", fontsize=20)
plt.ylabel("ΔΦ", fontsize=20)
plt.title("Pt Density Quarks Jets", fontsize=20)


# In[21]:


np.save('./Quark_images_NEW.npy', im_to_plt_q)


# # ELSE

# In[13]:


def Max_pt_prop(jets_event):
    max_pt = [jet.PT for jet in jets_event]
    ind_max = max_pt.index(max(max_pt))
    
    massimo_jet = jets_event[ind_max]
    
    return massimo_jet.NCharged, massimo_jet.NNeutrals, massimo_jet.Flavor


# In[19]:


signal = []
# Loop over all events
#for entry in range(0, numberOfEntries):
for entry in range(0, 20):
    
    event = np.array([])
    # Load selected branches with data from specified event
    #jet analysis
    treeReader.ReadEntry(entry)
    number_of_jets = branchJet.GetEntries()
    
    if branchJet.GetEntries() > 0:
        jets_pt = []
        jets_mass = []
        jets_eta = []
        jets_phi = [] 
        all_jets = []
        
        for j in range(branchJet.GetEntries()):
            jets_pt.append(branchJet.At(j).PT)
            jets_mass.append(branchJet.At(j).Mass)
            jets_eta.append(branchJet.At(j).Eta)
            jets_phi.append(branchJet.At(j).Phi)
            all_jets.append(branchJet.At(j))
            
        NCh, NNeu, Flava = Max_pt_prop(all_jets)
        jets_pt = sorted(jets_pt, reverse = True)
        jets_mass = sorted(jets_mass, reverse = True)
        jets_eta = sorted(jets_eta, reverse = True)
        jets_phi = sorted(jets_phi, reverse = True)
        
        jet_max_pt = jets_pt[0]
        jets_max_mass = jets_mass[0]
        jets_max_eta = jets_eta[0]
        jets_max_phi = jets_phi[0]
        
    else:
        jet_max_pt = 0
        jets_max_mass = 0
        jets_max_eta = 0
        jets_max_phi = 0
     
    if branchElectron.GetEntries() > 0:
        electron_Pt = branchElectron.At(0).PT
    else:
        electron_Pt = 0
        
    if branchMuon.GetEntries() > 0:
        muon_Pt = branchMuon.At(0).PT
        muon_eta = branchMuon.At(0).Eta
        muon_phi = branchMuon.At(0).Phi
        muon_T = branchMuon.At(0).T
    else:
        muon_Pt = 0
        muon_eta = 0
        muon_phi = 0
        muon_T = 0
        
    if branchMET.GetEntries() > 0:
        met_pt = branchMET.At(0).MET
        met_m = branchMET.At(0).Eta
    else:
        met_pt = 0
        met_m = 0
        
    event = np.r_[event, (number_of_jets, jet_max_pt, jets_max_mass, jets_max_eta, jets_max_phi, electron_Pt, muon_Pt, muon_eta, muon_phi, muon_T, met_pt, met_m, NCh, NNeu, Flava)]
    
    
    #fatjets
    
    number_of_fat_jets = branchFatJet.GetEntries()
    if branchFatJet.GetEntries() > 0:
        fat_jets_pt = []
        fat_jets_mass = []
        
        for j in range(branchFatJet.GetEntries()):
            fat_jets_pt.append(branchFatJet.At(j).PT)
            fat_jets_mass.append(branchFatJet.At(j).Mass)
            
        fat_jets_pt = sorted(fat_jets_pt, reverse = True)
        fat_jets_mass = sorted(fat_jets_mass, reverse = True)
        
        fat_jet_max_pt = fat_jets_pt[0]
        fat_jets_max_mass = fat_jets_mass[0]
        
    else:
        fat_jet_max_pt = 0
        fat_jets_max_mass = 0
        
        
    #photon
    number_of_photon = branchPhoton.GetEntries()
    if branchPhoton.GetEntries() > 0:
        photon_pt = []
        photon_energy = []
        photon_Eta = []
        photon_Phi = []
        
        for p in range(branchPhoton.GetEntries()):
            photon_pt.append(branchPhoton.At(p).PT)
            photon_energy.append(branchPhoton.At(p).E)
            photon_Eta.append(branchPhoton.At(p).Eta)
            photon_Phi.append(branchPhoton.At(p).Phi)
            
        photon_pt = sorted(photon_pt, reverse = True)
        photon_energy = sorted(photon_energy, reverse = True)
        photon_Eta = sorted(photon_Eta, reverse = True)
        photon_Phi = sorted(photon_Phi, reverse = True)
        
        photon_max_pt = photon_pt[0]
        photon_max_energy = photon_energy[0]
        photon_max_eta = photon_Eta[0]
        photon_max_phi = photon_Phi[0]
    
    else:
        photon_max_pt = 0
        photon_max_energy = 0
        photon_max_eta = 0
        photon_max_phi = 0
            
        
    event = np.r_[event, (number_of_fat_jets, fat_jet_max_pt, fat_jets_max_mass,number_of_photon, photon_max_pt, photon_max_energy, photon_max_eta, photon_max_phi)]
    
    signal.append(list(event))


# In[21]:


len(signal[1])


# # BACKGROUND 

# In[14]:


data_path2 = "/Users/boldrinicoder/MADGRAPH5/MG5_aMC_v2_6_4/ttjj_semilep/Events/run_04/tag_1_delphes_events.root"

f_back = ROOT.TFile.Open(data_path2)
chain = ROOT.TChain("Delphes")
chain.Add(data_path2)

treeReader2 = ROOT.ExRootTreeReader(chain)
numberOfEntries_bkg = treeReader2.GetEntries()

branchJet = treeReader2.UseBranch("Jet")
branchFatJet = treeReader2.UseBranch("FatJet")
branchMuon = treeReader2.UseBranch("Muon")
branchElectron = treeReader2.UseBranch("Electron")
branchMET = treeReader2.UseBranch("MissingET")
branchPhoton = treeReader2.UseBranch("Photon")

bkg = []
# Loop over all events
for entry in range(0, numberOfEntries_bkg):
    
    event = np.array([])
    # Load selected branches with data from specified event
    #jet analy2sis
    treeReader2.ReadEntry(entry)
    number_of_jets = branchJet.GetEntries()
    
    if branchJet.GetEntries() > 0:
        jets_pt = []
        jets_mass = []
        jets_eta = []
        jets_phi = []
        all_jets = []
        
        for j in range(branchJet.GetEntries()):
            jets_pt.append(branchJet.At(j).PT)
            jets_mass.append(branchJet.At(j).Mass)
            jets_eta.append(branchJet.At(j).Eta)
            jets_phi.append(branchJet.At(j).Phi)
            all_jets.append(branchJet.At(j))
            
        NCh, NNeu, Flava = Max_pt_prop(all_jets)
        jets_pt = sorted(jets_pt, reverse = True)
        jets_mass = sorted(jets_mass, reverse = True)
        jets_eta = sorted(jets_eta, reverse = True)
        jets_phi = sorted(jets_phi, reverse = True)
        
        jet_max_pt = jets_pt[0]
        jets_max_mass = jets_mass[0]
        jets_max_eta = jets_eta[0]
        jets_max_phi = jets_phi[0]
        
    else:
        jet_max_pt = 0
        jets_max_mass = 0
        jets_max_eta = 0
        jets_max_phi = 0
     
    if branchElectron.GetEntries() > 0:
        electron_Pt = branchElectron.At(0).PT
    else:
        electron_Pt = 0
        
    if branchMuon.GetEntries() > 0:
        muon_Pt = branchMuon.At(0).PT
        muon_eta = branchMuon.At(0).Eta
        muon_phi = branchMuon.At(0).Phi
        muon_T = branchMuon.At(0).T
    else:
        muon_Pt = 0
        muon_eta = 0
        muon_phi = 0
        muon_T = 0
        
    if branchMET.GetEntries() > 0:
        met_pt = branchMET.At(0).MET
        met_m = branchMET.At(0).Eta
    else:
        met_pt = 0
        met_m = 0
        
    event = np.r_[event, (number_of_jets, jet_max_pt, jets_max_mass, jets_max_eta, jets_max_phi, electron_Pt, muon_Pt,muon_eta, muon_phi, muon_T, met_pt, met_m, NCh, NNeu, Flava)]
    
    
    #fatjets
    
    number_of_fat_jets = branchFatJet.GetEntries()
    if branchFatJet.GetEntries() > 0:
        fat_jets_pt = []
        fat_jets_mass = []
        
        for j in range(branchFatJet.GetEntries()):
            fat_jets_pt.append(branchFatJet.At(j).PT)
            fat_jets_mass.append(branchFatJet.At(j).Mass)
            
        fat_jets_pt = sorted(fat_jets_pt, reverse = True)
        fat_jets_mass = sorted(fat_jets_mass, reverse = True)
        
        fat_jet_max_pt = fat_jets_pt[0]
        fat_jets_max_mass = fat_jets_mass[0]
        
    else:
        fat_jet_max_pt = 0
        fat_jets_max_mass = 0
        
    
        #photon
    number_of_photon = branchPhoton.GetEntries()
    if branchPhoton.GetEntries() > 0:
        photon_pt = []
        photon_energy = []
        photon_Eta = []
        photon_Phi = []
        
        for p in range(branchPhoton.GetEntries()):
            photon_pt.append(branchPhoton.At(p).PT)
            photon_energy.append(branchPhoton.At(p).E)
            photon_Eta.append(branchPhoton.At(p).Eta)
            photon_Phi.append(branchPhoton.At(p).Phi)
            
        photon_pt = sorted(photon_pt, reverse = True)
        photon_energy = sorted(photon_energy, reverse = True)
        photon_Eta = sorted(photon_Eta, reverse = True)
        photon_Phi = sorted(photon_Phi, reverse = True)
        
        photon_max_pt = photon_pt[0]
        photon_max_energy = photon_energy[0]
        photon_max_eta = photon_Eta[0]
        photon_max_phi = photon_Phi[0]
    
    else:
        photon_max_pt = 0
        photon_max_energy = 0
        photon_max_eta = 0
        photon_max_phi = 0
            
            
    event = np.r_[event, (number_of_fat_jets, fat_jet_max_pt, fat_jets_max_mass,number_of_photon, photon_max_pt, photon_max_energy, photon_max_eta, photon_max_phi)]
    
    bkg.append(list(event))
        


# In[15]:


len(bkg[0])


# # VISUALISE DISTRIBUTIONS

# In[16]:


dataset_bkg = np.array(bkg)
dataset_sig = np.array(signal)

TH = ROOT.TH1F("Fat_bkg", "Fat_bkg", 7, 0, 7)
TH.SetLineColor(ROOT.kRed)
TH.SetLineWidth(3)
TH1 = ROOT.TH1F("Fay_fig", "Fat_fig", 7, 0, 7)
TH1.SetLineColor(ROOT.kBlue)
TH1.SetLineWidth(3)
for b,s in zip(dataset_bkg, dataset_sig):
    TH.Fill(b[9])
    TH1.Fill(s[9])
    
c = ROOT.TCanvas("c", "c", 1000,1000, 1000,800)
TH.Draw("hist")
TH1.Draw("hist same")
c.Draw()
    


# In[17]:


TH = ROOT.TH1F("jet_bkg", "jet_bkg", 15, 0, 15)
TH.SetLineColor(ROOT.kRed)
TH.SetLineWidth(3)
TH1 = ROOT.TH1F("jet_fig", "jet_figa", 15, 0, 15)
TH1.SetLineColor(ROOT.kBlue)
TH1.SetLineWidth(3)
for b,s in zip(dataset_bkg, dataset_sig):
    TH.Fill(b[0])
    TH1.Fill(s[0])
    
c = ROOT.TCanvas("c", "c", 1000,1000, 1000,800)
TH.Draw("hist")
TH1.Draw("hist same")
c.Draw()


# # CREATING DATASET  TRAINING

# In[18]:


from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
import matplotlib.pyplot as plt
import math as mt

def plot_confusion_matrix(y_true, y_pred, classes, normalize=False, title=None, cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    #classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        
    fig, ax = plt.subplots(figsize=(5,5))
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax
    


# In[19]:


dataset_signal = np.c_[dataset_sig, (np.zeros(dataset_sig.shape[0]))]
dataset_background = np.c_[dataset_bkg, (np.ones(dataset_bkg.shape[0]))]


# In[20]:


dataset_signal.shape


# In[21]:


final_dataset = np.concatenate((dataset_signal, dataset_background), axis = 0)


# In[22]:


np.random.shuffle(final_dataset)


# In[33]:


import sklearn
from sklearn.model_selection import train_test_split

def target_split(dataset, target):
    y = dataset[:,target]
    X = np.delete(dataset, target, 1)
    
    return X, y


# In[34]:


X, y = target_split(final_dataset, -1)


# In[35]:


from sklearn.preprocessing import StandardScaler, LabelBinarizer
from sklearn.model_selection import StratifiedShuffleSplit
import keras
import tensorflow
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.callbacks import EarlyStopping

scaler = StandardScaler()
X = scaler.fit_transform(X)

lb = LabelBinarizer()
y = lb.fit_transform(y)

x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)


# In[36]:


y_test


# In[37]:


model = Sequential()
model.add(Dense(300, input_dim=x_train.shape[1], activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(units=200, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(units=100, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(units=30, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(units=20, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(units=10, activation='relu'))
model.add(Dense(units=y_train.shape[1], activation='sigmoid'))

early_stop = EarlyStopping(monitor='val_loss', min_delta=1e-5, patience=8, verbose=1, mode='auto', baseline=None)

model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

history = model.fit(x_train, y_train, epochs=1000, validation_data=(x_test, y_test), callbacks= [early_stop])


# In[55]:


fig = plt.figure(figsize=(12,7))
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Loss')
plt.legend(['Train','Test'])
plt.show()


# In[56]:


fig = plt.figure(figsize=(12,7))
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title("Accuracy")
plt.legend(['train', 'test'])
plt.show()


# In[57]:


from sklearn.metrics import accuracy_score

y_pred_proba = model.predict(x_test)


# In[58]:


for i in range(len(y_pred_proba)):
    if y_pred_proba[i] >= 0.5:
        y_pred_proba[i] = 1
    else:
        y_pred_proba[i] = 0


# In[59]:


print("accuracy score: {:.2f}".format(accuracy_score(y_test, y_pred_proba)))


# In[60]:


#We feed into the function the test, the prediction, the list of the class labels for talc variable and the title
plot_confusion_matrix(y_test, y_pred_proba, [0,1] , normalize=True, title = "Confusion matrix S/B = 0,1")


# In[61]:


y_pred_proba = model.predict(x_test)

legend = ROOT.TLegend(0.11, 0.8, 0.4,0.9)
histo_S = ROOT.TH1F("h1","h1", 100, 0, 1.01)
legend.AddEntry(histo_S, "Signal")
histo_B = ROOT.TH1F("h2","h2", 100, 0, 1.01)
legend.AddEntry(histo_B, "Background")
k1 = ROOT.TCanvas("c1","c1",50,50,1000,800)


for i in range(0,len(y_pred_proba)):
    if(y_test[i] == 1):
        histo_B.Fill(y_pred_proba[i])
    else:
        histo_S.Fill(y_pred_proba[i])

histo_S.SetFillColor(ROOT.kRed)
histo_S.SetLineColor(ROOT.kRed)
histo_B.SetLineColor(ROOT.kBlue)
histo_B.SetFillColor(ROOT.kBlue)
histo_S.SetFillStyle(3002)
histo_B.SetFillStyle(3004)
histo_S.Draw("")
histo_B.Draw("same")
k1.SetLogy()
ROOT.gStyle.SetOptStat(0)
legend.Draw()
k1.Draw()


# In[62]:


bins = histo_S.GetNbinsX()
Significance = ROOT.TH1F("sig","Significance",bins,0,1)
c_sig = ROOT.TCanvas("c2","c2",1000,1000,1000,1000)
#Significance.SetStats(0000)
Significance.SetLineWidth(4)
Significance.SetLineColor(ROOT.kRed)
for i in range(0,bins):
    int_S = histo_S.Integral(i,bins)
    int_B = histo_B.Integral(i,bins)
    if(int_B != 0):
        sig = int_S / mt.sqrt(int_B)
        Significance.SetBinContent(i,sig)
#Significance.SetXTitle("NN Response")
#Significance.SetYTitle("Significance")
#c_sig.SetTitle("Significance")
Significance.Draw("histo")
c_sig.Draw("")


# In[63]:


from sklearn.metrics import roc_auc_score, roc_curve, auc

def plot_roc_curve(y_test, pred):
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='b',
    label='0 pred power', alpha=.8)
    fp , tp, th = roc_curve(y_test, pred)
    roc = roc_auc_score(y_test, pred)
    plt.plot(fp, tp, 'r', label='ROC binary categorizzation (AUC = %0.3f)' %(roc))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic curve')
    plt.legend(loc="lower right")
    

plt.figure(figsize=(13,8))
plot_roc_curve(y_test, y_pred_proba)


# In[47]:


h_sig = ROOT.TH1F("t", "t", 100,0,1400)
h_sig.SetLineColor(ROOT.kRed)
h_sig.SetLineWidth(2)
h_sig.SetFillStyle(3003)
h_bkg = ROOT.TH1F("t1", "t1", 100,0,1400)
h_bkg.SetLineColor(ROOT.kBlue)
h_bkg.SetLineWidth(2)
h_bkg.SetFillStyle(3003)


for i, j in zip(dataset_signal, dataset_bkg):
    h_sig.Fill(i[1])
    h_bkg.Fill(j[1])
    
    
c = ROOT.TCanvas("c", "c", 1000,1000,1000,800)
h_sig.Draw("hist")
h_bkg.Draw("hist same")
c.Draw()

l = []
for bins in range(1400):
    sig_bin = h_sig.GetBinContent(bins)
    bkg_bin = h_bkg.GetBinContent(bins)
    
    if (sig_bin and bkg_bin) != 0:
        l.append((sig_bin-bkg_bin)**2/(sig_bin+bkg_bin))
    else:
        l.append(0)
sig = sum(l)
    


# In[48]:


sig/2


# In[32]:


import pandas as pd


# In[33]:


columns = ['number_of_jets', 'jet_max_pt', 'jets_max_mass', 'jets_max_eta', 'jets_max_phi', 'electron_Pt', 'muon_Pt', 'muon_eta', 'muon_phi', 'muon_T', 'met_pt', 'met_m', 'NCh', 'NNeu', 'Flava','number_of_fat_jets', 'fat_jet_max_pt', 'fat_jets_max_mass', 'number_of_photon', 'photon_max_pt', 'photon_max_energy', 'photon_max_eta', 'photon_max_phi', 'label']


# In[34]:


len(columns)


# In[35]:


dataset_signal.shape


# In[36]:


data_sig = pd.DataFrame(dataset_signal, columns = columns)


# In[37]:


data_sig.to_csv("./data_signal.csv", columns = columns, index=False)


# In[38]:


data_bkg = pd.DataFrame(dataset_background, columns = columns)
data_bkg.to_csv("./data_background.csv", columns = columns, index=False)


# In[39]:


f = pd.read_csv("./data_background.csv")


# In[ ]:





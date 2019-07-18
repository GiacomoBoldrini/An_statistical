#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random

#keras imports
from keras.models import Model
from keras.layers import Dense, Input, Conv2D, Dropout, Flatten
from keras.layers import MaxPooling2D, BatchNormalization, Activation
from keras.utils import plot_model
from keras import backend as K
from keras import metrics
from keras.callbacks import EarlyStopping, ReduceLROnPlateau, TerminateOnNaN
import tensorflow as tf 


# In[2]:


gluon_im = np.load("/Users/bcoder/Bonesini_qq_gg/image_creators/Gluons_images_NEW.npy")


# In[3]:


plt.figure(figsize=(10,10))
SUM_Image = np.sum(gluon_im, axis = 0)
plt.imshow(SUM_Image/float(gluon_im.shape[0]), origin='lower', norm=LogNorm(vmin=0.01))
plt.colorbar()
plt.xlabel("Δη", fontsize=20)
plt.ylabel("ΔΦ", fontsize=20)
plt.title("Gluons", fontsize=20)
#plt.savefig('../graphs/Gluons_old_jets.pdf')


# In[4]:


quark_im = np.load("/Users/bcoder/Bonesini_qq_gg/image_creators/Quark_images_NEW.npy")


# In[5]:


plt.figure(figsize=(10,10))
SUM_Image = np.sum(quark_im, axis = 0)
plt.imshow(SUM_Image/float(quark_im.shape[0]), origin='lower', norm=LogNorm(vmin=0.01))
plt.colorbar()
plt.xlabel("Δη", fontsize=20)
plt.ylabel("ΔΦ", fontsize=20)
plt.title("Quark", fontsize=20)
#plt.savefig('../graphs/Quarks_old_jets.pdf')


# # CNN TRAINING

# In[6]:


def Creating_train_test(X, y, size = 0.2):
    
    X = list(X)
    y = list(y)
    X_test = []
    y_test = []
    
    test_size = int(size*len(X))
    counting = 0
    while(counting != test_size):
        try:
            r = random.randint(0, len(X))
            x_i = X.pop(r)
            y_i = y.pop(r)
            X_test.append(x_i)
            y_test.append(y_i)
            
            counting += 1
        except:
            continue
    
    X_train = np.array(X)
    y_train = np.array(y)
    X_test = np.array(X_test)
    y_test = np.array(y_test)
    
    return X_train, X_test, y_train, y_test


# In[7]:


shape = min(quark_im.shape[0], gluon_im.shape[0])
quark_im = quark_im[:shape]
gluon_im = gluon_im[:shape]


# In[8]:


y_gluons = [1]*gluon_im.shape[0]
y_quarks = [0]*quark_im.shape[0]

y = y_gluons + y_quarks
y = np.array(y)

X = np.concatenate((gluon_im, quark_im), axis=0)

print(X.shape, y.shape)


# In[9]:


x_train, x_test, y_train, y_test = Creating_train_test(X, y, size = 0.1)
print(x_train.shape, x_test.shape, y_train.shape, y_test.shape)


# In[10]:


x_train = x_train.reshape((x_train.shape[0], x_train.shape[1], x_train.shape[2], 1))
x_test = x_test.reshape((x_test.shape[0], x_test.shape[1], x_test.shape[2], 1))


# In[11]:


img_rows = x_train.shape[1]
img_cols = x_train.shape[2]
dropoutRate = 0.1


# In[12]:


image_shape = (img_rows, img_cols, 1)
####
inputImage = Input(shape=(image_shape))
x = Conv2D(5, kernel_size=(5,5), data_format="channels_last", strides=(1, 1), padding="same")(inputImage)
#x = Conv2D(10, kernel_size=(5,5), data_format="channels_last", strides=(1, 1), padding="same")(inputImage)
x = BatchNormalization()(x)
x = Activation('relu')(x)
x = MaxPooling2D( pool_size = (5,5))(x)
x = Dropout(dropoutRate)(x)
#
x = Conv2D(3, kernel_size=(3,3), data_format="channels_last", strides=(1, 1), padding="same")(x)
x = BatchNormalization()(x)
x = Activation('relu')(x)
x = MaxPooling2D( pool_size = (3,3))(x)
x = Dropout(dropoutRate)(x)
#Here the output as a 1xN vector of convolution
x = Flatten()(x)
#Here enters DNN
x = Dense(50, activation='relu')(x)
x = Dropout(dropoutRate)(x)
x = Dense(20, activation='relu')(x)
x = Dropout(dropoutRate)(x)
x = Dense(5, activation='relu')(x)
x = Dropout(dropoutRate)(x)
#
output = Dense(1, activation='sigmoid')(x)
####
model = Model(inputs=inputImage, outputs=output)


# In[13]:


model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.summary()

tf.keras.utils.plot_model(model, to_file='../NN_images/model_CNN.png')


# In[14]:


from keras.utils import plot_model
plot_model(model, to_file='model_CNN.png')


# In[15]:


batch_size = 128
n_epochs = 20


# In[16]:


history = model.fit(x_train, y_train, epochs=n_epochs, batch_size=batch_size, verbose = 1,
                validation_data=(x_test, y_test),
                callbacks = [
                EarlyStopping(monitor='val_loss', patience=10, verbose=1),
                ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=2, verbose=1),
                TerminateOnNaN()])


# In[17]:


model_json = model.to_json()
with open("CNN.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model.save_weights("CNN.h5")
print("Saved model to disk")


# # Plotting metrics and results

# In[18]:


from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels


# In[27]:


plt.figure(figsize=(15,10))
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title('Model accuracy', size=18)
plt.ylabel('accuracy', size=18)
plt.xlabel('epoch', size=18)
plt.legend(['train', 'test'], loc='lower right', prop={'size': 16})
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNN_acc.pdf")


# In[28]:


plt.figure(figsize=(15,10))
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model Loss', size=18)
plt.ylabel('Loss value', size=18)
plt.xlabel('epoch', size=18)
plt.legend(['train', 'test'], loc='lower right', prop={'size': 16})
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNN_loss.pdf")


# In[21]:


y_pred = model.predict(x_test)
y_pred = [round(i[0]) for i in y_pred]


# In[22]:


from sklearn.metrics import roc_auc_score, roc_curve, auc

y_pred = model.predict(x_test)
def plot_roc_curve(y_test, pred):
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='b',
    label='0 pred power', alpha=.8)
    fp , tp, th = roc_curve(y_test, pred)
    roc = roc_auc_score(y_test, pred)
    plt.plot(fp, tp, 'r', label='ROC binary categorizzation (AUC = %0.3f)' %(roc))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic curve')
    plt.legend(loc="lower right", prop={'size': 16})
    

plt.figure(figsize=(13,8))
plot_roc_curve(y_test, y_pred)
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNN_ROC.pdf")


# In[26]:


y_pred_proba = model.predict(x_test)
y_pred_sig = []
y_pred_bkg = []
for i in range(0,len(y_pred_proba)):
    if(y_test[i] == 1):
        y_pred_bkg.append(y_pred_proba[i][0])
    else:
         y_pred_sig.append(y_pred_proba[i][0])

plt.figure(figsize=(20,10))
n, bins, _ = plt.hist(y_pred_sig, 100, histtype='step', fill=False, linewidth=2, label = "Quarks")
n1, bins1, _ = plt.hist(y_pred_bkg, 100, histtype='step', fill=False, linewidth=2, label = "Gluons")
#plt.yscale('log', nonposy='clip')
plt.legend(loc = "upper center",  borderpad=1, fontsize=25)
plt.xlabel("Prediction Value", size=15)
plt.ylabel("Counts", size=15)
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNN_prob_dist.pdf")

np.save("./CNN_predictions.npy", y_pred_proba)


# In[ ]:





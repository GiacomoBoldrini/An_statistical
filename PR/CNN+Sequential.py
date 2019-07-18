#!/usr/bin/env python
# coding: utf-8

from keras.layers import Conv2D, MaxPooling2D, Flatten
from keras.layers import Input, LSTM, Embedding, Dense
from keras.models import Model, Sequential
from keras.layers import Dense, Dropout
from keras.layers import BatchNormalization, Activation
from keras.utils import plot_model
from keras import backend as K
from keras import metrics
from keras.callbacks import EarlyStopping, ReduceLROnPlateau, TerminateOnNaN
import keras
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels


# # Loading both images and variables for each jet 



def makePlot(feature_name, signal_data, background_data, Save=False):
    #plt.subplots()
    # notice the use of numpy masking to select specific classes of jets
    sig = signal_data[feature_name]
    bkg = background_data[feature_name]
    plt.figure(figsize=(13,8))
    # then plot the right quantity for the reduced array
    plt.hist(sig, 50, density=True, histtype='step', fill=False, linewidth=1.5, label = "sig")
    plt.hist(bkg, 50, density=True, histtype='step', fill=False, linewidth=1.5, label = "bkg")
    plt.yscale('log', nonposy='clip')    
    plt.legend(fontsize=12, frameon=False)  
    plt.xlabel(feature_name, fontsize=15)
    plt.ylabel('Prob. Density (a.u.)', fontsize=15)
    #plt.show()
    if Save:
        plt.savefig("../graphs/{}.pdf".format(feature_name))
    #del fig, ax
    #return fig, ax


#loading attributes
quark_var = np.load("../FNN_creator/q_var_arr.npy")
gluon_var = np.load("../FNN_creator/g_var_arr.npy")


columns = ['jet_PT', 'jet_Eta', 'jet_Phi', 'jets_number', 'jet_Mass', 'jet_Mass_PT_ratio', 'jet_Charged_Multiplicity', 'jet_Neutral_Multiplicity', 'jet_ECAL_HCAL_ratio', 'MeanSqDeltaR', 'PTD', 'tau1', 'frac_tau2_tau1', 'tau2', 'frac_tau3_tau2','tau3', 'frac_tau4_tau3', 'tau4', 'frac_tau5_tau4']


quark_pd_var = pd.DataFrame(quark_var, columns=columns)
gluon_pd_var = pd.DataFrame(gluon_var, columns=columns)


#for i in range(len(columns)):
#    makePlot(columns[i], quark_pd_var, gluon_pd_var, Save=True)

#adding label
shape = min(quark_var.shape[0], gluon_var.shape[0])

quark_var = quark_var[:shape]
gluon_var = gluon_var[:shape]


final_var_data = np.concatenate((gluon_var, quark_var), axis = 0)
final_y_var = np.concatenate((np.ones(gluon_var.shape[0]), np.zeros(quark_var.shape[0])), axis = 0)


var_shape = (final_var_data.shape[0], )


print("Final var shape:{},{}".format(final_var_data.shape, final_y_var.shape))


#loading images

gluons_images = np.load("/Users/bcoder/Bonesini_qq_gg/image_creators/Gluons_images_NEW.npy")
quarks_images = np.load("/Users/bcoder/Bonesini_qq_gg/image_creators/Quark_images_NEW.npy")

#inferior count class
shape = min(gluons_images.shape[0], quarks_images.shape[0])

#equal size class

gluons_images = gluons_images[:shape]
quarks_images = quarks_images[:shape]


print(gluons_images.shape)
print(quarks_images.shape)

y_gluons = [1]*gluons_images.shape[0]
y_quarks = [0]*quarks_images.shape[0]

y_im = y_gluons + y_quarks
y_im = np.array(y_im)

X_im = np.concatenate((gluons_images, quarks_images), axis=0)

print(X_im.shape, y_im.shape)


def Creating_train_test(X_var, X_im,  y_var, y_im, size = 0.2):
    
    X_var = list(X_var)
    y_var = list(y_var)
    X_im = list(X_im)
    y_im = list(y_im)
    X_var_test = []
    y_var_test = []
    X_im_test = []
    y_im_test = []
    
    test_size = int(size*len(X_var))
    counting = 0
    while(counting != test_size):
        try:
            #index
            r = random.randint(0, len(X_var)-counting)
            x_var_i = X_var.pop(r)
            x_im_i = X_im.pop(r)
            y_var_i = y_var.pop(r)
            y_im_i = y_im.pop(r)
            X_var_test.append(x_var_i)
            y_var_test.append(y_var_i)
            X_im_test.append(x_im_i)
            y_im_test.append(y_im_i)
            
            counting += 1
        except:
            continue
    
    X_var_train = np.array(X_var)
    y_var_train = np.array(y_var)
    X_var_test = np.array(X_var_test)
    y_var_test = np.array(y_var_test)
    
    X_im_train = np.array(X_im)
    y_im_train = np.array(y_im)
    X_im_test = np.array(X_im_test)
    y_im_test = np.array(y_im_test)
    
    return X_var_train, X_var_test, y_var_train, y_var_test, X_im_train, X_im_test, y_im_train, y_im_test


X_var_train, X_var_test, y_var_train, y_var_test, X_im_train, X_im_test, y_im_train, y_im_test = Creating_train_test(final_var_data, X_im, final_y_var, y_im)


print("Var dataset train test: {}, {}, {}, {}".format(X_var_train.shape, X_var_test.shape, y_var_train.shape, y_var_test.shape))

print("Var dataset train test: {}, {}, {}, {}".format(X_im_train.shape, X_im_test.shape, y_im_train.shape, y_im_test.shape))

#reshaping images for the CNN
X_im_train = X_im_train.reshape((X_im_train.shape[0], X_im_train.shape[1], X_im_train.shape[2], 1))
X_im_test = X_im_test.reshape((X_im_test.shape[0], X_im_test.shape[1], X_im_test.shape[2], 1))

img_rows = X_im_train.shape[1]
img_cols = X_im_train.shape[2]
dropoutRate = 0.1


image_shape = (img_rows, img_cols, 1)
###


image_input = Input(shape=(image_shape), name = 'image_input')

x = Conv2D(10, kernel_size=(5,5), data_format="channels_last", strides=(1, 1), padding="same")(image_input)

x = BatchNormalization()(x)
x = Activation('relu')(x)
x = MaxPooling2D( pool_size = (5,5))(x)
x = Dropout(dropoutRate)(x)
#
x = Conv2D(5, kernel_size=(5,5), data_format="channels_last", strides=(1, 1), padding="same")(x)
x = BatchNormalization()(x)
x = Activation('relu')(x)
x = MaxPooling2D( pool_size = (5,5))(x)
x = Dropout(dropoutRate)(x)
#Here the output as a 1xN vector of convolution
x = Flatten()(x)

#concatenating kinematical variables
kin_var_input = Input(shape=(X_var_train.shape[1],), name = 'kin_var_input')
merged = keras.layers.concatenate([x, kin_var_input])

#Adding DNN
x = Dense(200, activation='relu')(merged)
x = Dropout(dropoutRate)(x)
x = Dense(100, activation='relu')(x)
x = Dropout(dropoutRate)(x)
x = Dense(50, activation='relu')(x)
x = Dropout(dropoutRate)(x)
#
output = Dense(1, activation='sigmoid', name='output')(x)

# This is our final model:
vqa_model = Model(inputs=[image_input, kin_var_input], outputs=output)

from keras.utils import plot_model
plot_model(vqa_model, to_file='model_FFNN+CNN.png')



vqa_model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])
vqa_model.summary()

print()
history = vqa_model.fit({'image_input': X_im_train, 'kin_var_input': X_var_train},
          {'output': y_im_train},
          epochs=20, batch_size=128, shuffle = True, validation_data=([X_im_test, X_var_test], y_im_test))

model_json = vqa_model.to_json()
with open("CNN+FFNN.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
vqa_model.save_weights("CNN+FFNN.h5")
print("Saved model to disk")



#saving accuracy and loss
plt.figure(figsize=(15,10))
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title('Model accuracy', size=18)
plt.ylabel('accuracy', size=18)
plt.xlabel('epoch', size=18)
plt.legend(['train', 'test'], loc='lower right', prop={'size': 16})
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNNFFNN_acc.pdf")

plt.figure(figsize=(15,10))
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model Loss', size=18)
plt.ylabel('Loss value', size=18)
plt.xlabel('epoch', size=18)
plt.legend(['train', 'test'], loc='lower right', prop={'size': 16})
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNNFFNN_loss.pdf")

from sklearn.metrics import roc_auc_score, roc_curve, auc

y_pred = vqa_model.predict([X_im_test, X_var_test])
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
plot_roc_curve(y_im_test, y_pred)
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNNFFNN_ROC.pdf")

y_pred = vqa_model.predict([X_im_test, X_var_test])
y_pred_sig = []
y_pred_bkg = []
for i in range(0,len(y_pred)):
    if(y_im_test[i] == 1):
        y_pred_bkg.append(y_pred[i][0])
    else:
        y_pred_sig.append(y_pred[i][0])

plt.figure(figsize=(20,10))
n, bins, _ = plt.hist(y_pred_sig, 100, histtype='step', fill=False, linewidth=2, label = "Quarks")
n1, bins1, _ = plt.hist(y_pred_bkg, 100, histtype='step', fill=False, linewidth=2, label = "Gluons")
#plt.yscale('log', nonposy='clip')
plt.legend(loc = "upper center",  borderpad=1, fontsize=25)
plt.xlabel("Prediction Value", size=15)
plt.ylabel("Counts", size=15)
plt.savefig("/Users/bcoder/Bonesini_qq_gg/NN_images/CNNFNN_prob_dist.pdf")

np.save("./CNNFFNN_predictions.npy", y_pred)


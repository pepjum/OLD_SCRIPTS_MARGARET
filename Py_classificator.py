### clasificador en python



f = open("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_EXPERIMENT_WITHOUT_AMBIGUOUS_TyD_PREPARED.txt", "r")

data = f.readlines()

file_sep = map(lambda x: x.split(), data)

#map(func, iteration)
#a = filter(lambda Y: Y[6] == 'T', file_sep)
labels = [item[6] for item in file_sep]
#
labels_rep = [w.replace('T', '1') for w in labels]
labels_rep = [w.replace('D', '0') for w in labels_rep]

class_names = [item[6] for item in file_sep]

for line in file_sep:
    if line[6] == 'T':
        line[6] = '1'
    else :
        line[6] = '0'

##### separar entre train y test y tener los datos preparados
from random import sample

import tensorflow as tf
from tensorflow import keras

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, Dropout, MaxPooling2D
from tensorflow.keras.preprocessing.image import ImageDataGenerator

import os
import numpy as np
import matplotlib.pyplot as plt


train_idx = sample(range(len(file_sep)), 756100)
test_idx   = list(set(range(len(file_sep))) - set(train_idx))
train_data = [file_sep[i] for i in train_idx]
test_data  = [file_sep[i] for i in test_idx]

train_label = [file_sep[i][6] for i in train_idx]
test_label  = [file_sep[i][6] for i in test_idx]

train_data_s = [item[1:6] for item in train_data]
test_data_s = item[1:6] for item in test_data]

####### configurar las capas

model = keras.Sequential([
    keras.layers.Flatten(input_shape=(756100, 5)),
    keras.layers.Dense(128, activation='relu'),
    keras.layers.Dense(10, activation='sigmoid')
])

# compilar modelo
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])

### entrenar modelo

model.fit(train_data_s, train_label, epochs=10)

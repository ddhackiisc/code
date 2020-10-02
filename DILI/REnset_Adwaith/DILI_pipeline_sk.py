# import necesary modules
import sys
import gc
from keras.layers.normalization import BatchNormalization
from keras.layers.convolutional import Conv1D
from keras.layers.convolutional import AveragePooling1D
from keras.layers.convolutional import MaxPooling1D
from keras.layers.convolutional import ZeroPadding1D
from keras.layers.core import Activation
from keras.layers.core import Dense
from keras.layers import Flatten
from keras.layers import Input, Reshape
from keras.models import Model
from keras.layers import add
from keras.regularizers import l2
from keras.losses import BinaryCrossentropy, CategoricalCrossentropy, MeanAbsoluteError, MeanSquaredError
from keras import backend as K
import tensorflow as tf
import pandas as pd
import numpy as np
from keras.optimizers import Adam, Nadam, RMSprop, SGD
from keras.metrics import AUC, Accuracy
import sklearn.model_selection as sk
from keras import backend as K
import csv
from multiprocessing import Process, Array

# Define the model inside the class Resnet


class Resnet:
    @staticmethod
    def residual_module(data, K, stride, activation, kernel_size, chanDim=-1, red=False, reg=0.0001, bnEps=2e-5, bnMom=0.9, i=-1, j=-1):
        # the shortcut branch of a residual module is initialized as the input
        shortcut = data

        # First block of resudual module
        bn1 = BatchNormalization(
            epsilon=bnEps, momentum=bnMom, name='bn1_'+str(i)+'_'+str(j))(data)
        act1 = Activation(
            activation, name='act1_'+str(i)+'_'+str(j))(bn1)
        conv1 = Conv1D(K, kernel_size, use_bias=False, kernel_regularizer=l2(
            reg), padding='same', name='conv1_'+str(i)+'_'+str(j))(act1)

        # 2nd block of residual module
        bn2 = BatchNormalization(
            epsilon=bnEps, momentum=bnMom, name='bn2_'+str(i)+'_'+str(j))(conv1)
        act2 = Activation(
            activation, name='act2_'+str(i)+'_'+str(j))(bn2)
        conv2 = Conv1D(K, kernel_size, use_bias=False, kernel_regularizer=l2(
            reg), padding='same', name='conv2_'+str(i)+'_'+str(j))(act2)

        # In the beginning of each stage of the resnet you will have adjust the shortcut dimensions
        if red:
            shortcut = Conv1D(K, kernel_size, use_bias=False, kernel_regularizer=l2(
                reg), padding='same', name='conv_srt_'+str(i)+'_'+str(j))(act1)

        # add the second residual moduleand the shortcut
        x = add([conv2, shortcut],
                name='add_'+str(i)+'_'+str(j))

        return x

    @staticmethod
    def build(input_shape, classes, stages, filters, activation, fc_layers, kernel_size=9, stride=1, reg=0.0001, bnEps=2e-5, bnMom=0.9):
        # create input and apply BN
        inputs = Input(shape=input_shape)
        x = BatchNormalization(
            epsilon=bnEps, momentum=bnMom, name='bn_input')(inputs)
        x = Reshape((881, 1), name='Reshape')(x)

        # Create the starting conv layer
        x = Conv1D(filters[0], kernel_size, use_bias=False, padding='same',
                   kernel_regularizer=l2(reg), name='conv_first_layer')(x)
        x = BatchNormalization(
            epsilon=bnEps, momentum=bnMom, name='bn_first_layer')(x)
        x = Activation(
            activation, name='act_first_layer')(x)

        # loop over the stages
        for i in range(0, len(stages)):
            # first layer in each stage should reduce the dimensions of input
            if i == 0 and filters[0] == filters[1]:
                x = Resnet.residual_module(
                    x, filters[i+1], stride, activation, kernel_size, bnEps=bnEps, bnMom=bnMom, i=i+1, j=1)
            else:
                x = Resnet.residual_module(
                    x, filters[i+1], stride, activation, kernel_size, red=True, bnEps=bnEps, bnMom=bnMom, i=i+1, j=1)

            for j in range(0, stages[i] - 1):
                # apply a renset module
                x = Resnet.residual_module(
                    x, filters[i+1], stride, activation, kernel_size, bnEps=bnEps, bnMom=bnMom, i=i+1, j=j+2)

        x = BatchNormalization(
            epsilon=bnEps, momentum=bnMom, name='bn_last_layer')(x)
        x = Activation(activation)(x)

        # Final dense layer applied with softmax activation
        x = Flatten(name='flatten')(x)
        x = Dense(fc_layers, kernel_regularizer=l2(
            reg), name='out_dense_layer')(x)
        x = Dense(classes, activation='sigmoid',
                  name='last_layer', kernel_regularizer=l2(reg))(x)
        # x = Activation('softmax')(x)
        # x = Activation('sigmoid', name='out_act')(x)

        # create the model
        model = Model(inputs, x, name='resnet')

        # return the model
        return model

# Define functions for the metrics


def sensitivity(y_true, y_pred):
    true_positives = K.sum(
        K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(
        K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())


def specificity(y_true, y_pred):
    true_negatives = K.sum(
        K.round(K.clip((1-y_true) * (1-y_pred), 0, 1)))
    possible_negatives = K.sum(
        K.round(K.clip(1-y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())


def f1(y_true, y_pred):  # taken from old keras source code
    true_positives = K.sum(
        K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(
        K.round(K.clip(y_true, 0, 1)))
    predicted_positives = K.sum(
        K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / \
        (predicted_positives + K.epsilon())
    recall = true_positives / \
        (possible_positives + K.epsilon())
    f1_val = 2*(precision*recall) / \
        (precision+recall+K.epsilon())
    return f1_val


def matthews_correlation(y_true, y_pred):
    y_pred_pos = K.round(K.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos
    y_pos = K.round(K.clip(y_true, 0, 1))
    y_neg = 1 - y_pos
    tp = K.sum(y_pos * y_pred_pos)
    tn = K.sum(y_neg * y_pred_neg)
    fp = K.sum(y_neg * y_pred_pos)
    fn = K.sum(y_pos * y_pred_neg)
    numerator = (tp * tn - fp * fn)
    denominator = K.sqrt(
        (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return numerator / (denominator + K.epsilon())


# define a function to train, cross-validate and predict on test a model with given hyperparameters
def train_predict(X_train, X_test, y_train, y_test, model_no, lossf, kernel_size, stages, filters, stride, fc_layers, activation, optimizer, reg, batch_size, epochs):
    # Initialize arrays for metrics
    loss = []
    acc = []
    auc = []
    f1s = []
    mcc = []
    sens = []
    spec = []

    # # Define K-fold cross validator
    # kfold = sk.KFold(n_splits=10, shuffle=True)
    #
    # fold_no = 1
    #
    # for train, val in kfold.split(X_train, y_train):
    #     # with mirrored_strategy.scope():
    #     ret_obj = Array('d',  7, lock=False)
    #
    #     def run(ret_obj):
    #         # gpu_devices = tf.config.experimental.list_physical_devices(
    #         #     'GPU')
    #         # tf.config.experimental.set_memory_growth(
    #         #     gpu_devices[0], True)
    #         # tf.config.experimental.set_memory_growth(
    #         #     gpu_devices[0], True)
    #         my_model = Resnet.build((881,), 2, stages, filters, activation,
    #                                 fc_layers, kernel_size=kernel_size, stride=stride, reg=reg)
    #         my_model.compile(loss=lossf, optimizer=optimizer, metrics=[
    #                          'accuracy', AUC(), f1, matthews_correlation, sensitivity, specificity])
    #         print(
    #             "Training for fold {}....".format(fold_no))
    #
    #         my_model.fit(X_train[train], y_train[train],
    #                      batch_size=batch_size, epochs=epochs, verbose=0)
    #
    #         result = my_model.evaluate(
    #             X_train[val], y_train[val], verbose=0)
    #
    #         for i in range(7):
    #             ret_obj[i] = result[i]
    #
    #     pro = Process(target=run, args=[ret_obj])
    #     pro.start()
    #     pro.join()
    #
    #     scores = ret_obj
    #     loss.append(scores[0])
    #     acc.append(scores[1]*100)
    #     auc.append(scores[2])
    #     f1s.append(scores[3])
    #     mcc.append(scores[4])
    #     sens.append(scores[5])
    #     spec.append(scores[6])
    #
    #     fold_no = fold_no + 1

    print("Validating on independent test data set")
    # with mirrored_strategy.scope():

    ret_obj = Array('d',  7, lock=False)

    def run2(ret_obj):
        # gpu_devices = tf.config.experimental.list_physical_devices(
        #     'GPU')
        # tf.config.experimental.set_memory_growth(
        #     gpu_devices[0], True)
        # tf.config.experimental.set_memory_growth(
        #     gpu_devices[0], True)

        my_model = Resnet.build((881,), 2, stages, filters, activation,
                                fc_layers, kernel_size=kernel_size, stride=stride, reg=reg)
        my_model.compile(loss=lossf, optimizer=optimizer, metrics=[
                         'accuracy', AUC(), f1, matthews_correlation, sensitivity, specificity])
        my_model.fit(X_train, y_train,
                     batch_size=batch_size, epochs=epochs, verbose=0)
        result = my_model.evaluate(
            X_test, y_test, verbose=0)

        for i in range(7):
            ret_obj[i] = result[i]

    pro = Process(target=run2, args=[ret_obj])
    pro.start()
    pro.join()

    scores = ret_obj

    return [model_no, lossf, kernel_size, stages, filters, stride, fc_layers, activation, optimizer, reg, batch_size, epochs, np.mean(loss), np.mean(acc), np.mean(auc), np.mean(f1s), np.mean(mcc), np.mean(sens), np.mean(spec), scores[0], scores[1]*100, scores[2], scores[3], scores[4], scores[5], scores[6]]


# Load the data and split it into trianing and independent test set
data = pd.read_csv('Merged_input_with_labels.csv',
                   header=None).to_numpy(dtype=np.int8)
input = data[:, :881]
label = data[:, 881:883]
X_train, X_test, y_train, y_test = sk.train_test_split(
    input, label, test_size=0.1, random_state=42)


lossf = [MeanSquaredError()]
kernel_size = [9]
stages_filters = [[[3], [9, 9]], [[3], [9, 18]], [[3], [18, 18]], [[3], [18, 36]], [[3, 3], [9, 9, 18]], [[3, 3], [9, 18, 36]], [[3, 3], [18, 18, 36]], [[3, 3], [18, 36, 72]], [[3, 3, 3], [9, 9, 18, 36]], [[3, 3, 3], [
    9, 18, 36, 72]], [[3, 3, 3], [18, 18, 36, 72]], [[3, 3, 3], [18, 36, 72, 144]], [[3, 3, 3, 3], [9, 9, 18, 36, 72]], [[3, 3, 3, 3], [9, 18, 36, 72, 144]], [[3, 3, 3, 3], [18, 18, 36, 72, 144]], [[3, 3, 3, 3], [18, 36, 72, 144, 288]]]
stride = [1]
fc_layers = [970]
activation = ['relu', 'sigmoid']
reg = [0.0001]
optimizer = [Adam(learning_rate=0.00001)]
batch_size = [47]
epochs = [50, 1500]

Results = []
model_no = 1
print('\n---------------------------------------------------------------------\nTotal Models to evaluate: {}\n---------------------------------------------------------------------\n'.format(
    len(lossf)*len(kernel_size)*len(stages_filters)*len(stride)*len(fc_layers)*len(activation)*len(reg)*len(optimizer)*len(batch_size)*len(epochs)))
for l in lossf:
    for ks in kernel_size:
        for s, f in stages_filters:
            for stri in stride:
                for fc in fc_layers:
                    for ac in activation:
                        for r in reg:
                            for o in optimizer:
                                for b in batch_size:
                                    for e in epochs:
                                        print('\n-----------\n Model {} \n-----------\n'.format(
                                            model_no))
                                        print('Parameters : \n          Loss fn: {}\n          Kernel size: {}\n          Stages: {}\n          Filters: {}\n          Strides: {}\n          Full_connected layers: {}\n          Activation function: {}\n          Reg: {}\n          Optimizer: {}\n          Batch_size: {}\n          Epochs: {}\n'.format(
                                            l, ks, s, f, stri, fc, ac, r, o, b, e))
                                        score = train_predict(
                                            X_train, X_test, y_train, y_test, model_no, l, ks, s, f, stri, fc, ac, o, r, b, e)
                                        print('\nScores : \n       Val loss: {}\n       Val Accuracy(%): {}\n       Val AUC: {}\n       Val F1: {}\n       Val MCC: {}\n       Val Sensitivity: {}\n       Val Specifity: {}\n       Test loss: {}\n       Test Accuracy(%): {}\n       Test AUC: {}\n       Test F1: {}\n       Test MCC: {}\n       Test Sensitivity: {}\n       Test Specifity: {}\n       '.format(
                                            score[12], score[13], score[14], score[15], score[16], score[17], score[18], score[19], score[20], score[21], score[22], score[23], score[24], score[25]))
                                        model_no = model_no + 1
                                        Results.append(
                                            score)

with open('Results.csv'.format(n), 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Model no.', 'lossf', 'kernel_size', 'stages', 'filters', 'stride', 'fc_layers', 'activation', 'optimizer', 'reg', 'batch_size', 'epochs', 'val_loss',
                        'val_acc(%)', 'val_auc', 'val_f1s', 'val_mcc', 'val_sens', 'val_spec', 'test_loss', 'test_acc(%)', 'test_auc', 'test_f1s', 'test_mcc', 'test_sens', 'test_spec'])
    csvwriter.writerows(Results)

print("Completed")

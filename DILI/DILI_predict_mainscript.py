import tensorflow as tf
import pandas as pd
import numpy as np
import csv

import pickle

#   Hello, add the necessary packages



def process_smiles_input(file_location):

    #input code to read, process file into fingerprints
    #this will probably be a text file with smiles strings
    #we can ask them to keep it as some specific format
    #maybe \n seperated of \t or ',' separated text file
    #choose which one is most convenient

    #returning the processed numpy array

    return fingerprint_input



def load_model(file_location):

    file = open(file_location, 'rb')

    model = pickle.load(file)

    file.close()

    return model




model_file_location = input("Please enter the pickle file name:- ")

#CHANGE THIS TEXT TO YOUR PREFERRED FORMAT
smiles_file_location = input("Please enter the smiles intput file name. Ensure that the input is a txt file with tab seperated canonical SMILES:-")



model = load_model(model_file_location)
training_data = process_smiles_input(smiles_file_location)


######### code to predit and give results. 0 for non tox, 1 for tox######################






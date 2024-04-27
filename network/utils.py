import json
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

def load_data_from_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract features and label
    features = []
    label = None
    for key, value in data.items():
        if key == 'label':
            label = value
        else:
            features.append(value)
    
    # Convert to numpy array
    features = np.array(features)
    label = np.array(label)
    
    return features, label

def read_csv_to_np_arrays(file_path):
    df = pd.read_csv(file_path)
    df_features_10 = df.iloc[:, :10]
    df_features_3 = df.iloc[:, 10:]
    features_10 = df_features_10.to_numpy()
    features_3 = df_features_3.to_numpy()
    return features_10, features_3


def k_fold_cross_validation(k, data, labels):
    kf = KFold(n_splits=k)
    fold_data = []

    for train_index, test_index in kf.split(data):
        X_train, X_test = data[train_index], data[test_index]
        y_train, y_test = labels[train_index], labels[test_index]
        fold_data.append((X_train, y_train, X_test, y_test))

    return fold_data
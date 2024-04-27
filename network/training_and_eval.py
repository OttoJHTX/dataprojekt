from utils import load_data_from_json, read_csv_to_np_arrays, k_fold_cross_validation
from network_funcs import define_model, train_model, save_model
from tensorflow.keras.models import load_model
import os
import numpy as np


file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(file_dir)


##################################################################################################################################
############################################# LOAD PRE-EXISTING MODEL OR DEFINE NEW ##############################################
##################################################################################################################################

model_name = "testModel"
load_existing_model = True
batch_size = 32
epochs = 10
final_test_N = 10

model = load_model(f'models/{model_name}.h5') if load_existing_model else None
print("Loaded model") if model else None

model = define_model(learning_rate = 0.001, num_layers=3, num_nodes=[10, 20, 20], activation='relu', output_nodes=3) if not model else model

##################################################################################################################################
############################################# LOADING OF DATA & BATCHING (K-FOLD CV) #############################################
##################################################################################################################################

data_points, labels = read_csv_to_np_arrays("data/data.csv")

random_indices = np.random.choice(len(data_points), final_test_N, replace=False)
final_test_data = data_points[random_indices]
final_test_labels = labels[random_indices]

k = int(np.ceil(len(data_points) / batch_size))
fold_data = k_fold_cross_validation(k, data_points, labels)

##################################################################################################################################
############################################## K-FOLD CROSS VALIDATION TRAINING ##################################################
##################################################################################################################################


for fold_idx, (train_data, train_labels, test_data, test_labels) in enumerate(fold_data):
    print(f"Training fold {fold_idx + 1}/{k}...")
    
    model = train_model(model, train_data, train_labels, epochs=epochs, batch_size=batch_size)

    test_loss, test_accuracy = model.evaluate(test_data, test_labels)
    print(f"Fold {fold_idx + 1}/{k}: Test loss: {test_loss}, Test accuracy: {test_accuracy}")

##################################################################################################################################
######################################################### FINAL TEST #############################################################
##################################################################################################################################

final_test_loss, final_test_accuracy = model.evaluate(final_test_data, final_test_labels)
print(f"Final Test Loss: {final_test_loss}, Final Test Accuracy: {final_test_accuracy}")

save_model(model, f"models/{model_name}.h5")
print(f"Saved {model_name}.h5")

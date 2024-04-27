from network_funcs import define_model, train_model, save_model
from tensorflow.keras.models import load_model
import numpy as np
import os

file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(file_dir)

model_name = "testModel"
load_model = False

loaded_model = load_model('models/testModel.h5') if load_model else None
print("Loaded model") if loaded_model else None

num_rows = 100000
num_cols = 10

# Generate random data with a normal distribution
data = np.random.randn(num_rows, num_cols)
labels = np.zeros((num_rows, 3))

for i in range(len(labels)):
    index = np.random.choice([0, 1, 2], p=[1/3, 1/3, 1/3])
    labels[i, index] = 1

network = define_model()
epochs, batch_size = 1, 32
network = train_model(network, data, labels, epochs, batch_size)

new_data = np.random.randn(1, num_cols)
prediction = network.predict(new_data)

print(prediction)

save_model(network, f"models/{model_name}.h5")

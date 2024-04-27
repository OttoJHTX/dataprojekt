import numpy as np
from tensorflow import keras
from keras import layers



# def define_model():
#     # Define the model architecture
#     model = keras.Sequential([
#         layers.Dense(10, activation='relu', input_shape=(10,), name='layer1'),
#         layers.Dense(20, activation='relu', name='layer2'),
#         layers.Dense(20, activation='relu', name='layer3'),
#         layers.Dense(3, activation='softmax', name='layer4')
#     ])
#     optimizer = keras.optimizers.Adam(learning_rate=0.001)  # Set the learning rate here

#     model.compile(
#         optimizer=optimizer,
#         loss='categorical_crossentropy',
#         metrics=['accuracy']
#     )
    
#     return model

def define_model(learning_rate = 0.001, num_layers=3, num_nodes=[10, 20, 20], activation='relu', output_nodes=3):
    # Define the model architecture
    model = keras.Sequential()
    
    # Add input layer
    model.add(layers.Dense(num_nodes[0], activation=activation, input_shape=(num_nodes[0],), name='layer1'))
    
    # Add hidden layers
    for i in range(1, num_layers):
        model.add(layers.Dense(num_nodes[i], activation=activation, name='layer{}'.format(i+1)))
    
    # Add output layer
    model.add(layers.Dense(output_nodes, activation='softmax', name='output_layer'))
    
    optimizer = keras.optimizers.Adam(learning_rate=learning_rate)  # Set the learning rate here

    model.compile(
        optimizer=optimizer,
        loss='categorical_crossentropy',
        metrics=['accuracy']
    )
    
    return model


def train_model(model, data, labels, epochs, batch_size):
    model.fit(
        data,
        labels,
        epochs=epochs,
        batch_size=batch_size
    )
    
    return model

def save_model(model, file_path):
    model.save(file_path)
    print("Model saved successfully at:", file_path)


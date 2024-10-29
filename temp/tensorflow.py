import numpy as np
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Reshape, Flatten, Dropout, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import LeakyReLU
import pickle
import os

# Mechanistic simulation function (MDASimulation) code goes here...
# You can place the full code from the mechanistic simulation shared above here, if needed.

####################################################################################

# GAN Configuration
latent_dim = 100  # Size of latent space (random noise vector)
adam = Adam(learning_rate=0.0002, beta_1=0.5)

# GAN Generator
def build_generator(latent_dim):
    model = Sequential([
        Dense(256, input_dim=latent_dim),
        LeakyReLU(alpha=0.2),
        Dense(512),
        LeakyReLU(alpha=0.2),
        Dense(1024),
        LeakyReLU(alpha=0.2),
        Dense(5000, activation='tanh')  # Assuming 5000 features in simulated data
    ])
    return model

# GAN Discriminator
def build_discriminator(input_shape=(5000,)):  # Adjust input shape based on data
    model = Sequential([
        Dense(512, input_shape=input_shape),
        LeakyReLU(alpha=0.2),
        Dropout(0.4),
        Dense(256),
        LeakyReLU(alpha=0.2),
        Dropout(0.4),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
    return model

# GAN Model (Discriminator not trainable in GAN model)
def build_gan(generator, discriminator):
    discriminator.trainable = False
    gan_input = Input(shape=(latent_dim,))
    x = generator(gan_input)
    gan_output = discriminator(x)
    model = Model(gan_input, gan_output)
    model.compile(optimizer=adam, loss='binary_crossentropy')
    return model

# Load Real Dataset (single-cell data)
def load_real_data(filepath):
    with open(filepath, 'rb') as f:
        real_data = pickle.load(f)
    # Normalize and shape the data if necessary
    return real_data

# GAN Training Function
def train_gan(generator, discriminator, gan, epochs=10000, batch_size=64, real_data=None, sim_data_path=None):
    half_batch = batch_size // 2

    for epoch in range(epochs):
        # 1. Train the Discriminator
        # Generate "real" examples (from the mechanistic simulation or real data)
        real_samples = real_data[np.random.randint(0, real_data.shape[0], half_batch)]
        
        # Generate "fake" examples
        noise = np.random.normal(0, 1, (half_batch, latent_dim))
        generated_samples = generator.predict(noise)

        # Train on real and fake data
        d_loss_real = discriminator.train_on_batch(real_samples, np.ones((half_batch, 1)))
        d_loss_fake = discriminator.train_on_batch(generated_samples, np.zeros((half_batch, 1)))
        d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)

        # 2. Train the Generator (via GAN model, where discriminator weights are frozen)
        noise = np.random.normal(0, 1, (batch_size, latent_dim))
        g_loss = gan.train_on_batch(noise, np.ones((batch_size, 1)))

        # Print progress
        if epoch % 100 == 0:
            print(f"{epoch} [D loss: {d_loss[0]:.4f}, acc.: {100 * d_loss[1]:.2f}%] [G loss: {g_loss:.4f}]")
            
    return generator

# Main Script to Run Mechanistic Simulation and Train GAN
def main(sim_data_output_folder, real_data_path, epochs=10000):
    # Generate simulated data with mechanistic simulation
    patSeq = "ATCG" * 10000  # Placeholder sequence
    matSeq = "CGTA" * 10000  # Placeholder sequence
    simulated_data = MDASimulation(
        patSeq, matSeq, output_folder=sim_data_output_folder
    )

    # Load real single-cell data
    real_data = load_real_data(real_data_path)

    # Initialize GAN
    generator = build_generator(latent_dim)
    discriminator = build_discriminator()
    gan = build_gan(generator, discriminator)

    # Train GAN with simulated and real data
    trained_generator = train_gan(
        generator, discriminator, gan,
        epochs=epochs, real_data=real_data
    )

    # Save the trained generator
    generator.save(os.path.join(sim_data_output_folder, "trained_generator.h5"))

# Specify paths for output and real data
output_folder = "output"  # Directory to save generated data
real_data_file = "real_single_cell_data.pkl"  # Path to real single-cell data

# Run main function to simulate data and train GAN
main(output_folder, real_data_file, epochs=10000)

# Siamese Net

This repository provides a TensorFlow/Keras implementation of a Siamese neural network for comparing paired feature representations from the same neurons across different data modalities. The project is designed for experiments on single-cell in-vitro recordings, where feature vectors extracted from two related views are compared through a shared embedding space.

## Overview

The main idea is to learn embeddings that make corresponding samples from two modalities more similar, while pushing unrelated pairs apart. The implementation uses a Siamese architecture with:

- a configurable base network for each input branch,
- a cosine similarity head for comparing embeddings,
- optional dropout and regularization,
- custom loss functions for weighted regression tasks.

## Project structure

- [src/Siamese_net/model.py](src/Siamese_net/model.py): implementation of the Siamese model class.
- [src/Siamese_net/loss.py](src/Siamese_net/loss.py): weighted loss functions such as MAE, MSE, RMSE, and risk-aware variants.
- [src/Siamese_net/utils.py](src/Siamese_net/utils.py): small utility helpers.
- [experiments](experiments): Jupyter notebooks demonstrating usage and analysis workflows.
- [data](data): data directory used by the experiments.

## Installation

1. Clone the repository:

   ```bash
   git clone <repository-url>
   cd Siamese_net
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

   On Windows PowerShell:

   ```powershell
   .venv\Scripts\Activate.ps1
   ```

3. Install the package and its dependencies:

   ```bash
   pip install -e .
   ```

4. Make sure TensorFlow and supporting packages are available:

   ```bash
   pip install tensorflow h5py numpy
   ```

## Quick start

The core model is exposed through the SiameseModel class.

```python
from Siamese_net.model import SiameseModel

input_dim1 = (600,)
input_dim2 = (600,)

model = SiameseModel(
    input_dim1=input_dim1,
    input_dim2=input_dim2,
    base_dims=(256, 128),
    embedding_dim=64,
    dropout_rate=0.2,
)

model.compile(optimizer="adam", loss="mse")
model.summary()
```

You can then train the model on paired feature arrays for your dataset.

## Usage notes

- The model expects paired inputs representing two related views of the same sample.
- Input dimensions should match the feature dimensionality of each modality.
- The notebooks in [experiments](experiments) are a good starting point for adapting the workflow to your own data.
- The losses in [src/Siamese_net/loss.py](src/Siamese_net/loss.py) can be used when you want weighted or risk-aware objectives.

## License

This project is distributed under the MIT license, as declared in the package metadata.


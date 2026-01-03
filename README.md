# Disentangling and Interpreting Internal Representations of ProtT5 Using Sparse Crosscoders

This repository contains the code accompanying the paper:

**Disentangling and Interpreting Internal Representations of ProtT5 Using Sparse Crosscoders**

---

## Overview

Protein language models (pLMs) such as **ProtT5** achieve state-of-the-art performance across a wide range of biological tasks, yet their internal representations remain difficult to interpret. While ProtT5 embeddings encode rich information about protein structure and function, the biological concepts captured within the model’s hidden layers are largely implicit.

This project applies **sparse crosscoders** to the internal activations of the ProtT5 encoder to enable **mechanistic interpretability across network depth**. By learning a **single shared dictionary of sparse features across all layers**, the model disentangles biologically meaningful concepts and tracks how they are represented and transformed throughout the network.

Unlike layer-wise sparse autoencoders, sparse crosscoders:

* Capture **cross-layer features**
* Reduce redundant representations
* Allow meaningful connections to “jump” across identity transformations
* Enable more robust feature-level interpretation

---

## Method Summary

* **Model**: ProtT5-XL-U50 encoder
* **Input**: Hidden activations from all encoder layers (excluding the final layer)
* **Architecture**: Sparse crosscoder with a shared feature dictionary
* **Loss**: L2-of-norms sparsity objective
* **Dataset**: Millions of protein sequences sampled primarily from UniRef50 (optionally augmented with BFD)

---

## Analysis and Evaluation

The repository supports three main analysis pipelines:

### 1. Biological Validation

Learned features are evaluated using a Swiss-Prot–based concept annotation pipeline. Protein-level annotations are converted to amino-acid–level concepts, and alignment between crosscoder features and known biological annotations (e.g. active sites, binding motifs, domains) is quantified using F1 scores.

### 2. Automated Feature Interpretation

Features are automatically annotated using large language models (LLMs) based on sequences that maximally activate each feature, following prior work in mechanistic interpretability.

### 3. Application: Interpretable Directed Evolution

As an exploratory use case, the sparse crosscoder is used as an **interpretable signal** for guiding mutations in an active learning-assisted directed evolution (ALDE) setting. Instead of predicting fitness directly, mutations are selected to maximize activation of features associated with desirable biological properties.

---

## Repository Structure (Planned / Expected)

```text
.
├── data/                # Dataset preparation and preprocessing
├── protT5/              # ProtT5 inference and activation extraction
├── crosscoder/          # Sparse crosscoder architecture and training
├── analysis/            # Feature analysis and biological validation
├── interpretation/      # Automated feature interpretation with LLMs
├── evolution/           # ALDE experiments and mutation steering
└── scripts/             # Training and evaluation scripts
```

---

## Goals

* Disentangle ProtT5’s internal representations into human-interpretable features
* Identify biologically meaningful cross-layer feature circuits
* Enable mechanistic insight into how protein properties are encoded
* Explore interpretable steering of protein sequences for downstream design tasks

---

## Citation

If you use this code, please cite the associated paper:

```bibtex
@article{tawana2025prott5crosscoders,
  title   = {Disentangling and Interpreting Internal Representations of ProtT5 Using Sparse Crosscoders},
  author  = {Tawana, Sohrab},
  year    = {2025}
}
```

---

## Acknowledgements

This work builds on prior research in protein language models and mechanistic interpretability, including ProtT5, sparse autoencoders for pLMs, and sparse crosscoders originally proposed for large language models.

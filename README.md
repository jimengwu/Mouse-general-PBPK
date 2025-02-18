
# MLR Model for Nanoparticle Biodistribution

## Overview

This project introduces a predictive framework for nanoparticle biodistribution using **Physiologically Based Pharmacokinetic (PBPK) Modeling**, **Quantitative Structure-Activity Relationship (QSAR)** principles, and **Multivariate Linear Regression (MLR)**. The model aims to predict nanoparticle biodistribution based solely on their **physicochemical properties**, offering an alternative to traditional animal-based studies and enhancing the Safe and Sustainable by Design (SSbD) frameworks for nanoparticle risk assessment.

The approach integrates **Bayesian analysis** with **Markov Chain Monte Carlo (MCMC)** simulations to fit PBPK models and generate kinetic parameters, providing a robust, non-animal alternative for early-stage nanoparticle evaluation and design.

## Model Description

- **PBPK Modeling**: Physiologically based pharmacokinetic (PBPK) models describe how a nanoparticle distributes and behaves within a biological system. These models rely on physiological parameters (e.g., blood flow, organ size) and physicochemical properties of nanoparticles to predict their pharmacokinetic profile.

- **Bayesian Analysis with MCMC**: Bayesian inference and MCMC simulations were employed to fit the PBPK model to experimental data, providing parameter estimates that capture uncertainty and variability in the model.

- **Multivariate Linear Regression (MLR)**: Used to establish relationships between nanoparticle properties (such as size, zeta potential, coating) and their biodistribution. MLR helps in identifying the most influential predictors for nanoparticle pharmacokinetics.


## Installation & Usage

1. Clone this repository to your local machine.
2. Install required R packages, indicating in the model files.   
3. Load the model script into your R environment.


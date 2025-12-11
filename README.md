[![DOI](https://zenodo.org/badge/340074964.svg)](https://zenodo.org/badge/latestdoi/340074964)

  
# PINF

## Learning Dynamical Transitions with Physics-Informed Normal-Form Priors 


This repository contains the PINF method applied to the different dynamical systems models.

* [Summary](#Summary)
* [Scripts](#Scripts)
* [Working example: saddle node](#Working-example-saddle-node)
* [Deep Learning for Normal Form Identification](#deep-learning-for-normal-form-identification)
* [Scripts ML](#scripts-for-machine-learning)
* [Contributing](#Contributing)
* [License](#License)
* [Requirements](#Requirements)



## Summary

#### Authors
Juan Pablo Munoz Diaz,  Juan P. Bernal, Alberto Maillo, Ali Balubaid, Subash Balsamy, Sumeer Khan, Vincenzo Lagani, David Gomez Cabrero, Narsis Kiani, Jesper Tegner

#### Abstract

Complex systems, including biological networks, ecosystems, climate dynamics, neural populations, and synthetic gene circuits, often undergo transitions between distinct dynamical states. Modeling these shifts is challenging when governing equations are unknown or data are sparse. Here, we introduce Physics-Informed Normal Forms (PINF), motivated by converging evidence for low-dimensional manifolds in empirical data, nonlinear dynamical models, and deep learning representations. PINF embeds canonical bifurcation structures as inductive biases to infer hidden control parameters directly from time-series observations. By leveraging the manifold hypothesis, PINF constrains learning to physically meaningful dynamics, enhancing interpretability and predictive robustness. We validate PINF on standard bifurcation models and demonstrate generalization to nonlinear systems without access to underlying equations. PINF further classifies composite bifurcations and infers parameters of interacting hybrid normal forms through differentiable simulation. PINF establishes a robust, interpretable, and generalizable tool for data-driven discovery of effective low-dimensional models of dynamical transitions, bridging nonlinear dynamics and machine learning to advance physics-grounded modeling of complex systems.

For additional information, please refer to the [Publication](https://arxiv.org/abs/2304.02443)

#### PINF Workflow
![Figure 1](https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/PINFWorkflow.png)


## Scripts
  
PINF is entirely developed in MATLAB, tested version 2020. There is no limitation in OS and no specific computational requirement.

First, clone this repository and set the working directory:

```
#clone repository
git clone https://github.com/munozdjp/PINF.git

# set working dir
cd PINF_Code/<script>
```

Each script generates the results for:
  * Prediction of learned variable: a figure including the learned variable and the ground truth. 
  * Noise analysis: a figure with the performance of PINF under the influence of noise with 5 different variances. 

The list of the scripts that you can find on this GitHub:

1- [Saddle-node bifurcation](https://github.com/munozdjp/PINF/blob/main/PINF_Code/Scripts_Models/Saddle_V4.m): prediction for switching fix point. 

2- [Pitchfork-bifurcation](https://github.com/munozdjp/PINF/blob/main/PINF_Code/Scripts_Models/Pitchfork_V4.m): prediction for growing fix point dynamic.

3- [Hopf-Bifurcation](https://github.com/munozdjp/PINF/blob/main/PINF_Code/Scripts_Models/Hopf_V4.m): learned hidden variable on oscillator to fixpoint dynamic. 

4- [Lorentz](https://github.com/munozdjp/PINF/blob/main/PINF_Code/Scripts_Models/Lorentz_V4.m): learned hidden variable on different regimes of chaotic system. 

5- [FitzHugh-Nagumo](https://github.com/munozdjp/PINF/blob/main/PINF_Code/Scripts_Models/FitzHugh_V4.m): to run the script, please download this [dataset](https://figshare.com/s/a1e42815cf89b4eff381) and save it in the same directory as the script.

The extra scripts are complementary validation of PINF method. 
  
Working example: saddle node 
----------------------------

The saddle node bifurcation is a dynamical system that shows a switching point dynamic, from a lower state to a high upper state. 

Here you can find the complete script of [saddle node](https://github.com/munozdjp/PINF/blob/main/PINF_Code/SaddleNodeLeft2Rigth.m) 
  
First, we need to define the static variables. We recommend using the default values:

```
%timestep
dt=0.1; 
%final time point
Initinterval = 16;
%initial time point
neg_inter = 0;
```  

The output is the following 5 figures:

**Fig 1:** shows the ground truth of the alpha polynomial dynamic.
  
<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig1_saddlenode.png" alt="My Image" style="width: 399px; height: 300px;">
    
**Fig 2:** state variable when a hidden variable changes over time.
  
<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig2_saddlenode.png" alt="My Image" style="width: 399px; height: 300px;">
  
**Fig 3:** state variable vs hidden variable, PINF prediction. 
  
<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig3_saddlenode.png" alt="My Image" style="width: 646px; height: 300px;">

**Fig 4:** shows how PINF works when there are different kinds of noise disturbances.
  
<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig4_saddlenode.png" alt="My Image" style="width: 695px; height: 300px;">
  
**Fig 5:** shows how one sample of noise was analyzed.    

<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig5_saddlenode.png" alt="My Image" style="width: 399px; height: 300px;">
  
**Fig 6:** error distribution for different noise samples.        

<img src="https://github.com/munozdjp/PINF/blob/main/figures%20saddle%20node/fig6_saddlenode.png" alt="My Image" style="width: 417px; height: 300px;">

## Deep Learning for Normal Form identification

### Setup Instructions

To access the [Python scripts](https://github.com/munozdjp/PINF/tree/main/PINF_Code/DeepLearningNormalForms), modify your current working directory by executing the following command:

```
# Set working directory
cd PINF_Code/DeepLearningNormalForms
```
### Create Conda Environment

Initialize the environment using the provided environment.yml file:

```sh
conda env create -f environment.yml
```

To activate the newly created environment named "SVM_RF_MLP", use the following command:

```
conda activate SVM_RF_MLP
```

### Download Data Set

Obtain the [data set of the normal forms](https://doi.org/10.6084/m9.figshare.25664514.v1) from Figshare and ensure it is placed in the current directory or same directory of the scripts.
 
### Scripts for Machine Learning

In the directory of [Deep Learning Normal Form](https://github.com/munozdjp/PINF/tree/main/PINF_Code/DeepLearningNormalForms), you will find various scripts that generate their respective results:

Run each script with the following command:

```
python <scriptname>
```

1- [Classification Comparison: MLP, RF, KNN](https://github.com/munozdjp/PINF/blob/main/PINF_Code/DeepLearningNormalForms/ClassificationComparisonRF_MLP_KNN.py): generate a comparison of accuracy for a library of four normal forms and compare the accuracy between Multilayer Perceptron (MLP), Random Forest (RF), and K-Nearest Neighbor (KNN).

2- [Noise robustness for Random Forest on Normal Forms](https://github.com/munozdjp/PINF/blob/main/PINF_Code/DeepLearningNormalForms/RForest_Noise.py): generate a comparison of classification using different noise additions for a library of four bifurcation normal forms: Saddle, Pitchfork, Transcritical, and Hopf, plus a combination of two normal forms, Hopf-Saddle and Hopf-Pitchfork.

3- [Downsampling Accuracy](https://github.com/munozdjp/PINF/blob/main/PINF_Code/DeepLearningNormalForms/DownsamplingComparison.py): generate a downsampling size of the training data from a reduction of 2% to 0.5% of the original size.

Results from these scripts are stored in a new directory titled "figures results deep learning"

**Working example ML**
----------------------

From the [Comparison Classifier Script](https://github.com/munozdjp/PINF/blob/main/PINF_Code/DeepLearningNormalForms/ClassificationComparisonRF_MLP_KNN.py) you will obtain the following plot 

**Fig 1:** show a comparison of accuracy and noise robusteness for three algorithm MLP, RF and KNN. 

<img src="https://github.com/munozdjp/PINF/blob/main/figures%20classification%20ML_RF_KNN/accuracy_vs_noise_level.png" alt="My Image" style="width: 700px; height: 420px;">


## Contributing
Simply fork the repository and make the appropriate changes to contribute to PINF. Once complete, send a pull request and we will evaluate your contributions.

## License
PINF is licensed under the GNU License. Refer to LICENSE for additional details.
  
## Requirements

Before running the Matlab Scripts, please ensure you have the following MATLAB toolboxes installed:

- Optimization Toolbox
- SimBiology
- Statistics and Machine Learning Toolbox

These toolboxes are essential for the execution of the PINF method scripts. They provide functions and features that are not available in the standard MATLAB environment.

To install these toolboxes, you can use the Add-On Explorer in MATLAB. Simply go to the Home tab, and in the Environment section, click Add-Ons > Get Add-Ons. Search for the required toolboxes and install them.


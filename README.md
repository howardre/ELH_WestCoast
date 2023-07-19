# Projecting Species Distributions Using Size-Structured and Size-Aggregated Models

The is part of an NSF funded project that will create species distribution maps using all life stages of West Coast fishes. Portions of this work will comprise a chapter in my PhD dissertation.

### Data
Multiple data sources have potential applications to this project. Currently we are using:
- The NOAA **Rockfish Recruitment and Ecosystem Assessment Survey**, data for which are available upon request. This survey largely target pre-recruits but it initially started with a focus on rockfish. The survey operates in the late spring/early summer off the California Coast.
- The NOAA **Pre-recruit Survey**, data for which are available upon request. This survey follows the same protocols as the Rockfish Recruitment and Ecosystem Assessment Survey, but operates off the Oregon and Washington Coasts. 


### Methods
#### Variable Coefficient Generalized Additive Models (VGAMs)
VGAMs allow a coefficient to vary with change in another variable. In this application, the coefficient is location of a species which varies with a given environmental variable. Our initial model use climate indices as a stand-in for this variable while future models will likely use sea surface temperature.
- [VGAMs](code/RREAS_model_exploration.RMD/) applied to the RREAS data

#### Predictive-Process Generalized Linear Mixed Models (GLMMs)
This is a work in progress and will include future applications of the [sdmTMB package](https://github.com/pbs-assess/sdmTMB).
- Code for these models is found [here](code/sdmTMB_models.R/)

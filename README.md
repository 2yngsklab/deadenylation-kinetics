# deadenylation-kinetics
Tools related to mathematical models and their applications for estimating the kinetics of a single deadenylation events.

- Example gel image ([gel.tif](gel.tif)) and pre-processed text ([converted.txt](converted.txt)). Briefly, a synthetic RNA sequence of seven nucleotides (5′-UCUACAU-3′) and a 20 adenosine stretch (A20) was used as a substrate to assess the deadenylation reaction of the full human CNOT complex. Other synthetic RNAs such as A1, A0, and UCU were used as size markers. The first two lanes are the marker and control lanes used during the image processing steps. The neighboring gel lanes correspond to deadenylation reactions after predefined time points such as 2 minutes, 4 minutes, 6 minutes, and 8 minutes.

- Example script for data analysis ([script.R](script.R))

- Example script for computer simulation ([simulate.R](simulate.R))

## Custom R functions ([deadenylation.R](deadenylation.R))

### Parameter estimation
`rxnrate`: The `rxnrate` (Reaction Rate) function defines the mathematical model for a sequence of deadenylation events and is used as input for the ODE solver. As input, this function takes three parameters: 
1. list of the time points of interest,
2. concentration of each RNA species, and 
3. the rate parameters for each deadenylation event.

`simulation`: The `simulation` function computes the concentration of each RNA species for a given time point and solves the ODE defined by the `rxnrate` function. As input, this function takes three parameters: 
1. the rate parameters for each deadenylation event, 
2. the length of synthetic RNA (i.e., the number of RNA species), and 
3. the list of time points of interest.

`ssq`: The `ssq` (Sum of Squared Residuals) function calculates the difference between the predicted and observed concentrations of each RNA species and is key for the LM routine in the `minpack.lm` R package. As input, this function takes the preprocessed image data and the rate parameters for each deadenylation event. The function generates the predicted concentration by the simulation function.

### Data visualization
`GetMetaPlot`: The `GetMetaPlot` function returns a single line plot that represents the cumulative intensity value for each pixel row of the text file (.txt) converted from the image processing software. The vertical lines in red indicate the local minima computed based on the meta line plot which may include false positives. This metaplot and the following lane plots are used to curate the local minima manually and include false negatives.

`GetLanePlot`: The `GetLanePlot` function returns a multi-panel plot that explores the intensity values for each lane. Each lane may be selected by the user. The vertical lines in red indicate the manually curated local minima used to bin individual bands and digitize the preprocessed text file (.txt). 

`GetHeatMap`: The `GetHeatMap` function returns a heatmap that represents the digital version of the original gel image to improve data visualization and interpretation. Analogous to the *in vitro* experiment, the reaction time points are in the x-axis, and the individual RNA species are shown in the y-axis. This digitization corrects the nonlinear effects of electrophoresis and enables linear interpretation of the deadenylation reaction. The heatmap values are the normalized intensity values. 

`GetStepPlot`: The `GetStepPlot` returns a step plot that represents the parameter estimate of deadenylation at each nucleotide position. The x-axis of the step plot indicates the nucleotide position, and the y-axis indicates the estimated rate of deadenylation. The error bars represent the standard error of the parameter estimation.

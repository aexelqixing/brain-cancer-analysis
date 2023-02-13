# Simulations

I made multiple different simulations, the first ones in Java because Java was more familiar to me at the time, and the more complex ones in NetLogo. 

# Java Models

These were created with the intent of having multiple Cell classes in a single Environment, and are explained in more detail in this folder. 

# oxygen-model.nlogo

This model simply shows how oxygen is introduced in the environment, with some patches randomly absorbing oxygen and others diffusing it. It shows the first component of all the other models and how the differential equation works. 

# SPARC vs TGFBI

These models showcase how the expressions of SPARC and TGFBI affect the movement and division rate of cells. These rules were made by reading past literature, such as [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6859668/), [this paper](https://pubmed.ncbi.nlm.nih.gov/22911700/) and [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5294964/). The basics of EMT were found in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2689101/). 

# spreading cancer cells

This is the last model made for this project, and took all of the findings from previous models and incorporated it into one about spreading from the primary site. The basics of EMT and hypoxia data sets were utilized for this model. 

# CreatingGraphs.ipynb

Since I had multiple csv files, I turned them all into multiple dataframes and created multiple graphs in python using the package matplotlib. 

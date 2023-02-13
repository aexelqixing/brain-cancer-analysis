# Spreading Cancer Cell Model

This file goes more into detail of how the spreading cell models work. 

## healthy-cell-model.nlogo

This only creates a simulation of healthy cells, adding onto the oxygen-model in the previous folder. These were heavily influenced by [this paper](https://scholarscompass.vcu.edu/cgi/viewcontent.cgi?article=3830&context=etd).

## spreading-cancer-cell.nlogo

Using the simulation of healthy cells, I only modified it a little bit to include what happens when cancer cells spread throughout. The code itself has more details on how it was implemented. 

## CSCs & CSNCs

50 simulations were ran for the spreading cancer cell model, and different csvs were created for the cell cycle distribution for cancer stem cells and cancer non stem cells. These were then averaged together to make a stacked area chart. 
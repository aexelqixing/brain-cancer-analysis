# SPARC vs. TGFBI

This README focuses more on the individual files themselves.

## SPARC-factor.nlogo

SPARC treatment in previous experiments mostly affected cell proliferation and migration, so in this file, the SPARC expression mostly affected how much cells could divide (they divided more if their SPARC expression was greater) and how much they could move (more if SPARC was greater). The SPARC expression changed randomly at every timestamp to make use of random mutations, but not enough to change the actual behavior of the cell. 

## TGFBI-factor.nlogo

TGFBI factors are associated with cancer progression and promote cancer growth within the cell, so in this file, TGFBI affected how much oxygen a cell could consume and how fast it could move. Because a cell could grow quicker the more TGFBI was expressed, the faster it could divide, but it was not a direct influence. 

## control, SPARC and TGFBI folders

The individual folders contain csv files of 50 simulations each of the population of cancer non-stem cells and cancer stem cells. These were then averaged over all 50 simulations and turned into a graph. 
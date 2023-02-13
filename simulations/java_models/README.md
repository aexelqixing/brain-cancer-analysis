# Java Models

Here I will explain the model in much more detail. 

## Cell.java

This class only contains the multiple attributes for each Cell, which can be seen in the actual file itself. 

## Environment.java

This class actually makes the environment for all the cells, and starts with all of the Cells in the middle. Then, they can divide randomly and also migrate to other places, but only if those spaces are empty. This model was influenced heavily by [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6587968/). However, it was ultimately decided to discard this project in Java and move on to a more niche programming language for agent-based modeling - NetLogo. 
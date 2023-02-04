# Python Notebook Scripts

This was where I created all of the support vector machines and compared them with other machine learning algorithms. 

All files were the same, except they compared the different subtypes of brain cancer with normal tissues. For example, the file 'shap_with_astrocytoma.ipynb' was an astrocytoma classifier, and 'shap_with_GBM.ipynb' was a GBM classifier. 

## Reading in the Data

The [pandas](https://pandas.pydata.org) library was used in order to read in the csv files and wrangle them in a way suitable for a machine learning algorithm to learn from them. 

## Splitting the Data and Creating Algorithms

The [scikit-learn](https://scikit-learn.org/stable/) library was then used to do a 70-30 train-test-split with the data and create many machine learning models, such as:
- support vector machines
- random forest classifier
- logistic regression
- k nearest neighbors
- decision tree classifier
- multilayer perceptron neural network

This library was also used to measure the performance of each algorithm using metrics such as 
- F1 score: harmonic mean of precision and recall, and therefore provides a balance between them
- Accuracy: the fraction of correct predictions made by the model over all predictions
- Recall: the fraction of cancer samples that are correctly detected by the model
- Precision: the fraction of correct cancer predictions made by the model out of total cancer predictions
- ROC AUC score (Receiver Operating Characteristic, Area Under the Curve): plots the true positive rate against false positive rate, and the area under the curve represents the model's performance

## Selecting the Best Parameters

The GridSearchCV function from scikit-learn is then used to select the best parameters for each model. 

## Comparing Performances

The [matplotlib](https://matplotlib.org) library is then used to plot all the results in a readable table than can be easily understood. 

## References

These notebooks were influenced by the tutorial notebooks on how to use SHAP, which was not utilized in this project. 
- https://shap.readthedocs.io/en/latest/example_notebooks/tabular_examples/model_agnostic/Iris%20classification%20with%20scikit-learn.html
- https://www.kaggle.com/code/joshuaswords/predicting-a-stroke-shap-lime-explainer-eli5
---
layout: default
title: The Model
---

From the Utah State University Wetland Ecology Lab.

Coming Soon!

## Model Framework

The framework for this model is essentially an integrated population model (IPM) characterizing plant life stage transitions. Each life stage transition is modeled using a Bayesian hierarchical model of logistic regression. The response variable for each hierarchical model is a matrix of transitions (represented as binary - 0s and 1s) for each individual plant at each timestep.

For more details, see the following contents:

- [model parameterization]({{ site.baseurl }}/model/ModelSpecs)
- [model validation]({{ site.baseurl }}/model/ModelValidation)
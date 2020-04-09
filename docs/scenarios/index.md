---
title: "Model Scenarios"
---

## The Data

The seed-based wetland restoration model is based on plants grown under controlled experimental conditions. These conditions include 2 moisture treatments to mimic the following sets of field conditions:

- dry early seeding
- wet early seeding
- dry late seeding
- wet late seeding
- dry climate change scenario
- wet climate change scenario

Data on seedling traits, weekly survival, and end-of-season cover were then collected for all combinations of species and sources grown under the above conditions. For more details on how this data was collected, [see here](link to page with more details on data collection).

## Model Basics

The above data were used to build a model which simulates weekly rates of seed germination, seedling survival, and adult survival throughout a 12-week growing period after seeds are sown. These transition rates vary by our *primary predictors*: seeded species, source population, moisture level, and temperature. For more details on the model, [see here](link describing model in full detail, including equations). The model was validated using a statistical method called cross-validation, where a subset of the original data is held out of the analysis, then predictions from the model are compared to these data. The model was further field-verified using data from plants grown in both greenhouse and field experiments. More details on model validation can be found [here](link to detailed model validation page).

## Model Outputs

For each scenario, the model outputs include estimated end-of-season cover for each species-source combination and a graph of population composition at each week for each species-source combination. These outputs were derived by running the model to predict transition rates and end-of-season cover based on each species, source, moisture level, and temperature combination.

Outputs for each scenario can be found at the links below:

- [dry early seeding]({{ site.baseurl }}/SBWR/scenarios/DryEarlySeeding)
- [wet early seeding]({{ site.baseurl }}/SBWR/scenarios/WetEarlySeeding)
- [dry late seeding]({{ site.baseurl }}/SBWR/scenarios/DryLateSeeding)
- [wet late seeding]({{ site.baseurl }}/SBWR/scenarios/WetLateSeeding)
- [dry climate change scenario]({{ site.baseurl }}/SBWR/scenarios/DryClimateScenario)
- [wet climate change scenario]({{ site.baseurl }}/SBWR/scenarios/WetClimateScenario)

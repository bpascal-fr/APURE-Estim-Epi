# Covid19 new infection counts time series

## Johns Hopkins University repository for worldwide data


Reported new cases collected by National Health Authorities of 200+ countries made available by Johns Hopkins University [JHU](https://coronavirus.jhu.edu/).

> *The last release has been done on October 3, 2023.*

`COVID-19-JH_confirmed-worldwide_stored.csv` downloaded from [JHU repository](https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv) on *March 27, 2024*.


## Santé Publique France data

Reported new cases in all French departments collected by Regional Health Agencies (*Agences Régionales de Santé*, ARS) made available by the French government through the SI-DEP framework on the portal [data.gouv.fr/](https://www.data.gouv.fr/fr/datasets/donnees-de-laboratoires-pour-le-depistage-a-compter-du-18-05-2022-si-dep/).

> *Daily new cases from May 13, 2020 to June 27, 2023.*

`Covid_French_Departments_stored.csv` downloaded from [Santé Publique France repository](https://www.data.gouv.fr/fr/datasets/r/426bab53-e3f5-4c6a-9d54-dba4442b3dbc) on *March 31, 2024*.

## French departments connectivity

The multivariate, spatially regularized, reproduction number estimators rely on a prior spatial connectivity pattern. 
In the present work, French departments are considered connected if they share a terrestrial border, which has been encoded manually by the contributors in

> `French_Contiguous_Departments.mat`: adjacency matrix encoding the contiguity of the 96 French metropolitan departments.

Overseas French territories are considered not connected to any other French department.

## Synthetic example

Synthetic piecewise linear reproduction coefficient of length T = 300 used as ground truth in Section IV of the paper:
> Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a Nonstationary Autoregressive Model for Data-Driven Reproduction Number Estimation. *arXiv preprint*. [arXiv:]().
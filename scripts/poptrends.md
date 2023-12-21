---
title: "Population Trends Wadden Sea"
author: "Helmut Hillebrand"
date: '2023-12-21'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

The Wadden Sea is changing. Direct local human impacts such as dredging, fisheries and tourism coincide with regional stressors such as climate change and eutrophication. To understand the response of biota to these changes many projects and monitoring programs have gathered time series of communities. These data are mainly used to assess community wide metrics of biodiversity or to analyse population trends for certain species of interest. However, biodiversity metrics offer only a limited view on changes, as some only capture net changes (i.e. more immigrations than extinctions reflected by species richness), whereas others only reflect responses of dominant species. If, for example, all species would decline at the same rate without already going extinct, most biodiversity assessments would not pick that up. On the other hand, focusing on a subset of species may not capture the full range of population trends. There is a risk we are mainly focusing on extremes, i.e. species that are rapidly declining or have become dominant or even invasive. 

In order to gain a more comprehensive view of the shifts in Wadden Sea biodiversity, we need a systematic and quantitative overview of individual species trends. These can be used to identify winners and loosers under current and past environmental change by comparing the trend per species across locations. By aggregating at higher taxonomic levels such as orders, classes or phyla, we can ask whether certain groups are especially endangered profiting, which also can enable addressing vulnerability. In the first approach, we will focus on simply analysing the trends per species, but in later iterations of the synthesis we can also compare this to potential drivers and co-dependencies between species (e.g. in a trophic context). 

## Workflow

Population trends can be linear or non-linear of different polynomial orders. As we want to synthesize across species and phyla, I suggest to limit the polynomial to second order. 


![Fig.1: Workflow of the PopTrend Synthesis project](figures/workflow.png)

### Install libraries
```{r, echo = TRUE, message=FALSE}
library(tidyverse)
library(vegan)
library(reshape2)
library(lubridate)
```


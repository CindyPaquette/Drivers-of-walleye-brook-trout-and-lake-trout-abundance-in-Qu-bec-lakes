# drivers of fish abundance in Quebec lakes

This repository contains code for analyzing drivers of walleye, brook trout and lake trout across Quebec lakes using Random Forest models.

---

##  Project Structure

├── Code_fishabundanceV2/ # R scripts for modeling and plotting

├── README.md # This file

---

##  Summary

This project evaluates the environmental drivers of fish abundance in Quebec lakes, focusing on:
- Map visualizations in Mercator projection (WGS 84)
- PCoA to represent variation in fish communities across lakes 
- Species relative abundance modeling using Random Forests (RF)

---

## Dependencies

This project was developed in R (version 4.4.1). Key packages used:

library(party) 
library(caret) 
library(tidyverse)
library(permimp)
library(pdp)
library(gridExtra)
library(sf)
library(ggspatial)
library(geojsonio)
library(readxl)
library(vegan)
library(BiodiversityR)
library(sp)
library(corrplot)

See Code_fishabundanceV2.R for full details.
---
## Data
The fish data are not openly available for resource/fisheries management purposes. Data can be obtained by contacting the Québec Ministère de l’Environnement, de la Lutte contre les Changements Climatiques, de la Faune et des Parcs. 

---
## Citation
If you use this code or data, please cite:
Paquette, C. et al. (2025). Drivers of walleye, brook trout, and lake trout abundance in Québec lakes. Canadian journal of fisheries and aquatic sciences. 
---
## License 
This project is licensed under the CC-BY License.

---

## Contact
Cindy Paquette 
Postdoctoral Researcher, Université du Québec à Trois-Rivières
Cindy.Paquette@uqtr.ca



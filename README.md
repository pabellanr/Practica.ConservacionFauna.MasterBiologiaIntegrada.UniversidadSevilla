# Practica.ConservacionFauna.MasterBiologiaIntegrada.UniversidadSevilla
Scripts para las prácticas de la asignatura "Gestión y Conservación de Fauna Terrestre y Marina" del Máster Universitario en Biología Avanzada: Investigación y Aplicación de la Universidad de Sevilla

1. Assessing the performance of protected area networks in representing a regional species pool

It computes the number of taxa represented in a reserve network for a given conservation target and assesses whether this level of representation is significantly lower or greater than expected by chance. An analysis for multiple thresholds to consider a grid cell as protected is undertaken. The inputs are: (1) matrix of presence/absence data by species and area; (2) a table with area codes and the percentage of overlap with protected areas for each area.

gap_analysis.r                                 matrix.txt                                cells.txt

Citation and details:

Abellán P, Sánchez-Fernández D. 2015. A gap analysis comparing the effectiveness of Natura 2000 and national protected area networks in representing European amphibians and reptiles. Biodiversity and Conservation 24: 1377-1390.

2. Assessing species’ representation in protected area networks

For each species, the level of representativeness in a protected area network is computed as the mean percentage of spatial overlap between those planning units in which the species occurs in the study area and the protected areas (grid cells). The inputs are: (1) matrix of presence/absence data by species and area; (2) a table with area codes and the percentage of overlap with protected areas for each area.

mean_percentage_overlap.r                                      matrix.txt                           cells.txt

Citation and details:

Sánchez-Fernández D, Abellán P. 2015. Using null models to identify under-represented species in protected areas: A case study using European amphibians and reptiles. Biological Conservation 184: 290–299.

3. Tutorial – Manipulating and analyzing spatial data in R
This tutorial provides an introduction to the manage and analysis of spatial data in R. It requires basic knowledge of R.

Tutorial.pdf                                                   data.zip

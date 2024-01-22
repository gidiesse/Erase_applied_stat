# ERASE 
Erase was a project developed in collaboration with Human Technopole for the Applied Statistics course (A.Y. 2022/2023) at Politecnico di Milano. The project was supervised by Doctor Michela Massi and Doctor Alessandro Ravenae. 

## Project Aim:
Recent research in archaegenomics (the combination of genetics and archeology) has revealed that modern-day Europeans carry genetic contributions from five main ancestral sources: Western European Hunter-Gatherers (WHG), Eastern European Hunter-Gatherers (EHG), Caucauses Hunter-Gatheres (CHG), Anatolian Neolithic farmers (Anatolia_N) and, in some cases, North African populations (North_Africa) (predominantly found in Southern Europe). The objectives of this study were two-fold:  
1. To develop a novel statistical approach for analyzing the evolution of these ancestral contributions over time starting from a PC analysis. Our aim is to explain how the mean composition of individuals varied through time and space.
2. To investigate the CHG ancestry in order to better understand its evolution from ancient to modern times. Caucasus Hunter Gatherer is an Ancestry characteristic of certain ancient populations such as Neolithic populations from Iran. Specifically, we are interested in investigating the spread of CHG across Eurasia to gain a clearer picture of migration patterns at the time.

## Dataset: 
The dataset is composed of two parts:  
1. **metadata file (data/erase_metadata.csv)** : this file was composed of 15 columns and 6580 rows, they most important are:
   - Macro_Period: the period the individual belonged to. We had 8 options: Neolithic Age (10000 - 2000 BC), Stone Age (4000 - 2000 BC), Bronze Age (3300 - 1200 BC), Iron Age (1200 - 550 BC), Classical Age (8 BC - 5 AD), Middle Ages (500 - 1500 AD), Modern Age (1500 - 1900 AD) and Present Day (now). 
   - Main_Ancestry: the 5 main ancestral sources: WHG, EHG, CHG, Anatolia_N and North_Africa.
   - Macroarea: the 12 main areas the samples came from: Middle East, Balkans, Northern Europe, Central Europe, Western Europe, Eastern Europe, Penninsular Italy, Continental Italy, Sardinia, Sicily, Northern Africa and Asia. 
   - long_num: the longitudianl coordinate
   - lat_num: the latitudinal coordinate
   - date: the date of the sample
2. **snp files (data/snp)**: these are a type of file to store genetic information, they were used to generate the PCAs.
3. *PCA files (data/PCA)**: these contains the first 100 PCA components

## Methods and Scripts
1. **Exploratory Analysis**: this is the preliminary investigation into the Metadata. We were particularly interested in figuring out where our NAs lay and if there was any way to interpolate them. We noticed that the features date and Macro_Period were linked and also that the features Macroarea and lat_num, lat_long were linked so we thought that for individuals where only one of the features in the pair was missing we could use the other to interpolate it.
Files used here: ```data/erase_metadata.csv```

2. **genome_eda**: in this script we began to explore the PCA scores. In particular we explored the relations between the first three scores, also differentiating between modern and ancient samples and by main ancestry. 
Files used here: ```/data/snp/adna_pca_newiid_2convertV3.bed```, ```data/PCA/adna_pca_newiid_2convertV3_100.evec```, ```data/PCA/adna_pca_newiid_2convertV3_100.eval```

3. **individuals_deconvolution**: In this script we explored different ways of creating a representation for the individuals based on their genetic composition. By exploiting the individuals with an assigned ancestry, we built a compositional representation for each individuals. We tried using Gaussian Mixture Models and a C-means clustering. Both of these methods did not give us good results and they were quickly abandoned. 
The method we used is the following: we computed the mean of the pca scores for each of the 5 ancestries, based on the labelled individuals, resulting in 5 individuals that were the "representatives" of their ancestry. These would be our centroids. 
After that, we took the remaining 6513 individuals without an ancestry and computed the cosine similarity
between their PCA scores and the PCA scores of each of the 5 "centroid" individuals. We then used a 
softmax to normalize the similarity scores so that now each individual had 5 scores all between 0 and 1
which summed to 1 - essentially each score represented how much the individual belonged to each ancestry. 
In this way we performed a deconvolution of the individuals based on the ancestry, obtaining compositional data for each individual.
Files used here: ```data/PCA/adna_pca_newiid_2convertV3_100.evec```, ```data/erase_metadata_interpolated.csv```
Output: ```./corr_convolution.csv``` this is a file where we take ```data/erase_metadata_interpolated.csv``` and add 6 additional columns: one for each ancestry containing a value between 0 and 1 denoting how much that individual "belongs" to the specific ancestry and one last one containing a number from 1 to 5 which says which ancestry we've assigned the individual to based on the highest score of the composition (1 - WHG, 2-EHG, 3-Anatolia_N, 4-CHG, 5-North_Africa). 

4. ***geo_comp***: This is the file where we fit spatial statistic models. The aim of our project is understanding how main ancestry evolved in the different macroareas over time. The compositional data we obtained belongs to the Aitchison Simplex, a constrained space with **n - 1 ** degrees of freedom (where **n** is the number of components in the mixture). 
We performed a CLR transformation to map the data back into a Euclidean space, then a PCA to find an orthonormal basis. In this space we carried out the geospatial analysis and finally used the inverse ILR transformation to go back to the  Aitchison Simplex, where it's possible to interpret the data as a composition. 
Files used here: ```./corr_convolution.csv```
Output: ```data/comp_time_space.csv``` this is a file that contains the result of the spatial statistics analysis; for each Macro_Area and Macro_Period the mean coefficients of the linear geospatial model is reported, resulting in a symbolic composition for an individual that lived in a specific Macro_Area and Macro_Period.

5. **many_plots**: This is the file where we made most of the plots used to represent our data on the poster
Files used here: ```./corr_convolution.csv```

6. **world_plots**: This is the file we used to make the plots where we plotted our observations on top of maps from google maps
Files used here: ```./corr_convolution.csv```

## Conclusion:
We came to 2 main conclusions: 
1. We noticed that as time increased the genetic composition of the individuals followed a general trend of flattening - with the predominance of different ancestries gradually decreasing. This was particularly prominent in areas with strong migration pattens eg Western Europe and much less prominent in isolated areas eg Sicily or Sardinia. We believe that this is due to an elevated mixture among individuals from different ancestries that crossed paths at different points in history.
2. We chose to focus specifically on the CHG (Caucasus Hunter Gatherer) ancestry, which was predominant in the region sorrounding the Caucasus. In the Neolithic, this ancestry was not very present in Eurasia and in fact was the perdominant ancestry in just one macroarea (Eastern Europe, where the Caucasus region is situated geographically). There was a very significant shift in the CHG composition during the Bronze Age, suggesting a high level of mixture between the CHG poplation and the other populations present in Europe. By the Classical Age, CHG was the prominent ancestry in six of our macroareas (out of twelve), showing a very significant expansion. Notably, even in the islands, which we saw were the most resistant to change in composition, CHG increased from 10-15% to 20-25% by the modern age.

## Final notes: 
- In the Plots folder, the plots from each script (except for scripts 5 and 6 which were exclusively made to generate plots and so it would have been counter-intuitive to save each plot) have been saved for easier access
- Two animations were made to showcase:
     1. Which was the predominant ancestry in each Macro_Area over time [youtube link](https://www.youtube.com/watch?v=xGQ-XL5xd68)
     2. How CHG composition changed in each Macro_Area over time [youtube link](https://www.youtube.com/watch?v=VSL53vpeXYM)
- The link to the pdf of our [poster](https://www.dropbox.com/scl/fi/0gbapieyx1j0i9q8etmwu/APPLIED-STATISTICS-FINAL-84.1-118.9cm.pdf?rlkey=316mg3jez6vws31b9vejmb3xv&dl=0)






# ERASE 
Erase was a project developed in collaboration with Human Technopole for the Applied Statistics course (A.Y. 2022/2023) at Politecnico di Milano. The project was supervised by doctor Michela Massi and Doctor Alessandro Ravenae. 

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

3. **individuals_deconvolution**: In this script we explored different ways of assigning ancestries to individuals who didn't have one (6513 of the 6580 individuals in our dataset do not have an ancestry assigned). We tried three different approaches. The first 2 were Gaussian Mixture Models and a C-means clustering (having such few samples with labeled ancestries - only 67 out of 6580 - we needed to use unsupervised methods). The code for this is at the end of the file but since both of these had disastrous results, they were quickly abandoned. 
Our last idea was to do a sort of deconvolution. Essentially we divided the 67 individuals into their 
own respective ancestries. We then computed the mean of the pca scores for each ancestry resulting 
into 5 individuals that were the "representatives" of their ancestry. These would be our centroids. 
After that, we took the remaining 6513 individuals without an ancestry and computed the cosine similarity
between their PCA scores and the PCA scores of each of the 5 "centroid" individuals. We then used a 
softmax to normalize the similarity scores so that now each individual had 5 scores all between 0 and 1
which summed to 1 - essentially each score represented how much the individual belonged to each ancestry. 
We then assigned each individual to the ancestry that had the biggest score. This approach worked rather 
well. In fact of the 67 individuals that already had ancestry labels, using this method 66 were assigned 
corrrectly.
Files used here: "data/PCA/adna_pca_newiid_2convertV3_100.evec", "data/erase_metadata_interpolated.csv"
Output: "./corr_convolution.csv" this is a file where we take "data/erase_metadata_interpolated.csv" and add 6 additional columns: one for each ancestry containing a value between 0 and 1 denoting how much that individual "belongs" to the specific ancestry and one last one containing a number from 1 to 5 which says which ancestry we've assigned the individual to (1 - WHG, 2-EHG, 3-Anatolia_N, 4-CHG, 5-North_Africa). 

4. geo_comp: This is the file where we use spatial statistics.  The crux of our project is understanding how main ancestry evolved in the different macroareas over time. The challenge was that during the deconvolution (previous R script) we generated compositional data, by adding the constraint that the values had to sum to 1 we lost a degree of freedom and as such are no longer in a Euclidean space and couldn't use any of the statistical methods we saw in class as in class we always oeprated under the assumption of being in a Euclidean space.
As a result, we performed an ILR transformation of our compositional data and that brought us back into a Euclidean space. We then performed another PCA to find a 4-dimensional orthonormal basis. We then performed our geospatial analysis and finally used the inverse ILR transformation to go back to our starting space where we could then interpret the results.
Files used here: "./corr_convolution.csv"
Output: "data/comp_time_space.csv" this is a file that contains the result of the spatial statistics analysis for each Macro_Area and Macro_Period which is calculated by taking the mean for each individual that belongs to that Macro_Area and Macro_Period (NOT A GOOD EXPLANATION).

5. many_plots: This is the file where we made most of the plots used to represent our data on the poster
Files used here: "./corr_convolution.csv"

6. world_plots: This is the file we used to make the plots where we plotted our observations on top of maps from google maps
Files used here: "./corr_convolution.csv"

Conclusion:
We came to 2 main conslusions: 
1) ERASE or how individuals composition flattened over time: we noticed that as time increased the genetic composition of the individuals followed a general trend of flattening - with the predominance of different ancestries gradually decreasing. This was particularly prominent in areas with strong migration pattens eg Western Europe and much less prominent in isolated areas eg Sicily or Sardinia. We believe that this is due to an elevated mixture among individuals from different ancestries that crossed paths at different points in history.
2) MIGRATION - CHG takes over Eurasia: We chose to focus specifically on the CHG (Caucasus Hunter Gatherer) ancestry, which was predominant in the region sorrounding the Caucasus. In the Neolithic, this ancestry was not very present in Eurasia and in fact was the perdominant ancestry in just one macroarea (Eastern Europe, where the Caucasus region is situated geographically). There was a very significant shift in the CHG composition during the Bronze Age, suggesting a high level of mixture between the CHG poplation and the other populations present in Europe. By the Classical Age, CHG was the prominent ancestry in six of our macroareas (out of twelve), showing a very significant expansion. Notably, even in the islands, which we saw were the most resistant to change in composition, CHG increased from 10-15% to 20-25% by the modern age.

Final notes: 
- in the Plots folder, the plots from each script (except for scripts 5 and 6 which were exclusively made to generate plots and so it would have been counter-intuitive to save each plot) have been saved for easier access
- Two animations were made to showcase 1) which was the predominant ancestry in each Macro_Area in each era (https://www.youtube.com/watch?v=xGQ-XL5xd68) and 2) how CHG composition changed in each Macro_Area over the different eras (https://www.youtube.com/watch?v=VSL53vpeXYM)
- The link to the pdf of our poster is: https://drive.google.com/drive/folders/1iPBIb2cv1xdwvPw2e0haM4TqVpg7U6SS







Title:Ericoid mycorrhizal growth response is influenced by host plant phylogeny

Authors: Alice S. Neves, Laura G. van Galen, Martin Vohník, Martina Peter, Elena Martino, Thomas W. Crowther, Camille S. Delavaux

=================================================================================
1. Data processing.R

- In this R script the dry root biomass is calculated using a conversion ratio obtained with the dry root biomass data available from the samples that we were able to weight using the scale (wet weight : dry weight ratio).

- The total dry biomass is calculated as well as the mycorrhizal growth response (MGR) for above, belowground and total biomass. The shoot : root dry weight ratio and the difference between shoot : root dry weight ratio between treated and control seedlings (here as MGR ratio)  was also calcuated. These values are added to the master dataframe that is exported as a csv table and will be further used in the next R scripts.

---------------------------------------------------------------------------------
2. differences_estimates.R

- the difference between the real dry biomass values and the estimates originated using the wet weight : dry weight ratio are also plotted with the R2 values to check if the ratio is aproppriate to be used

---------------------------------------------------------------------------------
3. colonization_assessment.R

- In this R script the proportion and percentage of colonized roots is calculated for fungal coils, hyphae and conidiophores (spikes)

- Proportion of colonization is added to the master dataframe

- The efficiency of the inoculation is tested by analysing if the colonization of treated samples is significantly higher than the controls
	- we were only able to do so for hyphal coils
	- hyphae and conidiophores have raw data plotted

- It was tested if mycorrhizal growth response (MGR) can be predicted by proportion of colonization

---------------------------------------------------------------------------------
4. plant_models_MGR.R

- In this R script the MGR is modelled using linear mixed effects models
- marginal square means and pairwaise comparisons between plant-fungi combinations

- model for shoot : root ratio MGR

--------------------------------------------------------------------------------- 
5. response_specificity_plots.R

- generate plots of response specificy for fungi and plant regarding above, belowground and total MGR

---------------------------------------------------------------------------------
6. ericoid_phylogeny.R

- In this R script a phylogenetic tree for the plant species used in the experiment is obtained from the GBOTB
- The MGR is ploted in a heath map and the existence of phylogenetic signal for the MGR in the plant phylogeny is tested (using Blomberg's K and Pagel's lambda)

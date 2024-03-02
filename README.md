**Explanation of the scripts and spreadsheets of each analysis - with their respective variables**
---------------------------------------------------------------------------------------------------

**_Impact of the EFN traits on the potential indirect effects between plants via dominant and subordinate ant species:_**

**Script File: “JEcol - indirect_effects_plants_traits.R”**


This script describes our analyses to examine the influence of two extrafloral nectar (EFN) traits – the average number of active nodes and secretory area of EFNs – on potential indirect effects among Bignonieae plant species. 
Promoter degree quantifies the extent to which a particular plant species affects another, while receptor degree measures the degree to which one plant species is influenced by another.


**Data File: “dataset_muller_planttraits.txt”:**

Column A (Species) = Author-assigned plant species identification number.
Column B (Species_name) = Name of each plant species.
Column C (average_activenodes_perplot) = Average number of active nodes per plot for each species.
Column D (Secretoryarea_perplantnode_mm2) = Average size (in mm²) of secretory area per node for each species.
Column E (promoter_all_ants) = Promoter degree of each plant species when sharing the entire ant assemblage.
Column F (receptor_all_ants) = Receptor degree of each plant species when sharing the entire ant assemblage.
Column G (promoter_dominant) = Promoter degree of each plant species when sharing only dominant ant species.
Column H (receptor_dominant) = Receptor degree of each plant species when sharing only dominant ant species.
Column I (promoter_subordinate) = Promoter degree of each plant species when sharing only subordinate ant species.
Column J (receptor_subordinate) = Receptor degree of each plant species when sharing only subordinate ant species.




**_Impact of indirect effects on leaf herbivory patterns:_**


**Script File: “JEcol - indirect_effects_herbivory.R”**

This script describes our analyses to evaluate whether indirect effects of promoter plant species influence, positively or negatively, leaf herbivory patterns of receptor plant species.



**Data File: “herb_muller.txt”:**

Column A (Species) = Author-assigned plant species identification number.
Column B (Species_name) = Name of each plant species.
Column C (mean_herb_sp) = Average leaf herbivory per species.
Column D (total_leaflets_sampled_sp) = Total number of leaflets sampled per plant species.
Column E (mean_prop_leaflets_herb_sp) = Average proportion of leaflets damaged by herbivory in each plant species.
Column F (promoter_all_ants) = Promoter degree of each plant species when sharing the entire ant assemblage.
Column G (receptor_all_ants) = Receptor degree of each plant species when sharing the entire ant assemblage.
Column H (promoter_dominant) = Promoter degree of each plant species when sharing only dominant ant species.
Column I (receptor_dominant) = Receptor degree of each plant species when sharing only dominant ant species.
Column J (promoter_subordinate) = Promoter degree of each plant species when sharing only subordinate ant species.
Column K (receptor_subordinate) = Receptor degree of each plant species when sharing only subordinate ant species.





**_Local impact of the promoter plant species on the ant visitation pattern and plant herbivory on neighbour and non-neighbour plants:_**

**Script File: “JEcol - Indirect_effects_local_impact_ants_herbivory.R”**

This script explores competition and facilitation hypotheses at a local scale (within our plots), considering the neighborhood of the most promoter plant per plot.

**Data File: “plots_promoter_ants_herb.txt”**

Column A (plot) = Plot identification.
Column B (specie) = Name of each plant species.
Column C (id) = Author's control number for each plant.
Column D (plant) = Classification of the plant in the plot as: promoter of indirect effects (1promoter), neighbor (2neighbor - plant neighboring the most promoting plant within a diameter of 10 m per plot) and non-neighbor (3n_neighbor - plant outside the diameter of 10 m regarding to the promoter plant).
Column E (sum_dom) = The total sum of dominant ant workers that the plant interacted with per plot.
Column F (occ_dom) = Occurrence of dominant ants on the plant per plot (0 = non-occurrence; 1 = occurrence).
Column G (sum_sub) = The total sum of subordinate ant workers that the plant interacted with per plot.
Column H (occ_sub) = Occurrence of subordinate ants on the plant per plot (= non-occurrence; 1 = occurrence).
Column I (mean_herb) = Average number per plot of leaflets damaged by herbivore of each plant.

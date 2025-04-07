README
The following readme describes the structure of the R code used to explore how the the New Zealand Seafloor Community Classification (NZSCC) and associated biodiversity layers could be used in Systematic Conservation Planning. The work builds on code and analysis from:
1. Stephenson et al. (2021) Species composition and turnover models provide robust approximations of biodiversity in marine conservation planning. Ocean and Coastal Management.
2. Tablada, J.,Geange, S. & Stephenson, F. (in review) Evaluation of current spatial management areas using the New Zealand Seafloor Community Classification. Aquatic Conservation: Marine And Freshwater Ecosystems
3. Tablada, J.,Geange, S. Hiddink, JG. & Stephenson, F. (draft). Enhancing representativity of protected seafloor communities in Aotearoa New Zealand

R CODE:
1. Generating_Inter_intra_group_similarity.R		---		
Description: Output spatial layers of within (intra) and between (inter) group similarity for the NZSCC.
2. 2. Using the NZSCC for MSP.R		---		
Description: Calculate the representativity of Spatial Management Areas as assessed by NZSCC group extent, within and between group proportion and species richness

DATA:
DF.source - data for demersal fish (for years 2070 - 2022)
Pred_1km.CMB.source Ð environmental data used at a 1km grid resolution saved as a dataframe

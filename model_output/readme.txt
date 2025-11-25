# README FILE FOR LANDLAB MODEL OUTPUTS
#
# PURPOSE: Contains model outputs reported in Rossi et al. (in review)
#
# FILE NAMING STRUCTURE:
# 'k1*.***L***bh*Ui*.**' prefix for main model (slope dependent weathering)
# 'noslopek1*.***' prefix for alternative model (no slope dependent weathering)
# '_initial.csv', '_final.csv', or '_series.csv' suffix
#
# DESCRIPTION:
# These are model outputs used create Figures 5, S9-10, S12-14 and Table S3. The prefix for
# the main model contains the value for k1 [-] in the Strudley et al (2008) soil production 
# function, the total hillslope length [m], the initial bump heights [m], and the initial
# uplift rates as a fraction of the new uplft rate [-]. Other parameters for soil creep 
# (d1=0.5 m; Sc=0.75) and soil production (d2=0.5 m; d3=0.1 m; w0={192,240,288,240} m/Ma for 
# left-, top-, right-, and bottom-facing slopes) are held constant. The new uplift rate is
# 75 m/Ma. Grid spacing is 2 m, and the weathering time step is 1 year. The 'noslope' scenarios
# use a 400-m swath, bump heights of 3 m, and 0.25 initial uplift rates. 
#
# REFERENCES:
# Rossi, M.W., Tucker, G.E., Anderson, S.P., Anderson, R.S., and McGlinchy, J., in review, 
# Coevolving topography, patchy soils, and forest structure. 
#
# Strudley, M. W., Murray, A. B., & Haff, P. K. (2006). Emergence of pediments, tors, and 
# piedmont junctions from a bedrock weathering–regolith thickness feedback. Geology, 34(10), 
# 805-808.

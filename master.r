rm(list = ls())
MS.PATH = paste(getwd(),"/Data and Code/",sep="")

# Re-produce Figure 1
source(paste(MS.PATH, 'code/Fig1.r', sep=''))

# Re-produce Figure 2
source(paste(MS.PATH, 'code/Fig2.r', sep=''))


# ------------------------------------------------------------------
# Bioassay Data Analysis
# The following scripts analyze the Bioassay data to reproduce:
# 1. Table 1
# 2. Figure 3
# 3. Supplementary Table D.1
# ------------------------------------------------------------------
source(paste(MS.PATH, 'code/Table1.r', sep=''))  # Table 1
source(paste(MS.PATH, 'code/Fig3.r', sep=''))    # Figure 3
source(paste(MS.PATH, 'code/TableD1.r', sep='')) # Supplementary Table D.1



# ------------------------------------------------------------------
# La Paloma Stream Data Analysis
# The following scripts analyze the La Paloma Stream data to reproduce:
# 1. Table 2
# 2. Figure 4
# 3. Supplementary Table D.2
# ------------------------------------------------------------------
source(paste(MS.PATH, 'code/Table2.r', sep=''))  # Table 2
source(paste(MS.PATH, 'code/Fig4.r', sep=''))    # Figure 4
source(paste(MS.PATH, 'code/TableD2.r', sep='')) # Supplementary Table D.2



# ------------------------------------------------------------------
# Lake Michigan Water Chemistry Data Analysis
# The following scripts analyze the Lake Michigan Water Chemistry Data to reproduce:
# 1. Table 3
# 2. Supplementary Table D.3
# 2. Figure 5
# ------------------------------------------------------------------
source(paste(MS.PATH, 'code/Table3.r', sep='')) # Table 3
source(paste(MS.PATH, 'code/TableD3.r', sep='')) # Supplementary Table D.3
source(paste(MS.PATH, 'code/Fig5.r', sep=''))   # Figure 5



# Re-produce Table 4
source(paste(MS.PATH, 'code/Table4.r', sep=''))

# Re-produce Table 5
source(paste(MS.PATH, 'code/Table5.r', sep=''))

# Re-produce Supplementary Table E.1
source(paste(MS.PATH, 'code/TableE1.r', sep=''))

# Re-produce Supplementary Table E.2
source(paste(MS.PATH, 'code/TableE2.r', sep=''))



# Re-produce Figure 6 
# Note: This script generates two separate output files ('fig6_nu4' and 'fig6_nu10'),
# which are combined into the single Figure 6.
source(paste(MS.PATH, 'code/Fig6.r', sep=''))

# Re-produce Figure 7
# Note: This script generates two separate output files ('fig7_1' and 'fig7_2'),
# which are combined into the single Figure 7.
source(paste(MS.PATH, 'code/Fig7.r', sep=''))

# Re-produce Figure 8
# Note: This script generates two separate output files ('fig8_1' and 'fig8_2'),
# which are combined into the single Figure 8 in the supplementary manuscript.
source(paste(MS.PATH, 'code/Fig8.r', sep=''))




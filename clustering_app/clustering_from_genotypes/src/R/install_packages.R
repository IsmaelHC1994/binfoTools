requiredPackages = c('argparse', 'ggplot2', 'dplyr', 'factoextra', 'NbClust')
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p, Ncpus=4)
}

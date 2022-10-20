# Transmission-of-trained-immunity
Re-sequenced data for the reviewers from the paper "Transmission of trained immunity and heterologous resistance to infections across generations."

The script "1. Analysis of merged data using original script" re-capitulates the original anlysis based on the merged data (original data + additionally sequenced data).

The script "2.1 DESeq Analysis cMOPs based on merged data" was used for the detailed analysis of the merged dataset, testing different pre-filtering cutoffs on the subset of cMOP samples only. 

The script "2.2 DESeq Analysis Ly6Ch based on merged data" was used for the detailed analysis of the merged dataset, testing different pre-filtering cutoffs on the subset of Ly6Ch monocyte samples only. 

All analyses were performed using the docker image "jsschrepping/r_docker:jss_R412_S41".

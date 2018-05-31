## Prepare deployment data

library(tidyverse)

# read in deployment csv
d <- read_csv("lnfs_deployment_data.csv")
d <- d %>% dplyr::select(animal_id_1, recovery_date_win, deploy_length_cm, deploy_girth_cm, deploy_mass_kg) %>% filter(!is.na(recovery_date_win))
colnames(d) <- c('id', 'recovery_date', 'length', 'girth', 'mass')

save(d, file = 'lnfs_recovered_deployments_2016-17.Rdata')

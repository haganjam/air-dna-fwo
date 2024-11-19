
# download belgian breeding bird survey

# https://www.gbif.org/dataset/81c5a091-6e94-40db-a2a4-48f4de42d410
# 81c5a091-6e94-40db-a2a4-48f4de42d410

# load rgbif library
library(rgbif)

# extract the florabank1 data
occ_download(
  pred_default(), 
  pred("datasetKey", "81c5a091-6e94-40db-a2a4-48f4de42d410"), 
  format = "SIMPLE_CSV"
)

# check the status
occ_download_wait('0021646-241107131044228')

# download the data locally
occ_download_get('0021646-241107131044228',
                 path = "02-wp-1-2/raw-data/breeding-bird-survey/") 

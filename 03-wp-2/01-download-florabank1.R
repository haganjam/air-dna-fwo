
# download florabank1 from gbif

# https://doi.org/10.3897/phytokeys.12.2849
# 271c444f-f8d8-4986-b748-e7367755c0c1

# load rgbif library
library(rgbif)

# extract the florabank1 data
occ_download(
  pred_default(), 
  pred("datasetKey", "271c444f-f8d8-4986-b748-e7367755c0c1"), 
  format = "SIMPLE_CSV"
)

# check the status
occ_download_wait('0014629-241107131044228')

# download the data locally
occ_download_get('0014629-241107131044228',
                 path = "02-wp-1-2/raw-data/florabank1/")

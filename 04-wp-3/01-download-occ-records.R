
# download occurrence records

# set-up gbif credentials
usethis::edit_r_environ()

# load the rgbif package
library(rgbif)

# use occ_search() for the Eurasian lynx
lynx_name <- name_backbone("Lynx lynx")
head(lynx_name)

# extract the lynx data
occ_download(
  pred_default(), 
  pred("taxonKey", lynx_name$usageKey), 
  format = "SIMPLE_CSV"
  )

# check status
occ_download_wait('0014538-241107131044228')

# save the data locally
occ_download_get('0014538-241107131044228',
                 path = "03-wp-3/raw-data/occurrence-data")

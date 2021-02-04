This is a companion piece to **GBIF filtering checklist** the blog post here[](). 

Here I will take you through some additional filters that you might want to add. The rest of these filters are a little bit more difficult and might involve more judgment calls. So I put them in this repository. 

## Filter **metagenomics**

**Metagenomics** (see previous [blog post](https://data-blog.gbif.org/post/gbif-molecular-data-quality/)) is a new publishing area for GBIF. Without going into too many details, metagenomics samples the environment for DNA and then matches them against an existing reference database. Especially with non-microorganisms these matches can often be incorrect or suspicious. 

Currently, there is not a great way for filtering for only metagenomics datasets. There are some [dicussions](https://github.com/gbif/registry/issues/247)
 about additional dataset categories but these have not be implemented yet. 

As a researcher you might want to check records with **material sample** or those published by [MGnify](https://www.gbif.org/publisher/ab733144-7043-4e88-bd4f-fca7bf858880) 

```R
# here is a script that will remove *most* metagenomics records 
gbif_download %>%
filter(!basisOfRecord == "MATERIAL_SAMPLE") %>%
filter(!publishingOrgKey == "ab733144-7043-4e88-bd4f-fca7bf858880")
```

## Filter **automated-ids** 

Some datasets use images to automatically identify 


## Filter not in **IUCN range**

**CoordinateCleaner** also includes a function for filtering based on expert distribution polygons. 

Here you must [download](https://www.iucnredlist.org/resources/spatial-data-download) your polygons first. Not all taxa will have a reliable range polygon.  

In my experience, unless your group is well is well-studied (Mammals, Birds, Reptiles see list [here](https://www.iucnredlist.org/resources/spatial-data-download)) it might be hard to get complete enough coverage for this to be worthwhile, but the coverage is always improving. 

Here is an example of filtering using the range from the damselfly **Calopteryx xanthostoma**. Here you can [download](/post/2020-12-08-typical-user-gbif-data-cleaning_files/Calopteryx xanthostoma/) the range shapefile for this species for the example below. The [gbif_download](/post/2020-12-08-typical-user-gbif-data-cleaning_files/Calopteryx xanthostoma.csv) used for this example. 


```R 
library(dplyr)
library(CoordinateCleaner) 

# read in directory where you downloaded the shapefile
range_shp = sf::st_read("Calopteryx xanthostoma") 

gbif_download = readr::read_tsv("Calopteryx xanthostoma.csv") %>%
filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # need to remove missing coordiantes

spdf = range_shp %>% as("Spatial")
class(spdf) # should be SpatialPolygonsDataFrame

gbif_download %>%
cc_iucn(
range = spdf,
lon = "decimalLongitude",
lat = "decimalLatitude",
species = "species",
buffer = 5, # buffer in decimal degrees
value = "clean"
)
```

It is usually good to add a buffer to the polygons, to catch any occurrences that might be just outside of the range. I added a 5 decimal degree buffer in this example. 

Also beware some IUCN polygons are **very large** and might run slow or crash your session. 

## Filter **gridded datasets**

Rasterized or gridded datasets are common on GBIF. These are datasets where location information is pinned to a low-resolution grid. 

GBIF has an experimental tag for datasests which exhibit a certain about of "griddyness". You can read more [here](). To remove gridded datasets that might be above the acceptable resolution. 

Most publishers actually fill in  

## Filter **spatial outliers** 

Sometimes a range polygon is not avaiable or you want a more general purpose way to flag suspicous points. 

I have found **DBSCAN** to be an effective way to remove **spatial outliers** in patchy GBIF data. You can read more about DBSCAN from a [previous post](https://data-blog.gbif.org/post/outlier-detection-using-dbscan/). 

The following could be put at the end of a pipeline. Note that you would need to split by species, if you had multiple species. 

```R 
library(dplyr)
library(dbscan)

gbif_download = readr::read_tsv("Calopteryx xanthostoma.csv") %>%
filter(!is.na(decimalLongitude) | !is.na(decimalLatitude))

gbif_download %>% 
mutate(cluster = 
dbscan::dbscan(
as.matrix(.[,c("decimalLatitude","decimalLongitude")]), 
eps = 15, minPts = 3)$cluster %>%
as.factor() 
) %>%
mutate(dbscan_outlier = ifelse(cluster == 0,TRUE,FALSE)) %>%
filter(!dbscan_outlier)
```

## Filter **environmental outliers**

Removing **environmental outliers** can also be important for certain applications. This can be done with using **reverse jackknifing**. In this case you must also download the environmental data you are interested in using. I will be using [bioclim](https://www.worldclim.org/data/bioclim.html) and the R package [biogeo](https://cran.r-project.org/web/packages/biogeo/biogeo.pdf) for the `biogeo::rjack()` function.  

```R 
# The following could be put at the end of a pipeline. 
# Note that you would need to split by species, if you had multiple species.

library(sp)
library(raster)
library(dplyr)
library(purrr)

path = "" # where you want to save the raster data
r = raster::getData('worldclim', var='bio',res=10,path=path) 

gbif_download = readr::read_tsv("Calopteryx xanthostoma.csv") %>%
filter(!is.na(decimalLongitude) | !is.na(decimalLatitude)) 

bioclim_data = sf::st_as_sf(
gbif_download,
coords = c("decimalLongitude", "decimalLatitude"),
crs=4326
) %>%
raster::extract(r,.) %>% # extract values from points
as.data.frame() 

gbif_bioclim = cbind(gbif_download,bioclim_data) %>% 
mutate(row_number = row_number()) # use to match index

gbif_download_with_outliers = bioclim_data %>% 
na.omit() %>% # remove missing climate values
map(~ biogeo::rjack(.x)) %>% # run the reverse jackknife outlier search
compact() %>% # remove columns with no outliers 
map(~gbif_bioclim$row_number %in% .x) %>% 
bind_rows() %>% 
setNames(paste0("rjack_outlier_", names(.))) %>%
mutate(number_of_outliers = rowSums(.)) %>% 
cbind(gbif_bioclim,.) %>%
glimpse() 

# should find two environmental outliers 
gbif_download_with_outliers %>% 
filter(number_of_outliers > 5)
```

You will want to do any **outlier detection at the end** of any cleaning pipeline, since noise can reduce the effectiveness of outlier detection.


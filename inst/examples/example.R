# the purpose of this script is to demonstrate how to
# install and and use the `cloud2trees` package from GitHub
# https://github.com/georgewoolsey/cloud2trees

# this script was tested on a Windows machine using
# a fresh install of R version 4.4.2
# that is, no packages had been installed or used in R prior

# install pkgbuild
install.packages("pkgbuild")

# check for Rtools which is required to build packages
pkgbuild::check_build_tools(debug = TRUE)
### ... there will be a message:
### ... "building r package from source requires installation of additional build tools"
### ... select "yes"
### ... an error will be issued "Error: Could not find tools necessary to compile a package"
### ... but the build tools should begin installing (minimize R window to see if its back there)
### ... after the download completes, try again...
# check for Rtools which is required to build packages
pkgbuild::check_build_tools(debug = TRUE)
### ... response should be: "Your system is ready to build packages!"

# install remotes package
install.packages("remotes")

# install tidyverse package...not required but tidyverse is {insert fire emoji}
install.packages("tidyverse")

# install lasR for point cloud processing
install.packages("lasR", repos = 'https://r-lidar.r-universe.dev')

# install github package from "tiagodc/TreeLS"
remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)

# install github package from "olgaviedma/LadderFuelsR"
remotes::install_github(repo = "olgaviedma/LadderFuelsR", upgrade = F)

# install github package from "DRAAlmeida/leafR"
remotes::install_github(repo = "DRAAlmeida/leafR", upgrade = F)

# install github package from "georgewoolsey/cloud2trees"
remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)

#############################################
## run some tests
#############################################
library(tidyverse)
library(cloud2trees)
# path to las data
# a test las file but this could also be a directory path with >1 .las|.laz files
## ... notice, we didn't directly install "lidR" above
## ... it was installed with the other packages as a dependency
i <- system.file(package="lidR", "extdata", "MixedConifer.laz")
# run it
cloud2trees_ans <- cloud2trees::cloud2trees(
  output_dir = "c:/Users/gwoolsey/Downloads/"
  , input_las_dir = i
)
# did it do it?
cloud2trees_ans %>% names()
# is there a DTM?
## ... notice, we didn't directly install "terra"
cloud2trees_ans$dtm_rast %>%
  terra::plot()
# is there a CHM?
cloud2trees_ans$chm_rast %>%
  terra::plot()
# are there trees?
cloud2trees_ans$treetops_sf %>% dplyr::glimpse()
# would a CHM with tree crowns overlaid look neat?
cloud2trees_ans$chm_rast %>%
  terra::as.data.frame(xy=T) %>%
  dplyr::rename(f=3) %>%
  ggplot2::ggplot() +
  ggplot2::geom_tile(
    mapping = ggplot2::aes(x=x,y=y,fill=f)
  ) +
  ggplot2::geom_sf(
    data = cloud2trees_ans$crowns_sf
    , color = "gray33"
    , fill = NA
    , lwd = 1
  ) +
  ggplot2::scale_fill_viridis_c(option = "plasma") +
  ggplot2::labs(fill="CHM ht. (m)") +
  ggplot2::theme_void()
# neat-o!

# what if we try to estimate tree dbh?
# add_dbh_crowns_sf <- cloud2trees::trees_dbh(tree_list = cloud2trees_ans$crowns_sf)

# ### ... ERROR:
# # Error in cloud2trees::trees_dbh(tree_list = cloud2trees_ans$crowns_sf) :
# #   Treemap data has not been downloaded to package contents. Use `get_treemap()` first.
# # If you supplied a value to the `input_treemap_dir` parameter check that directory for data.

# the message tells us to use get_treemap() first
cloud2trees::get_treemap()
# this took ~3 mins for me with a fast internet connection (900mbps+)

# try to estimate tree dbh now!
add_dbh_crowns_sf <- cloud2trees::trees_dbh(tree_list = cloud2trees_ans$crowns_sf)

# check for dbh
add_dbh_crowns_sf %>% dplyr::glimpse()

# plot it?
add_dbh_crowns_sf %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = tree_height_m, y = dbh_cm)) +
  ggplot2::geom_point(color = "navy", alpha = 0.6) +
  ggplot2::labs(x = "tree ht. (m)", y = "tree DBH (cm)") +
  ggplot2::scale_x_continuous(limits = c(0,NA)) +
  ggplot2::scale_y_continuous(limits = c(0,NA)) +
  ggplot2::theme_light()

# works.

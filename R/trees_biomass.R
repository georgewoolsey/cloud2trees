#' @title Estimate tree biomass for a tree list
#'
#' @description
#' `trees_biomass()` uses the input tree list (e.g. as exported by [raster2trees()]) with the columns
#' `treeID`, `tree_x`, `tree_y` to attach species information using USDA Forest Inventory and Analysis (FIA) codes.
#' If a spatial data frame of points is the input tree list, then the columns `tree_x`, `tree_y` are not required.
#'
#' FIA Forest Type Group Code is attached to each tree in the tree list based on the spatial overlap with the Forest Type Groups of the Continental United States dataset [Wilson 2023](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d).
#'
#' The simplified process for attaching forest type group to a tree is:
#'
#' * Forest type group 30-m raster (Wilson 2023) was aggregated to 90-m to make the data more accessible over the entire continental US
#' * Nearest neighbor imputation is used to fill forest type data if a tree falls inside a no-forest cell in the original data
#' * The FIA forest type group is applied to a tree based on spatial overlap
#'
#' @param tree_list data.frame. A data frame with the columns `treeID`, `tree_x`, `tree_y`, and `tree_height_m`.
#' If an `sf` class object with POINT geometry (see [sf::st_geometry_type()]), the program will use the data "as-is" and only require the `treeID` column.
#' @param crs string. A crs string as returned from [sf::st_crs()] or the EPSG code of the x,y coordinates.
#' Defaults to the crs of the `tree_list` data if of class "sf".
#' @param study_boundary sf. The boundary of the study are to define the area of the regional model.
#' If no boundary given, regional model will be built from location of trees in the tree list.
#' @param input_foresttype_dir directory where Forest Type Groups data exists. Use [get_foresttype()] first.
#'
#' @references
#' * [Forest Type Groups of the Continental United States](https://www.arcgis.com/home/item.html?id=10760c83b9e44923bd3c18efdaa7319d)
#' Wilson, B.T. (2023). Forest Type Groups of the Continental United States.
#'
#' @return Returns a list of objects: tree_list = spatial data frame of individual trees; foresttype_rast = raster of forest types in the area.
#'
#' @examples
#'  \dontrun{
#'  # example tree list
#'  tl <- dplyr::tibble(
#'      treeID = c(1:21)
#'      , tree_x = rnorm(n=21, mean = 458064, sd = 11)
#'      , tree_y = rnorm(n=21, mean = 4450074, sd = 11)
#'    )
#'  # call the function
#'  tl_type <- trees_biomass(tree_list = tl, crs = "32613")
#'  # what?
#'  tl_type %>% class()
#'  # a list, but what is in it?
#'  tl_type %>% names()
#'  # plot the tree_list spatial points
#'  tl_type$tree_list %>% ggplot2::ggplot() + ggplot2::geom_sf(ggplot2::aes(color=forest_type_group))
#'  # plot the foresttype_rast raster
#'  tl_type$foresttype_rast %>% terra::plot()
#'  }
#' @export
#'
trees_biomass <- function(
  tree_list
  , crs = NA
  , study_boundary = NA
  , input_foresttype_dir = NULL
){
  ####################################################################
  # check external data
  ####################################################################
    # find external data
    find_ext_data_ans <- find_ext_data(
      input_foresttype_dir = input_foresttype_dir
    )
    # if can't find external foresttype data
    if(is.null(find_ext_data_ans$foresttype_dir)){
      stop(paste0(
        "Forest Type Group data has not been downloaded to package contents. Use `get_foresttype()` first."
        , "\nIf you supplied a value to the `input_foresttype_dir` parameter check that directory for data."
      ))
    }
  ##################################
  # convert to spatial points data
  ##################################
  tree_tops <- check_spatial_points(tree_list, crs)

  # get rid of columns we'll create
    tree_tops <- tree_tops %>%
      # throw in hey_xxxxxxxxxx to test it works if we include non-existent columns
      dplyr::select( -dplyr::any_of(c(
        "hey_xxxxxxxxxx"
        , "cell"
        , "crown_dia_m"
        , "crown_length_m"
        , "crown_volume_m3"
        , "tree_kg_per_m3"
        , "cruz_biomass_kg"
      )))

  ##################################
  # load foresttype data. see get_foresttype()
  ##################################
    # get the foresttype data
    foresttype <- terra::rast(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype.tif")
      )
    # get the lookup
    foresttype_lookup <- readr::read_csv(
        file.path(find_ext_data_ans$foresttype_dir, "foresttype_lookup.csv")
        , progress = F
        , show_col_types = F
      ) %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything()
        , as.character
      ))

  ##################################
  # create extent if empty
  ##################################
    # set study boundary to tree extent if missing
    if(
      !(
        c(
          inherits(study_boundary, "sf")
          , inherits(study_boundary, "sfc")
          , inherits(study_boundary, "SpatVector")
        ) %>% any()
      )
    ){
      study_boundary <- tree_tops %>%
        sf::st_bbox() %>%
        sf::st_as_sfc()
    }

  ################################################################
  # calc_rast_cell_trees
    # function to aggregate tree list to the raster cell level
    # and join to the raster cell overlap data generated via
    # calc_rast_cell_overlap() within the function
  ################################################################
    calc_rast_cell_trees_ans <- calc_rast_cell_trees(
      rast = foresttype
      , tree_list = tree_tops
      , poly_extent = study_boundary
      , calc_tree_level_cols = T
    )

  ################################################################
  # distribute_stand_fuel_load
    # use our `get_cruz_stand_kg_per_m3()` function to calculate
    # the stand level CBH in kilograms per cubed meter
    # and distribute this across the tree list
  ################################################################
    distribute_stand_fuel_load_ans <- distribute_stand_fuel_load(
      cell_df = calc_rast_cell_trees_ans$cell_df
      , tree_list = calc_rast_cell_trees_ans$tree_list
    )

  # return
  return(
    distribute_stand_fuel_load_ans$tree_list
  )
}
#######################################################
# intermediate function 1
#######################################################
  # [Cruz et al. (2003)](https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5)
  # developed models to predict canopy fuel stratum at the stand level for
  # four coniferous forest types common in the western US:
  #   Douglas-fir, ponderosa pine, lodgepole pine, and mixed conifer
  #   map these to the FIA forest type groups and apply the model estimates
  get_cruz_stand_kg_per_m3 <- function(forest_type_group_code, basal_area_m2_per_ha, trees_per_ha){
    forest_type_group_code <- dplyr::coalesce(as.numeric(forest_type_group_code), 0)
    # Cruz et al. (2003)
    # https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5
    # Page 46, Table 4
    if(!is.na(forest_type_group_code) && forest_type_group_code == 200){
      #Douglas-Fir Group
      b0 = -7.380
      b1 = 0.479
      b2 = 0.625
    }else if(!is.na(forest_type_group_code) && forest_type_group_code == 220){
      #Ponderosa Pine Group
      b0 = -6.649
      b1 = 0.435
      b2 = 0.579
    }else if(!is.na(forest_type_group_code) && forest_type_group_code == 280){
      #Lodgepole Pine Group
      b0 = -7.852
      b1 = 0.349
      b2 = 0.711
    }else if(!is.na(forest_type_group_code) && forest_type_group_code %in% c(120,260,320) ){
      #Mixed Conifer Group
      b0 = -8.445
      b1 = 0.319
      b2 = 0.859
    }else{
      #No Cruz et al. formulas for these YET!!!
      b0 = as.numeric(NA)
      b1 = as.numeric(NA)
      b2 = as.numeric(NA)
    }
    #Apply Cruz et al. if species found
    if(!is.na(b0)){
      return(
        exp(b0 + b1 * log(basal_area_m2_per_ha) + b2 * log(trees_per_ha))
      )
    }else{
      return(as.numeric(NA))
    }
  }
#######################################################
# intermediate function 2
#######################################################
  # function to ingest a raster and a polygon and calculate
  # the area of each raster cell that overlaps with the polygon
  calc_rast_cell_overlap <- function(rast, poly, buff = 100) {
    if(!inherits(rast, "SpatRaster")){
      stop("must pass a SpatRaster object to `rast`")
    }
    # convert to terra vector with same projection
    if(
      !inherits(poly, "SpatVector")
      && (inherits(poly, "sf") || inherits(poly, "sfc"))
    ){
      poly_vect <- poly %>%
        sf::st_make_valid() %>%
        sf::st_union() %>%
        terra::vect() %>%
        terra::project(terra::crs(rast))
    }else if(inherits(poly, "SpatVector")){
      poly_vect <- poly %>%
        terra::union() %>%
        terra::project(terra::crs(rast))
    }else{
      stop("must pass a spatial SpatVector or sf object to `poly`")
    }

    # crop the raster to our extent with a buffer and change the cell values to the total cell area
    r_crop <- rast %>%
      terra::crop(poly_vect %>% terra::buffer(width = buff)) %>%
      terra::cellSize(transform = F) # converts cell value to cell area

    # get a data frame of the cell numbers with the area of the raster cells that overlap with the poly extent
    overlap_df_temp <- terra::rasterize(
      x = poly_vect
      , y = r_crop
      , field = c(1)
      , cover = T
    ) %>%
    terra::as.data.frame(cells = T) %>%
    dplyr::rename(pct_overlap = layer)

    # create a data frame of the cropped raster and join our overlap data to create a "stand" data frame
    stand_df <-
      r_crop %>%
      terra::as.data.frame(xy = T, cells = T, na.rm = F) %>%
      # join on pct overlap with las_ctg
      dplyr::left_join(
        overlap_df_temp
        , by = "cell"
      ) %>%
      dplyr::mutate(
        overlap_area_m2 = ( area*dplyr::coalesce(pct_overlap, 0) )
        , overlap_area_ha = overlap_area_m2 / 10000
      )

    return(list(
      df = stand_df
      , rast = r_crop
    ))
  }
#######################################################
# intermediate function 3
#######################################################
  # calculate tree-level values we'll need for the back transformation
  # from stand-level metrics to tree level metrics such as crown volume
  # in kilograms per cubed meter so that the mean stand kg/m3 * tree m3 = tree kg
  calc_tree_level_cols <- function(df){
    # # proportion of crown length to have max crown width at
    # # 0.5 is a perfect ellipsoid, 0 is more conical, 1 looks like an icecream cone
    # # link to forest type
    # ht_to_max = 0.5

    ## check for cols
    nms <- names(df) %>% dplyr::coalesce("")
    has_cols <- c("crown_area_m2", "tree_height_m", "tree_cbh_m") %>%
      purrr::map(function(x){
        stringr::str_equal(tolower(nms), x) %>%
        max() # do any columns match, T=1
      }) %>%
      unlist() %>%
      min()
    if(has_cols==0){
      stop("the `df` data does not contain the columns `crown_area_m2`, `tree_height_m`, and `tree_cbh_m`, ensure columns exist")
    }

    ## basal area
    if(
      !(stringr::str_detect(nms, "basal_area_m2") %>% any())
      && (stringr::str_detect(nms, "dbh_cm") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          dbh_cm = as.numeric(dbh_cm)
          , basal_area_m2 = pi * (((dbh_cm/100)/2)^2)
        )
    }else if(
      !(stringr::str_detect(nms, "basal_area_m2") %>% any())
      && (stringr::str_detect(nms, "dbh_m") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          dbh_m = as.numeric(dbh_m)
          , basal_area_m2 = pi * ((dbh_m/2)^2)
        )
    }else if(
      (stringr::str_detect(nms, "basal_area_m2") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          basal_area_m2 = as.numeric(basal_area_m2)
        )
    }else{
      stop(paste0(
        "the `df` data does not contain the columns `basal_area_m2`, `dbh_cm`, or `dbh_m`"
        , "\n .... at least one of these columns must be present"
      ))
    }

    ## apply the crown calculations
    r_df <- df %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        crown_area_m2 = as.numeric(crown_area_m2)
        , tree_height_m = as.numeric(tree_height_m)
        , tree_cbh_m = as.numeric(tree_cbh_m)
      ) %>%
      # apply the calculations
      dplyr::mutate(
        #calculate crown diameter (m) = SQRT(area/pi) * 2
        crown_dia_m = sqrt(crown_area_m2 / pi) * 2
        #calculate crown length (m)
        , crown_length_m = tree_height_m - tree_cbh_m
        # calculate crown volume [4/3 * pi * a * b * c]
        # which is [4/3 * pi * (crownlength/2) * (max crown radius) * (max crown radius)]
        ## !!!!!!!!!!!!!!!!!!11 source??????????????????????????????????????????????????????????????
        , crown_volume_m3 = (4/3) * pi * ((crown_length_m/2)) * ((crown_dia_m/2)^2)
        # #calculate height to max crown
        # , height_to_max = (crown_length_m * ht_to_max) + tree_cbh_m
      )
    return(r_df)
  }
#######################################################
# intermediate function 4
#######################################################
  # function to aggregate tree list to the raster cell level
  # and join to the raster cell overlap data by calling
  # `calc_rast_cell_overlap()` within the function
  calc_rast_cell_trees <- function(rast, tree_list, poly_extent, buffer = 100, calc_tree_level_cols=T){
    if(!inherits(tree_list, "sf")){
      stop("must pass a spatial sf object to `tree_list`")
    }
    # check if not points
    if( min(sf::st_is(tree_list, type = c("POINT", "MULTIPOINT"))) == 0 ){
      stop(paste0(
        "data passed to `tree_list` is not point or multipoint data"
        , "\n see sf::st_geometry_type"
      ))
    }

    # calc_rast_cell_overlap
    overlap_ans <- calc_rast_cell_overlap(rast = rast, poly = poly_extent, buff = buffer)

    # attach cell to trees
    tree_list$cell <-
      # use terra::extract to get the cell id
      terra::extract(
        x = overlap_ans$rast
        , y = tree_list %>%
          terra::vect() %>%
          terra::project(terra::crs(overlap_ans$rast))
        , cells = T # cell numbers are also returned
      ) %>%
      dplyr::pull(cell)

    # check calc_tree_level_cols
    if(calc_tree_level_cols==T){
       tree_list <- tree_list %>%
         calc_tree_level_cols()
    }

    # aggregate tree list to cell
    nms <- names(tree_list) %>% dplyr::coalesce("")
    if(stringr::str_detect(nms, "forest_type_group_code") %>% any()){
      trees_agg <- tree_list %>%
        sf::st_drop_geometry() %>%
        # aggregate to stand level
        dplyr::group_by(cell, forest_type_group_code)
    }else{
      trees_agg <- tree_list %>%
        sf::st_drop_geometry() %>%
        # aggregate to stand level
        dplyr::group_by(cell)
    }

    # summarize
    if(calc_tree_level_cols==T){
      trees_agg <- trees_agg %>%
        dplyr::summarise(
          trees = dplyr::n()
          , basal_area_m2 = sum(basal_area_m2, na.rm = T)
          , mean_crown_length_m = mean(crown_length_m, na.rm = T)
          , sum_crown_volume_m3 = sum(crown_volume_m3, na.rm = T)
        ) %>%
        dplyr::ungroup()
    }else if(stringr::str_detect(nms, "basal_area_m2") %>% any()){
      trees_agg <- trees_agg %>%
        dplyr::summarise(
          trees = dplyr::n()
          , basal_area_m2 = sum(basal_area_m2, na.rm = T)
        ) %>%
        dplyr::ungroup()
    }else{
      trees_agg <- trees_agg %>%
        # aggregate to stand level
        dplyr::group_by(cell) %>%
        dplyr::summarise(
          trees = dplyr::n()
        ) %>%
        dplyr::mutate(basal_area_m2 = as.numeric(NA)) %>%
        dplyr::ungroup()
    }

    # join to raster area data
    r_df <- overlap_ans$df %>%
      dplyr::left_join(trees_agg, by = "cell") %>%
      dplyr::mutate(
        basal_area_m2_per_ha = basal_area_m2/overlap_area_ha
        , trees_per_ha = trees/overlap_area_ha
      )

    #return
    return(list(
      cell_df = r_df
      , tree_list = tree_list
      , rast = overlap_ans$rast
    ))
  }
#######################################################
# intermediate function 5
#######################################################
  # use our `get_cruz_stand_kg_per_m3()` function to calculate
  # the stand level CBH in kilograms per cubed meter
  distribute_stand_fuel_load <- function(cell_df, tree_list) {
    # calculate fuel loading at the stand levesl
    cell_df <- cell_df %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>% # this is key
      dplyr::mutate(
        kg_per_m3 = get_cruz_stand_kg_per_m3(
          forest_type_group_code = forest_type_group_code
          , basal_area_m2_per_ha = basal_area_m2_per_ha
          , trees_per_ha = trees_per_ha
        )
      ) %>%
      dplyr::ungroup() %>%
      # tertiary columns
      dplyr::mutate(
        # get cfl in kg/m2
        kg_per_m2 = mean_crown_length_m * kg_per_m3
        # get stand biomass in kg at the stand level
        , biomass_kg = kg_per_m2 * overlap_area_m2
        # single tree CBD in kg/m3 will be constant by stand/cell
        , tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3
      )

    # apply stand fuel load to trees
    tree_list <- tree_list %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
        cell_df %>%
          dplyr::select(cell,tree_kg_per_m3)
        , by = "cell"
      ) %>%
      dplyr::mutate(cruz_biomass_kg = tree_kg_per_m3*crown_volume_m3)

    #return
    return(list(
      cell_df = cell_df
      , tree_list = tree_list
    ))
  }

#' @title internal functions to estimate tree biomass
#'
#' @description
#' internal functions to extract raster values at point locations used by `trees_biomass*()` functions
#' for the most part, these functions are used to distribute raster cell (i.e. "stand") level estimates of fuel load
#' to individual trees
#'
#' @param forest_type_group_code numeric. as extracted by [trees_type()]
#' @param basal_area_m2_per_ha numeric.
#' @param trees_per_ha numeric.
#'
#' @keywords internal
#'
get_cruz_stand_kg_per_m3 <- function(forest_type_group_code, basal_area_m2_per_ha, trees_per_ha){
  #######################################################
  # intermediate function 1
  #######################################################
    # [Cruz et al. (2003)](https://scholar.google.com/scholar?cluster=316241498622221569&oi=gsb&hl=en&as_sdt=0,5)
    # developed models to predict canopy fuel stratum at the stand level for
    # four coniferous forest types common in the western US:
    #   Douglas-fir, ponderosa pine, lodgepole pine, and mixed conifer
    #   map these to the FIA forest type groups and apply the model estimates
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
    if(!inherits(poly, "SpatVector") &&
       ( inherits(poly, "sf") || inherits(poly, "sfc") )
    ){
      poly_vect <- poly %>%
        sf::st_union() %>%
        sf::st_transform(terra::crs(rast)) %>%
        terra::vect()
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
        , rast_epsg_code = terra::crs(rast, describe=T)$code[1]
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
  # this function will throw an error if the columns "crown_area_m2", "tree_height_m", "tree_cbh_m"
    # (and "dbh_cm" or "basal_area_m2") are not in the data or have all missing values
  calc_tree_level_cols <- function(df){
    # # proportion of crown length to have max crown width at
    # # 0.5 is a perfect ellipsoid, 0 is more conical, 1 looks like an icecream cone
    # # link to forest type
    # ht_to_max = 0.5

    ## check for cols
    nms <- names(df) %>% dplyr::coalesce("")
    check_df_cols_all_missing(
      df
      , col_names = c("crown_area_m2", "tree_height_m", "tree_cbh_m")
      , all_numeric = T
    )
    ## basal area
    if(
      !(stringr::str_equal(nms, "basal_area_m2") %>% any())
      && (stringr::str_equal(nms, "dbh_cm") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          dbh_cm = as.numeric(dbh_cm)
          , basal_area_m2 = pi * (((dbh_cm/100)/2)^2)
        )
    }else if(
      !(stringr::str_equal(nms, "basal_area_m2") %>% any())
      && (stringr::str_equal(nms, "dbh_m") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          dbh_m = as.numeric(dbh_m)
          , basal_area_m2 = pi * ((dbh_m/2)^2)
        )
    }else if(
      (stringr::str_equal(nms, "basal_area_m2") %>% any())
    ){
      df <- df %>%
        dplyr::mutate(
          basal_area_m2 = as.numeric(basal_area_m2)
        )
    }else{
      stop(paste0(
        "the data does not contain the columns `basal_area_m2`, `dbh_cm`, or `dbh_m`"
        , "\n .... at least one of these columns must be present"
      ))
    }

    check_df_cols_all_missing(
      df
      , col_names = c("basal_area_m2")
      , all_numeric = T
    )

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
        , crown_length_m = ifelse(tree_height_m<tree_cbh_m, tree_height_m*0.5, tree_height_m - tree_cbh_m)
        # calculate crown volume [4/3 * pi * a * b * c]
        # which is [4/3 * pi * (crownlength/2) * (max crown radius) * (max crown radius)]
        ## !!!!!!!!!!!!!!!!!!11 source??????????????????????????????????????????????????????????????
        , crown_volume_m3 = (4/3) * pi * ((crown_length_m/2)) * ((crown_dia_m/2)^2)
        # #calculate height to max crown
        # , height_to_max = (crown_length_m * ht_to_max) + tree_cbh_m
      )
    # return
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
          sf::st_transform(terra::crs(overlap_ans$rast)) %>%
          terra::vect()
        , cells = T # cell numbers are also returned
      ) %>%
      dplyr::pull(cell)

    # check calc_tree_level_cols
    if(calc_tree_level_cols==T){
       tree_list <- tree_list %>%
         calc_tree_level_cols()
    }

    # check col names
    nms <- names(tree_list) %>% dplyr::coalesce("")

    # summarize
    if(calc_tree_level_cols==T){
      trees_agg <- tree_list %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(
          trees = dplyr::n()
          , basal_area_m2 = sum(basal_area_m2, na.rm = T)
          , mean_crown_length_m = mean(crown_length_m, na.rm = T)
          , mean_crown_dia_m = mean(crown_dia_m, na.rm = T)
          , sum_crown_volume_m3 = sum(crown_volume_m3, na.rm = T)
        ) %>%
        dplyr::ungroup()
    }else if(stringr::str_equal(nms, "basal_area_m2") %>% any()){
      trees_agg <- tree_list %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(
          trees = dplyr::n()
          , basal_area_m2 = sum(basal_area_m2, na.rm = T)
        ) %>%
        dplyr::ungroup()
    }else{
      trees_agg <- tree_list %>%
        sf::st_drop_geometry() %>%
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
  # use our calculate stand-level CBD for cruz (lf already has CBD)
  # the stand level CBD in kilograms per cubed meter is distributed to trees
  distribute_stand_fuel_load <- function(cell_df, tree_list, cbd_method = "cruz", max_crown_kg_per_m3 = 2) {
    # check the method
    cbd_method <- dplyr::coalesce(cbd_method[1], "") %>% tolower() %>% stringr::str_squish()
    if( !(cbd_method %in% c("cruz", "landfire")) ){
      stop("`cbd_method` parameter must be one of \"cruz\" or \"landfire\"")
    }

    # for cruz, need to get the forest type group code
    # we'll do this based on the most common code in the tree data
    nms <- names(tree_list) %>% dplyr::coalesce("")
    if(
      !(stringr::str_equal(nms, "cell") %>% any())
    ){
      stop("tree_list data must have the column `cell`...should be in the output of calc_rast_cell_trees()?")
    }
    if(
      cbd_method == "cruz"
      && !(stringr::str_equal(nms, "forest_type_group_code") %>% any())
    ){
      stop(paste0(
        "cbd_method set to `cruz` but the `forest_type_group_code` does not exist in the tree list data"
        , "\n try running cloud2trees::trees_type() first"
      ))
    }
    if(
      cbd_method == "landfire"
      && !(stringr::str_equal(nms, "landfire_cell_kg_per_m3") %>% any())
    ){
      stop(paste0(
        "cbd_method set to `landfire` but the `landfire_cell_kg_per_m3` does not exist in the tree list data"
        , "\n try running cloud2trees::trees_landfire_cbd() first"
      ))
    }
    #############################################
    # CRUZ
    #############################################
    if(cbd_method == "cruz"){
      # select most common forest type group in a cell based on tree list
      # this should be extracted by trees_type()
      cell_ft_temp <- tree_list %>%
        sf::st_drop_geometry() %>%
        dplyr::mutate(
          is_na_ft = ifelse(
            is.na(as.numeric(forest_type_group_code))
            , 1
            , 0
          )
        ) %>%
        dplyr::group_by(cell, forest_type_group_code, is_na_ft) %>%
        dplyr::summarise(trees = dplyr::n()) %>%
        dplyr::group_by(cell) %>%
        dplyr::arrange(cell, is_na_ft, desc(trees)) %>% # sorts non-na first even if more trees in na
        dplyr::filter(dplyr::row_number()==1) %>%
        dplyr::ungroup()

      # calculate fuel loading at the stand level
      cell_df <- cell_df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(cell_ft_temp %>% dplyr::select(cell,forest_type_group_code), by = "cell") %>%
        dplyr::rowwise() %>% # this is key
        dplyr::mutate(
          cruz_stand_kg_per_m3 = get_cruz_stand_kg_per_m3(
            forest_type_group_code = forest_type_group_code
            , basal_area_m2_per_ha = basal_area_m2_per_ha
            , trees_per_ha = trees_per_ha
          )
        ) %>%
        dplyr::ungroup() %>%
        # tertiary columns
        dplyr::mutate(
          # get cfl in kg/m2
          kg_per_m2 = mean_crown_length_m * cruz_stand_kg_per_m3
          # get stand biomass in kg at the stand level
          , biomass_kg = kg_per_m2 * overlap_area_m2
          # single tree CBD in kg/m3 will be constant by stand/cell
          , cruz_tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3
        )

      # check max_crown_kg_per_m3
      max_crown_kg_per_m3 <-
        ifelse(is.null(max_crown_kg_per_m3),NA,max_crown_kg_per_m3) %>%
        as.numeric() %>%
        dplyr::coalesce(1e10) # set really high so no replacement if missing
      # do we need to do it?
      if(
        # are there records?
        !is.na(max(cell_df$cruz_tree_kg_per_m3, na.rm = T))
        # are there records over the max?
        && max(cell_df$cruz_tree_kg_per_m3, na.rm = T) > max_crown_kg_per_m3
        # are there records under the max to generate replacement?
        && ( cell_df %>%
          dplyr::filter(
            cruz_tree_kg_per_m3<max_crown_kg_per_m3
            & !is.na(cruz_tree_kg_per_m3)
          ) %>%
          nrow() ) > 0
      ){
        # replacement value is median of cells with < max
        new_tree_kg_per_m3 <- cell_df %>%
          dplyr::filter(cruz_tree_kg_per_m3<max_crown_kg_per_m3) %>%
          dplyr::pull(cruz_tree_kg_per_m3) %>%
          stats::median(na.rm = T)
        # replace value in data
        cell_df <- cell_df %>%
          dplyr::mutate(
            cruz_tree_kg_per_m3 = ifelse(
              !is.na(cruz_tree_kg_per_m3) & cruz_tree_kg_per_m3>max_crown_kg_per_m3
              , dplyr::coalesce(new_tree_kg_per_m3, cruz_tree_kg_per_m3)
              , cruz_tree_kg_per_m3
            )
          )
      }

      # apply stand fuel load to trees
      tree_list <- tree_list %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
          cell_df %>%
            dplyr::select(cell, cruz_tree_kg_per_m3, cruz_stand_kg_per_m3)
          , by = "cell"
        ) %>%
        dplyr::mutate(cruz_crown_biomass_kg = cruz_tree_kg_per_m3*crown_volume_m3) %>%
        dplyr::rename(cruz_stand_id = cell)
      # rename cell data
      cell_df <- cell_df %>% dplyr::rename(cruz_stand_id = cell)
    }

    #############################################
    # LANDFIRE
    #############################################
    if(cbd_method == "landfire"){
      # select median landfire_cell_kg_per_m3 in a cell based on tree list
      # this should be extracted by trees_landfire_cbd()
      cell_ft_temp <- tree_list %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(
          landfire_stand_kg_per_m3 = median(as.numeric(landfire_cell_kg_per_m3), na.rm=T)
        ) %>%
        dplyr::ungroup()

      # calculate fuel loading at the stand level
      cell_df <- cell_df %>%
        dplyr::ungroup() %>%
        dplyr::left_join(cell_ft_temp %>% dplyr::select(cell,landfire_stand_kg_per_m3), by = "cell") %>%
        # tertiary columns
        dplyr::mutate(
          # get cfl in kg/m2
          kg_per_m2 = mean_crown_length_m * landfire_stand_kg_per_m3
          # get stand biomass in kg at the stand level
          , biomass_kg = kg_per_m2 * overlap_area_m2
          # single tree CBD in kg/m3 will be constant by stand/cell
          , landfire_tree_kg_per_m3 = biomass_kg / sum_crown_volume_m3
        )

      # check max_crown_kg_per_m3
      max_crown_kg_per_m3 <- as.numeric(max_crown_kg_per_m3) %>%
        dplyr::coalesce(1e10) # set really high so no replacement if missing
      # do we need to do it?
      if(
        # are there records?
        !is.na(max(cell_df$landfire_tree_kg_per_m3, na.rm = T))
        # are there records over the max?
        && max(cell_df$landfire_tree_kg_per_m3, na.rm = T) > max_crown_kg_per_m3
        # are there records under the max to generate replacement?
        && ( cell_df %>%
          dplyr::filter(
            landfire_tree_kg_per_m3<max_crown_kg_per_m3
            & !is.na(landfire_tree_kg_per_m3)
          ) %>%
          nrow() ) > 0
      ){
        # replacement value is median of cells with < max
        new_tree_kg_per_m3 <- cell_df %>%
          dplyr::filter(landfire_tree_kg_per_m3<max_crown_kg_per_m3) %>%
          dplyr::pull(landfire_tree_kg_per_m3) %>%
          stats::median(na.rm = T)
        # replace value in data
        cell_df <- cell_df %>%
          dplyr::mutate(
            landfire_tree_kg_per_m3 = ifelse(
              !is.na(landfire_tree_kg_per_m3) & landfire_tree_kg_per_m3>max_crown_kg_per_m3
              , dplyr::coalesce(new_tree_kg_per_m3, landfire_tree_kg_per_m3)
              , landfire_tree_kg_per_m3
            )
          )
      }

      # apply stand fuel load to trees
      tree_list <- tree_list %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
          cell_df %>%
            dplyr::select(cell, landfire_tree_kg_per_m3, landfire_stand_kg_per_m3)
          , by = "cell"
        ) %>%
        dplyr::mutate(landfire_crown_biomass_kg = landfire_tree_kg_per_m3*crown_volume_m3) %>%
        dplyr::rename(landfire_stand_id = cell) %>%
        # drop landfire_cell_kg_per_m3 so we don't get confused
        dplyr::select( -dplyr::any_of(c(
          "hey_xxxxxxxxxx"
          , "landfire_cell_kg_per_m3"
        )))
      # rename cell data
      cell_df <- cell_df %>% dplyr::rename(landfire_stand_id = cell)
    }

    #############################################
    # return data
    #############################################
    #return
    return(list(
      cell_df = cell_df
      , tree_list = tree_list
    ))
  }
#######################################################
# intermediate function 14
#######################################################
  # check vector for forest_type_group_code available in cruz
  has_cruz_forest_type_group_code <- function(x) {
    # check if data.frame
    if(inherits(x, "data.frame")){
      nms <- x %>% names() %>% dplyr::coalesce("")
      if(
        !(stringr::str_equal(nms, "forest_type_group_code") %>% any())
      ){
        stop(paste0(
          "data must contain `forest_type_group_code` column"
        ))
      }
      x <- x$forest_type_group_code
    }
    # cruz codes
    cruz_codes <- c(
      200 #Douglas-Fir Group
      , 220 #Ponderosa Pine Group
      , 280 #Lodgepole Pine Group
      , c(120,260,320) #Mixed Conifer Group
    )
    x <- as.numeric(x)
    # check it
    n_cruz_codes <- x[x %in% cruz_codes] %>% length()
    # return
    return(
      dplyr::coalesce(n_cruz_codes, 0) > 0
    )
  }
#######################################################
# intermediate function 14.4
#######################################################
  # check the `method` argument of the trees_biomass() function
  check_biomass_method <- function(method) {
      # clean method
      method <- dplyr::coalesce(method, "") %>%
        tolower() %>%
        stringr::str_squish()
      # potential methods
      pot_methods <- c("cruz", "landfire") %>% unique()
      find_method <- paste(pot_methods, collapse="|")
      # can i find one?
      which_methods <- stringr::str_extract_all(string = method, pattern = find_method) %>%
        unlist() %>%
        unique()
      # make sure at least one is selected
      # get_list_diff() from get_url_data.R
      n_methods_not <- get_list_diff(pot_methods, which_methods) %>% length()
      if(n_methods_not>=length(pot_methods)){
        stop(paste0(
          "`method` parameter must be one or multiple of:\n"
          , "    "
          , paste(pot_methods, collapse=", ")
        ))
      }else{
        return(which_methods)
      }
  }
#######################################################
# intermediate function 15
#######################################################
  # check data for column missing in data or if all records missing
  # will throw an error if either condition
  check_df_cols_all_missing <- function(df, col_names, all_numeric=T, check_vals_missing=T) {
    # check if data.frame
    if(!inherits(df, "data.frame")){
      stop(paste0(
        "`df` must be a data.frame"
      ))
    }
    ######################################
    # ensure all columns exist
    ######################################
      nms <- names(df) %>% dplyr::coalesce("")
      has_cols <- col_names %>%
        purrr::map(function(x){
          stringr::str_equal(nms, x) %>%
          any() # do any columns match
        }) %>%
        unlist()

      if(all(has_cols)==F){
        stop(paste0(
          "the data does not contain the columns: "
          , paste(col_names[!has_cols], collapse = ", ")
          , "\n this data must exist"
        ))
      }
    ######################################
    # ensure all columns aren't missing all data
    ######################################
    if(check_vals_missing==F){return(T)}
    if(all_numeric == T){
      all_missing_cols <- df %>%
        sf::st_drop_geometry() %>%
        dplyr::select(dplyr::all_of(col_names)) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(as.numeric(.))))) %>%
        tidyr::pivot_longer(dplyr::everything()) %>%
        dplyr::mutate(all_missing = value==nrow(df)) %>%
        dplyr::filter(all_missing==T) %>%
        dplyr::pull(name)
    }else{
      all_missing_cols <- df %>%
        sf::st_drop_geometry() %>%
        dplyr::select(dplyr::all_of(col_names)) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.)))) %>%
        tidyr::pivot_longer(dplyr::everything()) %>%
        dplyr::mutate(all_missing = value==nrow(df)) %>%
        dplyr::filter(all_missing==T) %>%
        dplyr::pull(name)
    }

    if(length(all_missing_cols)>0){
      stop(paste0(
        "the columns listed below have all missing data:\n   "
        , paste(all_missing_cols, collapse = ", ")
      ))
    }

    # return
    return(T)
  }


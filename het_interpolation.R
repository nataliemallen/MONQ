# Set working directory
setwd("/Users/natal/Documents/Purdue/MONQ/shortread_analysis")

# read your csv
genetic_data <- read_csv("prelim_het.csv")

# Spatial Interpolation of Heterozygosity Values
# Libraries needed for modern GIS analysis and visualization
library(sf)         # Simple Features for spatial vector data
library(terra)      # Modern replacement for raster package
library(gstat)      # Geostatistical modeling
library(stars)      # Spatiotemporal arrays
library(ggplot2)    # Plotting
library(readr)      # CSV reading
library(maps)       # Map data
library(mapdata)    # Additional map data

# Check the data structure
head(genetic_data)

# 2. Convert the data to modern spatial objects using sf and terra
# Create an sf object for the point data - here coordinates are already in x,y columns
sf_data <- st_as_sf(genetic_data, 
                    coords = c("x", "y"), 
                    crs = 4326)  # WGS84 (standard for lon/lat)

# Create a SpatVector from terra package (needed for some operations)
# Use sf directly to avoid coordinate issues
terra_points <- vect(sf_data)

# 3. Create a grid for interpolation using terra
# Define the boundaries (adjust based on your data coverage)
# These coordinates roughly cover Mexico, Texas, and New Mexico
ext_buffer <- 1  # buffer around the data extent

# Carefully check and build the extent
x_min <- min(genetic_data$x, na.rm = TRUE) - ext_buffer
x_max <- max(genetic_data$x, na.rm = TRUE) + ext_buffer
y_min <- min(genetic_data$y, na.rm = TRUE) - ext_buffer
y_max <- max(genetic_data$y, na.rm = TRUE) + ext_buffer

# Print the extents to verify them
cat("X extent:", x_min, "to", x_max, "\n")
cat("Y extent:", y_min, "to", y_max, "\n")

# Create the extent object safely
data_extent <- ext(x_min, x_max, y_min, y_max)

# Create a raster template for interpolation with explicit resolution
x_res <- (x_max - x_min) / 100
y_res <- (y_max - y_min) / 100

rast_template <- rast(data_extent, 
                      resolution = c(x_res, y_res),
                      crs = "EPSG:4326")

# 4. Prepare data for gstat (still needed for interpolation)
# For gstat compatibility, we need to use sp objects but we can create them properly
library(sp)  # We still need sp for gstat compatibility

# Create spatial points from data frame
sf_df <- data.frame(x = genetic_data$x, 
                    y = genetic_data$y, 
                    heterozygosity = genetic_data$Heterozygosity)

# Create spatial points data frame - use x,y as coordinates directly
coordinates(sf_df) <- c("x", "y")

# Set the CRS properly - sp still needs CRS object
sp::proj4string(sf_df) <- sp::CRS("+proj=longlat +datum=WGS84")

# Create a grid for gstat - properly setting up coordinates
grid_df <- as.data.frame(xyFromCell(rast_template, 1:ncell(rast_template)))
names(grid_df) <- c("x", "y")
coordinates(grid_df) <- c("x", "y")
sp::proj4string(grid_df) <- sp::CRS("+proj=longlat +datum=WGS84")
gridded(grid_df) <- TRUE

# Perform IDW interpolation
idw_model <- gstat(formula = heterozygosity ~ 1, 
                   locations = sf_df,
                   nmax = 12, # Maximum number of neighboring points to consider
                   set = list(idp = 2)) # Power parameter for IDW

# Apply interpolation to the grid
idw_interpolation <- predict(idw_model, grid_df)

# Convert to SpatRaster for visualization
idw_raster <- rast(idw_interpolation)

# 5. Create a nice visualization
# First, get map data for context
usa <- map_data("usa")
mexico <- map_data("worldHires", region = "Mexico")

# Convert idw_interpolation to a data frame for ggplot
idw_df <- as.data.frame(idw_interpolation)
# No need to rename, the column names are already correct from coordinates() above

# Plot the results
ggplot() +
  # Add the interpolated surface
  geom_tile(data = idw_df, 
            aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis_c(name = "Heterozygosity", 
                       option = "plasma",
                       limits = c(min(genetic_data$Heterozygosity), 
                                  max(genetic_data$Heterozygosity))) +
  # Add state/country borders for context
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 0.5) +
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 0.5) +
  # Add the original data points
  geom_point(data = genetic_data, aes(x = x, y = y), 
             color = "black", size = 2) +
  # Customize the plot
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Interpolation of Heterozygosity Values",
       subtitle = "Inverse Distance Weighting (IDW) Method",
       x = "Longitude", 
       y = "Latitude")

# 6. Alternative: Kriging interpolation (generally more accurate than IDW)
# First fit a variogram model to the data
v <- variogram(heterozygosity ~ 1, sf_df)
v_fit <- fit.variogram(v, model = vgm(psill = var(genetic_data$Heterozygosity)/2, 
                                      "Sph", range = mean(diff(range(genetic_data$x)))/3, 
                                      nugget = var(genetic_data$Heterozygosity)/10))
plot(v, v_fit) # Check how well the variogram model fits

# Create kriging model
kriging_model <- gstat(formula = heterozygosity ~ 1, 
                       locations = sf_df,
                       model = v_fit)

# Apply kriging to the grid
kriging_interpolation <- predict(kriging_model, grid_df)

# Convert kriging_interpolation to a data frame for ggplot
# Use proper column names based on the actual data structure
kriging_df <- as.data.frame(kriging_interpolation)

# Create a nice visualization for kriging results
ggplot() +
  # Add the interpolated surface
  geom_tile(data = kriging_df, 
            aes(x = x, y = y, fill = var1.pred)) +
  scale_fill_viridis_c(name = "Heterozygosity", 
                       option = "viridis",
                       limits = c(min(genetic_data$Heterozygosity), 
                                  max(genetic_data$Heterozygosity))) +
  # Add state/country borders for context
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 0.5) +
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black", size = 0.5) +
  # Add the original data points  
  geom_point(data = genetic_data, aes(x = x, y = y), 
             color = "black", size = 2) +
  # Customize the plot
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Interpolation of Heterozygosity Values",
       subtitle = "Ordinary Kriging Method",
       x = "Longitude", 
       y = "Latitude")

#create R object for country boundarues
require("maptools")
CountryPolySP <- readShapePoly("/Net/Groups/BSY/people/tkoch/gadm_database/world_country_admin_boundary_shapefile_with_fips_codes.shp")
save("CountryPolySP",file="CountryPolySP.rdata")

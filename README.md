## Install packages

```R
install.packages("lidR")
install.packages("raster")
install.packages("sf")
install.packages("terra")
install.packages("moments")
```

## Import libraries

```R
library(lidR)
library(raster)
library(sf)
library(terra)
library(moments)
```

## Read LAS data

```R
las = lidR:readLAS("<path>/<to>/file.las")
```

## Preprocessing LAS data

```R
# Terrain normalize
las = normalize_height(las, knnidw())
```

## Create Canopy Height Model

```R
# Create CHM map with 0.3 m resolution
chm = rasterize_canopy(las, res=0.3, dsmtin())

# Smoothen CHM map with 3 window kernel
chm = terra::focal(chm, matrix(1,3,3), fun = median)
```

## Segment trees

```R
# li2012 algorithm
seg_las = segment_trees(las, li2012())


# Dalponte2016
treetops = locate_trees(las, lmf(5))
seg_las = segment_trees(las, dalponte2016(chm, treetops))

```

## Calculate tree metrics

```R

  getmode = function(v){
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  getAllMetrics = function(Z, Intensity, ReturnNumber, NumberOfReturns, above=2.0, minht=1.37){
    
    
    Intensity = as.double(Intensity[which(Z>=minht)])
    ReturnNumber = as.double(ReturnNumber[which(Z>=minht)])
    NumberOfReturns = as.double(NumberOfReturns[which(Z>=minht)])
    Z = as.double(Z[which(Z>=minht)])

    metrics = list(
     Total_all_return_count = length(Z),
     Total_first_return_count = length(Z[which(ReturnNumber==1)]),
     Total_all_return_count_above_minht = length(Z[which(Z>above)]),
     Return_1_count_above_minht = length(Z[which(Z>above & ReturnNumber==1)]),
     Return_2_count_above_minht = length(Z[which(Z>above & ReturnNumber==2)]),
     Return_3_count_above_minht = length(Z[which(Z>above & ReturnNumber==3)]),
     Return_4_count_above_minht = length(Z[which(Z>above & ReturnNumber==4)]),
     Return_5_count_above_minht = length(Z[which(Z>above & ReturnNumber==5)]),
     Return_6_count_above_minht = length(Z[which(Z>above & ReturnNumber==6)]),
     Return_7_count_above_minht = length(Z[which(Z>above & ReturnNumber==7)]),
     Return_8_count_above_minht = length(Z[which(Z>above & ReturnNumber==8)]),
     Return_9_count_above_minht = length(Z[which(Z>above & ReturnNumber==9)]),
     HMIN = min(Z),
     HMAX = max(Z),
     HMEAN = mean(Z),
     HMODE = getmode(Z),
     HMEDIAN = as.double(median(Z)),
     HSD = sd(Z),
     HVAR = var(Z),
     HCV = as.numeric(sd(Z) / mean(Z)),
     HKUR = kurtosis(Z),
     HSKE = skewness(Z),
     H01TH = quantile(Z, 0.01),
     H05TH = quantile(Z, 0.05),
     H10TH = quantile(Z, 0.10),
     H15TH = quantile(Z, 0.15),
     H20TH = quantile(Z, 0.20),
     H25TH = quantile(Z, 0.25),
     H30TH = quantile(Z, 0.30),
     H35TH = quantile(Z, 0.35),
     H40TH = quantile(Z, 0.40),
     H45TH = quantile(Z, 0.45),
     H50TH = quantile(Z, 0.50),
     H55TH = quantile(Z, 0.55),
     H60TH = quantile(Z, 0.60),
     H65TH = quantile(Z, 0.65),
     H70TH = quantile(Z, 0.70),
     H75TH = quantile(Z, 0.75),
     H80TH = quantile(Z, 0.80),
     H90TH = quantile(Z, 0.90),
     H95TH = quantile(Z, 0.95),
     H99TH = quantile(Z, 0.99),
     Canopy_relief_ratio = (mean(Z)-min(Z))/(max(Z)-min(Z)),
     IMIN = min(Intensity),
     IMAX = max(Intensity),
     IMEAN = mean(Intensity),
     IMODE = getmode(Intensity),
     IMEDIAN = as.double(median(Intensity)),
     ISD = sd(Intensity),
     IVAR = var(Intensity),
     ICV = as.numeric(sd(Intensity) / mean(Intensity)),
     IKUR = kurtosis(Intensity),
     ISKE = skewness(Intensity),
     I01TH = quantile(Intensity, 0.01),
     I05TH = quantile(Intensity, 0.05),
     I10TH = quantile(Intensity, 0.10),
     I15TH = quantile(Intensity, 0.15),
     I20TH = quantile(Intensity, 0.20),
     I25TH = quantile(Intensity, 0.25),
     I30TH = quantile(Intensity, 0.30),
     I35TH = quantile(Intensity, 0.35),
     I40TH = quantile(Intensity, 0.40),
     I45TH = quantile(Intensity, 0.45),
     I50TH = quantile(Intensity, 0.50),
     I55TH = quantile(Intensity, 0.55),
     I60TH = quantile(Intensity, 0.60),
     I65TH = quantile(Intensity, 0.65),
     I70TH = quantile(Intensity, 0.70),
     I75TH = quantile(Intensity, 0.75),
     I80TH = quantile(Intensity, 0.80),
     I90TH = quantile(Intensity, 0.90),
     I95TH = quantile(Intensity, 0.95),
     I99TH = quantile(Intensity, 0.99),
     Pentage_first_returns_above = length(Z[which(ReturnNumber==1&Z>above)])/length(Z)*100,
     Percentage_all_returns_above = length(Z[which(Z>above)])/length(Z)*100,
     X_All_returns_above_Total_first_returns_100 = length(Z[which(Z>above)])/length(Z[which(ReturnNumber==1)])*100,
     First_returns_above_above = length(Z[which(Z>above & ReturnNumber==1)]),
     All_returns_above_above = length(Z[which(Z>above)]),
     Percentage_first_returns_above_mean = length(Z[which(ReturnNumber==1&Z>mean(Z))])/length(Z)*100,
     Percentage_first_returns_above_mode = length(Z[which(ReturnNumber==1&Z>getmode(Z))])/length(Z)*100,
     Percentage_all_returns_above_mean = length(Z[which(Z>mean(Z))])/length(Z)*100,
     Percentage_all_returns_above_mode = length(Z[which(Z>getmode(Z))])/length(Z)*100,
     X_All_returns_above_mean_Total_first_returns_100 = length(Z[which(Z>mean(Z))])/length(Z[which(ReturnNumber==1)])*100,
     X_All_returns_above_mode_Total_first_returns_100 = length(Z[which(Z>getmode(Z))])/length(Z[which(ReturnNumber==1)])*100,
     First_returns_above_mean = length(Z[which(ReturnNumber==1 & Z > mean(Z))]),
     First_returns_above_mode = length(Z[which(ReturnNumber==1 & Z > getmode(Z))]),
     All_returns_above_mean = length(Z[which(Z>mean(Z))]),
     All_returns_above_mode = length(Z[which(Z>getmode(Z))])
    )
    
    return(metrics)
  }

#tree metrics calculation

# Calculate metrics by canopy convex shape
 metrics_tree = crown_metrics(seg_las, func = ~getAllMetrics(Z, Intensity, ReturnNumber, NumberOfReturns, above = 2.0, minht = 1.37), geom="convex")

# Calculate metrics by canopy locations
 metrics_tree = crown_metrics(seg_las, func = ~getAllMetrics(Z, Intensity, ReturnNumber, NumberOfReturns, above = 2.0, minht = 1.37))

  
```

## Export tree crowns with metrics

```R

# save metrics as canopy boundary
sf::st_write(sf::st_as_sf(metrics_tree), "<path>/<to>/file.geojson")

# Save tree tops locations
sf::st_write(sf::st_as_sf(treetops), "<path>/<to>/file.geojson")


writeRaster(chm, "<path>/<to>/file.tif")
```
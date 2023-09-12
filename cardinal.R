library(Cardinal)
library(glue)

## --- CROPPING OUT BACKGROUND -----------------

imzmlFile <- "Thyroid/PTC_MAD12-372-1_20220823_CLMC_Centroid.raw/PTC_MAD12-372-1_20220823_CLMC_Centroid_1"
#imzmlFile <- "Thyroid/PTC_MAD12-372-1_20220823_CLMC.raw/PTC_MAD12-372-1_20220823_CLMC" ## takes long time to process because not centroided?

imzml_data <- readImzML(imzmlFile)

## View entire image
#image(imzml_data, mz=885.55)
image(imzml_data, mz=514.2670)

## K-MEANS CLUSTER K=2, REMOVE BACKGROUND CLUSTER
## cluster with k=2 to remove background
## need to do peak picking beforehand?

#plot(imzml_data, pixel=775:785, xlim=c(880,890),key=FALSE,superpose=TRUE)
## peaks are already aligned, whether in profile or centroided
## already no noise?

data_mean <- summarizeFeatures(imzml_data, "mean") ## summarizes an MSImagingExperiment by feature, calculates mean spectrum
#plot(data_mean)

data_ref <- data_mean %>% ##  local maxima of the mean spectrum are used as the reference
  peakPick(SNR=3) %>%
  peakAlign(ref="mean",
            tolerance=0.5,
            units="mz") %>%
  peakFilter(freq.min = 0.1) %>%
  process()
## keeping 1422 peaks

data_peaks <- imzml_data %>%
  normalize(method="tic") %>%
  peakBin(ref=mz(data_ref),
          tolerance=0.5,
          units="mz") %>%
  process()


#plot(data_peaks, pixel=775:785, xlim=c(880,890),key=FALSE,superpose=TRUE)

## use spatial shrunken centroids to segment the dataset.
set.seed(1)

## method --> method to use to calculate the spatial smoothing weights. 
## 'gaussian' method refers to spatially-aware (SA) weights
##'adaptive' refers to spatially-aware structurally-adaptive (SASA) weights.
## s --> The sparsity threshold parameter by which to shrink the t-statistics.

data_segm <- spatialShrunkenCentroids(data_peaks, method="adaptive", r=2, s=c(5,10,15,20), k=2)
image(data_segm)

## Keep sample region
cluster_indices <- c(1)
pixels_in_cluster <- data_segm@resultData@listData[[1]]$class%in%cluster_indices
data_tissue <- data_peaks[,pixels_in_cluster]

## --- CLUSTER REGIONS -------------------------

#data_ssc_1 <- spatialShrunkenCentroids(data_tissue, method="adaptive", r=2, s=c(10,15,20), k=c(2,3,5,10))
#image(data_ssc_1, col=hcl.colors(10,"viridis"))

#data_ssc_2 <- spatialShrunkenCentroids(data_tissue, method="adaptive", r=2, s=c(5,10,15), k=c(3,4,5,6))
#image(data_ssc_2, col=hcl.colors(6,"viridis"))

data_ssc <- spatialShrunkenCentroids(data_tissue, method="adaptive", r=2, s=5, k=4)
image(data_ssc,col=hcl.colors(4,"Viridis"))

## --- m/z VALUES ASSOCIATED WITH CLUSTERS -------
## find the m/z values associated with different regions of the tissue

## Class 1
#topFeatures(data_ssc, class==1)
#image(data_tissue, mz=810.5211)

## Loop through classes and top feature mz values and save images

## helper function
find_top_mz <- function(file,region_class){
  features <- topFeatures(file,class == region_class)
  mz_values <- features@listData$mz
  #print(mz_values)
}

## TO VIEW IN R
for (r in 1:length(ROI_classes)){
  #print(r)
  mz_values <- find_top_mz(data_ssc,r)
  for (m in 1:length(mz_values)){
    #print(mz_values[[m]])
    print(image(data_tissue, mz=mz_values[[m]]))
    title(glue("Class {r} Feature"), adj=0)
  }
}

## TO PRINT
for (r in 1:length(ROI_classes)){
  #print(r)
  mz_values <- find_top_mz(data_ssc,r)
  for (m in 1:length(mz_values)){
    #print(mz_values[[m]])
    jpeg(file=glue("class-{r}_mz-{format(round(mz_values[[m]],4))}.jpeg"))
    print(image(data_tissue, mz=mz_values[[m]]))
    title(glue("Class {r} Feature"), adj=0)
    dev.off()
  }
}



## --- CONVERT PROCESSED TO CONTINUOUS IMZML --------------

imzmlFile <- "Thyroid/PTC_MAD12-372-1_20220823_CLMC_Centroid.raw/PTC_MAD12-372-1_20220823_CLMC_Centroid_1"

imzml_data <- readImzML(imzmlFile) ## load imzml and ibd

data_cont <- as(imzml_data,"MSContinuousImagingExperiment") ## taking a long time, warning message that result is no longer an s4 object

writeImzML(data_cont, "data_test",outformat = "imzML", mz.type = "32-bit float",intensity.type = "32-bit float")

## --- CONVERT BETWEEN PROCESSED AND CONTINUOUS

#image(data_tissue, mz=514.2670)
#data_tissue_proc <- as(data_tissue,"MSProcessedImagingExperiment")

## export ROI as imzML
#writeMSIData(data_tissue_proc,"test",outformat = "imzML")
#writeImzML(data_tissue_proc, "test4")

## View ROI
#image(data_tissue, mz=885.55)
image(data_tissue, mz=514.2670)



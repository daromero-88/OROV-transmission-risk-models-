#####OROV TRANSMISSION RISK MODELING######

--------------------------------------------------------
  #'  Author: 
  #'  Daniel Romero-Alvarez
  #'  Luis Osorio-Olvera 
  #'  
  #' The following functions and scripts will allow the user to replicate 
  #' the models presented in the manuscript: 
  #' "Transmission Risk of Oropouche fever in the Americas"
  #' 
  #' The models include: 
  #' Hypervolumes
  #' Convex Hulls
  #' Jacknife for Hypervolumes
  #' Jacknife for Convex hulls
  #' 
  #' The script also includes the filtered points used for creating the models 
  #' and the environmental predictors (raw and principal components used)
  #' 
  #' The script replicates the results using the entire America continent 
  #' for creating the models and estimating the occurrence contribution.
  #' If needed, environmental predictors can be changed in order to obtain 
  #' models for each studied area presented (see main text and supplementary material)
  
  --------------------------------------------------------
  
#LIBRARIES-----------------------------------------------

#Spatial analysis
library(raster)
library(sp)
library(rgeos)
library (rgdal)
library (maptools)
library(mapdata)
library(spThin)

#Spetial shape objects for plotting
library (sf)

#Hypervolumes
library (Rcpp)
library (hypervolume)

#To use the check.in.range function, presence of values in one dimension
library (spatstat.utils)

#For functions provided by Luis Osorio-Olvera: 
library(rgl)
library(geometry)
library(ptinpoly)

#For evaluation in pair-wise dimensions: 
library(concaveman)
library(sp)

#lineplots: 
library (ggplot2)


#WORKING DIRECTORY-----------------------------

#change working directory accordingly to your needs: 

setwd ('/Users/daniel/Documents/OROPUCHE/PREDICTION/NEW_PROJ_2/_MANUSCRIPT1/DATA_GITHUB')


#DATA----------------------------------------------------

#OCURRENCES: 
#Curated points:
oro_thin = read.csv ('oro2_thin_thin1.csv')


#SHAPES: 
#WorldShape 
data("wrld_simpl", package = "maptools")
WGS84 = crs(wrld_simpl) # geographic projection

#unique regions of the world 
unique (wrld_simpl@data$REGION)

#subsetting to the region of the Americas: 
ss = subset (wrld_simpl, REGION == '19') 

#specific shapes for the Americas: 
#Americas downloaded from 
#https://tapiquen-sig.jimdofree.com/english-version/free-downloads/americas/
amr = readOGR ('./SHAPES/Americas/Americas.shp')
amr = buffer(amr, width = 0.1, dissolve = T)

#America provinces: 
provs_wr = readOGR ('./SHAPES/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp')

provs_df = crop (provs_wr, amr) #an error will show up but should be ignored
class(provs_df)


#Complete South America 
#https://tapiquen-sig.jimdofree.com/english-version/free-downloads/south-america/
sa1 = readOGR('./SHAPES/south_america/South_America.shp')
sa1 = buffer(sa1, width = 0.1, dissolve = T) #adding some degrees of buffer to capture all the pixels

#Impartial bufer: 
buff = readOGR ('./SHAPES/buffer/buffer1.shp')


#visualization of occurrences: 
plot (ss, col = '#dbdbdb')
zoom(ss, col = '#dbdbdb', ext=drawExtent(), new=TRUE, useRaster=FALSE)
points (oro_thin[,3:4], col = c('gold', 'red'), pch = c(16, 3), cex = 1.2)



#ENVIRONMENTAL PREDICTORS: 
#Raw climatic predictors: 

#buffer area
buf_land_ras = stack (list.files ('./raw_predictors/buf_land_ENVS', full.names = T, pattern = '.asc'))
names (buf_land_ras)

#South America
sa1 = stack (list.files ('./raw_predictors/sa1', full.names = T, pattern = '.asc'))
names (sa1)

#Americas continent 
amr1 = stack (list.files('./raw_predictors/amr1', full.names = T, pattern = '.asc'))
names (amr1)


#Principal components: 

#Buffer area: 
pca_buff_res = read.table ('./PCs/1_pca_buff/Initial/pca_results.txt', header = T, sep = '\t') #PCA output
pca_buff = stack(list.files('./PCs/1_pca_buff/Initial', full.names = T, pattern = 'asc')) #reading five PCs
pca_buf2 = pca_buff[[c(1, 2, 3)]] #selecting three first PCs
names (pca_buf2)

#South America: 
pca_sa_res = read.table ('./PCs/1_pca_sa/Initial/pca_results.txt', header = T, sep = '\t')
pca_sa = stack(list.files('./PCs/1_pca_sa/Initial', full.names = T, pattern = 'asc'))
pca_sa2 = pca_sa[[c(1, 2, 3)]]
names (pca_sa2)

#The Americas continent: 
pca_amr_res = read.table ('./PCs/1_pca_amr/Initial/pca_results.txt', header = T, sep = '\t')
pca_amr = stack(list.files('./PCs/1_pca_amr/Initial', full.names = T, pattern = 'asc'))
pca_amr2 = pca_amr[[c(1, 2, 3)]]
names (pca_amr2)


#FUNCTIONS-----------------------------------------------

#FUNCTION TO DEVELOP THREE DIMENSIONAL CONVEX HULLS--------------------------- 

#'Provided by Luis Osorio-Olvera,
#'Modified by Daniel Romero-Alvarez,

convexhull3d <- function(data,test_points=NULL,plot=F,centroid=NULL, background = NULL){
  if(dim(data)[2]==3){
    # Code for creating the convex hull 
    ch_data = convhulln(data, output.options = 'FA')  #recover indexes that define the convex hull, and with FA, it recovers the volume and other attributes
    ts.surf1 <- t(ch_data$hull) 
    
    convexhull_faces <- cbind(data[ts.surf1,1], #obtaining the vertices of the hull in the parameter space
                              data[ts.surf1,2],
                              data[ts.surf1,3]) 
    
    #this ifelse portion is defined because the centroid should be defined in the function! 
    if(is.null(centroid)){
      centroide <- colSums(convexhull_faces)/(dim(convexhull_faces)[1])
    }else{
      centroide <- centroid
    }
    
    #Plotting:   
    if(plot){
      convex_mesh <- tmesh3d(
        vertices = t(data),
        indices = ts.surf1,
        homogeneous = FALSE)
      
      plot3d(data,size=5)
      plot3d(convexhull_faces,type="p", size=10,col="blue",add=T) # plotting the vertices
 
      wire3d(convex_mesh, col = 4)
     }
    
  }
  
  #Are evaluating testing points inside the hull?: 
  if(dim(test_points)[1]>=1 && dim(test_points)[2]==3){
    in_convex <- pip3d(Vertices = as.matrix(data),
                       Faces = t(ts.surf1),
                       Queries = as.matrix(test_points))
    in_convex <- sapply(in_convex,function(x){
      if(x==1) return(1)
      if(x==0) return(1)
      if(x==-1) return(0)
    })
    dist_eucli <- unlist(lapply(X = 1:dim(test_points)[1],
                                FUN=function(x) euc.dist(centroide,
                                                         test_points[x,])))
    dconvex <- data.frame(dist_eucli,in_convex)
  }
  
  #Where are the background points in the hull?: 
  if(dim(background)[1]>1 && dim(background)[2]==3){
    in_convex_bc <- pip3d(Vertices = as.matrix(data),
                          Faces = t(ts.surf1),
                          Queries = as.matrix(background))
    in_convex_bc <- sapply(in_convex_bc,function(x){
      if(x==1) return(1)
      if(x==0) return(1)
      if(x==-1) return(0)
    })
    dist_eucli_bc <- unlist(lapply(X = 1:dim(background)[1],
                                   FUN=function(x) euc.dist(centroide,
                                                            background[x,])))
    dconvex_bc <- data.frame(dist_eucli_bc,in_convex_bc)
    
    #Collecting results: the centroid of the 3d hull, the list of evaluating points, the volume of the hull, the list of background points
    results <- list(centroid=centroide,
                    testing = dconvex, 
                    vol = ch_data[3], 
                    backgroound = dconvex_bc)
  }
  
  return(results)
}

#measuring Euclidean distances: 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) 


#SCRIPTS-------------------------------------------------


#HYPERVOLUME: OBSERVED MODEL, CALIBRATION LOOP, EVALUATION STATS, MODELS---------------

#PCs in the Americas: 
pca_amr2

#Original points
oro_thin 

#Extracting values considering rasters provided, raw data for the hypervolume function: 
clusters8 = data.frame(cbind (oro_thin[,1:4], extract(pca_amr2, oro_thin[,3:4])))

#hypervolume (svm) with provided information --> values of each environmental layer at each point: 
svm_amr = hypervolume_svm(data = clusters8[,5:7])
plot (svm_amr)

#total volume: 
svm_amr_tvol = svm_amr@Volume #obtaining the volumne of the model 

#Projecting to geography: 
svm_amr_g = hypervolume_project (hv = svm_amr,rasters = pca_amr2, type = "inclusion")
plot (svm_amr_g) #plotting projected model 

#total number of pixels: 
svm_amr_pxs = rasterToPoints(svm_amr_g) #transforming the projected model to points 
dim(svm_amr_pxs) #3438508 #total pixels 

svm_amr_tpxs = nrow(svm_amr_pxs[which (svm_amr_pxs[,3]==1),]) #colleting only those pixels with suitable value = 1

#writing results for the observed model, environmental dataframe:
#environmental values of the model
write.csv (svm_amr@RandomPoints, paste ('./HV_results/hvamr_pca_tot_e','.csv', sep = ''), row.names = F)

#writing results for the observed model, raster file:
writeRaster(svm_amr_g, filename = paste ('./HV_results/hvamr_pca_tot_g', sep = ''), format = 'ascii',
            bylayer = T, suffix = paste ('hvamr_clim_tot_g', j, sep = ''), overwrite = T,
            NAflag = -9999)


#Model calibration: 

eval_in_e4 = list () #collecting bins in E
eval_in_e4_df = NULL #dataframe of results in E
eval_in_g4 = list () #collecting bins in G 
eval_in_g4_df = NULL #dataframe of results in G

a3 = stack () #collecting the stack of raster files 

#CALIBRATION LOOP (50 boostrap):

for (j in 1:50){
  #subsetting database in 70 calibration, 30 evaluation: 
  eval30 = clusters8 [sample(nrow(clusters8),round(nrow(clusters8)*0.30)),]
  #dim(eval30)
  cal70 = clusters8 [!(clusters8$ID %in% eval30$ID),]
  #dim(cal70)
  
  cal70_val = extract (pca_amr2, cal70[,3:4]) #values for model calibration
  eval30_val = extract (pca_amr2, eval30[,3:4]) #values for model evaluation
  
  svm_e = hypervolume_svm(data = cal70_val) #calibrating the model, only three dimensions used here
  
  #Evaluation in the environmental space
  #Are the evaluation points contained in the environmental cloud developed with the calibration points? 
  #Approach by Luis Osorio-Olvera: 
  
  #Dimension 1: 
  polygons1 = concaveman(as.matrix(svm_e@RandomPoints[,1:2]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in1 = point.in.polygon(eval30_val[,1],eval30_val[,2],polygons1[,1],polygons1[,2]) 
  
  #Dimension 2: 
  polygons2 <- concaveman(as.matrix(svm_e@RandomPoints[,2:3]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in2 <- point.in.polygon(eval30_val[,2],eval30_val[,3],polygons2[,1],polygons2[,2]) 
  
  #Dimension 3: 
  polygons3 <- concaveman(as.matrix(svm_e@RandomPoints[,c(1,3)]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in3 <- point.in.polygon(eval30_val[,1],eval30_val[,3],polygons3[,1],polygons3[,2]) 
  
  #plotting if needed: 
  #plot (polygons2, xlim = c(50,500), ylim = c(600,3000))
  #lines(polygons2)
  #points (eval30_val[,c(2,3)], col = 'red')
  
  red1_res = mean(c(mean (ids_in1), mean (ids_in2), mean (ids_in3))) #obtaing the mean across dimensions, meassure of performance in E
  
  eval_in_e4[[length(eval_in_e4)+1]] = red1_res #collecting the mean as a list 
  
  #Dataframe of replicates data in the environmental space (E): 
  df_e = cbind (replicate = j, #replicate number 
                perform = red1_res, #collecting the evaluation statistic of the replicate
                vol = svm_e@Volume, #collecting volume of the replicate 
                proportion = svm_e@Volume/svm_amr_tvol) #obtaining the proportion of volume as a function of the full model 
  
  eval_in_e4_df = rbind (eval_in_e4_df,df_e) #adding info to the results dataframe 
  
  ##If environmental points for each replicate are needed, the following lines write that data:
  ##writing the random points around the occurrences (the environmental buffer)
  #write.csv (svm_e@RandomPoints, paste ('./HV_results/env_sa_', j, '.csv', sep = ''), row.names = F)
  
  #Projection to geographic space: 
  svm_g = hypervolume_project (hv = svm_e, rasters = pca_amr2, type = "inclusion") #using the corresponding rasters! 
  
  ##If needed, you can write the geographic rasters for each replicate
  #writeRaster(svm_g, filename = paste ('./HV_results/rep_sa_', j, sep = ''), format = 'ascii',
  #            bylayer = T, suffix = paste ('rep_sa_', j, sep = ''), overwrite = T,
  #            NAflag = -9999)
  
  #Evaluation in the Geographic space (G):
  #Are the evaluation points contained in the raster obtained with the calibration points? 
  
  r_t_p = rasterToPoints(svm_g) #transforming raster to points 
  prop_pxs = nrow(r_t_p[which(r_t_p[,3]==1),])/svm_amr_tpxs #obtaining the proportion of suitable pixels per replicate as function of the full raster 
  
  tst = mean(extract (svm_g, eval30[,3:4])) #the mean of values across the evaluation points 
  eval_in_g4[[length(eval_in_g4)+1]] = tst #collecting the mean as a list 
 
  #Dataframe of replicates data in the geographic space (G): 
  df_g = cbind (replicate = j, #replicate number
                perform = tst, #collecting the evaluation statistic of the replicate
                tot_suit_pxs = nrow(r_t_p[which(r_t_p[,3] ==1),]), #collecting number of suitable pixels of the replicate 
                prop_pxs = prop_pxs) #collecting the proportion of suitable pixels per replicate 
  
  eval_in_g4_df = rbind (eval_in_g4_df,df_g) #adding info to the results dataframe 
  
  #rasters#
  a3 = stack (a3, svm_g) #collecting each raster per replicate 
}

#POST-LOOP: 

#Checking frequencies of evaluation stat in E and G: 
par(mfrow = c(2,1))
plot (eval_in_e4_df[,2], type = 'h')
plot (eval_in_g4_df[,2], type = 'h')

#Checking the correlation between volume and performance in E:
plot (eval_in_e4_df[,3], eval_in_e4_df[,2])

#Checking the correlation between suitable pixles and performance in G:
plot (eval_in_g4_df[,3], eval_in_g4_df[,2])

core4_e = cor.test(eval_in_e4_df[,3], eval_in_e4_df[,2], method = 'pearson') #is this relationship statistically significant in E
core4_g = cor.test(eval_in_g4_df[,3], eval_in_g4_df[,2], method = 'pearson') #is this relationship statistically significant in G

#Dataframe with performance statistics in E: 
re_svm_e4 = c(mean = mean (unlist (eval_in_e4)), #evaluation stat
              sd = sd (unlist (eval_in_e4)), #standard deviation of evaluation stat
              vol_tot = svm_amr_tvol, #volume of full model
              vol_median3 = median (eval_in_e4_df[,3]), #median of the replicates volume 
              vol_max3 = quantile (eval_in_e4_df[,3], 0.975),  #volume 97.5 percentile of the replicates
              vol_min3 = quantile (eval_in_e4_df[,3], 0.025), #volume 2.5 percentile of the replicates
              cor_perf_vol = core4_e$estimate, #correlation value between volume and performance
              cor_pe = core4_e$p.value) #p-value 

#Dataframe with performance statistics in G:
re_svm_g4 = c(mean = mean (unlist (eval_in_g4)), #evaluation stat
              sd = sd (unlist (eval_in_g4)), #standard deviation of evaluation stat
              pxls_tot = svm_amr_tpxs, #pixels of full model
              pxls_median3 = median (eval_in_g4_df[,3]), #median of the replicates suitable pixels 
              pxls_max3 = quantile (eval_in_g4_df[,3], 0.975), #pixels 97.5 percentile of the replicates
              pxls_min3 = quantile (eval_in_g4_df[,3], 0.025), #pixels 2.5 percentile of the replicates
              cor_perf_vol = core4_g$estimate, #correlation value between pixels and performance
              cor_pe = core4_g$p.value) #performance for niche: 

#Writing statistics dataframe
write.csv (re_svm_g4, './HV_results/stat_svm_g4.csv', row.names = T)
write.csv (re_svm_e4, './HV_results/stat_svm_e4.csv', row.names = T)

#Writing replicates dataframe
write.csv (eval_in_e4_df, './HV_results/svm_e4_df.csv', row.names = T)
write.csv (eval_in_g4_df, './HV_results/svm_g4_df.csv', row.names = T)

#Writing rasters: 

a_median3 = overlay(a3, fun = median, na.rm = T) #median raster across replicates
a_low3 = calc(a3, fun = function(x) {quantile(x,probs = c(0.025),na.rm=TRUE)}) #lower percentile raster across replicates
a_high3 = calc(a3, fun = function(x) {quantile(x,probs = c(0.975),na.rm=TRUE)}) #higher percentile raster across replicates
a_dif3 = a_high3-a_low3 #difference between percentiles
a_sum3 = a_high3+a_low3 #sum between percentiles 

a_all3= stack (a_median3, a_low3, a_high3, a_dif3, a_sum3) #stack of developed models 
nm5 = c('median', 'low', 'high', 'range', 'sum') #names for the stack 

#Wrting the rasters in a loop: 

for (jj in 1:5){
  writeRaster(a_all3[[jj]], filename = paste ('./HV_results/ras_output/', nm5[jj], sep = ''), format = 'ascii',
              bylayer = T, suffix = paste ('hvamrclim_', j, sep = ''), overwrite = T,
              NAflag = -9999)
}


#CONVEX HULLS: OBSERVED MODEL, CALIBRATION LOOP, EVALUATION STATS, MODELS---------------

#environments: 
rp_pca_amr2= rasterToPoints(pca_amr2)

#Model using all the points: 
clusnich_amr = cbind (oro_thin[,1:4],extract(pca_amr2, oro_thin[,3:4])) #adding all the data together

chlfull_amr=  convexhull3d(data = clusnich_amr[,5:7], 
                           test_points= clusnich_amr[,5:7], #same points as in previous line 
                           plot = T, #when T, it will depict the model in 3d 
                           background = rp_pca_amr2[,3:5]) #add the background points 

#checking points inside the hull: 
points3d(rp_pca_amr2[,3:5][chlfull_amr[[4]]$in_convex_bc==1,],
         col="red",size=5)

#checking background points outside the hull: 
points3d(rp_pca_amr2[,3:5][chlfull_amr[[4]]$in_convex_bc==0,],
         col="grey",size=1)

chl_amr_vol = chlfull_amr[[3]][[1]] #collecting the volume of the convex hull 

#Creating a database with all the data (important step, otherwise some data is going to be lost): 
nicha_amr_ex = data.frame(rp_pca_amr2[chlfull_amr[[4]]$in_convex_bc==0,]) #points outside the hull 
nicha_amr_in = data.frame(rp_pca_amr2[chlfull_amr[[4]]$in_convex_bc==1,]) #points inside the hull 

nicha_amr_ex$ind = rep (0, length(nicha_amr_ex[,1])) #adding value of 0 to the points outside the hull 
nicha_amr_in$ind = rep (1, length(nicha_amr_in[,1])) #adding value of 1 to the points inside the hull
nicha_amr_df = rbind (nicha_amr_ex, nicha_amr_in) #combining the dataframe 

#writing the dataframe complete dataframe: points and values
#writing results for the observed model, environmental dataframe:
write.csv (nicha_amr_df,'./CH_results/chamr_pca_tot_e.csv', row.names = F)

#Projection to geography: 
amr_ras = rasterize (nicha_amr_df[,1:2], pca_amr2[[1]], field = nicha_amr_df[,6])
amr_spx = dim(nicha_amr_in)[1] #total number of suitable pixels

#writing results for the observed model, raster file:
writeRaster(amr_ras, filename = './CH_results/chamr_pca_tot_g', format = 'ascii',
            bylayer = T, suffix = 'chamr_pca_tot_g', overwrite = T,
            NAflag = -9999)

#CALIBRATION LOOP (50 boostrap): 

eval_ch_e3 = list () #collecting bins in E
eval_ch_e3_df = NULL #dataframe of results in E
eval_ch_g3 = list () #collecting bins in G 
eval_ch_g3_df = NULL #dataframe of results in G

b3= stack () #stack of raster files 

for (j in 1:50){
  #subsetting database in 70 calibration, 30 evaluation: 
  eval30 = clusnich_amr [sample(nrow(clusnich_amr),round(nrow(clusnich_amr)*0.30)),]
  #dim(eval30)
  cal70 = clusnich_amr [!(clusnich_amr$ID %in% eval30$ID),]
  #dim(cal70)
  
  #Convex hull in environmental space: 
  ee = convexhull3d(data = cal70[,5:7], 
                    test_points= eval30[,5:7], 
                    plot = F , 
                    background = rp_pca_amr2[,3:5])
  
  per_ee = mean (ee$testing[,2]) #collecting evaluation statistic
  
  #checking points inside/outside convex hull: 
  # points3d(eval30[,4:6][ee[[2]]$in_convex==1,],
  #          col="red",size=5)
  # 
  # points3d(eval30[,4:6][ee[[2]]$in_convex==0,],
  #          col="darkgreen",size=5)
  
  vlm = ee[[3]][[1]] #collecting volume of convex hull for each replicate  
  
  eval_ch_e3[[length(eval_ch_e3)+1]]= per_ee #collecting statistic per replicate 
  
  #Dataframe of replicates in the environmental space (E): 
  res_rw = cbind (repl = j, #replicate number 
                  perform = per_ee, #performance statistic
                  vol = vlm,  #volume of the replicate convex hull
                  proportion = vlm/chl_amr_vol) #proportion of replicates' volume as function of the convex hull volumne using all the points 
  
  eval_ch_e3_df = rbind (eval_ch_e3_df, res_rw) #collecting information for the replicates dataframe
  
  #Projecting to geography: 
  ee_ex = data.frame(rp_pca_amr2[ee[[4]]$in_convex_bc==0,]) #points outside CH
  ee_in = data.frame(rp_pca_amr2[ee[[4]]$in_convex_bc==1,]) #points inside CH
  ee_ex$ind = rep (0, length(ee_ex[,1])) #adding value = 0, points outside
  ee_in$ind = rep (1, length(ee_in[,1])) #adding value = 1, points inside
  ee_env = rbind (ee_ex, ee_in) #dataframe with values and coordinates 
  
  st_pxs = dim(ee_in)[1] #collecting total suitable pixels per replicate 
  
  #Converting dataframe to raster 
  amr_ras = rasterize (ee_env[,1:2], pca_amr2[[1]], field = ee_env[,6])
  
  #plot(amr_ras)
  #points(eval30[,2:3], col ='red')

  #Evaluation statistics in geography: 
  tst2 = mean (extract (amr_ras, eval30[,3:4])) #collecting statistic in geography
  eval_ch_g3[[length(eval_ch_g3)+1]] = tst2 #collecting statistic as a list 
  
  #Dataframe of replicates in the geography (G): 
  res_geo = cbind (repl = j, #replicate number 
                   perform = tst2, #performance statistic
                   suit_pxls = st_pxs, #pixels of the replicate
                   proportion = st_pxs/amr_spx) #proportion of replicates' pixels as function of the convex hull total pixels with all the points 
  
  eval_ch_g3_df = rbind (eval_ch_g3_df, res_geo) #collecting dataframe 
  
  #Collecting replicate rasters 
  b3 = stack (b3, amr_ras) 
}

#POST-LOOP: 

#Checking frequencies of evaluation stat in E and G: 
par(mfrow = c(2,1))
plot (eval_ch_e3_df[,2], type = 'h')
plot (eval_ch_g3_df[,2], type = 'h')

#Checking the correlation between volume and performance in E:
plot (eval_ch_e3_df[,3], eval_ch_e3_df[,2])

#Checking the correlation between suitable pixles and performance in G:
plot (eval_ch_g3_df[,3], eval_ch_g3_df[,2])

core3_che = cor.test(eval_ch_e3_df[,3], eval_ch_e3_df[,2], method = 'pearson') #is this relationship statistically significant in E
core3_chg = cor.test(eval_ch_g3_df[,3], eval_ch_g3_df[,2], method = 'pearson') #is this relationship statistically significant in G

#Dataframe with performance statistics in E: 
re_chl_e3 = c(mean = mean (unlist (eval_ch_e3)), #evaluation statistic
              sd = sd (unlist (eval_ch_e3)), #standard deviation 
              vol_tot = chl_amr_vol, #volume using all the points 
              vol_median8 = median (eval_ch_e3_df[,3]), #median of the volume from replicates
              vol_max8 = quantile (eval_ch_e3_df[,3], 0.975), #higher quantile across replicates
              vol_min8 = quantile (eval_ch_e3_df[,3], 0.025), #lower quantile across replicates
              cor_perf_vol = core3_che$estimate, #correlation estimate 
              cor_pe = core3_che$p.value) #p value 

#Dataframe with performance statistics in G:
re_chl_g3 = c(mean = mean (unlist (eval_ch_g3)), 
              sd = sd (unlist (eval_ch_g3)), 
              pxls_tot = amr_spx, 
              pxls_median8 = median (eval_ch_g3_df[,3]),
              pxls_max8 = quantile (eval_ch_g3_df[,3], 0.975), 
              pxls_min8 = quantile (eval_ch_g3_df[,3], 0.025),
              cor_perf_vol = core3_chg$estimate, 
              cor_pe = core3_chg$p.value) 

#Writing statistics dataframe
write.csv (re_chl_e3, './CH_results/stat_chl_e3.csv', row.names = T)
write.csv (re_chl_g3, './CH_results/stat_chl_g3.csv', row.names = T)

#Writing replicates dataframes
write.csv (eval_ch_e3_df, './CH_results/chl_e3_df.csv', row.names = T)
write.csv (eval_ch_g3_df, './CH_results/chl_g3_df.csv', row.names = T)

#Definitive rasters: 

b_median3 = overlay(b3, fun = median, na.rm = T) #median across replicates 
b_low3 = calc(b3, fun = function(x) {quantile(x,probs = c(0.025),na.rm=TRUE)}) #higher percentile raster across replicates
b_high3 = calc(b3, fun = function(x) {quantile(x,probs = c(0.975),na.rm=TRUE)}) #lower percentile raster across replicates
b_dif3 = b_high3-b_low3 #difference between percentiles 
b_sum3 = b_high3+b_low3 #sum between percentiles 

b_all3= stack (b_median3, b_low3, b_high3, b_dif3, b_sum3) #stack of rasters 
nm5 = c('median', 'low', 'high', 'range', 'sum') #names for the rasters 

#Wrting the rasters in a loop: 

for (jj in 1:5){
  writeRaster(b_all3[[jj]], filename = paste ('./CH_results/ras_output/', nm5[jj], sep = ''), format = 'ascii',
              bylayer = T, suffix = paste ('chamr_pca_', j, sep = ''), overwrite = T,
              NAflag = -9999)
}


#OCCURRENCE CONTRIBUTION: JACKKNIFE IN CONVEX HULLS----------------------------------

#Select the raster for the full model: 
fm3 = raster ('./CH_results/chamr_pca_tot_g.asc')

#Raster stacks for calibration and projection 
pca_amr2
rp_pca_amr2 = rasterToPoints(pca_amr2)

#Total volumne obtained wih previous codes 
chpca_amr_tot_vol = 15.32279 #from results table 1, main text

#Total pixels obtained wih previous codes 
chpca_amr_tot_pxls = 98754 #from results table 1, main text

#Database with occurrences and climatic values 
jk_hvpca_amr = cbind (oro_thin[,1:4],extract(pca_amr2, oro_thin[,3:4]))

eval_in_e_jk3 = list () #collecting bins in E
eval_in_e_jk3_df = NULL #dataframe of results in E
eval_in_g_jk3 = list () #collecting bins in G 
eval_in_g_jk3_df = NULL #dataframe of results in G

for (j in 1:length(jk_hvpca_amr[,1])){
  #subsetting database in n-1 bins, Jackknife approach: 
  eval_pt = jk_hvpca_amr [j,]
  #dim(eval_pt)
  cal_pts = jk_hvpca_amr [-j,]
  #dim(cal_pts)
  
  ee = convexhull3d(data = cal_pts[,5:7], 
                    test_points= eval_pt[,5:7], 
                    plot = F, 
                    background = rp_pca_amr2[,3:5])
  
  ##absent
  #points3d(eval_pt[,4:6][ee[[2]]$in_convex==0,],
  #        col="red",size=5)
  ##present
  #points3d(eval_pt[,4:6][ee[[2]]$in_convex==1,],
  #        col="red",size=5)
  
  perf = ee$testing[,2] #point predicted/not predicted
  vlm = ee[[3]][[1]] #volume
  
  eval_in_e_jk3 [[length(eval_in_e_jk3)+1]] = perf
  
  res_rw = cbind (repl = j, 
                  ev_long = eval_pt[,2], 
                  ev_lat = eval_pt[,3],
                  perform = perf,
                  vol = vlm, 
                  proportion = vlm/chpca_amr_tot_vol)
  
  eval_in_e_jk3_df = rbind (eval_in_e_jk3_df, res_rw)
  
  #in geography: 
  ee_ex = data.frame(rp_pca_amr2[ee[[4]]$in_convex_bc==0,])
  ee_in = data.frame(rp_pca_amr2[ee[[4]]$in_convex_bc==1,])
  ee_ex$ind = rep (0, length(ee_ex[,1]))
  ee_in$ind = rep (1, length(ee_in[,1]))
  ee_env = rbind (ee_ex, ee_in)
  
  #write.csv (ee_env, paste ('./JACKKNIFE_CH/ch_pca_amr_', j, '.csv', sep = ''), row.names = F)
  
  st_pxs = dim(ee_in)[1] #environmental rows 
  
  #geographical projection and raster: 
  amr_ras = rasterize (ee_env[,1:2], pca_amr2[[1]], field = ee_env[,6])
  
  ##writing the raster files (svm models):
  #writeRaster(amr_ras, filename = paste ('./JACKKNIFE_CH/repch_pca_amr_', j, sep = ''), format = 'ascii',
  #            bylayer = T, suffix = paste ('repch_pca_amr_', j, sep = ''), overwrite = T,
  #            NAflag = -9999)
  
  # plot(amr_ras)
  # points(eval_pt[,2:3], col ='red')
  
  #evaluation
  tst2 = extract (amr_ras, eval_pt[,3:4]) #extracting raster value, if it is inside = 1, otherwise = 0, the average will be the evaluation statistic
  eval_in_g_jk3[[length(eval_in_g_jk3)+1]] = tst2 #collecting if point was predicted or not
  
  res_geo = cbind (repl = j, 
                   ev_long = eval_pt[,3], 
                   ev_lat = eval_pt[,4], 
                   perform = tst2, 
                   suit_pxls = st_pxs, 
                   proportion = st_pxs/chpca_amr_tot_pxls)
  
  eval_in_g_jk3_df = rbind (eval_in_g_jk3_df, res_geo)
}


#checking e vs g: 
#frequencies
par(mfrow = c(2,1))
plot (eval_in_e_jk3_df[,4], type = 'h')
points (eval_in_e_jk3_df[,6], col = 'blue', pch = 15)

plot (eval_in_g_jk3_df[,4], type = 'h')
points (eval_in_g_jk3_df[,6], col = 'blue', pch = 15)

#dividing the database per contribution occs: 
test3 = data.frame(eval_in_g_jk3_df)
class(test3)

vc3 = list()

for (jj in test3[,6]){
  if (jj>=1){
    vc3[[length(vc3)+1]] = 1
  }else{
    if (jj<1 & jj>= 0.90){
      vc3[[length(vc3)+1]] = 2
    }else {
      if (jj<0.90){
        vc3[[length(vc3)+1]] = 3
      }
    }
  }
}

test3$cat = unlist(vc3)
test_ord3 = test3[order(test3[,6]),]

#Plotting the contribution of occurrences and mapping: 
par(mfrow = c(1,2))
plot (test_ord3[,6], col = c('black', 'blue','red')[test_ord3[,7]], pch = 16)
plot (fm3)
plot (wrld_simpl, add = T)
points (test_ord3[,2:3], col= c('black', 'blue','red')[test_ord3[,7]], pch = 16, cex = 0.5)

#correlation volume performance/pixels performance
plot (eval_in_e_jk3_df[,5], eval_in_e_jk3_df[,4])
plot (eval_in_g_jk3_df[,5], eval_in_g_jk3_df[,4])

core_e_jk3 = cor.test(eval_in_e_jk3_df[,5], eval_in_e_jk3_df[,4], method = 'pearson') #statistically significant
core_g_jk3 = cor.test(eval_in_g_jk3_df[,5], eval_in_g_jk3_df[,4], method = 'pearson') #statistically significant

#performance stat: 
re_svm_e_jk3= c(mean = mean (unlist (eval_in_e_jk3)), 
                sd = sd (unlist (eval_in_e_jk3)), 
                vol = chpca_amr_tot_vol, 
                cor_perf_vol = core_e_jk1$estimate, 
                cor_pe = core_e_jk1$p.value) #performance for niche: 

re_svm_g_jk3 = c(mean = mean (unlist (eval_in_e_jk3)), 
                 sd = sd (unlist (eval_in_g_jk3)), 
                 pxls = chpca_amr_tot_pxls, 
                 cor_perf_vol = core_g_jk1$estimate, 
                 cor_pe = core_g_jk1$p.value) #performance for niche: 

#re_svm_past_g = mean (unlist (eval_in_g_past))
#re_svm_past_e = mean (unlist (eval_in_e_past))

write.csv (re_svm_g_jk3, './JACKKNIFE_CH/stat_svm_g_jk3.csv', row.names = T)
write.csv (re_svm_e_jk3, './JACKKNIFE_CH/stat_svm_e_jk3.csv', row.names = T)

write.csv (eval_in_e_jk3_df, './JACKKNIFE_CH/svm_e_jk3_df.csv', row.names = T)
write.csv (eval_in_g_jk3_df, './JACKKNIFE_CH/svm_g_jk3_df.csv', row.names = T)


#OCCURRENCE CONTRIBUTION: JACKKNIFE IN HYPERVOLUMES----------------------------------

#Example with Hypervolume definitive model: 

#full model
fm1 = raster('./HV_results/hvamr_pca_tot_g.asc')

#envs: 
pca_amr2

#tot_vol
pca_amr_tot_vol = 20.15296 #from results table

#tot_pxls
pca_amr_tot_pxls = 141431 #from results table

#database for the current excercise 
jk_amr = cbind (oro_thin[,1:4],extract(pca_amr2, oro_thin[,3:4]))

eval_in_e_jk1 = list () #collecting bins in E
eval_in_e_jk1_df = NULL #dataframe of results in E
eval_in_g_jk1 = list () #collecting bins in G 
eval_in_g_jk1_df = NULL #dataframe of results in G


for (j in 1:length(jk_amr[,1])){ #jackknife will loop for the total number of occurrences
  #subsetting database in n-1 for calibration and evaluation: 
  eval_pt = jk_amr [j,]
  #dim(eval_pt)
  cal_pts = jk_amr [-j,]
  #dim(cal_pts)
  
  svm_e = hypervolume_svm(data = cal_pts[,5:7]) #calibrating the model, only three dimensions
  
  #Luis approach, for details check index: APPROACH LUIS OSORIO
  #dimen1: 
  polygons1 <- concaveman(as.matrix(svm_e@RandomPoints[,1:2]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in1 <- point.in.polygon(eval_pt[,5],eval_pt[,6],polygons1[,1],polygons1[,2]) 
  
  #dimen2: 
  polygons2 <- concaveman(as.matrix(svm_e@RandomPoints[,2:3]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in2 <- point.in.polygon(eval_pt[,6],eval_pt[,7],polygons2[,1],polygons2[,2]) 
  
  #dimen3: 
  polygons3 <- concaveman(as.matrix(svm_e@RandomPoints[,c(1,3)]),concavity = 2) #creates a concave matrix, higher parameters less overfitted
  ids_in3 <- point.in.polygon(eval_pt[,5],eval_pt[,7],polygons3[,1],polygons3[,2]) 
  
  # plot (polygons1, xlim = c(0, 10), ylim = c(-5, 5)) 
  # lines(polygons1)
  # points (eval_pt[,c(4,5)], col = 'red', cex= 0.5, pch = 16)
  
  red1_res = mean(c(mean (ids_in1), mean (ids_in2), mean (ids_in3)))
  
  eval_in_e_jk1[[length(eval_in_e_jk1)+1]] = red1_res #collecting when the algorithm is collecting points in E
  
  #dataframe of results
  df_e = cbind (replicate = j, #coordinate excluded from the analysis
                ev_long = eval_pt[,3],
                ev_lat = eval_pt[,4],
                perform = red1_res, #addding if point excluded was obtained
                #past = mn2, #adding values from the points of the past 
                vol = svm_e@Volume, #collecting volume of replicate niche
                proportion = svm_e@Volume/pca_amr_tot_vol) #calculating proportion vs volume all points
  
  eval_in_e_jk1_df = rbind (eval_in_e_jk1_df,df_e) #adding info to the dataframe
  
  ##writing the random points around the occurrences (the environmental buffer)
  #write.csv (svm_e@RandomPoints, paste ('./JACKKNIFE_HV/env_amr_', j, '.csv', sep = ''), row.names = F)
  
  #projection to G
  svm_g = hypervolume_project (hv = svm_e, rasters = pca_amr2, type = "inclusion") #using the first three PCA rasters! 
  
  ##writing the raster files (svm models):
  #writeRaster(svm_g, filename = paste ('./JACKKNIFE_HV/rep_amr_', j, sep = ''), format = 'ascii',
  #            bylayer = T, suffix = paste ('rep_amr_', j, sep = ''), overwrite = T,
  #            NAflag = -9999)
  
  #evaluation in G jackknife present points
  r_t_p = rasterToPoints(svm_g) #matrix counting 0s and 1s, 
  prop_pxs = nrow(r_t_p[which(r_t_p[,3] ==1),])/pca_amr_tot_pxls #total pixels calculated as suitable with n-1 niche divided total pixels of niche using all points
  
  tst = mean(extract (svm_g, eval_pt[,3:4])) #extracting raster value, if it is inside = 1, otherwise = 0, the average will be the jackknife evaluation statistic
  eval_in_g_jk1[[length(eval_in_g_jk1)+1]] = tst #collecting if point was predicted or not
  
  # #eval in G with points of the past 
  # m1 = mean (extract (svm_g, eval_p[,2:3])) #0 = absent, 1 = present 
  # eval_in_g_past[[length(eval_in_g_past)+1]] = m1
  # 
  #results dataframe
  df_g = cbind (replicate = j, #coordinate excluded from the analysis
                ev_long = eval_pt[,3],
                ev_lat = eval_pt[,4],
                perform = tst, #if point excluded was correctly predicted by the n-1 model
                #past = m1, #adding if the points of the past were obtained in each iteration
                tot_suit_pxs = nrow(r_t_p[which(r_t_p[,3] ==1),]), #number of suitable pixels
                prop_pxs = prop_pxs) #proportion of suitable pixels as compared with the model with all the points
  
  eval_in_g_jk1_df = rbind (eval_in_g_jk1_df,df_g)
}


#checking e vs g: 
#frequencies
par(mfrow = c(2,1))
plot (eval_in_e_jk1_df[,4], type = 'h')
points (eval_in_e_jk1_df[,6], col = 'blue', pch = 15)

plot (eval_in_g_jk1_df[,4], type = 'h')
points (eval_in_g_jk1_df[,6], col = 'blue', pch = 15)

#dividing the database per contribution occs: 
test = data.frame(eval_in_g_jk1_df)

#reading and visualization of this table in the next section... 

#correlation volume performance/pixels performance
plot (eval_in_e_jk1_df[,5], eval_in_e_jk1_df[,4])
plot (eval_in_g_jk1_df[,5], eval_in_g_jk1_df[,4])

core_e_jk1 = cor.test(eval_in_e_jk1_df[,5], eval_in_e_jk1_df[,4], method = 'pearson') #statistically significant
core_g_jk1 = cor.test(eval_in_g_jk1_df[,5], eval_in_g_jk1_df[,4], method = 'pearson') #statistically significant

#performance stat: 
re_svm_e_jk1= c(mean = mean (unlist (eval_in_e_jk1)), 
                sd = sd (unlist (eval_in_e_jk1)), 
                vol = pca_amr_tot_vol, 
                cor_perf_vol = core_e_jk1$estimate, 
                cor_pe = core_e_jk1$p.value) #performance for niche: 

re_svm_g_jk1 = c(mean = mean (unlist (eval_in_e_jk1)), 
                 sd = sd (unlist (eval_in_g_jk1)), 
                 pxls = pca_amr_tot_pxls, 
                 cor_perf_vol = core_g_jk1$estimate, 
                 cor_pe = core_g_jk1$p.value) #performance for niche: 

#re_svm_past_g = mean (unlist (eval_in_g_past))
#re_svm_past_e = mean (unlist (eval_in_e_past))

write.csv (re_svm_g_jk1, './JACKKNIFE_HV/stat_svm_g_jk1.csv', row.names = T)
write.csv (re_svm_e_jk1, './JACKKNIFE_HV/stat_svm_e_jk1.csv', row.names = T)

write.csv (eval_in_e_jk1_df, './JACKKNIFE_HV/svm_e_jk1_df.csv', row.names = T)
write.csv (eval_in_g_jk1_df, './JACKKNIFE_HV/svm_g_jk1_df.csv', row.names = T)


#plotting and visualization: 
test = read.csv ('./JACKKNIFE_HV/svm_g_jk1_df.csv', header = T)

class(test)

vc = list()

for (jj in test[,7]){
  if (jj>=1){
    vc[[length(vc)+1]] = 1
  }else{
    if (jj<1 & jj>= 0.90){
      vc[[length(vc)+1]] = 2
    }else {
      if (jj<0.90){
        vc[[length(vc)+1]] = 3
      }
    }
  }
}

unlist(vc)
test$cat = unlist(vc)
test2 = test[order(test[,7]),]

par(mfrow = c(1,2))
plot (test2[,7], col = c('black', 'blue','red')[test2[,8]], pch = 16, 
      ylab = 'Proportion of suitable pixels', xlab = 'Points',
      main = 'Occurrence contribution')
plot (fm1)
plot (wrld_simpl, add = T)
points (test[,3:4], col= c('black', 'blue','red')[test[,8]], pch = 16, cex = 1)


#POST MODELLING STEPS-------------------------------------------------

#ENVIRONMENTAL VALUES IN OUTBREAK LOCALITIES--------------------------

#DATA: 

#points with localities with decreasing more than 10% the geographial output: 
test2

#environmental predictors of the Americas: 
amr1

#APPLICATION: 
#transform in spatial point dataframe: 
test3 = SpatialPointsDataFrame(coords = test2[,3:4], 
                                data = test2,
                               proj4string = WGS84)


#creating buffer: 
dist1 = 100000 #100 km = 100,000 meters #distant buffer

buff_gen = buffer(test3, width = dist1, dissolve = F) 
plot (buff_gen)


#splitting buffers to obtain data: 
all1 = buff_gen[-c(1:2),] #all points except the identified points (<0.90)

north_occ = buff_gen[1,] #northern site 

sout_occ = buff_gen[2,] #southern site

#visualization
plot (all1)
plot (north_occ, col = 'blue', add =T)
plot (sout_occ, col = 'red', add =T)


#STATISTICS: 

#Northern site: 
t1 = crop (amr1, north_occ)
t1 = mask (t1, north_occ)
north_m = cellStats (t1, mean) #mean 
north_sd =  cellStats (t1, sd) #standard deviation
north_counts = dim(rasterToPoints(t1))[1] #counts
catn = rep ('north',2) #category 

north_res = data.frame (rbind (mean = north_m, sd = north_sd))

#fixing for MERRACLIM corrections: 
north_res$X5m_mean_00s_bio1 = north_res$X5m_mean_00s_bio1/10
north_res$X5m_mean_00s_bio7 = north_res$X5m_mean_00s_bio7/10
north_res$X5m_mean_00s_bio12 = north_res$X5m_mean_00s_bio12/100000

#adding category 
north_res$cat = catn


#Southern site: 
t2 = crop (amr1, sout_occ)
t2 = mask (t2, sout_occ)
south_m = cellStats (t2, mean) #mean 
south_sd =  cellStats (t2, sd) #shiva
south_counts = dim(rasterToPoints(t2))[1] #counts
cats = rep ('south',2) #category

south_res = data.frame (rbind (mean = south_m, sd = south_sd))

#fixing for MERRACLIM corrections: 
south_res$X5m_mean_00s_bio1 = south_res$X5m_mean_00s_bio1/10
south_res$X5m_mean_00s_bio7 = south_res$X5m_mean_00s_bio7/10
south_res$X5m_mean_00s_bio12 = south_res$X5m_mean_00s_bio12/100000

#addding category 
south_res$cat = cats



#Averages: 
t3 = crop (amr1, all1)
t3 = mask (t3, all1)
avgs_m = cellStats (t3, mean) #mean 
avgs_sd =  cellStats (t3, sd) #sd
avgs_counts = dim(rasterToPoints(t3))[1] #counts
cata = rep ('average',2) #category 

avgs_res = data.frame (rbind (mean = avgs_m, sd = avgs_sd))

#fixing for MERRACLIM corrections: 
avgs_res$X5m_mean_00s_bio1 = avgs_res$X5m_mean_00s_bio1/10
avgs_res$X5m_mean_00s_bio7 = avgs_res$X5m_mean_00s_bio7/10
avgs_res$X5m_mean_00s_bio12 = avgs_res$X5m_mean_00s_bio12/100000

#addng category
avgs_res$cat = cata


#writing values: 
write.csv (north_res, './env_values/north_res.csv', row.names = T)
write.csv (south_res, './env_values/south_res.csv', row.names = T)
write.csv (avgs_res, './env_values/avgs_res.csv', row.names = T)

#VISUALIZATION OF POINT PLOTS WITH GGPLOT2-----------------------------------
#plotting stuff including errors bars: 
#overall dataset: 
overdf = rbind (t(avgs_res), t(north_res), t(south_res))

#eliminate category: 
overdf = data.frame (overdf[-c(4,8,12),])

#vector of categories for the new dataframe: 
cat = c('Avg.', 'North','South')

#BIO1 dataframe: 
bio1= overdf[c(1,4,7),]
bio1$cat = cat
bio1$mean = as.numeric (as.character(bio1$mean))
bio1$sd = as.numeric (as.character(bio1$sd))

#BIO7 dataframe: 
bio7= overdf[c(3,6,9),]
bio7$cat = cat
bio7$mean = as.numeric (as.character(bio7$mean))
bio7$sd = as.numeric (as.character(bio7$sd))

#BIO12 dataframe: 
bio12= overdf[c(2,5,8),]
bio12$cat = cat
bio12$mean = as.numeric (as.character(bio12$mean))
bio12$sd = as.numeric (as.character(bio12$sd))


#ABOUT GGPLOT: 
#' you have to FOR SURE understand the syntax of the package: 
#' http://www.sthda.com/english/wiki/be-awesome-in-ggplot2-a-practical-guide-to-be-highly-effective-r-software-and-data-visualization#geom_errorbar-error-bars
#' https://www.geeksforgeeks.org/plot-mean-and-standard-deviation-using-ggplot2-in-r/
#' Other important points at the end of this section

#' Ggplot is super conversome if you 
#' ignore the syntax... thus here are details to get this plot together
#' 

#Mian plot, this is a dot plot so we define x and y, and also the categories: 
bio1_graph = ggplot(bio1, aes(x=cat, y=mean, color = cat)) + 
      geom_point()+
      #manually changing the colors and modifying the legend: 
      scale_color_manual(values = c('gold', 'darkred', 'brown1'), 
                         guide = guide_legend(title= 'Sites', 
                                              override.aes = list(size = 0.5)))+
      #creating the error bars, I am using the variables from the database 
      #notice with = 0 so no horizontal ticks appear in the error bars 
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = cat), width=0,
                    position=position_dodge(0.05))+
      #adding names to each axis and to the title 
      xlab('Occurrences')+
      ylab('Means +/- SD')+
      ggtitle('BIO1')+
      #eliminating the awful ggplot background and controling 
      #the size of the labels and fonts and the position of tht title
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            plot.title = element_text (size = 14, hjust = 0.5), 
            axis.title=element_text(size=14), 
            axis.text=element_text(size=11), 
            legend.title = element_text(size=14),
            legend.text = element_text(size=14))

#now the other plots: 

bio7_graph = ggplot(bio7, aes(x=cat, y=mean, color = cat)) + 
  geom_point()+
  #manually changing the colors and modifying the legend: 
  scale_color_manual(values = c('gold', 'darkred', 'brown1'), 
                     guide = guide_legend(title= 'Sites', 
                                          override.aes = list(size = 0.5)))+
  #creating the error bars, I am using the variables from the database 
  #notice with = 0 so no horizontal ticks appear in the error bars 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = cat), width=0,
                position=position_dodge(0.05))+
  #adding names to each axis and to the title 
  xlab('Occurrences')+
  ylab('Means +/- SD')+
  ggtitle('BIO7')+
  #eliminating the awful ggplot background and controling 
  #the size of the labels and fonts and the position of tht title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text (size = 14, hjust = 0.5), 
        axis.title=element_text(size=14), 
        axis.text=element_text(size=11), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))


bio12_graph = ggplot(bio12, aes(x=cat, y=mean, color = cat)) + 
  geom_point()+
  #manually changing the colors and modifying the legend: 
  scale_color_manual(values = c('gold', 'darkred', 'brown1'), 
                     guide = guide_legend(title= 'Sites', 
                                          override.aes = list(size = 0.5)))+
  #creating the error bars, I am using the variables from the database 
  #notice with = 0 so no horizontal ticks appear in the error bars 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = cat), width=0,
                position=position_dodge(0.05))+
  #adding names to each axis and to the title 
  xlab('Occurrences')+
  ylab('Means +/- SD')+
  ggtitle('BIO12')+
  #eliminating the awful ggplot background and controling 
  #the size of the labels and fonts and the position of tht title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text (size = 14, hjust = 0.5), 
        axis.title=element_text(size=14), 
        axis.text=element_text(size=11), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))
  

#EXPORTING PDFs: 

pdf (paste ('./env_values/bio1_graph', '.pdf', sep = ''), width = 10.5, height = 8, paper = 'letter')
bio1_graph
dev.off()


pdf (paste ('./env_values/bio7_graph', '.pdf', sep = ''), width = 10.5, height = 8, paper = 'letter')
bio7_graph
dev.off()

pdf (paste ('./env_values/bio12_graph', '.pdf', sep = ''), width = 10.5, height = 8, paper = 'letter')
bio12_graph
dev.off()

#' Legend can be controlled from the theme argument
#' https://www.statology.org/ggplot2-legend-size/
#' Eliminating background depends highly on theme again where you can control 
#' the size of different parts of the plot 
#' 


#RANDOMIZATION TEST - EVI--------------------------------------

#DATA PREPARATION: 

#' #rasters for 2019 and 2003 obtained through Google Earth: 

evi2019 = raster ('./vegetation_cover/EVI/mean2019_evi.tif')
evi2003 = raster ('./vegetation_cover/EVI/mean2003_evi.tif')

#obtaining the difference: 
dif1= evi2019 - evi2003

#resample to experiment resolution:
#' The data obtained is a 500 meters resolution, we need to change it to the 
#' 5 arc/sec 
#' resampling is done using one of the environmental predictors 
#' because we want the resolution of the variables used in the 
#' modeling to develop this experiment. 

dif2= resample (dif1, amr1[[1]], method = 'bilinear') #raster ready for cropping: 


#reading final model 
hv1 = raster ('./DEF_MODEL/both_predictors.asc')

#visualization: 
par (mfrow = c(1,2))
plot (hv1)
plot (dif1)


#raster cropping rasters function: 
#' https://github.com/daromero-88/Gadgets-and-scripts-for-niche-modeling-/blob/main/mask_ras.R
#' Data needed: raster 1, cropping raster, projection
#' Here the WGS84 has the projection of from the wrld_spl build in the maptools package

r2 = mask_ras(dif2, hv1, projection_def = WGS84)
plot (r2)

#writing difference raster: 
writeRaster(r2, filename = './vegetation_cover_results/dif_evi_cropped', format = 'ascii', 
            bylayer = T, suffix = './vegetation_cover_results/dif_evi_cropped', 
            overwrite = T, NAflag = -9999)


#RANDOMIZATION TEST ITSELF: 

#obtaining observed statistic: 
#' Change statistic accordingly
#' You can select the code until 'CHANGE UP HERE' and replace
#' median for mean, or other statistic

evi_median = median(extract(r2, oro_thin[,3:4])) #sensitivity, number of points correctly predicted as present
draws = dim(oro_thin)[1]

#obtaining the values from the raster: 
ex1 = data.frame(rasterToPoints(r2))
#aa = is.na(ex1$layer) == TRUE #are any NA values? 
#ex1$layer[aa] #no, there are no NAs


#generating null models, n = 1000:

evi_null = NULL #dataframe to be filled 

#jj = 1 #debugging

for (jj in 1:1000){
  samp = sample(ex1$layer, draws, replace = TRUE)
  stat1 = median(samp)
  pre_df = cbind (rep = jj, stat = stat1) #joining a dataframe with the replicate and statistic
  evi_null = rbind (evi_null, pre_df) #writing the results
}


#xx = range(c(max(evi_null[,2]), min(evi_null[,2]), evi_median)) #in case x limits are needed

#dev.new()
pdf (paste ('./vegetation_cover_results/evi_median', '.pdf', sep = ''), width = 10.5, height = 8, paper = 'letter')
par(mar = c(10,4,3,3))
hist (evi_null[,2], main = 'Randomization test \n\ EVI-median', 
      xlab = 'Values of the median', col = 'grey')
abline (v = evi_median, col = 'red')
qn975 = quantile (evi_null[,2], probs = 0.975)
qn025 = quantile (evi_null[,2], probs = 0.025)
abline (v = qn975, lty = 3, col = 'black') 
abline (v = qn025, lty = 3, col = 'black')
dev.off()
#model is statistically significant, observation is uncontained in the null distribution 


#writing rasters and appropriate files: 

#null distribution: 
write.csv (evi_null, './vegetation_cover_results/evi_median_null.csv', row.names = F)

'CHANGE UP HERE'


#RANDOMIZATION TEST - NDVI--------------------------------------

#DATA PREPARATION: 

#' #rasters for 2019 and 2003 obtained through Google Earth: 

ndvi2019 = raster ('./vegetation_cover/ndvi/mean2019_ndvi.tif')
ndvi2003 = raster ('./vegetation_cover/ndvi/mean2003_ndvi.tif')

#obtaining the difference: 
dif_ndvi = ndvi2019 - ndvi2003

#resample to experiment resolution:
#' The data obtained is a 500 meters resolution, we need to change it to the 
#' 5 arc/sec 
#' resampling is done using one of the environmental predictors 
#' because we want the resolution of the variables used in the 
#' modeling to develop this experiment. 

dif_ndvi2 = resample (dif_ndvi, amr1[[1]], method = 'bilinear') #raster ready for cropping: 

#reading final model 
hv1 = raster ('./DEF_MODEL/both_predictors.asc')

#visualization: 
par (mfrow = c(1,2))
plot (hv1, main = 'model')
plot (dif_ndvi2, main = 'ndvi')


#raster cropping rasters function: 
#' https://github.com/daromero-88/Gadgets-and-scripts-for-niche-modeling-/blob/main/mask_ras.R
#' Data needed: raster 1, cropping raster, projection
#' Here the WGS84 has the projection of from the wrld_spl build in the maptools package

r3 = mask_ras(dif_ndvi2, hv1, projection_def = WGS84)
plot (r3)

#writing difference raster: 
writeRaster(r3, filename = './vegetation_cover_results/dif_ndvi_cropped', format = 'ascii', 
            bylayer = T, suffix = './vegetation_cover_results/dif_ndvi_cropped', 
            overwrite = T, NAflag = -9999)


#RANDOMIZATION TEST ITSELF: 

#obtaining observed statistic: 
#' Change statistic accordingly
#' You can select the code until 'CHANGE UP HERE' and replace
#' mean for mean, or other statistic

ndvi_mean = mean(extract(r3, oro_thin[,3:4])) #sensitivity, number of points correctly predicted as present
draws = dim(oro_thin)[1]

#obtaining the values from the raster: 
ex1 = data.frame(rasterToPoints(r3))
#aa = is.na(ex1$layer) == TRUE #are any NA values? 
#ex1$layer[aa] #no, there are no NAs


#generating null models, n = 1000:

ndvi_null = NULL #dataframe to be filled 

#jj = 1 #debugging

for (jj in 1:1000){
  samp = sample(ex1$layer, draws, replace = TRUE)
  stat1 = mean(samp)
  pre_df = cbind (rep = jj, stat = stat1) #joining a dataframe with the replicate and statistic
  ndvi_null = rbind (ndvi_null, pre_df) #writing the results
}


#xx = range(c(max(ndvi_null[,2]), min(ndvi_null[,2]), ndvi_mean)) #in case x limits are needed

#dev.new()
pdf (paste ('./vegetation_cover_results/ndvi_mean', '.pdf', sep = ''), width = 10.5, height = 8, paper = 'letter')
par(mar = c(10,4,3,3))
hist (ndvi_null[,2], main = 'Randomization test \n\ ndvi-mean', 
      xlab = 'Values of the mean', col = 'grey')
abline (v = ndvi_mean, col = 'red')
qn975 = quantile (ndvi_null[,2], probs = 0.975)
qn025 = quantile (ndvi_null[,2], probs = 0.025)
abline (v = qn975, lty = 3, col = 'black') 
abline (v = qn025, lty = 3, col = 'black')
dev.off()
#model is statistically significant if the observation is uncontained in the null distribution 


#writing rasters and appropriate files: 

#null distribution: 
write.csv (ndvi_null, './vegetation_cover_results/ndvi_mean_null.csv', row.names = F)

'CHANGE UP HERE'


#POPULATION IN RISK - FULL DISTRIBUTION-----------------------

#data: 
#selected model: 
hv1 = raster ('./DEF_MODEL/both_predictors.asc')

#population data
#obtained from: https://hub.worldpop.org/geodata/summary?id=24777

#pops = raster ('./pop_in_risk/ppp_2020_1km_Aggregated.tif') #**
#**please download from the corresponding link since such big files cannot be uploaded to GITHUB

#cropping to the studied area

pops_amr = crop(pops, amr)
pops_amr2 = mask (pops_amr, amr) #changing the name prevent an important error

plot (pops_amr2)

#matching the resolution of the model: 
pops_res = resample (pops_amr, amr1[[1]], method = 'bilinear')
plot (pops_res)

#cropping to the selected model: 
#mask_rar function --> https://github.com/daromero-88/Gadgets-and-scripts-for-niche-modeling-/blob/main/mask_ras.R

pop_def = mask_ras (pops_res, hv1, WGS84)
plot (pop_def)

writeRaster(pop_def, filename = './pop_in_risk/pop_cropped', format = 'ascii', 
            bylayer = T, suffix = './pop_in_risk/pop_cropped', 
            overwrite = T, NAflag = -9999)

#obtaining global statistics: 
pop_count = cellStats(pop_def, sum) #4944341


#' PROXY OF LOCAL INDCIDENCE: 
#' Numerator= population pixels in OROV province/state
#' Denomiantor= population pixels in province/state
#' Close to 1 --> more risk 
#' Close to 0 --> less risk 

#developing a chorolopeth map: 
provs_df #shapefile with information of states/provinces in the Americas 

#' checking the looping variable: 
#' This variable is selected because it has values across the entire 
#' dataframe because subet function does to accept NA as value descriptors

colnames(provs_df@data)
length (provs_df@data$woe_name) #669 provinces/states in the continent 

#how many NAs are in the variable name: 
aa = provs_df@data$adm1_code[is.na(provs_df@data$adm1_code)]
length(aa) #0

#' if needed, it is possible to eliminate the NAs in the shapefile instantly using 
#' the following approach: there where 4 NAs in the name variable thus, a total of 665 rows 
#' are recovered for the object prov3: 

#prov3 = provs_df[!is.na(provs_df$name),]

#I can create a data base with #665 entries and then add it to the shape file: 

sts2 = NULL

#w2 = provs_df@data$name[1] #debugging

for (w2 in provs_df@data$adm1_code){
  ss1 = subset (provs_df, adm1_code == w2)
  #numerator: 
  orov_pxl = crop (pop_def, ss1) #cropping from population in OROV distribution 
  orov_pxl = mask (orov_pxl, ss1)
  num1 = cellStats(orov_pxl, sum)
  #denominator: 
  tot_pxls = crop (pops_res, ss1) #cropping from overall population 
  tot_pxls = mask (tot_pxls, ss1)
  dem1 = cellStats(tot_pxls, sum)
  #local incidence: 
  li = num1/dem1
  ct = as.character(ss1@data$geonunit)
  nm = as.character(ss1@data$name)
  code = as.character (ss1@data$adm1_code)
  pre_df = cbind (cod = code, country = ct, second = nm, 
                  pop_risk = num1, pop_tot = dem1, incid= li)
  sts2 = rbind (sts2, pre_df)
}

#transforming to dataframe and fixing class 
sts2 = data.frame (sts2)
sts2$incid = as.numeric (as.character(sts2$inc)) #incidence as numeric 
sts2$incid [is.na (sts2$incid)] = 0 #replacing NaN with 0 (created whenever 0/0=NaN)

sts2$pop_risk = as.numeric (as.character(sts2$pop_risk)) #population in risk as numeric 
sts2$pop_tot = as.numeric (as.character(sts2$pop_tot)) #population total as numeric 

dim(sts2) #same length as state/province shape

#Both datasets end with the same codes 
tail (provs_df@data$adm1_code)
tail (sts2$cod)

#Adding data to shape: 
provs3_db = cbind (provs_df@data, sts2) #full database 
provs3 = SpatialPolygonsDataFrame(provs_df, provs3_db, match.ID = T) #transforming to shape file
class(provs3)

plot (provs3)
colnames(provs3@data) #only columns 60-65

#leaving only useful columns: 
provs4 = provs3[,c(60:65)]
colnames(provs4@data)

#writing datbase and shapefile: 

write.csv (provs4@data, './pop_in_risk/full_dist_database.csv', row.names = F) #database

writeOGR(provs4, layer= 'provs4', 
         dsn = './pop_in_risk/provs4', 
         driver = 'ESRI Shapefile') #shapefile 



#VISUALIZATION: 
#reading the object so we can establish an sf object: 

#defining a sf object for plotting: 
ch1 = st_read('./pop_in_risk/provs4', quiet = TRUE)

#defining a particular color palette 
con_pal1 = colorRampPalette(c('#e7e1ef', '#dd1c77'))

#checking the distribution of the variable I am looking for: 
hist (provs4@data$incid) #it is convenient to break for quantiles 


#plot: 
plot (ch1['incid'], 
      main = 'Local incidence', 
      breaks = c(0, 0.001, 0.1, 0.2, 0.5, 0.7, 1),  
      pal = con_pal1)



#POPULATION IN RISK - OVERLAPPING ONLY-----------------------

#data: 
#selected model: 
hv_overlap = raster ('./DEF_MODEL/overlap_only.asc')

#population data
#obtained from: https://hub.worldpop.org/geodata/summary?id=24777

pops = raster ('./pop_in_risk/ppp_2020_1km_Aggregated.tif')

#cropping to the studied area

pops_amr = crop(pops, amr)
pops_amr2 = mask (pops_amr, amr) #changing the name prevent an important error

plot (pops_amr2)

#matching the resolution of the model: 
pops_res = resample (pops_amr, amr1[[1]], method = 'bilinear')
plot (pops_res)

#cropping to the selected model: 
#mask_rar function --> https://github.com/daromero-88/Gadgets-and-scripts-for-niche-modeling-/blob/main/mask_ras.R

pop_overlap = mask_ras (pops_res, hv_overlap, WGS84)
plot (pop_overlap)

writeRaster(pop_overlap, filename = './pop_in_risk/pop_crop_overlap', format = 'ascii', 
            bylayer = T, suffix = './pop_in_risk/pop_crop_overlap', 
            overwrite = T, NAflag = -9999)

#obtaining global statistics: 
(pop_count = cellStats(pop_overlap, sum)) #2409182


#' PROXY OF LOCAL INDCIDENCE: 
#' Numerator= population pixels in OROV province/state
#' Denomiantor= population pixels in province/state
#' Close to 1 --> more risk 
#' Close to 0 --> less risk 

#developing a chorolopeth map: 
provs_df #shapefile with information of states/provinces in the Americas 

#' checking the looping variable: 
#' This variable is selected because it has values across the entire 
#' dataframe because subet function does to accept NA as value descriptors

colnames(provs_df@data)
length (provs_df@data$woe_name) #669 provinces/states in the continent 

#how many NAs are in the variable name: 
aa = provs_df@data$adm1_code[is.na(provs_df@data$adm1_code)]
length(aa) #0

#' if needed, it is possible to eliminate the NAs in the shapefile instantly using 
#' the following approach: there where 4 NAs in the name variable thus, a total of 665 rows 
#' are recovered for the object prov3: 

#prov3 = provs_df[!is.na(provs_df$name),]

#I can create a data base with #665 entries and then add it to the shape file: 

sts3 = NULL

#w2 = provs_df@data$name[1] #debugging

for (w2 in provs_df@data$adm1_code){
  ss1 = subset (provs_df, adm1_code == w2)
  #numerator: 
  orov_pxl = crop (pop_overlap, ss1) #cropping from population in OROV distribution 
  orov_pxl = mask (orov_pxl, ss1)
  num1 = cellStats(orov_pxl, sum) 
  #denominator: 
  tot_pxls = crop (pops_res, ss1) #cropping from overall population 
  tot_pxls = mask (tot_pxls, ss1)
  dem1 = cellStats(tot_pxls, sum)
  #local incidence: 
  li = num1/dem1 #local incidence numerator/denominator
  ct = as.character(ss1@data$geonunit) #country name 
  nm = as.character(ss1@data$name) #province/state name 
  code = as.character (ss1@data$adm1_code) #code for matching 
  pre_df = cbind (cod = code, country = ct, second = nm, 
                  pop_risk = num1, pop_tot = dem1, incid= li) 
  sts3 = rbind (sts3, pre_df)
}

#transforming to dataframe and fixing class 
sts3 = data.frame (sts3)
sts3$incid = as.numeric (as.character(sts3$inc)) #incidence as numeric 
sts3$incid [is.na (sts3$incid)] = 0 #replacing NaN with 0 (created whenever 0/0=NaN)

sts3$pop_risk = as.numeric (as.character(sts3$pop_risk)) #population in risk as numeric 
sts3$pop_tot = as.numeric (as.character(sts3$pop_tot)) #population total as numeric 

dim(sts3) #same length as state/province shape

#Both datasets end with the same codes 
tail (provs_df@data$adm1_code)
tail (sts3$cod)

#Adding data to shape: 
provs_overlap_dbs = cbind (provs_df@data, sts3) #full database 
provs_overlap = SpatialPolygonsDataFrame(provs_df, provs_overlap_dbs, match.ID = T) #transforming to shape file
class(provs_overlap)

plot (provs_overlap)
colnames(provs_overlap@data) #only columns 60-65

#leaving only useful columns: 
provs_overlap2 = provs_overlap[,c(60:65)]
colnames(provs_overlap2@data)

#writing datbase and shapefile: 

write.csv (provs_overlap2@data, './pop_in_risk/overlap_only_database.csv', row.names = F) #database

writeOGR(provs_overlap2, layer= 'provs_overlap2', 
         dsn = './pop_in_risk/provs_overlap2', 
         driver = 'ESRI Shapefile') #shapefile 


#VISUALIZATION: 
#reading the object so we can establish an sf object: 

#defining a sf object for plotting: 
ch2 = st_read('./pop_in_risk/provs_overlap2/', quiet = TRUE)

#defining a particular color palette 
con_pal1 = colorRampPalette(c('#e7e1ef', '#dd1c77'))

#checking the distribution of the variable I am looking for: 

hist (provs_overlap2@data$incid) #it is convenient to break for quantiles 


#visualization: 

plot (ch2['incid'], 
      main = 'Local incidence', 
      breaks = c(0, 0.001, 0.1, 0.2, 0.5, 0.7, 1),  
      pal = con_pal1)


#AND THIS IS IT :) 


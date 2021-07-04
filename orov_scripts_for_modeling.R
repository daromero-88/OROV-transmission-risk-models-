#####OROV TRANSMISSION RISK MODELING######

--------------------------------------------------------
  #'  Author: 
  #'  Daniel Romero-Alvarez
  #'  Luis Osorio-Olvera 
  #'  
  #' The following functions and scripts will allow the user to replicate 
  #' the models and simulations presented in the manuscript: 
  #' "Transmission Risk of Oropouche fever in the Americas"
  #' 
  
  --------------------------------------------------------
  
#LIBRARIES-----------------------------------------------

#Spatial analysis and niche models
library(raster)
library(sp)
library(rgeos)
library (rgdal)
library (maptools)
library(mapdata)
library (kuenm)
library (dismo)
library(spThin)
library (dbscan)

#Hypervolumes
#library (rgl)
library (Rcpp)
library (hypervolume)

#to use the check.in.range function, presence of values in one dimension
library (spatstat.utils)

#for functions provided by Luis Osorio-Olvera: 
library(rgl)
library(geometry)
library(ptinpoly)

#for evaluation in pair-wise dimensions: 
library(concaveman)
library(sp)

#graphics: 
library(ggplot2)


#DATA----------------------------------------------------


#FUNCTIONS-----------------------------------------------

#FUNCTION TO DEVELOP THREE DIMENSIONAL CONVEX HULLS--------------------------- 

####NEEDED REVIEW######

#'Provided by Luis Osorio-Olvera,
#'Modified by Daniel Romero-Alvarez,

convexhull3d <- function(data,test_points=NULL,plot=F,centroid=NULL, background = NULL){
  if(dim(data)[2]==3){
    # Codigo para graficar el convexhull 
    ch_data = convhulln(data, output.options = 'FA')  #recover indexes that define the convex hull, and with FA, it recovers the volume and other attributes
    ts.surf1 <- t(ch_data$hull) 
    
    convexhull_faces <- cbind(data[ts.surf1,1],
                              data[ts.surf1,2],
                              data[ts.surf1,3]) #obtaining the vertices of the hull in the parameter space
    
    #alphashp <- ashape3d(xyzmatrix(data),pert = T,
    #                     alpha = 130)
    
    #this ifelse portion is defined because the centroid can be defined in the function! 
    if(is.null(centroid)){
      centroide <- colSums(convexhull_faces)/(dim(convexhull_faces)[1])
    }else{
      centroide <- centroid
    }
    
    #start plotting if the function recommends:   
    if(plot){
      #convex_mesh <- as.mesh3d(alphashp)
      convex_mesh <- tmesh3d(
        vertices = t(data),
        indices = ts.surf1,
        homogeneous = FALSE)
      
      plot3d(data,size=5)
      plot3d(convexhull_faces,type="p", size=10,col="blue",add=T) #the vertices
      #plot3d(convexhull_faces,type="s",alpha = 0.3,size=0.6,col="blue",add=T) #the vertices, original by LOO
      #convex1 <-  rgl.triangles(convexhull_faces,
      #                          alpha=0.1,col=rgb(1,0,0.3))
      
      wire3d(convex_mesh, col = 4)
      # spheres3d(x=centroide[1], #theoretically this creates spheres but I am not sure
      #           y = centroide[2],
      #           z = centroide[3],
      #           col="blue",
      #           radius = 0.085,add=T)
    }
    
  }
  
  #Are evaluating testing points in the hull: 
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
  #Are background points in the hull: 
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
    
    results <- list(centroid=centroide,
                    testing = dconvex, 
                    vol = ch_data[3], 
                    backgroound = dconvex_bc)
    #the centroid of the 3d hull, the list of evaluating points, the volume of the hull, the list of background points
  }
  
  return(results)
}

#measuring Euclidean distance: 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) #we might eliminate that 'sum'


#SCRIPTS-------------------------------------------------


#HYPERVOLUME: OBSERVED MODEL, CALIBRATION LOOP, EVALUATION STATS, MODELS---------------

#climates in buffer: 
names (amr1)
amr2 = amr1[[c(1, 15, 4)]]
plot(amr2[[1]])
plot(amr2[[2]])
plot(amr2[[3]])

#original points
clusters7 = cbind (ID1 = 1:nrow (clusters), clusters) #adding ID for tracking points 

#extracting values considering rasters provided, raw data for the hypervolume function: 
clusters8 = data.frame(cbind (clusters7[,1:3], extract(amr2, clusters7[,2:3])))

#hypervolume (svm) with provided information --> values of each environmental layer at each point: 
svm_amr = hypervolume_svm(data = clusters8[,4:6])
plot (svm_amr)

#total volume: 
svm_amr_tvol = svm_amr@Volume #obtaining the volumne of the model 

#Projecting to geography: 
svm_amr_g = hypervolume_project (hv = svm_amr,rasters = amr2, type = "inclusion")
plot (svm_amr_g) #plotting projected model 

#total number of pixels: 
svm_amr_pxs = rasterToPoints(svm_amr_g) #transforming the projected model to points 
dim(svm_amr_pxs) #3438508 #total pixels 

svm_amr_tpxs = nrow(svm_amr_pxs[which (svm_amr_pxs[,3]==1),]) #colleting only those pixels with suitable value = 1

#writing results for the observed model, environmental dataframe:
#environmental values of the model
write.csv (svm_amr@RandomPoints, paste ('./hypervolume_all_points/jknife/hvamr_clim_tot_e','.csv', sep = ''), row.names = F)

#writing results for the observed model, raster file:
writeRaster(svm_amr_g, filename = paste ('./hypervolume_all_points/jknife/hvamr_clim_tot_g', sep = ''), format = 'ascii',
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
  cal70 = clusters8 [!(clusters8$ID1 %in% eval30$ID1),]
  #dim(cal70)
  
  cal70_val = extract (amr2, cal70[,2:3]) #values for model calibration
  eval30_val = extract (amr2, eval30[,2:3]) #values for model evaluation
  
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
  
  red1_res = mean(c(mean (ids_in1), mean (ids_in2), mean (ids_in3))) #obtaing the mean across dimensions 
  
  eval_in_e4[[length(eval_in_e4)+1]] = red1_res #collecting the mean as a list 
  
  #Dataframe of replicates data in the environmental space (E): 
  df_e = cbind (replicate = j, #replicate number 
                perform = red1_res, #collecting the evaluation statistic of the replicate
                vol = svm_e@Volume, #collecting volume of the replicate 
                proportion = svm_e@Volume/svm_amr_tvol) #obtaining the proportion of volume as a function of the full model 
  
  eval_in_e4_df = rbind (eval_in_e4_df,df_e) #adding info to the results dataframe 
  
  ##If environmental points for each replicate are needed, the following lines write that data:
  ##writing the random points around the occurrences (the environmental buffer)
  #write.csv (svm_e@RandomPoints, paste ('./hypervolume_all_points/jknife/env_sa_', j, '.csv', sep = ''), row.names = F)
  
  #Projection to geographic space: 
  svm_g = hypervolume_project (hv = svm_e, rasters = amr2, type = "inclusion") #using the corresponding rasters! 
  
  ##If needed, you can write the geographic rasters for each replicate
  #writeRaster(svm_g, filename = paste ('./hypervolume_all_points/jknife/rep_sa_', j, sep = ''), format = 'ascii',
  #            bylayer = T, suffix = paste ('rep_sa_', j, sep = ''), overwrite = T,
  #            NAflag = -9999)
  
  #Evaluation in the Geographic space (G):
  #Are the evaluation points contained in the raster obtained with the calibration points? 
  
  r_t_p = rasterToPoints(svm_g) #transforming raster to points 
  prop_pxs = nrow(r_t_p[which(r_t_p[,3]==1),])/svm_amr_tpxs #obtaining the proportion of suitable pixels per replicate as function of the full raster 
  
  tst = mean(extract (svm_g, eval30[,2:3])) #the mean of values across the evaluation points 
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
write.csv (re_svm_g4, './hypervolume_all_points/jknife/stat_svm_g4.csv', row.names = T)
write.csv (re_svm_e4, './hypervolume_all_points/jknife/stat_svm_e4.csv', row.names = T)

#Writing replicates dataframe
write.csv (eval_in_e4_df, './hypervolume_all_points/jknife/svm_e4_df.csv', row.names = T)
write.csv (eval_in_g4_df, './hypervolume_all_points/jknife/svm_g4_df.csv', row.names = T)

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
  writeRaster(a_all3[[jj]], filename = paste ('./hypervolume_all_points/jknife/hvamrclim_', nm5[jj], sep = ''), format = 'ascii',
              bylayer = T, suffix = paste ('hvamrclim_', j, sep = ''), overwrite = T,
              NAflag = -9999)
}


#CONVEX HULLS: OBSERVED MODEL, CALIBRATION LOOP, EVALUATION STATS, MODELS---------------

#environments: 
rp_pca_amr2= rasterToPoints(pca_amr2)

#Model using all the points: 
clusnich_amr = cbind (clusters2[,1:3],extract(pca_amr2, clusters2[,2:3])) #adding all the data together

chlfull_amr=  convexhull3d(data = clusnich_amr[,4:6], 
                           test_points= clusnich_amr[,4:6], #same points as in previous line 
                           plot = T, 
                           background = rp_pca_amr2[,3:5]) #add the background points 

#checking background points inside the hull: 

points3d(rp_pca_amr2[,3:5][chlfull_amr[[4]]$in_convex_bc==1,],
         col="red",size=5)

#checking background points outside the hull: 

points3d(rp_pca_amr2[,3:5][chlfull_amr[[4]]$in_convex_bc==0,],
         col="grey",size=1)

chl_amr_vol = chlfull_amr[[3]][[1]] #collecing the volume of the convex hull 

#Creating a database with all the data (important step, otherwise some data is going to be lost): 
nicha_amr_ex = data.frame(rp_pca_amr2[chlfull_amr[[4]]$in_convex_bc==0,]) #points outside the hull 
nicha_amr_in = data.frame(rp_pca_amr2[chlfull_amr[[4]]$in_convex_bc==1,]) #points inside the hull 

nicha_amr_ex$ind = rep (0, length(nicha_amr_ex[,1])) #adding value of 0 to the points outside the hull 
nicha_amr_in$ind = rep (1, length(nicha_amr_in[,1])) #adding value of 1 to the points inside the hull
nicha_amr_df = rbind (nicha_amr_ex, nicha_amr_in) #combining the dataframe 

#writing the dataframe complete dataframe: points and values
#writing results for the observed model, environmental dataframe:
write.csv (nicha_amr_df,'./hypervolume_all_points/jknife_nicheA2/chamr_pca_tot_e.csv', row.names = F)

#Projection to geography: 
amr_ras = rasterize (nicha_amr_df[,1:2], pca_amr2[[1]], field = nicha_amr_df[,6])
amr_spx = dim(nicha_amr_in)[1] #total number of suitable pixels

#writing results for the observed model, raster file:
writeRaster(amr_ras, filename = './hypervolume_all_points/jknife_nicheA2/chamr_pca_tot_g', format = 'ascii',
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
  cal70 = clusnich_amr [!(clusnich_amr$ID1 %in% eval30$ID1),]
  #dim(cal70)
  
  #Convex hull in environmental space: 
  
  ee = convexhull3d(data = cal70[,4:6], 
                    test_points= eval30[,4:6], #I should be able to not add this ####REVIEW THIS####
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
  tst2 = mean (extract (amr_ras, eval30[,2:3])) #collecting statistic in geography
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
write.csv (re_chl_e3, './hypervolume_all_points/jknife_nicheA2/stat_chl_e3.csv', row.names = T)
write.csv (re_chl_g3, './hypervolume_all_points/jknife_nicheA2/stat_chl_g3.csv', row.names = T)

#Writing replicates dataframes
write.csv (eval_ch_e3_df, './hypervolume_all_points/jknife_nicheA2/chl_e3_df.csv', row.names = T)
write.csv (eval_ch_g3_df, './hypervolume_all_points/jknife_nicheA2/chl_g3_df.csv', row.names = T)

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
  writeRaster(b_all3[[jj]], filename = paste ('./hypervolume_all_points/jknife_nicheA2/chamr_pca_', nm5[jj], sep = ''), format = 'ascii',
              bylayer = T, suffix = paste ('chamr_pca_', j, sep = ''), overwrite = T,
              NAflag = -9999)
}


#OCCURRENCE CONTRIBUTION: JACKKNIFE----------------------------------

#Select the raster for the full model: 
fm3 = raster ('./hypervolume_all_points/jknife_nicheA2/resuilts_amr_PCA/niche_amr_ras_g.asc')

#Raster stacks for calibration and projection 
pca_amr2
rp_pca_amr2

#Total volumne obtained wih previous codes 
chpca_amr_tot_vol = 15.32279 #from results table 1, main text

#Total pixels obtained wih previous codes 
chpca_amr_tot_pxls = 98754 #from results table 1, main text

#Database with occurrences and climatic values 
jk_hvpca_amr = cbind (clusters2[,1:3],extract(pca_amr2, clusters2[,2:3]))

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
  
  ee = convexhull3d(data = cal_pts[,4:6], 
                    test_points= eval_pt[,4:6], 
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
  
  write.csv (ee_env, paste ('./hypervolume_all_points/jknife_nicheA2/ch_pca_amr_', j, '.csv', sep = ''), row.names = F)
  
  st_pxs = dim(ee_in)[1] #environmental rows 
  
  #geographical projection and raster: 
  amr_ras = rasterize (ee_env[,1:2], pca_amr2[[1]], field = ee_env[,6])
  
  ##writing the raster files (svm models):
  writeRaster(amr_ras, filename = paste ('./hypervolume_all_points/jknife_nicheA2/repch_pca_amr_', j, sep = ''), format = 'ascii',
              bylayer = T, suffix = paste ('repch_pca_amr_', j, sep = ''), overwrite = T,
              NAflag = -9999)
  
  # plot(amr_ras)
  # points(eval_pt[,2:3], col ='red')
  
  #evaluation
  tst2 = extract (amr_ras, eval_pt[,2:3]) #extracting raster value, if it is inside = 1, otherwise = 0, the average will be the evaluation statistic
  eval_in_g_jk3[[length(eval_in_g_jk3)+1]] = tst2 #collecting if point was predicted or not
  
  res_geo = cbind (repl = j, 
                   ev_long = eval_pt[,2], 
                   ev_lat = eval_pt[,3], 
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

par(mfrow = c(1,2))
plot (test_ord3[,6], col = c('black', 'blue','red')[test_ord3[,7]], pch = 16)
plot (fm3)
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

write.csv (re_svm_g_jk3, './hypervolume_all_points/jknife_nicheA2/stat_svm_g_jk3.csv', row.names = T)
write.csv (re_svm_e_jk3, './hypervolume_all_points/jknife_nicheA2/stat_svm_e_jk3.csv', row.names = T)

write.csv (eval_in_e_jk3_df, './hypervolume_all_points/jknife_nicheA2/svm_e_jk3_df.csv', row.names = T)
write.csv (eval_in_g_jk3_df, './hypervolume_all_points/jknife_nicheA2/svm_g_jk3_df.csv', row.names = T)





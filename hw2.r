#HW2 JINGSHAN HUANG
#jhuang456@wisc.edu
rm(list=ls())
library('FITSio')
########data processing method######################################
zhishu = function(y,a){
  n = length(y) 
  
  y1 = numeric(n)
  y1[1] = y[1]
  for(i in 2:n){                
    y1[i] = a*y[i]+(1-a)*y1[i-1]
  }
  return(y1)
}

chafen = function(y){
  n = length(y)
  target = c(0,y[2:n] - y[1:(n-1)])
  return(target)
}

zengjian = function(diff_y){
  return(sign(diff_y))
}

bogu_omit_close = function(flux_sign,flux,closed = 5){
  n = length(flux_sign)
  mean_value = quantile(flux,0.3)
  last = -closed
  is_bogu = vector(mode = 'logical',length =n)
  for(i in (closed+1):(n-1)){
    if(!any(is_bogu[(i-closed):(i-1)]) &flux[i]<mean_value &flux_sign[i]==-1 & flux_sign[i+1]==1){
      is_bogu[i] = TRUE
    }
  }
  return(is_bogu)
}

process_mask = function(mask){
  n = length(mask)
  target = numeric(n)
  target[1] = mask[1] -1
  for(i in 2:n){
    j=i
    while(mask[j]-mask[j-1]==1){
      j = j-1
      if(j==1){
        break
      }
    }
    target[i] = mask[j] - 1
  }
  return(target)
}

process_data = function(dataset,islyman = 0){
  if(islyman == 1){
    dataset['zhishu'] = zhishu(dataset$FLUX,0.05)
  }else{
    #dataset = dataset[dataset$and_mask==0,]
    if(any(dataset$and_mask!=0)){
      #cat('yes,i have processed')
      mask = which(dataset$and_mask!=0)
      pmask = process_mask(mask)
      dataset$flux[mask] = dataset$flux[pmask]
    }
    n = length(dataset$flux)
    flux_sorted = sort(dataset$flux)
    flux_max5 = flux_sorted[n-4]
    flux_min5 = max(flux_sorted[5],-2)
    dataset$flux[dataset$flux>flux_max5] = flux_max5
    dataset$flux[dataset$flux<flux_min5] = flux_min5
    
    dataset['zhishu'] = zhishu(dataset$flux,0.05)
    #dataset['zhishu2'] = zhishu(dataset$flux,0.05)
  }
  dataset['diff'] = chafen(dataset$zhishu)
  dataset['sign_diff'] = zengjian(dataset$diff)
  dataset['isbogu'] = bogu_omit_close(dataset$sign_diff,dataset$zhishu,closed = 10)
  return(dataset)
}

lyman = readFrameFromFITS('cB58_Lyman_break.fit')
len_lyman = length(lyman$FLUX)
lyman = process_data(lyman,islyman = 1)

####for loop###############################################################
loop_function = function(dataset,corr_method = 'pearson'){
  dataset = process_data(dataset)
  len = length(dataset$flux)
  bogu = which(dataset$isbogu)
  loop_points = bogu[bogu>245 & bogu<(len-1934)]
  max_corr = 0
  max_point = 0
  for(point in loop_points){
    selected_data = as.numeric(dataset[(point - 245):(point + 1935),]$zhishu)
    current_corr = cor(selected_data,as.numeric(lyman$zhishu),method = corr_method)
    #cat(current_corr,end=' | ')
    if(current_corr > max_corr){
      max_point = point-245
      max_corr = current_corr
    }
  }
  return(c(max_corr,max_point))
}


path = './data'
filenames = dir(path)
score_list = numeric(100)
point_list = numeric(100)
i = 1
for(file in filenames){
  tryCatch({
    dataset = readFrameFromFITS(paste('data/',file,sep=''));
    dataset = process_data(dataset);
    result = loop_function(dataset);
    score = result[1];
    start_point = as.integer(result[2]);
    cat(sprintf('%3d : |score: %8.7f |start_point: %4d |filename: %s',i,score,start_point,file),end = '\n');
    score_list[i] = score;
    point_list[i] = start_point
  },error = function(e){cat(i,': error! filename is:',file,end = '\n')})
  
  i = i+1
}

df_result = data.frame(distance =score_list,spectrumID = filenames,i = point_list)
df_result = df_result[order(df_result$distance),]
write.csv(df_result,file = 'hw2.csv',row.names = FALSE)

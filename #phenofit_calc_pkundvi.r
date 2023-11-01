#PKU NDVI3g 物候提取

#参考范围用于crop
refer<-rast(res=1/12,xmin=-180,xmax=180,ymin=0,ymax=90,crs="EPSG:4326")

#确认变量与文件夹
VIs<-"NDVI"
data_folder<-'E:/huahao/pc_data/Rproject_PKU_NDVI/PKU_NDVI_NH_v2'

#确认年份
for(year_str in 1982:2020){
  #year_str<-1982
  
  #读取该年植被指数rast
  VI_rast <- list.files(paste0(data_folder,"/",year_str), pattern = ".tif$", full.names = T) %>% rast() %>% terra::crop(refer)
  raw_rast <- VI_rast

  phenofit_hh_cycles<-function(temp){
    #temp<-raw_rast[1218560]
    
    #phenofit全局参数
    nptperyear     <- 24
    wmin           <- 0.2
    TRS            <- 0.5 #c(0.1, 0.2, 0.5)
    result_names <- c("TRS5.sos1","TRS5.sos2","TRS5.eos1","TRS5.eos2",
                      "DER.sos1","DER.sos2","DER.pos1","DER.pos2","DER.eos1","DER.eos2",
                      "Greenup1","Greenup2","Maturity1","Maturity2","Senescence1","Senescence2","Dormancy1","Dormancy2")
    NA_result <- rep(NA,length(result_names))
    names(NA_result) <- result_names
    
    #
    GIMMS_DOY <- seq(1,360,15)
    raw_image_date <- data.frame(date=base::as.Date(as.integer(GIMMS_DOY), origin = "1982-12-31"))
    #
    image_date <- data.frame(date=base::as.Date(as.integer(GIMMS_DOY), origin = "1982-12-31"))
    
    #
    temp<-as.numeric(temp)#list向量化
    temp[is.infinite(temp)|is.nan(temp)] <- NA#把Inf更为NA
    
    #
    QC_vec <- rep(0, length(temp))
    QC_vec[temp == 0] <- 2
    #
    temp <- c(temp,QC_vec)
    
    #df
    df<- as.data.frame(matrix(temp,ncol = 2,byrow = F))
    names(df)<-(c("y","QC"))
    df<-base::merge(raw_image_date,cbind("date"=image_date,df),by="date",all.x=T)
    
    nrow1<-base::nrow(df[stats::complete.cases(df),])
    
    #一年中的有效观测数大于70%的观测总数
    if(nrow1 < (nptperyear*0.7)) return(NA_result)
    
    #插值
    df$y<-imputeTS::na_interpolation(df$y) 
    df<-cbind(df,base::as.data.frame(phenofit::qc_NDVI3g(df$QC)))
    df$t <- df$date+7
    #
    d <- data.table::data.table(base::subset(df,select=c("date","t","y","QC_flag","w")))
    
    #
    #d<-add_HeadTail(d, south=F, nptperyear = 23) 
    INPUT <- phenofit::check_input(d$t, d$y, d$w,QC_flag = d$QC_flag,
                                   nptperyear = nptperyear, south = FALSE,
                                   maxgap = base::round(nptperyear / 4), wmin = 0.2, na.rm=F,
                                   wsnow=0.8,ymin=0.08,alpha=0.02)#set ymin=0.08 for NDVI, ymin=0.05 for EVI, ymin=0.5 gC m-2 s-1 for GPP
    #plot_input(INPUT)
    
    #试验
    re1 = base::tryCatch({
      #
      brks <- phenofit::season_mov(INPUT,
                                   list(rFUN = "smooth_wSG", wFUN = "wTSM",#"smooth_wSG"smooth_wHANTS"smooth_wWHIT"
                                        #lambda = NULL,
                                        frame = floor(nptperyear/7) * 2 + 1, 
                                        maxExtendMonth = 3,
                                        wmin = wmin, r_min = 0.05,
                                        verbose=F
                                   ))
      
    }, warning = function(w) {
      "w"
    }, error = function(e) {
      "e"
    })
    if (any(re1 %in% c("w", "e"))) return(NA_result)
    
    #
    #lambda <- phenofit::init_lambda(INPUT$y)#只有16天的数据适用，其他情况应在season_mov中指定lambda = NULL，以使用V-curve理论寻找最优参数 
    brks <- phenofit::season_mov(INPUT,
                                 list(rFUN = "smooth_wSG", wFUN = "wTSM",#"smooth_wSG"smooth_wHANTS"smooth_wWHIT"
                                      #lambda = NULL,
                                      #frame = floor(nptperyear/7) * 2 + 1, 
                                      #maxExtendMonth = 3,
                                      wmin = wmin, r_min = 0.05,
                                      verbose=F
                                 ))
    # plot_season(INPUT, brks)
    
    if(all(is.null(brks))) return(NA_result)
    #brks结果存在
    if(is.null(brks$dt)) return(NA_result)#
    if(nrow(brks$dt)<=0) return(NA_result)#
    
    #
    fit  <- phenofit::curvefits(INPUT, brks,
                                options = list(
                                  methods = "Elmore",#c("AG", "Zhang", "Beck", "Elmore"), #,"klos",, 'Gu'
                                  wFUN = "wTSM",iters = 2,wmin = 0.2
                                ))
    
    ## check the curve fitting parameters
    # l_param <- get_param(fit)
    # dfit <- get_fitting(fit)
    # g<-plot_curvefits(fit, brks);grid::grid.newpage(); grid::grid.draw(g)
    
    ##Extract phenology
    #试验
    re1 = base::tryCatch({
      l_pheno <- phenofit::get_pheno(fit, TRS = TRS,IsPlot = FALSE) # %>% map(~melt_list(., "meth"))
    }, warning = function(w) {
      "w"
    }, error = function(e) {
      "e"
    })
    if (any(re1 %in% c("w", "e"))) return(NA_result)
    
    l_pheno <- phenofit::get_pheno(fit, TRS = TRS,IsPlot = FALSE) # %>% map(~melt_list(., "meth"))
    
    result_df<-base::as.data.frame(l_pheno$doy$Elmore)[,c("flag","TRS5.sos","TRS5.eos","DER.sos", "DER.pos", "DER.eos", "Greenup", "Maturity","Senescence", "Dormancy")]
    #year_str<-as.integer(strsplit(result_df$flag[1],"_")[[1]][1])
    result_df<-base::merge(data.frame(flag=paste0("1983_",c(1,2))),result_df,by="flag",all.x=T)
    
    return(base::unlist(result_df[,-1]))
    #c("TRS5.sos1","TRS5.sos2","TRS5.eos1","TRS5.eos2","DER.sos1","DER.sos2","DER.pos1","DER.pos2","DER.eos1","DER.eos2","Greenup1","Greenup2","Maturity1","Maturity2","Senescence1","Senescence2","Dormancy1","Dormancy2")
    
  }
  
  Sys.time()
  result_rast <- terra::app(raw_rast,phenofit_hh_cycles,cores=18)
  #result_rast2<-terra::app(raw_rast2,phenofit_hh_cycles,cores=12)
  Sys.time()
  #names(result_rast)<-c("TRS5.sos1","TRS5.sos2","TRS5.eos1","TRS5.eos2","DER.sos1","DER.sos2","DER.pos1","DER.pos2","DER.eos1","DER.eos2","Greenup1","Greenup2","Maturity1","Maturity2","Senescence1","Senescence2","Dormancy1","Dormancy2")
  
  result_folder<-paste0("PKU_NDVI_phenology_NH_v2/",year_str)
  dir.create(result_folder)
  writeRaster(result_rast,filename = paste0(result_folder,"/",result_rast@ptr$names,".tif"),overwrite=T)
  
  print(paste0(year_str," is done! - ",Sys.time()))
  gc()
}

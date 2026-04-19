begin


;;;;;;;;;;;;;;;;;; Data file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
psst_1co2 = "1xCO2_cam5_PSST_plev_file.nc"
som_1co2 = "1xCO2_cam5_SOM_plev_file.nc"

psst_2co2 = "2xCO2_cam5_PSST_plev_file.nc"
som_2co2 = "2xCO2_cam5_SOM_plev_file.nc"

psst_20BC = "20xBC_cam5_PSST_plev_file.nc"
som_20BC = "20xBC_cam5_SOM_plev_file.nc"

fp1 = addfile(psst_1co2,"r")
fs1 = addfile(som_1co2,"r")

fp2 = addfile(psst_2co2,"r")
fs2 = addfile ( som_2co2,"r")

fp3 = addfile(psst_20BC,"r")
fs3 = addfile ( som_20BC,"r")

;;;;;;;;;;;;;;;;;; Feedback Data file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f1 = addfile("2xCO2_cam5_feedback_values.nc","r")
f2 = addfile("20xBC_cam5_feedback_values.nc","r")

lat = f1->lat
lon = f1->lon
lev = f1->lev_p
lev@units = "mb"
lat@units = "degree_north"

;********************************** Gaussian weights for Area averages

  jlat  = dimsizes( lat )
  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  dlon   = abs(lon(2)-lon(1))*rr
  dx    = dlon*cos(lat*rad)
  dy     = new ( jlat, typeof(dx))
  dy(0)  = abs(lat(2)-lat(1))*rr
  dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5
  dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr
  area   = dx*dy
  clat   = cos(lat*rad)


ERF_2xCO2 								= f1->ERF
ERF_20xBC 								= f2->ERF
ERF_diff								= ERF_20xBC - ERF_2xCO2

LW_allsky_2xCO2							= f1->LW_all_sky_psst
SW_allsky_2xCO2							= f1->SW_all_sky_psst
LW_clrsky_2xCO2							= f1->LW_clear_sky_psst
SW_clrsky_2xCO2							= f1->SW_clear_sky_psst

LW_allsky_20xBC							= f2->LW_all_sky_psst
SW_allsky_20xBC							= f2->SW_all_sky_psst
LW_clrsky_20xBC							= f2->LW_clear_sky_psst
SW_clrsky_20xBC							= f2->SW_clear_sky_psst

Surf_T_tot_20xBC						= f2->TS
Surf_T_fst_20xBC						= f2->TP
Surf_T_slw_20xBC						= f2->surf_ts

Surf_T_tot_2xCO2						= f1->TS
Surf_T_fst_2xCO2						= f1->TP
Surf_T_slw_2xCO2						= f1->surf_ts


Ts_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(Surf_T_slw_2xCO2,0),1)
Ts_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(Surf_T_slw_20xBC,0),1)

Ts_total1 = dim_avg_n_Wrap(dim_avg_n_Wrap(Surf_T_tot_2xCO2,0),1)
Ts_total2 = dim_avg_n_Wrap(dim_avg_n_Wrap(Surf_T_tot_20xBC,0),1)

Ts_fast1 = wgt_areaave_Wrap(dim_avg_n_Wrap(Surf_T_fst_2xCO2,0),area,1,0)
Ts_fast2 = wgt_areaave_Wrap(dim_avg_n_Wrap(Surf_T_fst_20xBC,0),area,1,0)


drf1 = dim_avg_n_Wrap(dim_avg_n_Wrap((ERF_2xCO2/4.11),0),1)
drf2 = dim_avg_n_Wrap(dim_avg_n_Wrap((ERF_20xBC/2.59),0),1)


dLW_A1 = dim_avg_n_Wrap(dim_avg_n_Wrap((LW_allsky_2xCO2/4.11),0),1)
dLW_A2 = dim_avg_n_Wrap(dim_avg_n_Wrap((LW_allsky_20xBC/2.59),0),1)
dLW_A_diff = dLW_A2 - dLW_A1

dSW_A1 = dim_avg_n_Wrap(dim_avg_n_Wrap((SW_allsky_2xCO2/4.11),0),1)
dSW_A2 = dim_avg_n_Wrap(dim_avg_n_Wrap((SW_allsky_20xBC/2.59),0),1)
dSW_A_diff = dSW_A2 - dSW_A1

dLW_C1 = dim_avg_n_Wrap(dim_avg_n_Wrap((LW_clrsky_2xCO2/4.11),0),1)
dLW_C2 = dim_avg_n_Wrap(dim_avg_n_Wrap((LW_clrsky_20xBC/2.59),0),1)
dLW_C_diff = dLW_C2 - dLW_C1

dSW_C1 = dim_avg_n_Wrap(dim_avg_n_Wrap((SW_clrsky_2xCO2/4.11),0),1)
dSW_C2 = dim_avg_n_Wrap(dim_avg_n_Wrap((SW_clrsky_20xBC/2.59),0),1)
dSW_C_diff = dSW_C2 - dSW_C1

;;;;;;********************************** HGH CLOUDS ***************************************

High_cld_fst_20xBC						= fp3->CLDHGH  - fp1->CLDHGH
High_cld_fst_2xCO2						= fp2->CLDHGH  - fp1->CLDHGH

 
;;;;****************************************** HGH CLOUD fraction normalized with global mean surface temperature change***********************************************************

high_cld_fast1 = dim_avg_n_Wrap(dim_avg_n_Wrap(High_cld_fst_2xCO2*100/4.11,0),1)
high_cld_fast2 = dim_avg_n_Wrap(dim_avg_n_Wrap(High_cld_fst_20xBC*100/2.59,0),1)


high_cld_diff_20xBC_fast = (/high_cld_fast2 - high_cld_fast1/)

;;;;;;********************************** LOW CLOUDS ***************************************
  
Low_cld_fst_20xBC						= fp3->CLDLOW  - fp1->CLDLOW
Low_cld_fst_2xCO2						= fp2->CLDLOW  - fp1->CLDLOW


;;;;****************************************** LOW CLOUD fraction normalized with global mean surface temperature change***********************************************************


low_cld_fast1 = dim_avg_n_Wrap(dim_avg_n_Wrap(Low_cld_fst_2xCO2*100/4.11,0),1)
low_cld_fast2 = dim_avg_n_Wrap(dim_avg_n_Wrap(Low_cld_fst_20xBC*100/2.59,0),1)


low_cld_diff_20xBC_fast = (/low_cld_fast2 - low_cld_fast1/)

cloud_fst_20xBC						= fp3->CLOUD  - fp1->CLOUD
cloud_fst_2xCO2						= fp2->CLOUD  - fp1->CLOUD


cloud_fast1 = wgt_areaave_Wrap((dim_avg_n_Wrap(cloud_fst_2xCO2*100/4.11,0)),area,1,0)
cloud_fast2 = wgt_areaave_Wrap((dim_avg_n_Wrap(cloud_fst_20xBC*100/2.59,0)),area,1,0)

cloud_diff_20xBc_fast = (/cloud_fast2 - cloud_fast1/)

var_cloud = new((/3,17/), typeof(cloud_fst_2xCO2))
var_cloud(0,:) = cloud_fast1(:)
var_cloud(1,:) = cloud_fast2(:)
var_cloud(2,:) = cloud_diff_20xBc_fast(:)


;;;;******************************************* PLOTTING VARIABLES FOR HGH CLOUDS *****************************************************


var = new((/24,192/), typeof(ERF_2xCO2))
var(0,:) = drf1(:)
var(1,:) = drf2(:)
var(2,:) = (/drf2(:) - drf1(:)/)


var(3,:) = dLW_A1(:)
var(4,:) = dLW_A2(:) 
var(5,:) = dLW_A_diff(:) 


var(6,:) = dSW_A1(:) 
var(7,:) = dSW_A2(:) 
var(8,:) = dSW_A_diff(:) 

var(9,:) = dLW_C1(:) + dSW_C1(:)
var(10,:) = dLW_C2(:) + dSW_C2(:)
var(11,:) = ((dLW_C2(:) + dSW_C2(:)) - (dLW_C1(:) + dSW_C1(:)))


var(12,:) = dLW_C1(:)
var(13,:) = dLW_C2(:) 
var(14,:) = dLW_C_diff(:) 


var(15,:) = dSW_C1(:) 
var(16,:) = dSW_C2(:) 
var(17,:) = dSW_C_diff(:) 


var(18,:) = high_cld_fast1(:)
var(19,:) = high_cld_fast2(:) 
var(20,:) = high_cld_diff_20xBC_fast(:) 

var(21,:) = low_cld_fast1(:)
var(22,:) = low_cld_fast2(:)
var(23,:) = low_cld_diff_20xBC_fast(:)


var!1 = "lat"
var&lat = lat


values = decimalPlaces(dim_avg_n_Wrap(var,1),2,True)
print(values)	
;************************************

;***************************************
;::::::::::::::::::::::::::::::::::::::::::::: PLOTING START :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;*******************************plotting
 
pltName    =  "Figure2"  
pltType    =  "png"                                           
plot       =  new (9,"graphic")


wks_type = "png"
wks_type@wkWidth = 3000
wks_type@wkHeight = 3000

wks       =  gsn_open_wks(wks_type, pltName)
  FontH = 0.015
FontT = "helvetica"
FontTb = "helvetica-bold"

CentString = (/"   Effective Radiative Forcing~C~              (all-sky)","LW Radiative Forcing~C~          (all-sky)","SW Radiative Forcing~C~          (all-sky)","   Effective Radiative Forcing~C~              (clear-sky)","LW Radiative Forcing~C~         (clear-sky)" ,"SW Radiative Forcing~C~         (clear-sky)","    High Cloud fraction ~C~      (Fast response)","   Low cloud fraction ~C~      (Fast response)"," Clouds vertically ~C~    (Fast response)" /)
LeftString = (/"a)   ","b) ","c)  ", "d) ", "e) ", "f) ", "g) ", "h) ", "i) "/)
YAxisString = (/"Unit","Unit","Unit","Unit","Unit","Unit","% [ERF]~S~-1~N~ ","% [ERF]~S~-1~N~"/)     ;K~S~-1~N~
XAxisString = (/" ","","","","",""," "," "/)

;*********************************************************
 res = True
 FontT = 21
 FontTb = 22
 FontH = 0.015
 ;                                   XY axis
 ;************************************************************************************   
 res@gsnDraw          = False                   
 res@gsnFrame         = False         
 res@tiXAxisString    = "  "
 res@tiXAxisFont = FontT
 res@tiYAxisFont =FontT
 res@tiXAxisFontHeightF = FontH+0.02
 res@tiYAxisFontHeightF = FontH ;+0.02
 res@tiYAxisOffsetXF = 0.005
 
 res@tmXBLabelFont = FontT
 res@tmYLLabelFont = FontT
 res@tmXBLabelFontHeightF = FontH+0.015
 res@tmYLLabelFontHeightF = FontH+0.015	
 res@tmYROn = False          
 res@tmXTOn = False 
 res@tmXBLabelDeltaF = -0.5
 res@tmYLLabelDeltaF = -0.05
  ;******************************************************** grid resources
 res@tmXMajorGrid                = True          
 res@tmXMajorGridThicknessF      = 3.5        
 res@tmXMajorGridLineDashPattern = 2            
 res@tmXMajorGridLineColor      ="grey"
 res@tmYMajorGrid                = True         
 res@tmYMajorGridThicknessF      = 3.5          
 res@tmYMajorGridLineDashPattern = 2            
 res@tmYMajorGridLineColor      ="grey"
 thicknesses                 = (/  9.00,   9.0, 7.0/)     
 colors                      = (/"black", "red","black" /) 
 labels                      = (/"  2XCO~B~2~N~","  20XBC","  20XBC - 2XCO~B~2~N~"/)           
 xyDashPatterns = (/0,0,1/)
 ;***********************************************************************************    
 res@xyExplicitLegendLabels =  labels           
 res@xyLineThicknesses = thicknesses   
 res@xyDashPatterns  = xyDashPatterns
 res@xyLineColors           = colors
 res@gsnStringFont = FontT
 res@gsnCenterStringFontHeightF = FontH +0.02
 res@gsnLeftStringFontHeightF = FontH +0.02
 res@gsnRightStringFontHeightF = FontH +0.02
 res@gsnCenterString            = "High-lev Cloud change"
  res@tiYAxisFont = FontT
  res@tiYAxisFontHeightF = FontH+0.02
  res@tiMainFont = "Helvetica" 
  res@tmXBLabelFont = "Helvetica"
  res@tmYLLabelFont = "Helvetica" 
  res@tmYRLabelFont = "Helvetica" 
  res@gsnMaximize = True
                                            ;LEGENDS
 ;************************************************************************************   
  res@pmLegendSide           = "Top"            
  res@pmLegendParallelPosF   = .225               
  res@pmLegendOrthogonalPosF = -0.85              
  res@pmLegendWidthF         = 0.15              
  res@pmLegendHeightF        = 0.10             
  res@lgPerimOn              = False             
  res@lgLabelFontHeightF     = .015              
 res@xyExplicitLegendLabels =  labels          
 
j = 0
   do i = 0,7
    res@gsnCenterString = CentString(i)
    res@gsnLeftString = LeftString(i)
    res@tiYAxisString = YAxisString(i)
    res@tiXAxisString = XAxisString(i)
    plot(i) = gsn_csm_xy(wks, lat, var(j:j+2,:), res)
	j = j+3

  end do

XAxisString =(/"% [ERF]~S~-1~N~"/) 
YAxisString =(/"hPa"/) 
j = 0

    res@gsnCenterString = CentString(8)
    res@gsnLeftString = LeftString(8)
    res@tiXAxisString = XAxisString(0)
    res@tiYAxisString = YAxisString(0)
		 res@trYReverse        = True                       
    plot(8) = gsn_csm_xy(wks, var_cloud(j:j+2,:),lev, res)

lgres                    = True
 lgres@vpWidthF           = 0.14                  
 lgres@vpHeightF          = 0.10                 
 lgres@lgPerimOn          = False                
 lgres@lgDashIndexes      = xyDashPatterns         
lgres@lgLineLabelStrings = labels    
 lgres@lgLabelFontHeightF = 0.1                
 lgres@lgLineColors       = colors
 lgres@lgPerimOn          = False 
  lgres@lgLineThicknessF  =  thicknesses
  lgres@lgLabelFont = "Helvetica-bold"
 gsn_legend_ndc(wks,3,labels,0.10,0.61,lgres)
 
panelres = True
gsn_panel(wks,plot(0:8),(/3,3/),panelres) 


  
exit()
    

end

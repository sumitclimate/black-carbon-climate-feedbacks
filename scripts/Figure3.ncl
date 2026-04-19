begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Load DATA Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
psst_1co2 = "1xCO2_cam5_PSST_plev_file.nc"
som_1co2 = "1xCO2_cam5_SOM_plev_file.nc"

psst_2co2 = "2xCO2_cam5_PSST_plev_file.nc"
som_2co2 = "2xCO2_cam5_SOM_plev_file.nc"

psst_20BC = "20xBC_cam5_PSST_plev_file.nc"
som_20BC = "20xBC_cam5_SOM_plev_file.nc"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("surafce_temp_alpha_values.nc","r")

lat = psst_1co2->lat
lon = psst_1co2->lon
lev = psst_1co2->lev_p
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
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  global mean surface temperature change

ts1 = som_1co2->TS
ts2= som_2co2->TS

dts = ts2 - ts1
copy_VarMeta(ts2,dts)

surf_t_tot = dim_avg_n_Wrap(dts,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)

print(dsurf_t_tot)


tp1 = psst_1co2->TS
tp2 = psst_2co2->TS

dtp = tp2 - tp1
copy_VarMeta(tp2,dtp)

surf_t_fst = dim_avg_n_Wrap(dtp,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)

print(dsurf_t_fst)


dts_slow_co2 = (dts - dtp)/3.48
copy_VarMeta(tp2,dts_slow_co2)



surf_t_slow_co2 = dim_avg_n_Wrap(dts_slow_co2,0)
dsurf_t_slow = wgt_areaave_Wrap(surf_t_slow_co2,area,1,0)


print(dsurf_t_slow)



ts1 = som_1co2->TS
ts2= som_20BC->TS

dts = ts2 - ts1
copy_VarMeta(ts2,dts)

surf_t_tot = dim_avg_n_Wrap(dts,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)

print(dsurf_t_tot)


tp1 = psst_1co2->TS
tp2 = psst_20BC->TS

dtp = tp2 - tp1
copy_VarMeta(tp2,dtp)

surf_t_fst = dim_avg_n_Wrap(dtp,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)

print(dsurf_t_fst)


dts_slow_BC = (dts - dtp)/1.65
copy_VarMeta(tp2,dts_slow_BC)



surf_t_slow_BC = dim_avg_n_Wrap(dts_slow_BC,0)
dsurf_t_slow = wgt_areaave_Wrap(surf_t_slow_BC,area,1,0)


print(dsurf_t_slow)
;;;####################################################################################

alpha = ttest_file->alpha
 
alpha!0 = "lat"
alpha!1 = "lon"

alpha&lat = lat
alpha&lon = lon

dt_slow_diff = (/dts_slow_BC-dts_slow_co2/)
;;;################################################################

Ts_co2 = dim_avg_n_Wrap(dts_slow_co2,0)
Ts_BC = dim_avg_n_Wrap(dts_slow_BC,0)
Ts_diff = dim_avg_n_Wrap((dts_slow_BC-dts_slow_co2),0)

zonal_Ts_co2 = dim_avg_n_Wrap(Ts_co2,1)
zonal_Ts_BC = dim_avg_n_Wrap(Ts_BC,1)
zonal_Ts_diff = dim_avg_n_Wrap((Ts_diff),1)

Ts_co2!0 = "lat"
Ts_co2!1 = "lon"

Ts_co2&lat = lat
Ts_co2&lon = lon

copy_VarMeta(Ts_co2,Ts_BC)
copy_VarMeta(Ts_co2,Ts_diff)

var = new((/3,192,288/), typeof(Ts_co2))
var2 = new((/3,192/), typeof(Ts_co2))

var(0,:,:) = Ts_co2
var(1,:,:) = Ts_BC
var(2,:,:) = Ts_diff 

var2(0,:) = zonal_Ts_co2
var2(1,:) = zonal_Ts_BC
var2(2,:) = (/zonal_Ts_diff/) 

values = decimalPlaces(wgt_areaave_Wrap(var,area,1,0),2,True)
print(values)


;***************************************
FontT = "helvetica"
FontH = 0.025
FontTb = "helvetica-bold"
;*******************************plotting
pltName    =  "Figure3"  
pltType    =  "png"      
cmap = read_colormap_file("BlueWhiteOrangeRed")                                     
plot       =  new (3,"graphic")
wks       =  gsn_open_wks(pltType, pltName)
;********************************************************* resources

 res                 = True                   
 res@gsnDraw         = False                  
 res@gsnFrame        = False                   
 res@tiMainOn        = False             
 res@gsnLeftStringFontHeightF = FontH + 0.005
 res@gsnRightStringFontHeightF = FontH + 0.005
 res@gsnStringFont = FontT
 res@gsnAddCyclic = True
 res@lbLabelBarOn    = True                  
 res@lbTitleFont     = "helvetica"    
 res@lbLabelFont     = "helvetica"    
 res@lbLabelFontHeightF    = FontH+0.005
 res@lbTitleFontHeightF    = FontH +0.005
 res@lbBoxLinesOn  = True
 res@lbBoxSeparatorLinesOn = True
 res@lbLabelAutoStride = False
 res@lbAutoManage       = False             
 res@lbBoxEndCapStyle   = "TriangleBothEnds"
 res@lbTitlePosition    = "Bottom"      
 res@lbLabelStride      = 1 
 res@lbBoxLineThicknessF   = 0.1
 res@lbTitleString = "K K~S~-1~N~"

 res@cnLinesOn       = True  
 res@cnLineThicknessF = 0.4          
 res@cnLineLabelsOn  = False                   
 res@cnLevelSelectionMode = "ManualLevels"                   
 res@cnFillOn         = True                 
 res@cnFillPalette    = cmap(135::,:)           
 res@cnInfoLabelOn = False		
 res@mpProjection = "Robinson"
 res@mpGeophysicalLineThicknessF = 0.80
 res@mpGridLatSpacingF = 30.              
 res@mpGridLonSpacingF = 60.             
 res@mpGridAndLimbOn        = True
 res@mpGridLineDashPattern  = 2
 res@mpGridLatSpacingF = 30
 res@mpGridLineColor = "black"
 res@mpOutlineBoundarySets = "Geophysical"
 res@mpGridLineThicknessF = 2
 res@mpPerimLineThicknessF = 3
 res@gsnStringFontHeightF = FontH +0.01
 res@pmLabelBarHeightF = 0.1
 res@pmLabelBarWidthF = 0.59
 res@pmLabelBarOrthogonalPosF = 0.040
 res@gsnMaximize            = True             
 res@cnMinLevelValF   = 0.0
 res@cnMaxLevelValF   = 3.2
 res@cnLevelSpacingF = 0.4
 res@lbLabelStride = 2
 res@lbBoxLinesOn  = True
 res@lbBoxLineThicknessF = 2 

;***************************************
res@gsnLeftString = "a) 2XCO~B~2~N"  
res@gsnRightString = "1.0" 
plot(0) = gsn_csm_contour_map(wks,var(0,:,:), res)

res@gsnLeftString = "b)  20XBC"
res@gsnRightString = "1.0"  
plot(1) = gsn_csm_contour_map(wks,var(1,:,:), res)

delete(res@cnFillPalette)
delete(cmap)
cmap = read_colormap_file("BlueYellowRed") 


 res@lbLabelBarOn    = True
 res@cnMinLevelValF   = -1.2
 res@cnMaxLevelValF   = 1.2
 res@cnLevelSpacingF = 0.2
 res@lbLabelStride = 4
 res@cnFillPalette = cmap
 res@cnLinesOn       = True
 res@cnLineThicknessF = 0.2
res@gsnLeftString = "c) 20XBC~N~-~N~2XCO~B~2~N"   
res@gsnRightString = "0." 

plot(2) = gsn_csm_contour_map(wks,var(2,:,:), res)
 ;***********************************************************************************
  sgres                      = True		
  sgres@gsnDraw              = False		
  sgres@gsnFrame             = False		
  sgres@cnInfoLabelOn        = False		
  sgres@cnLinesOn            = False		
  sgres@cnLineLabelsOn       = False		
  sgres@cnFillScaleF         = 0.95		
  sgres@cnFillDotSizeF       = 0.0020
  sgres@gsnAddCyclic         = False
  sgres@cnFillOn             = True
  sgres@cnFillColors         = (/"transparent","transparent"/) 
  sgres@cnLevelSelectionMode = "ExplicitLevels"	         
  sgres@cnLevels             = 99.0                       
  sgres@lbLabelBarOn         = False

  sgres@tiMainString         = ""    
  sgres@gsnCenterString      = ""  
  sgres@gsnLeftString        = ""    
  sgres@gsnRightString       = ""   

  sig_plot = gsn_csm_contour(wks,alpha,sgres)


  opt                  = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17

  sig_plot = gsn_contour_shade(sig_plot,-999.,99.0,opt)

  overlay(plot(2),sig_plot)
 
 ;***********************************************************************************

txres = True
txres@txFontHeightF = 0.022
txres@txFont        = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks,"Surface temperature response ",0.53,0.77,txres)

txres1 = txres
txres1@txFontHeightF = 0.034
txres1@txFont        = "helvetica"
label = (/"60~S~o~N~S","30~S~o~N~S","Eq","30~S~o~N~N","60~S~o~N~N"/)

text = gsn_add_text(wks,plot(0),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(1),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(2),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
 

; 
;***************************** XY Plot resources *************************************
 res2                        = True              
 wks2       =  gsn_open_wks("png", "check")
 thicknesses                 = (/4.5/)     
 colors                      = (/ "black"/) 
 xyDashPatterns = (/0/)
 font = 21
 xyplot       =  new (3,"graphic")

 ;************************************************************************************   
 res2@tiYAxisString    = "Response"
 res2@tiXAxisString    = ""
 res2@tiXAxisFont = FontT
 res2@tiYAxisFont =FontT
 res2@tiXAxisFontHeightF = FontH+0.005
 res2@tiYAxisFontHeightF = FontH+0.005
 res2@tiMainFontHeightF  = FontH
 res2@tiYAxisFontColor = "black" 
 res2@tmXBLabelFont = FontTb
 res2@tmYLLabelFont = FontT
 res2@tmXBLabelFontHeightF = FontH+0.025
 res2@tmYLLabelFontHeightF = 0.0 
 res2@tmYROn = True           
 res2@tmXTOn = False          
 res2@tmXBFormat = ".1f"
 res2@gsnDraw          = False                   
 res2@gsnFrame         = False                  
  ;******************************************************** grid resources
 res2@tmXMajorGrid                = True          
 res2@tmXMajorGridThicknessF      = 2.5          
 res2@tmXMajorGridLineDashPattern = 2           
 res2@tmXMajorGridLineColor      ="grey"
 
 res2@tmYMajorGrid                = True          
 res2@tmYMajorGridThicknessF      = 2.5          
 res2@tmYMajorGridLineDashPattern = 2            
 res2@tmYMajorGridLineColor      ="grey"

 res2@gsnStringFont = FontT
 res2@gsnCenterStringFontColor = "black"
 res2@gsnStringFontHeightF = 0.005 
 res@xyLineThicknessF   = 3.0
 res2@tiYAxisString    = " "
 res2@tiXAxisString    = "K K~S~-1~N~"
  ;                                   LEGENDS
 ;************************************************************************************   
 res2@gsnCenterString   = "  "
 res2@xyLineThicknesses = thicknesses(:)   
 res2@xyDashPatterns  = xyDashPatterns(:)
 res2@xyLineColors           = colors(:)
 res2@trXMinF     = 0.0
 res2@trXMaxF     = 4.0
   
 xyplot(0) = gsn_csm_xy(wks,var2(0,:),lat,res2) 
 xyplot(1) = gsn_csm_xy(wks,var2(1,:),lat,res2) 
  
 res2@gsnCenterString   = " "
delete(res2@trXMinF )
delete(res2@trXMaxF )
 res2@tiYAxisString    = "   "
 res2@trXMinF     = -2.
 res2@trXMaxF     =  2.	
 xyplot(2) = gsn_csm_xy(wks,var2(2,:),lat,res2) 


newplot1 = gsn_attach_plots(plot(0),xyplot(0),res,res2)
newplot2 = gsn_attach_plots(plot(1),xyplot(1),res,res2)
newplot3 = gsn_attach_plots(plot(2),xyplot(2),res,res2)
 

 panres                        = True
 panres@gsnMaximize            = True             
 panres@pmTickMarkDisplayMode   = "Always"       
 panres@gsnFrame         = False 
 panres@lbTitleString = res@lbTitleString 
 panres@lbTitleFont     = FontT    
 panres@lbLabelFont     = FontT   
 panres@lbLabelFontHeightF    = FontH+0.0013
 panres@lbTitleFontHeightF    = FontH +0.002  
 panres@lbLabelAutoStride = False
 panres@lbLabelStride = 3
 panres@lbTitlePosition = res@lbTitlePosition 
 panres@lbBoxLineThicknessF   = 2.0
 panres@lbAutoManage          = False          
 panres@cnInfoLabelOn = False		
 panres@gsnPanelLabelBar    = False           
 panres@lbOrientation            = "horizontal"
 panres@lbBoxEndCapStyle   = "TriangleBothEnds"
 panres@gsnPanelMainString = " "
 panres@gsnPanelMainFontHeightF = 0.018 
 panres@gsnMaximize = True
 panres@gsnPanelRowSpec = False
 panres1 =panres
 panres1@pmLabelBarOrthogonalPosF = 0.070
 panres1@lbTitleString = "[K]"

gsn_panel(wks,plot(0:2),(/2,2/),panres1)

 frame(wks)
 exit()
 ;***********************************************************************************
 ;*************************	**********************************************************

end

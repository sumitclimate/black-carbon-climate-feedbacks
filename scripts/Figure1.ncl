begin

;;;;;;;;;;;;;;;;;; Feedback Data file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f1 = addfile("2xCO2_cam5_feedback_values.nc","r")
f2 = addfile("20xBC_cam5_feedback_values.nc","r")

lat = f1->lat
lon = f1->lon
lat@units = "degree_north"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("ERF_LW_SW_alpha_values.nc","r")

alpha1 = ttest_file->alpha 
alpha2 = ttest_file->alpha_lw
alpha3 = ttest_file->alpha_sw 

alpha1!0 = "lat"
alpha1!1 = "lon"

alpha1&lat = lat
alpha1&lon = lon

copy_VarMeta(alpha1,alpha2)
copy_VarMeta(alpha1,alpha3)
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
;********************************** 
; Variables for CO2
T_fast_co2 = f1->TP
ERF_co2 = f1->ERF
LW_rf_clear_co2 = f1->LW_clear_sky_psst 
SW_rf_clear_co2 = f1->SW_clear_sky_psst
LW_rf_all_co2 = f1->LW_all_sky_psst
SW_rf_all_co2 = f1->SW_all_sky_psst

; Variables for BC
T_fast_BC = f2->TP
ERF_BC = f2->ERF
LW_rf_clear_BC = f2->LW_clear_sky_psst 
SW_rf_clear_BC = f2->SW_clear_sky_psst
LW_rf_all_BC = f2->LW_all_sky_psst
SW_rf_all_BC = f2->SW_all_sky_psst

;Normalized ERF, LW, and SW with global mean ERF and taking time average
drf1 = dim_avg_n_Wrap((ERF_co2/4.11),0)
drf2 = dim_avg_n_Wrap((ERF_BC/2.59),0)
diff_rf = drf2 - drf1

LW1 = dim_avg_n_Wrap((LW_rf_all_co2/4.11),0)
LW2 = dim_avg_n_Wrap((LW_rf_all_BC/2.59),0)
diff_LW = LW2 - LW1

SW1 = dim_avg_n_Wrap((SW_rf_all_co2/4.11),0)
SW2 = dim_avg_n_Wrap((SW_rf_all_BC/2.59),0)
diff_SW = SW2 - SW1


var = new((/9,192,288/), typeof(drf1))
var(0,:,:) = drf1(:,:)
var(1,:,:) = drf2(:,:)
var(2,:,:) = diff_rf 
var(3,:,:) = LW1 
var(4,:,:) = LW2 
var(5,:,:) = diff_LW 
var(6,:,:) = SW1
var(7,:,:) = SW2 
var(8,:,:) = diff_SW 

values = decimalPlaces(wgt_areaave_Wrap(var,area,1,0),2,True)
print(values)

var!1 = "lat"
var&lat = lat
var!2 = "lon"
var&lon = lon
;*************************************** Plloting
FontT = "helvetica"
FontH = 0.025
FontTb = "helvetica-bold"
;*******************************plotting
pltName    =  "Figure1"  
pltType    =  "png"      
cmap = read_colormap_file("BlueWhiteOrangeRed")
cmap2 = read_colormap_file("BlueYellowRed")
pltType@wkWidth = 3000
pltType@wkHeight = 3000
plot       =  new (9,"graphic")
wks       =  gsn_open_wks(pltType, pltName)

;********************************************************* resources

 res                 = True                   
 res@gsnDraw         = False                  
 res@gsnFrame        = False                   
 res@tiMainOn        = False             
 res@gsnLeftStringFontHeightF = 0.035
 res@gsnRightStringFontHeightF = 0.035
 res@gsnStringFont = FontT
 res@lbLabelBarOn    = False                  
 res@lbTitleFont     = "helvetica"    
 res@lbLabelFont     = "helvetica"    
 res@lbLabelFontHeightF    = 0.010
 res@lbTitleFontHeightF    = 0.035
 res@lbBoxLinesOn  = True
 res@lbLabelAutoStride = False
 res@lbAutoManage       = False             
 res@lbBoxEndCapStyle   = "TriangleBothEnds"
 res@lbTitlePosition    = "Bottom"    
 res@lbLabelStride      = 1 
 res@cnLinesOn       = True  
 res@cnLineThicknessF = 0.4              
 res@cnLineLabelsOn  = False                   
 res@cnLevelSelectionMode = "ManualLevels"                   
 res@cnFillOn         = True                 
 res@cnFillPalette    = cmap         
 res@cnInfoLabelOn = False		
 res@mpProjection = "Robinson"
 res@mpGeophysicalLineThicknessF = 2
 res@mpGridLatSpacingF = 30.              
 res@mpGridLonSpacingF = 60.             
 res@mpGridAndLimbOn        = True
 res@mpGridLineDashPattern  = 1
 res@mpGridLatSpacingF = 30
 res@mpGridLineColor = "black"
 res@mpOutlineBoundarySets = "Geophysical"
 res@mpGridLineThicknessF = 7
 res@mpLimbLineThicknessF = 10
 res@mpGeophysicalLineThicknessF = 5
 res@mpPerimLineThicknessF = 6
 res@pmLabelBarHeightF = 0.14
 res@pmLabelBarWidthF = 0.59
 res@pmLabelBarOrthogonalPosF = 0.060
 res@lbTitleString = " "
 res@gsnMaximize            = True        
 res@cnMinLevelValF   = -3.0
 res@cnMaxLevelValF   = 3.0
 res@cnLevelSpacingF = 0.30
 res@lbLabelStride = 5
 res@lbBoxLinesOn  = True
 res@lbBoxLineThicknessF = 5


 ;***********************************************************************************
 res@gsnLeftString = "a) 2XCO~B~2~N" 
 res@gsnRightString = "1.0"  
plot(0) = gsn_csm_contour_map(wks,var(0,:,:), res)
 
 res@gsnLeftString = "d) 2XCO~B~2~N"  
 res@gsnRightString = "0.82"   
plot(3) = gsn_csm_contour_map(wks,var(3,:,:), res)
 
 res@gsnLeftString = "g) 2XCO~B~2~N"  
 res@gsnRightString = "0.18"    
plot(6) = gsn_csm_contour_map(wks,var(6,:,:), res)

 delete(res@gsnLeftString)
 res@gsnLeftString = "b) 20XBC"
 res@gsnRightString = "1.0"    
plot(1) = gsn_csm_contour_map(wks,var(1,:,:), res)
 res@gsnLeftString = "e) 20XBC"
 res@gsnRightString = "-0.87"    
plot(4) = gsn_csm_contour_map(wks,var(4,:,:), res)
 res@gsnLeftString = "h) 20XBC"
 res@gsnRightString = "1.87"    
plot(7) = gsn_csm_contour_map(wks,var(7,:,:), res)

 delete(res@cnFillPalette )
 delete(res@cnMinLevelValF )
 delete(res@cnMaxLevelValF )
 delete(res@cnLevelSpacingF )
 delete(res@lbLabelStride )

 res2 = res
 res2@cnMinLevelValF   = -4.0
 res2@cnMaxLevelValF   = 4.0
 res2@cnLevelSpacingF = 0.50
 res2@lbLabelStride = 4
 res2@lbLabelBarOn    = False    
 res2@cnFillPalette   = cmap
 res2@gsnLeftString = "c) 20XBC~N~-~N~2XCO~B~2~N"   
 res2@gsnRightString = "0"   
plot(2) = gsn_csm_contour_map(wks,var(2,:,:), res2)

  res3 = res2
 res3@gsnLeftString = "f) 20XBC~N~-~N~2XCO~B~2~N" 
 res3@gsnRightString = "-1.7"   
plot(5) = gsn_csm_contour_map(wks,var(5,:,:), res3)


   res4 = res2
 res4@lbLabelBarOn    = True  
 res4@gsnLeftString = "i) 20XBC~N~-~N~2XCO~B~2~N" 
 res4@gsnRightString = "1.69"   
 res4@lbLabelFontHeightF    = 0.030
plot(8) = gsn_csm_contour_map(wks,var(8,:,:), res4)



 ;*************** Stiplling for significant point in the difference plots *************************************
 
  sgres                      = True		
  sgres@gsnDraw              = False		
  sgres@gsnFrame             = False		
  sgres@cnInfoLabelOn        = False		
  sgres@cnLinesOn            = False		
  sgres@cnLineLabelsOn       = False		
  sgres@cnFillScaleF         = 0.7		
  sgres@gsnAddCyclic         = False
  sgres@cnFillOn             = True
  sgres@cnFillColors         = (/"transparent","transparent"/) 
  sgres@cnLevelSelectionMode = "ExplicitLevels"	        
  sgres@cnLevels             = 95.0                         
  sgres@lbLabelBarOn         = False
  sgres@tiMainString         = ""    
  sgres@gsnCenterString      = ""  
  sgres@gsnLeftString        = ""   
  sgres@gsnRightString       = ""  

  sig_plot0 = gsn_csm_contour(wks,alpha1,sgres)
  sig_plot1 = gsn_csm_contour(wks,alpha2,sgres)
  sig_plot2 = gsn_csm_contour(wks,alpha3,sgres)


  opt                  = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17

  sig_plot0 = gsn_contour_shade(sig_plot0,-999.,95.0,opt)
  sig_plot1 = gsn_contour_shade(sig_plot1,-999.,95.0,opt)
  sig_plot2 = gsn_contour_shade(sig_plot2,-999.,95.0,opt)

  overlay(plot(2),sig_plot0)
  overlay(plot(5),sig_plot1)
  overlay(plot(8),sig_plot2)
 
 ;***********************************************************************************
txres = True
txres@txFontHeightF = 0.022
txres@txFont        = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks," Effective radiative forcing  ",0.53,0.95,txres)
gsn_text_ndc(wks," Longwave radiative forcing ",0.53,0.71,txres)
gsn_text_ndc(wks," Shortwave radiative forcing ",0.53,0.47,txres)

txres1 = txres
txres1@txFontHeightF = 0.034
txres1@txFont        = "helvetica"
label = (/"60~S~o~N~S","30~S~o~N~S","Eq","30~S~o~N~N","60~S~o~N~N"/)

text = gsn_add_text(wks,plot(0),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(1),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(2),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
 ;***********************************************************************************
unit = " Wm~S~-2~N~"
 panres                        = True
 panres@gsnMaximize            = True             

 panres@pmTickMarkDisplayMode   = "Always"       
 panres@gsnFrame         = False   
 panres@gsnPanelLabelBar    = False    
 panres@lbTitleFont                = res@lbTitleFont 
 panres@lbTitleString = " " 
 panres@lbTitlePosition = "Bottom"
 panres@lbOrientation            = "horizontal"
 panres@lbLabelFont                = "helvetica"    
 panres@lbLabelFontHeightF    = 0.015
 panres@lbTitleFontHeightF    = 0.017
 panres@lbLabelAutoStride = False
 panres@lbLabelStride = 4
 panres@lbAutoManage          = False           
 panres@lbBoxLinesOn  = True
 panres@lbBoxLineThicknessF = 5
 panres@pmLabelBarHeightF = 0.069
 panres@pmLabelBarWidthF = 0.5
 panres@gsnPanelYWhiteSpacePercent = 5
 panres@gsnPanelXWhiteSpacePercent = 5
 panres@pmLabelBarParallelPosF = -0.13
panres@gsnPanelTop              = 0.99              
panres@gsnPanelBottom             = 0.69
panres@lbTitleString = " "
gsn_panel(wks,plot(0:2),(/1,3/),panres)

panres1 = panres
panres1@gsnPanelLabelBar = False
text = gsn_add_text(wks,plot(3),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(4),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(5),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
panres1@gsnPanelTop              = 0.75            
panres1@gsnPanelBottom             = 0.45
panres1@lbTitleString = " "
gsn_panel(wks,plot(3:5),(/1,3/),panres1)

panres2 = panres
panres2@gsnPanelLabelBar = True
text = gsn_add_text(wks,plot(6),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(7),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)
text = gsn_add_text(wks,plot(8),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres1)

panres2@gsnPanelTop              = 0.61             
panres2@gsnPanelBottom             = 0.05
panres2@pmLabelBarOrthogonalPosF = 0.070
panres2@lbTitleString = " "
gsn_panel(wks,plot(6:8),(/1,3/),panres2)




frame(wks)


exit()

end

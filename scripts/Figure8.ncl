begin
;;;;;;;;;;;;;;;;;; Feedback Data file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f1 = addfile("2xCO2_cam5_feedback_values.nc","r")
f2 = addfile("20xBC_cam5_feedback_values.nc","r")

lat = f1->lat
lon = f1->lon
lev = f1->lev_p
lev@units = "mb"
lat@units = "degree_north"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("Water_vapor_feedback_alpha_values.nc","r")
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
  
;**************************************************************
A = f1->WV_feed
B = f2->WV_feed
C = B-A

copy_VarMeta(B,C)
printVarSummary(C)
print(wgt_areaave_Wrap((dim_sum_n((dim_avg_n_Wrap(A,0)),0)),area,1,0))
print(wgt_areaave_Wrap((dim_sum_n((dim_avg_n_Wrap(B,0)),0)),area,1,0))
print(wgt_areaave_Wrap((dim_sum_n((dim_avg_n_Wrap(C,0)),0)),area,1,0))

var1 = dim_avg_n_Wrap(dim_avg_n_Wrap(A,0),2) 
var2 = dim_avg_n_Wrap(dim_avg_n_Wrap(B,0),2)
var3 = dim_avg_n_Wrap(dim_avg_n_Wrap(C,0),2)

var1!0 = "lev"
var1&lev = lev

var2!1 = "lat"
var2&lat = lat
var2!0 = "lev"
var2&lev = lev

var3!1 = "lat"
var3&lat = lat
var3!0 = "lev"
var3&lev = lev

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; tropopause position
ptrop = (250.0 - 150.0*clat^2.0)

t1 = f1->atm_t_slow
tropo_p =  conform(t1,ptrop,2)
copy_VarMeta(t1,tropo_p)
p = conform(t1,lev,1)
copy_VarMeta(t1,p)

WV_FEED_2xCO2 = mask(A, p.lt.tropo_p,False)
WV_FEED_20xBC = mask(B, p.lt.tropo_p,False)
WV_FEED_dif = mask(C, p.lt.tropo_p,False)

WV_FEED1 = dim_avg_n_Wrap(dim_sum_n_Wrap((dim_avg_n_Wrap(WV_FEED_2xCO2,0)),0),1)
WV_FEED2 = dim_avg_n_Wrap(dim_sum_n_Wrap((dim_avg_n_Wrap(WV_FEED_20xBC,0)),0),1)
WV_FEED_diff = dim_avg_n_Wrap(dim_sum_n_Wrap((dim_avg_n_Wrap(WV_FEED_dif,0)),0),1)

var = new((/3,192/), typeof(WV_FEED1))
var(0,:) = WV_FEED1(:)
var(1,:) = WV_FEED2(:)
var(2,:) = WV_FEED_diff(:)


var!1 = "lat"
var&lat = lat
;*************************************************
alpha = ttest_file->alpha
copy_VarCoords(var1,alpha)
;************************************************************************* plotting
pltName    =  "Figure8"  
pltType    =  "png"                                       
plot       =  new (4,"graphic")

cmap = read_colormap_file("BlueWhiteOrangeRed")
cmap2 = read_colormap_file("BlueYellowRed")

colorPal =  cmap(120::,:)
cnSpace = 0.08
FontH = 0.020
FontT = "helvetica"
FontTb = "helvetica-bold"
wks   = gsn_open_wks (pltType, pltName)        

 res                     = True        
 res@cnLabelBarEndStyle = "IncludeOuterBoxes"
 res@cnLinesOn       = True
 res@cnLineThicknessF = 0.2               
 res@cnLineLabelsOn  = False                   
 res@cnLevelSelectionMode = "ManualLevels"                   
 res@cnFillOn         = True                 
 res@cnFillPalette    = colorPal           
 res@cnInfoLabelOn = False	
 res@gsnDraw         = False                  
 res@gsnFrame        = False    
 res@gsnYAxisIrregular2Linear = True
 res@gsnStringFontHeightF = FontH +0.015
 res@gsnStringFont = FontT
 res@lbLabelBarOn    = False                  
 res@lbTitleFont     = FontT   
 res@lbLabelFont     = FontT   
 res@lbLabelFontHeightF    = FontH+0.01
 res@lbTitleFontHeightF    = FontH +0.01
 res@lbLabelAutoStride = False
 res@lbAutoManage       = False             
 res@lbTitlePosition    = "Bottom"      
 res@lbTitleString    = "Wm~S~-2~N~K~S~-1~N~/100 hPa"
 res@lbBoxEndCapStyle   = "TriangleBothEnds"
 res@pmLabelBarHeightF = 0.12
 res@pmLabelBarOrthogonalPosF = -0.01
 res@pmLabelBarWidthF = 0.72
 res@tmYLLabelDeltaF = -0.5
 res@tmXBLabelDeltaF = -0.5
 res@tiMainOn        = False      
 res@tiYAxisString      = "Pressure (hPa)"
 res@tiYAxisFont = FontT
 res@tiYAxisFontHeightF = FontH+0.015
 res@tiXAxisFontHeightF = FontH+0.015
 res@tiXAxisOffsetYF  = 0
 res@tiYAxisOffsetXF = 0.0002
 res@tiDeltaF = 0.9
 res@pmTickMarkDisplayMode = "Never"   
 res@tmXTOn = False           
 res@tmYROn = False           
 res@tmYLMode        = "Explicit"	
 res@tmYLValues      = (/800,600,400,200/)
 res@tmYLLabels      = "" + res@tmYLValues 
 res@tmXBLabelFontHeightF = FontH+0.01
 res@tmYLLabelFontHeightF = res@tmXBLabelFontHeightF
 res@trYMinF = 100.
 res@trYMaxF = 930.
 res@trYReverse  = True 
 
;******************************
 res@gsnLeftStringFontHeightF = FontH + 0.015
 res@gsnRightStringFontHeightF = FontH + 0.015
 res@gsnStringFont = FontT
 res@gsnLeftString = " "
 res@gsnRightString = ""
 res@gsnCenterString            = "  "	
 res@gsnCenterStringOrthogonalPosF = 0.025
 res@gsnLeftStringOrthogonalPosF = 0.025
 ;********************************************************************************
 res@cnMinLevelValF   = 0.0
 res@cnMaxLevelValF   = 0.36
 res@cnLevelSpacingF =   0.03
 res@lbLabelStride      = 2
 res@gsnLeftString            = "a) 2xCO~B~2~N~"
 res@gsnRightString = "1.53"
plot(0) = gsn_csm_contour(wks,var1,res)      
 
 res@gsnLeftString            = "b) 20xBC" 
 res@gsnRightString = "1.73"
plot(1) = gsn_csm_contour(wks,var2, res)

delete(res@cnFillPalette)
 res@lbLabelBarOn    = True                  
 res@cnMinLevelValF   = -0.24
 res@cnMaxLevelValF   = 0.24
 res@cnLevelSpacingF =   .03
 res@lbLabelStride      = 4
 res@gsnLeftString            = "c) 20xBC~N~-2xCO~B~2~N~"
 res@cnFillPalette    = cmap2  
 res@gsnRightString = "0.20"
plot(2) =  gsn_csm_contour(wks,var3, res)
 ;***********************************************************************************
 
  sgres                      = True		
  sgres@gsnDraw              = False		
  sgres@gsnFrame             = False		
  sgres@cnInfoLabelOn        = False		
  sgres@cnLinesOn            = False		
  sgres@cnLineLabelsOn       = False	
  sgres@cnFillScaleF         = 0.6		
  sgres@cnFillDotSizeF       = 0.0029
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
  sig_plot = gsn_csm_contour(wks,alpha,sgres)

  opt                  = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17

  sig_plot = gsn_contour_shade(sig_plot,-999.,80.0,opt)

  overlay(plot(2),sig_plot)
;**************************************************************************************************
 res2                        = True              
 thicknesses                 = (/  5.00,   5.0, 4.0/)     
 colors                      = (/"black", "red","black" /) 
 labels                      = (/"  2XCO~B~2~N~","  20XBC","  20XBC - 2XCO~B~2~N~"/)           
 xyDashPatterns = (/0,0,1/)
 font = 21
 xyplot       =  new (1,"graphic")
 res2@gsnLeftStringFontHeightF = FontH + 0.015
 res2@gsnRightStringFontHeightF = FontH + 0.015
 res2@tiYAxisString    = "WV Feedback (W m~S~-2~N~ K~S~-1~N~)"  
 res2@tiXAxisString    = ""
 res2@tiXAxisFont = FontT
 res2@tiYAxisFont =FontT
 res2@tiXAxisFontHeightF = FontH+0.015
 res2@tiYAxisFontHeightF = FontH+0.015
 res2@tiMainFontHeightF  = FontH
 res2@tiYAxisFontColor = "black" 
 res2@tmXBLabelFont = FontT
 res2@tmYLLabelFont = FontT
 res2@tmXBLabelFontHeightF = FontH+0.01
 res2@tmYLLabelFontHeightF = FontH+0.01 
 res2@tmYROn = False         
 res2@tmYLOn = True           
 res2@tmXTOn = False            
 res2@tmXBFormat = ".1f"
 res@tmXBLabelDeltaF = -0.5
 res@tmYLLabelDeltaF = -0.05
 res2@gsnDraw          = False                 
 res2@gsnFrame         = False                 
 ;******************************************************** grid resources
 res2@tmXMajorGrid                = True          
 res2@tmXMajorGridThicknessF      = 2.5           
 res2@tmXMajorGridLineDashPattern = 2           
 res2@tmXMajorGridLineColor      ="grey"
 
 res2@tmYMajorGrid                = True         
 res2@tmYMajorGridThicknessF      = 2.0           
 res2@tmYMajorGridLineDashPattern = 2            
 res2@tmYMajorGridLineColor      ="grey"

 res2@gsnStringFont = FontT
 res2@gsnCenterStringFontColor = "black"
 res2@gsnStringFontHeightF = FontH +0.015
 res2@tiXAxisString    = " "
 
 res2@gsnCenterString   = "  "
 res2@gsnLeftString       = "d) Zonal mean WV feedback"
 res2@xyLineThicknesses = thicknesses(:)   
 res2@xyDashPatterns  = xyDashPatterns(:)
 res2@xyLineColors           = colors(:)
 plot(3) = gsn_csm_xy(wks,lat,var,res2) 
 
  ;                                   LEGENDS
 ;************************************************************************************  
 lgres                    = True
 lgres@vpWidthF           = 0.14                  
 lgres@vpHeightF          = 0.10                 
 lgres@lgPerimOn          = False                
 lgres@lgDashIndexes      = xyDashPatterns         
 lgres@lgLineLabelStrings = labels    
 lgres@lgLabelFontHeightF = 0.12                
 lgres@lgLineColors       = colors
 lgres@lgPerimOn          = False 
 lgres@lgLineThicknessF  =  thicknesses
 lgres@lgLabelFont = "Helvetica"
 gsn_legend_ndc(wks,3,labels,0.58,0.45,lgres)
 
 ;***********************************************************************************
txres = True
txres@txFontHeightF = 0.022
txres@txFont        = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks,"Water vapor feedback  ",0.53,0.94,txres)

 ;***********************************************************************************

panres                     = True
panres@gsnDraw         = True                 
panres@gsnFrame        = False                   
panres@gsnMaximize         = True           
panres@gsnPanelLabelBar    = True            
panres@gsnPaperOrientation = "landscape"       
panres@lbTitleString = res@lbTitleString
panres@lbTitlePosition =  res@lbTitlePosition
panres@lbTitleFont   = FontT  
panres@lbLabelFont   = FontT    
panres@lbAutoManage          = False            
panres@lbOrientation            = "horizontal"
panres@cnInfoLabelOn = False		
panres@lbLabelFontHeightF    = res@lbLabelFontHeightF/2
panres@lbTitleFontHeightF    =  res@lbTitleFontHeightF/2 
panres@gsnPaperOrientation = "Landscape"       
panres@gsnMaximize = True
panres@pmLabelBarOrthogonalPosF = -0.01         
panres@pmLabelBarParallelPosF  = 0.05
panres@pmLabelBarHeightF = 0.06
panres@pmLabelBarWidthF = 0.46
panres@gsnPanelLeft  = 0.01
panres@lbLabelStride = 2
panres@lbBoxLineThicknessF   = 0.5
panres@gsnPanelTop              = 0.92
panres@gsnPanelBottom = 0.52
gsn_panel(wks,plot(0:1),(/1,2/),panres) 

panres@gsnPanelLabelBar    = False            
txres               = True
txres@txFont   = FontTb  
txres@txFontHeightF = 0.026
panres@gsnPanelTop              = 0.50
panres@gsnPanelBottom = 0.1
gsn_panel(wks,plot(2:3),(/1,2/),panres) 

frame(wks)
 

end

 
 
 

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Load DATA Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("Atm_temperature_alpha_values.nc","r")

lat = fs2 ->lat
lon = fs2->lon
lev = fs2->lev_p
lev@units = "mb"
; ;********************************** Gaussian weights for Area averages

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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;         Slow global mean surface temperature change

TS1 = fs1->TS
TS2= fs2->TS

dTS = TS2 - TS1
copy_VarMeta(TS2,dTS)

surf_t_tot = dim_avg_n_Wrap(dTS,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)

print(dsurf_t_tot)


TP1 = fp1->TS
TP2 = fp2->TS

dTP = TP2 - TP1
copy_VarMeta(TP2,dTP)

surf_t_fst = dim_avg_n_Wrap(dTP,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)
print(dsurf_t_fst)

dTS_slow_co2 = dTS - dTP
copy_VarMeta(TP2,dTS_slow_co2)

surf_t_slow_co2 = dim_avg_n_Wrap(dTS_slow_co2,0)
dsurf_t_slow_co2 = wgt_areaave_Wrap(surf_t_slow_co2,area,1,0)
print(dsurf_t_slow_co2)

TS1 = fs1->TS
TS2= fs3->TS
dTS = TS2 - TS1
copy_VarMeta(TS2,dTS)

surf_t_tot = dim_avg_n_Wrap(dTS,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)
print(dsurf_t_tot)

TP1 = fp1->TS
TP2 = fp3->TS

dTP = TP2 - TP1
copy_VarMeta(TP2,dTP)

surf_t_fst = dim_avg_n_Wrap(dTP,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)

print(dsurf_t_fst)
dTS_slow_BC = dTS - dTP
copy_VarMeta(TP2,dTS_slow_BC)

surf_t_slow_BC = dim_avg_n_Wrap(dTS_slow_BC,0)
dsurf_t_slow_BC = wgt_areaave_Wrap(surf_t_slow_BC,area,1,0)


print(dsurf_t_slow_BC)

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       

ts1 = fs1->T
ts2 = fs2->T
ts3 = fs3->T
tp1= fp1->T
tp2=fp2->T
tp3=fp3->T

dt_tot_co2 = ts2-ts1 
dt_fst_co2 = (tp2-tp1)

dt_tot_BC = ts3-ts1 
dt_fst_BC =(tp3-tp1)

dt_slow_co2 = dt_tot_co2 - dt_fst_co2
dt_slow_BC = dt_tot_BC - dt_fst_BC



;************************************************************ loading file required variables;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
AtmstempF_co2 = dt_fst_co2
AtmstempT_co2 = dt_tot_co2
AtmstempS_co2 = dt_slow_co2/dsurf_t_slow_co2

AtmstempF_BC = dt_fst_BC
AtmstempT_BC = dt_tot_BC
AtmstempS_BC = dt_slow_BC/dsurf_t_slow_BC


AtmstempF_diff = AtmstempF_BC - AtmstempF_co2
AtmstempT_diff = AtmstempT_BC - AtmstempT_co2
AtmstempS_diff = AtmstempS_BC - AtmstempS_co2

meta = fp2->T

copy_VarMeta(meta,AtmstempF_BC)
copy_VarMeta(meta,AtmstempT_BC)
copy_VarMeta(meta,AtmstempS_BC)
copy_VarMeta(meta,AtmstempF_co2)
copy_VarMeta(meta,AtmstempT_co2)
copy_VarMeta(meta,AtmstempS_co2)
copy_VarMeta(meta,AtmstempF_diff)
copy_VarMeta(meta,AtmstempS_diff)
copy_VarMeta(meta,AtmstempT_diff)

A = AtmstempS_co2
B = AtmstempS_BC
C = AtmstempS_diff

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; tropopause position
ptrop = (250.0 - 150.0*clat^2.0)

t1 = fs1->T
tropo_p =  conform(t1,ptrop,2)
copy_VarMeta(t1,tropo_p)
p = conform(t1,lev,1)
copy_VarMeta(t1,p)

AtmstempS_co2 = mask(AtmstempS_co2, p.lt.tropo_p,False)
AtmstempS_BC = mask(AtmstempS_BC, p.lt.tropo_p,False)
AtmstempS_diff = mask(AtmstempS_diff, p.lt.tropo_p,False)


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

alpha = ttest_file->alpha
copy_VarCoords(var1,alpha)

;************************************************************************* plotting
pltName    =  "Figure7"  
pltType    =  "png"                                       
plot       =  new (3,"graphic")

cmap = read_colormap_file("BlueWhiteOrangeRed")
cmap2 = read_colormap_file("BlueYellowRed")

colorPal =  cmap
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
 res@lbTitleString    = "K K~S~-1~N~"
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
res@cnMaxLevelValF   = 2.4	
res@cnMinLevelValF   = -2.4
res@cnLevelSpacingF = 0.2
res@lbLabelStride      = 4 

res@gsnLeftString            = "a) 2xCO~B~2~N~"
 res@gsnRightString = "0.35"
plot(0) = gsn_csm_contour(wks,var1,res)      
 
res@gsnLeftString            = "b) 20xBC" 
 res@gsnRightString = "0.77"
plot(1) = gsn_csm_contour(wks,var2, res)

delete(res@cnFillPalette)
 res@lbLabelBarOn    = True                  
res@cnMaxLevelValF   = 1
res@cnMinLevelValF   = -1
res@cnLevelSpacingF = 0.2
res@lbLabelStride      = 2
res@gsnLeftString            = "c) 20xBC~N~-2xCO~B~2~N~"
res@cnFillPalette    = cmap2  
 res@gsnRightString = "0.42"
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
;***********************************************************************************
txres = True
txres@txFontHeightF = 0.022
txres@txFont        = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks,"Atmospheric temperature response ",0.53,0.94,txres)
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
gsn_panel(wks,plot(2),(/1,1/),panres) 

frame(wks)     
 exit


end

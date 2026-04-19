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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("Clouds_fraction_alpha_values.nc","r")

lat = fs3 ->lat
lon = fs2->lon
lev = fs2->lev_p
lev@units = "mb"
p=lev
fillVal = -999999999999.0 

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
  
  alpha = ttest_file->alpha_cld

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;         Slow global mean surface temperature change


TS1 = fs1->TS
TS2= fs2->TS

dTS = TS2 - TS1
copy_VarMeta(TS2,dTS)

surf_t_tot = dim_avg_n_Wrap(dTS,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)

TP1 = fp1->TS
TP2 = fp2->TS

dTP = TP2 - TP1
copy_VarMeta(TP2,dTP)

surf_t_fst = dim_avg_n_Wrap(dTP,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)

dTS_slow_co2 = dTS - dTP
copy_VarMeta(TP2,dTS_slow_co2)

surf_t_slow_co2 = dim_avg_n_Wrap(dTS_slow_co2,0)
dsurf_t_slow_co2 = wgt_areaave_Wrap(surf_t_slow_co2,area,1,0)

TS1 = fs1->TS
TS2= fs3->TS

dTS = TS2 - TS1
copy_VarMeta(TS2,dTS)

surf_t_tot = dim_avg_n_Wrap(dTS,0)
dsurf_t_tot = wgt_areaave_Wrap(surf_t_tot,area,1,0)

TP1 = fp1->TS
TP2 = fp3->TS

dTP = TP2 - TP1
copy_VarMeta(TP2,dTP)

surf_t_fst = dim_avg_n_Wrap(dTP,0)
dsurf_t_fst = wgt_areaave_Wrap(surf_t_fst,area,1,0)

dTS_slow_BC = dTS - dTP
copy_VarMeta(TP2,dTS_slow_BC)

surf_t_slow_BC = dim_avg_n_Wrap(dTS_slow_BC,0)
dsurf_t_slow_BC = wgt_areaave_Wrap(surf_t_slow_BC,area,1,0)

;##################################################################################

cloud_tot_20xBC						= fs3->CLOUD - fs1->CLOUD
cloud_fst_20xBC						= fp3->CLOUD  - fp1->CLOUD
cloud_slw_20xBC						= cloud_tot_20xBC - cloud_fst_20xBC

cloud_tot_2xCO2						= fs2->CLOUD - fs1->CLOUD
cloud_fst_2xCO2						= fp2->CLOUD  - fp1->CLOUD
cloud_slw_2xCO2						= cloud_tot_2xCO2 - cloud_fst_2xCO2


T_tot_20xBC						= fs3->T - fs1->T
T_fst_20xBC						= fp3->T  - fp1->T
T_slw_20xBC						= (/T_tot_20xBC - T_fst_20xBC/)

T_tot_2xCO2						= fs2->T - fs1->T
T_fst_2xCO2						= fp2->T  - fp1->T
T_slw_2xCO2						= (/T_tot_2xCO2 - T_fst_2xCO2/)


V_tot_20xBC						= fs3->V - fs1->V
V_fst_20xBC						= fp3->V  - fp1->V
V_slw_20xBC						= (/V_tot_20xBC - V_fst_20xBC/)

V_tot_2xCO2						= fs2->V - fs1->V
V_fst_2xCO2						= fp2->V  - fp1->V
V_slw_2xCO2						= (/V_tot_2xCO2 - V_fst_2xCO2/)

U_tot_20xBC						= fs3->U - fs1->U
U_fst_20xBC						= fp3->U  - fp1->U
U_slw_20xBC						= (/U_tot_20xBC - U_fst_20xBC/)

U_tot_2xCO2						= fs2->U - fs1->U
U_fst_2xCO2						= fp2->U  - fp1->U
U_slw_2xCO2						= (/U_tot_2xCO2 - U_fst_2xCO2/)


omega_tot_20xBC						= fs3->OMEGA - fs1->OMEGA
omega_fst_20xBC						= fp3->OMEGA  - fp1->OMEGA
omega_slw_20xBC						= (/omega_tot_20xBC - omega_fst_20xBC/)

omega_tot_2xCO2						= fs2->OMEGA - fs1->OMEGA
omega_fst_2xCO2						= fp2->OMEGA  - fp1->OMEGA
omega_slw_2xCO2						= (/omega_tot_2xCO2 - omega_fst_2xCO2/)

cloud_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(cloud_slw_2xCO2/dsurf_t_slow_co2,0),2) 
cloud_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(cloud_slw_20xBC/dsurf_t_slow_BC,0),2) 


U_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(U_slw_2xCO2,0),2) 
U_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(U_slw_20xBC,0),2) 
 

V_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(V_slw_2xCO2,0),2) 
V_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(V_slw_20xBC,0),2) 

omega_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(omega_slw_2xCO2,0),2) 
omega_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(omega_slw_20xBC,0),2) 

plev = 100*conform(omega_slw_2xCO2,p,1)
T_slw_2xCO2 = where(T_slw_2xCO2.ne.0,T_slw_2xCO2,T_slw_2xCO2@_FillValue)
T_slw_20xBC = where(T_slw_20xBC.ne.0,T_slw_20xBC,T_slw_20xBC@_FillValue)

W_slw_2xCO2 = omega_to_w(omega_slw_2xCO2, plev, T_slw_2xCO2) 
W_slw_20xBC = omega_to_w(omega_slw_20xBC, plev, T_slw_20xBC) 


W_slow1 = dim_avg_n_Wrap(dim_avg_n_Wrap(W_slw_2xCO2,0),2) 
W_slow2 = dim_avg_n_Wrap(dim_avg_n_Wrap(W_slw_20xBC,0),2) 


cloud_diff_20xBC_slow = (/cloud_slow2 - cloud_slow1/)

V_diff_20xBC_slow = (/V_slow2 - V_slow1/)

W_diff_20xBC_slow = (/W_slow2 - W_slow1/)
;***************************************

cloud_slow = new((/3,17,192/), typeof(cloud_tot_2xCO2))

cloud_slow(0,:,:) = cloud_slow1(:,:) *100
cloud_slow(1,:,:) = cloud_slow2(:,:) *100
cloud_slow(2,:,:) = (/cloud_diff_20xBC_slow/) *100 

cloud_slow!2 = "lat"
cloud_slow&lat = lat
cloud_slow!1 = "lev"
cloud_slow&lev = lev

print(dim_sum_n_Wrap(dim_avg_n_Wrap(cloud_slow,2),1))

v = new((/3,17,192/), typeof(V_tot_2xCO2))
w = new((/3,17,192/), typeof(V_tot_2xCO2))

 wAve1   = avg(W_slow2(:,:))           
 vAve1   = avg(V_slow2(:,:))
 scale1  = fabs(vAve1/wAve1)
 scale1  = 100000
 
v(0,:,:) = V_slow1(:,:)
v(1,:,:) = V_slow2(:,:) 
v(2,:,:) = (/V_diff_20xBC_slow/)  

w(0,:,:) = W_slow1(:,:)*scale1
w(1,:,:) = W_slow2(:,:)*scale1
w(2,:,:) = (/W_diff_20xBC_slow/)*scale1
 

copy_VarCoords(cloud_slow,v)
copy_VarCoords(cloud_slow,w)

v!2 = "lat"
v&lat = lat
v!1 = "lev"
v&lev = lev

w!2 = "lat"
w&lat = lat
w!1 = "lev"
w&lev = lev

;************************************************************************* plotting
pltName    =  "Figure9"  
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
 res@lbTitleString    = "% K~S~-1~N~"
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
 res@gsnLeftStringFontHeightF = FontH + 0.015
 res@gsnRightStringFontHeightF = FontH + 0.015
 res@gsnStringFont = FontT
 res@gsnLeftString = " "
 res@gsnRightString = ""
 res@gsnCenterString            = "  "	
 res@gsnCenterStringOrthogonalPosF = 0.025
 res@gsnLeftStringOrthogonalPosF = 0.025
 res@cnMaxLevelValF   = 1.6	
 res@cnMinLevelValF   = -1.6
 res@cnLevelSpacingF = 0.2
 res@lbLabelStride      = 4 
 
 res@gsnLeftString = "a)"
 res@gsnCenterString            = " 2xCO~B~2~N~"
plot(0) = gsn_csm_contour(wks,cloud_slow(0,:,:),res)      

 res@gsnLeftString = "b)"
 res@gsnCenterString            = " 20xBC" 
plot(1) = gsn_csm_contour(wks,cloud_slow(1,:,:), res)

 res@gsnLeftString = "c)"
 res@gsnCenterString            = " 20xBC~N~-2xCO~B~2~N~"
 res@cnFillPalette    = cmap2  
plot(2) =  gsn_csm_contour(wks,cloud_slow(2,:,:), res)

;**************** VECTOR ******************
 vcres                           = True          
 vcres@gsnDraw                   = False
 vcres@gsnFrame                  = False
 vcres@vpWidthF = 1.5
 vcres@vpHeightF = 1.5
 vcres@vpClipOn     = False
 vcres@vcRefAnnoFontQuality = "High"
 vcres@vcVectorDrawOrder       = "PostDraw"       
 vcres@gsnAddCyclic            = False
 vcres@vcFillArrowsOn    = True
 vcres@vcFillArrowEdgeThicknessF = 3.0
 vcres@gsnRightString            = ""                
 vcres@gsnLeftString             = "" 
 vcres@vcMapDirection = False
 vcres@vcLineArrowHeadMaxSizeF = 0.20
 vcres@vcMinDistanceF  = 0.09          
 vcres@vcRefMagnitudeF           = 0.1     
 vcres@vcRefLengthF              = 0.01        
 vcres@vcLineArrowThicknessF   = 2.0              
 vcres@vcGlyphStyle            = "CurlyVector"    
 vcres@vcFillArrowHeadYF = 0.05
 vcres@vcLineArrowColor          = "Black"  
 vcres@vcRefAnnoOn = False
 vcres@vcRefAnnoOrthogonalPosF   = 0.1           
 vcres@vcRefAnnoParallelPosF     = 0.97            
 vcres@vcRefAnnoFontHeightF      = 0.085
 vcres@vcRefAnnoExplicitMagnitudeF = 0.2
 vcres@vcRefAnnoArrowMinOffsetF = 0.14
 vcres@vcRefAnnoArrowSpaceF = 1
 vcres@vcRefAnnoString1 = "0.1"
 vcres@vcRefAnnoString2 = "Reference vector"
 vectors = new(3,graphic)
 do i=0,2
	vectors(i) = gsn_csm_vector(wks,v(i,:,:),w(i,:,:),vcres)
    overlay(plot(i),vectors(i))
 end do
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
gsn_text_ndc(wks,"Cloud fraction change",0.53,0.94,txres)
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
panres@lbLabelFontHeightF    =  res@lbLabelFontHeightF/2
panres@lbTitleFontHeightF    = res@lbTitleFontHeightF/2 
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

txres               = True
txres@txFont   = FontTb  
txres@txFontHeightF = 0.026
   
panres@gsnPanelTop              = 0.50
panres@gsnPanelBottom = 0.10
gsn_panel(wks,plot(2),(/1,1/),panres) 

 

frame(wks)
 
 exit

end

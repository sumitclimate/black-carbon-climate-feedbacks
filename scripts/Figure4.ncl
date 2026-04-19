begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;DATA Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
psst_1co2 = "1xCO2_cam5_PSST_plev_file.nc"
som_1co2 = "1xCO2_cam5_SOM_plev_file.nc"

psst_2co2 = "2xCO2_cam5_PSST_plev_file.nc"
som_2co2 = "2xCO2_cam5_SOM_plev_file.nc"

psst_20BC = "20xBC_cam5_PSST_plev_file.nc"
som_20BC = "20xBC_cam5_SOM_plev_file.nc"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Kernel DATA Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fq = addfile("q.kernel.plev.nc","r")
ft = addfile("t.kernel.plev.nc","r")
fts = addfile("ts.kernel.nc","r")
fal = addfile("alb.kernel.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Std Dev Files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttest_file = addfile("surface_albedo_feedback_alpha_values.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Load  Data;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fp1 = addfile(psst_1co2,"r")

fp2 = addfile(psst_2co2,"r")

fpp2 = addfile(psst_20BC,"r")

fs1 = addfile(som_1co2,"r")
	
fs2 = addfile ( som_2co2,"r")

fss2 = addfile ( som_20BC,"r")

lat = fp1->lat
lon = fp1->lon
lev = fp1->lev_p
meta = fp2->alb
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
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       global mean surface temperature change


TS1 = fs1->TS
TS2= fs2->TS

dTS_co2 = TS2 - TS1
copy_VarMeta(TS2,dTS_co2)

tp1 = fp1->TS
tp2 = fp2->TS

dtp = tp2 - tp1
copy_VarMeta(tp2,dtp)

dTS_slow_co2 = dTS_co2 - dtp
copy_VarMeta(tp2,dTS_slow_co2)

surf_t_slow_co2 = dim_avg_n_Wrap(dTS_slow_co2,0)
dsurf_t_slow_co2 = wgt_areaave_Wrap(surf_t_slow_co2,area,1,0)


print(dsurf_t_slow_co2)


TS1 = fs1->TS
TS2= fss2->TS

dTS_BC = TS2 - TS1
copy_VarMeta(TS2,dTS_BC)

tp1 = fp1->TS
tp2 = fpp2->TS

dtp = tp2 - tp1
copy_VarMeta(tp2,dtp)

dTS_slow_BC = dTS_BC - dtp
copy_VarMeta(tp2,dTS_slow_BC)

surf_t_slow_BC = dim_avg_n_Wrap(dTS_slow_BC,0)
dsurf_t_slow_BC = wgt_areaave_Wrap(surf_t_slow_BC,area,1,0)


print(dsurf_t_slow_BC)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;           Slow Albedo feedback



alb1 = fs1->alb
alb2 = fs2->alb
dalb_tot = alb2 - alb1

copy_VarMeta(meta,dalb_tot)

alb3 = fp1->alb
alb4 = fp2->alb
dalb_fst_co2 = alb4 - alb3

copy_VarMeta(meta,dalb_fst_co2)

dalb_slow_co2 = dalb_tot - dalb_fst_co2
copy_VarMeta(meta,dalb_slow_co2)

alb1 = fs1->alb
alb2 = fss2->alb
dalb_tot = alb2 - alb1

copy_VarMeta(meta,dalb_tot)

alb3 = fp1->alb
alb4 = fpp2->alb
dalb_fst_BC = alb4 - alb3

copy_VarMeta(meta,dalb_fst_BC)

dalb_slow_BC = dalb_tot - dalb_fst_BC
copy_VarMeta(meta,dalb_slow_BC)


alb_resp_co2 = (dalb_slow_co2*100)/ dsurf_t_slow_co2

copy_VarMeta(meta,alb_resp_co2)
 
alb_resp_BC = (dalb_slow_BC*100)/ dsurf_t_slow_BC 

copy_VarMeta(meta,alb_resp_BC)

alb_resp_diff = alb_resp_BC - alb_resp_co2

copy_VarMeta(meta,alb_resp_diff)

;**************************************************************

alb_kernel = fal->FSNT
alb_kernel_clearsky = fal->FSNTC

alb_feed_co2 = (((dalb_slow_co2*100)*(alb_kernel)))/ dsurf_t_slow_co2    

alb_feed_co2_avg = dim_avg_n_Wrap(alb_feed_co2,0)

Albedo_feedback_slow_co2 = wgt_areaave_Wrap(alb_feed_co2_avg,area,1,0)

print(Albedo_feedback_slow_co2)


alb_feed_BC = (((dalb_slow_BC*100)*(alb_kernel)))/ dsurf_t_slow_BC   

alb_feed_BC_avg = dim_avg_n_Wrap(alb_feed_BC,0)

Albedo_feedback_slow_BC = wgt_areaave_Wrap(alb_feed_BC_avg,area,1,0)

print(Albedo_feedback_slow_BC)

alb_feed_diff = alb_feed_BC - alb_feed_co2

var1 = dim_avg_n_Wrap(alb_feed_co2,0) 
var2 = dim_avg_n_Wrap(alb_feed_BC,0)
var3 = dim_avg_n_Wrap(alb_feed_diff,0)
var4 = dim_avg_n_Wrap(alb_resp_co2,0)
var5 = dim_avg_n_Wrap(alb_resp_BC,0)
var6 = dim_avg_n_Wrap(alb_resp_diff,0)


var1!0 = "lat"
var1&lat = lat
var1!1 = "lon"
var1&lon = lon


var2!0 = "lat"
var2&lat = lat
var2!1 = "lon"
var2&lon = lon

var3!0 = "lat"
var3&lat = lat
var3!1 = "lon"
var3&lon = lon


;*************************************************
printVarSummary(var1)
printVarSummary(var4)

varResp = new((/3,dimsizes(lat)/),float)

varResp(0,:) = dim_avg_n_Wrap(var4,1)
varResp(1,:) = dim_avg_n_Wrap(var5,1)
varResp(2,:) = dim_avg_n_Wrap(var6,1)

varFeed = new((/3,dimsizes(lat)/), typeof(var1))
varFeed(0,:) = dim_avg_n_Wrap(var1,1)
varFeed(1,:) = dim_avg_n_Wrap(var2,1)
varFeed(2,:) = dim_avg_n_Wrap(var3,1)

 ;************************************* Plotting 

pltName    =  "Figure4"  
pltType    =  "png"                                           
plot       =  new (4,"graphic")
wks       =  gsn_open_wks(pltType, pltName)
pltType@wkWidth = 3000
pltType@wkHeight = 3000
cmap = read_colormap_file("BlueYellowRed")
colorPal =  cmap
cnSpace = 0.04
FontH = 0.015
FontT = "helvetica"
FontTb = "helvetica-bold"
;****************************** Resources              
 res = True
 res@gsnDraw         = False                  
 res@gsnFrame        = False                   
 res@tiMainOn        = False             
 res@lbLabelBarOn    = True                  
 res@lbTitleFont     = FontT    
 res@lbLabelFont     = FontT   
 res@lbLabelFontHeightF    = FontH+0.0035
 res@lbTitleFontHeightF    = FontH +0.007 
 res@pmLabelBarOrthogonalPosF = 0.01
 res@lbAutoManage       = False             
 res@lbTitlePosition    = "Bottom"      
 res@lbLabelStride      = 2 
 res@lbBoxEndCapStyle = "TriangleBothEnds"
 res@cnLabelBarEndStyle = "IncludeOuterBoxes"
 res@cnLinesOn       = True                   
 res@cnLineLabelsOn  = False                   
 res@cnLevelSelectionMode = "ManualLevels"                   
 res@cnFillOn         = True                 
 res@cnFillPalette    = colorPal           
 res@cnInfoLabelOn = False		
 res@cnFillOpacityF =0.9
 res@cnLineThicknessF     =  0.023              
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
 res@mpPerimOn = True
 res@mpPerimLineThicknessF = 2  
 res@gsnStringFontHeightF = FontH +0.01
 res@gsnStringFont = FontT
 res@lbLabelAutoStride = False
 res@gsnCenterStringFontColor = "black" 
 
 res@gsnCenterString            = "a)   Albedo Response ~C~     (20xBC - 2xCO~B~2~N~)";" (Solar-CO~B~2~N~)"
 res@lbTitleString    = "( % K~S~-1~N~)"
 res@gsnLeftString = ""
 res@gsnRightString = ""

 res@cnMinLevelValF   = -2.0
 res@cnMaxLevelValF   = 2.0
 res@cnLevelSpacingF = 0.2
 res@lbLabelStride = 4

plot(0) = gsn_csm_contour_map(wks, var6, res)

 ;***********************************************************************************
 delete(res@cnFillPalette)

 res@mpPerimOn         = True 	
 res@lbTitleString    = "(Wm~S~-2~N~K~S~-1~N~)"
 res@cnFillPalette = colorPal;(::-1,:)

 res@cnMinLevelValF   = -2.0
 res@cnMaxLevelValF   = 2.0
 res@cnLevelSpacingF = 0.2
 res@lbLabelStride = 4
 
 res@gsnCenterString = "b)   Albedo Feedback ~C~     (20xBC - 2xCO~B~2~N~)  "
 res@gsnCenterStringOrthogonalPosF = 0.02
 res@mpPerimOn = True

plot(1) = gsn_csm_contour_map(wks,var3, res)
 ;***********************************************************************************
 alpha_res = ttest_file->alpha_res
 alpha_feed = ttest_file->alpha_feed
 
alpha_feed!0 = "lat"
alpha_feed!1 = "lon"

alpha_feed&lat = lat
alpha_feed&lon = lon
copy_VarMeta(alpha_feed,alpha_res)

sgres                      = True		
sgres@gsnDraw              = False		
sgres@gsnFrame             = False		
sgres@cnInfoLabelOn        = False		
sgres@cnLinesOn            = False		
sgres@cnLineLabelsOn       = False		
sgres@cnFillScaleF         = 0.7	
sgres@cnFillDotSizeF       = 0.0022
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

sig_plot1 = gsn_csm_contour(wks,alpha_res,sgres)
sig_plot2 = gsn_csm_contour(wks,alpha_feed,sgres)

opt                  = True
opt@gsnShadeFillType = "pattern"
opt@gsnShadeHigh     = 17

sig_plot1 = gsn_contour_shade(sig_plot1,-999.,95.0,opt)
sig_plot2 = gsn_contour_shade(sig_plot2,-999.,95.0,opt)

overlay(plot(0),sig_plot1)
overlay(plot(1),sig_plot2)

 ;***********************************************************************************

txres = True
txres@txFontHeightF = FontH+0.008
txres@txFont        = FontT
txres@txFontColor = "black"
txres@amOrthogonalPosF = -1.93 
txres@amParallelPosF = -1.93 
label = (/"60~S~o~N~S","30~S~o~N~S","Eq","30~S~o~N~N","60~S~o~N~N"/)
text = gsn_add_text(wks,plot(0),label,(/-180,-180,-180,-180,-180/),(/-60,-30,0,30,60/),txres)
 ;***********************************************************************************

 panres                        = True
 panres@pmTickMarkDisplayMode   = "Always"       
 panres@gsnFrame         = False    
 panres@gsnDraw         = True 
 panres@lbLabelBarOn  = False
 panres@lbTitleFont                = res@lbTitleFont   
 panres@lbLabelFont                = res@lbLabelFont    
 panres@lbLabelFontHeightF    =  res@lbLabelFontHeightF
 panres@lbBoxLinesOn  = True
 panres@lbLabelAutoStride = False
 panres@lbLabelStride = res@lbLabelStride
 panres@lbAutoManage          = False            
 panres@cnInfoLabelOn = False		
 panres@gsnPanelLabelBar    = True           
 panres@lbOrientation            = "horizontal"
 panres@pmLabelBarHeightF = 0.05
 panres@pmLabelBarWidthF = 0.5
 panres@pmLabelBarOrthogonalPosF = 0.01
 panres@pmLabelBarParallelPosF = -0.2
 panres@gsnPanelLeft  = 0.05
 panres@gsnPanelRight  = 0.97
 panres@gsnMaximize            = True 
 panres@gsnPanelBottom           = 0.45
gsn_panel(wks,plot(0:1),(/1,2/),panres) 

;**************************************************************************************************
  res2                        = True            
 thicknesses                 = (/  4.00, 4.00, 3.0/)     
 colors                      = (/ "black","red","black"/) 
 xyDashPatterns = (/0, 0, 1/)
 font = 21
 ;                              XY axis
 ;************************************************************************************   
 res2@tiYAxisString    = "Response"
 res2@tiXAxisString    = ""
 res2@tiXAxisFont = FontT
 res2@tiYAxisFont =FontT
 res2@tiXAxisFontHeightF = FontH+0.01
 res2@tiYAxisFontHeightF = FontH+0.01
 res2@tiMainFontHeightF  = FontH
 res2@tiYAxisFontColor = "black" 
 res2@tmXBLabelFont = FontT
 res2@tmYLLabelFont = FontT
 res2@tmXBLabelFontHeightF = FontH+0.005
 res2@tmYLLabelFontHeightF = FontH+0.005
 res2@tmYROn = False         
 res2@tmXTOn = False       
 res2@gsnDraw          = False                   
 res2@gsnFrame         = False  
 
;*************** grid resources
 res2@tmXMajorGrid                = True        
 res2@tmXMajorGridThicknessF      = 2.0          
 res2@tmXMajorGridLineDashPattern = 2           
 res2@tmXMajorGridLineColor      ="grey"
 
 res2@tmYMajorGrid                = True       
 res2@tmYMajorGridThicknessF      = 2.0        
 res2@tmYMajorGridLineDashPattern = 2           
 res2@tmYMajorGridLineColor      ="grey"
 res2@tmBorderThicknessF = 2.0

 res2@gsnStringFont = FontT
 res2@gsnCenterStringFontColor = "black"
 res2@gsnStringFontHeightF = FontH +0.012

;                                   ASPECT RATIO
 ;*********************************************************************************  

 res2@vpWidthF         = 0.5                   
 res2@vpHeightF        = 0.4
 res2@tiYAxisString    = "( % K~S~-1~N~)"
 res2@tiXAxisString    = ""
  
 res2@trYMinF     = -6.0
 res2@trYMaxF     =  6.0		
  ;                                   LEGENDS
 ;************************************************************************************   
 res2@gsnCenterString   = "c)  Zonal average albedo ~C~             response     "
 res2@xyLineThicknesses = thicknesses(::-1)   
 res2@xyDashPatterns  = xyDashPatterns(::-1)
 res2@xyLineColors           = colors(::-1)

 ;***********************************************************************************    
 plot(2) = gsn_csm_xy(wks,lat,varResp(::-1,:),res2) 

 res2@gsnCenterString   = "d)  Zonal average albedo ~C~             feedback  "

 res2@tiYAxisString    = "   (Wm~S~-2~N~K~S~-1~N~)"
 res2@tiXAxisString    = ""
 
 plot(3) = gsn_csm_xy(wks,lat,varFeed(::-1,:),res2) 
   
 labels                      = (/"  2xCO~B~2~N~","  20xBC","  20xBC - 2xCO~B~2~N~"/)           

 lgres2                    = True
 lgres2@vpWidthF           = 0.12                   
 lgres2@vpHeightF          = 0.08                  
 lgres2@lgPerimOn         = False 
 lgres2@lgDashIndexes      = xyDashPatterns(::-1)       
 lgres2@lgLineLabelStrings = labels(::-1)     
 lgres2@lgLineColors       = colors(::-1)
 lgres2@lgPerimOn          = False 
 lgres2@lgLineThicknessF  =  thicknesses(::-1)
 lgres2@lgLabelFontHeightF     = FontH+0.1               
  lgres2@lgLabelFont = FontTb
 gsn_legend_ndc(wks,3,labels(::-1),0.22,0.45,lgres2)
  
;                                   Panel Plot
 ;*********************************************************************************  

 pres2                        = True
 pres2@gsnMaximize            = True            
 pres2@pmTickMarkDisplayMode   = "Always"       
 pres2@gsnFrame         = False  
 pres2@gsnDraw         = True    
 pres2@lbBoxLinesOn  = False
 pres2@gsnPanelLeft  = 0.03
 pres2@gsnPanelRight  = 0.95
 pres2@gsnPanelBottom           = 0.03
 pres2@gsnPanelTop              = 0.7
 panres@gsnAttachPlotsYAxis = True
 pres = True
gsn_panel(wks,plot(2:3),(/1,2/),pres2) 

frame(wks)
end



;This script will plot some figures for the my ERC proposal

FILTER=True 
HIGH_PASS=False
RUN_AVG=True
PLOT_COMPONENTS=True
ROUGHPLOTS=False

if FILTER then
  ;set up 15-30 band-pass filter (for annual mean data)
  fc_15   = 1/10.                        ; frequency of 15 years (in months)
  fc_30   = 1/50.                       ; frequency of 30 years (in months)
  sigma = 1.0                            ; Lanczos sigma
  nwgt  = 41                             ; loose 10 years at each end
  if HIGH_PASS then
    ihp   = 1                             ; choose band_pass
    wgts  = filwgts_lanczos (nwgt, ihp, fc_30, -999., sigma )
  else
    ihp   = 2                             ; choose band_pass
    wgts  = filwgts_lanczos (nwgt, ihp, fc_30, fc_15, sigma )
  end if
end if

;Load in the last millennium reanalysis global means (w.r.t. 1951-1981)
lmr_f=addfile("/data/aod/obs/LastMillenniumReanalysis/gmt_MCruns_ensemble_full.nc", "r")
lmr_gmt_all=lmr_f->gmt;load in all global mean timeseries
lmr_gmt_raw=dim_avg_Wrap(lmr_gmt_all);construction means of MC runs
if FILTER then
  ;apply band-pass filtering to LMR
  lmr_gmt=wgt_runave_n_Wrap(lmr_gmt_raw,wgts,1,0)
  ;average across the 20 different Monte-Carlo runs to get the best estimates
  lmr_gmt_be=dim_avg_n_Wrap(lmr_gmt,1)
  else
  ;average across the 20 different Monte-Carlo runs to get the best estimates
  lmr_gmt_be=dim_avg_n_Wrap(lmr_gmt_raw,1)
end if
;create samples of averaged quantities at 15-30 years
lmr_samples=lmr_gmt_be(862::15)
do i=0,dimsizes(lmr_samples)-1
  these_yrs=ispan(850+i*15,850+i*15+14,1)
  ;print(" "+i+": "+min(these_yrs)+"-"+max(these_yrs))
  if any(ismissing(lmr_gmt_be(these_yrs))) then
    lmr_samples(i)=lmr_samples@_FillValue
  else
    lmr_samples(i)=(/avg(lmr_gmt_be(these_yrs))/)
  end if
end do    
if ROUGHPLOTS then
  wks = gsn_open_wks("x11","Foo")
  plot = gsn_csm_xy(wks,ispan(1,dimsizes(lmr_gmt_be),1),lmr_gmt_be,False)
  delete(wks)
end if

;Load in the last millennium simulations
cesm1_lm_f=addfile("data\CESM1_LastMillenniumEnsemble.001-010.tas_global_avg_mon.850-2005.nc", "r")
lm_gmt_mon=cesm1_lm_f->tas_global_avg_mon;load in all global mean timeseries
lm_gmt_raw=month_to_annual(lm_gmt_mon(time|:,record|:),1);construction annual means from monthly ones
;conver to anomalies w.r.t. 1951-1980
lm_gmt_raw=lm_gmt_raw-conform(lm_gmt_raw,dim_avg_n(lm_gmt_raw(1951-850:1980-850,:),0),1)
if FILTER then
  ;apply band-pass filtering to Last Millennium Runs
  lm_gmt=wgt_runave_n_Wrap(lm_gmt_raw,wgts,1,0)
 else
  lm_gmt=lm_gmt_raw
end if
lm_samples=lm_gmt(12::15,:)
do i=0,dimsizes(lm_samples(:,0))-1
  these_yrs=ispan(i*15,i*15+14,1)
  ;print(" "+i+": "+min(these_yrs)+"-"+max(these_yrs))
  if any(ismissing(lm_gmt(these_yrs,:))) then
    lm_samples(i,:)=lm_samples@_FillValue
  else
    lm_samples(i,:)=(/dim_avg_n(lm_gmt(these_yrs,:),0)/)
  end if
end do    
if ROUGHPLOTS then
  wks = gsn_open_wks("x11","Foo")
  plot = gsn_csm_xy(wks,ispan(1,dimsizes(lm_gmt),1),lm_gmt,False)
  delete(wks)
end if
if PLOT_COMPONENTS then
  wks = gsn_open_wks("pdf","past2soon.LastMillennium_timeseries")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@tiXAxisString="Year (C.E.)"
  res@tiYAxisString="Temperature Anomaly (K)"
  res@xyMonoLineThickness=True
  res@xyLineThickness=0.5
  res@xyMonoLineColors=True
  res@xyLineColor="gray"
  res@trXMaxF=2000.
  res@gsnYRefLine=0.
  res@xyMonoDashPattern = True
  res@xyDashPattern = 0.
  plot = gsn_csm_xy(wks,ispan(850,850+dimsizes(lm_gmt_raw(:,0))-1,1),(/lm_gmt_raw(record|:,year|:)/),res)
  res@xyLineThickness=2.
  res@xyLineColor="black"
  oplot = gsn_csm_xy(wks,ispan(0,dimsizes(dim_avg_n(lmr_gmt_raw,1))-1,1),dim_avg_n(lmr_gmt_raw,1),res)
  overlay(plot, oplot)
  draw(plot)
  frame(wks)
  delete([/wks,plot,oplot,res/])
end if

print("The variance of the LMR is "+stddev(lmr_samples(0:75))^2)
print("The variance of the LM is "+stddev(lm_samples(0:75,:))^2)
scaling_ratio=(stddev(lmr_samples)^2)/(stddev(lm_samples)^2)
print("Therefore scaling ratio of CESM * "+scaling_ratio)

if PLOT_COMPONENTS then
  array_lm=ndtooned(lm_samples(0:75,:))
  array_lmr=ndtooned(conform(lm_samples(0:75,:),lmr_samples(0:75),0))
  wks = gsn_open_wks("pdf","past2soon.histo")              ; send graphics to PNG file
  gsn_define_colormap(wks,"temp1")
  res                            = True          ; plot mods desired
  res@gsnHistogramCompare        = True
  res@gsnHistogramComputePercentages = True
  res@gsnHistogramNumberOfBins=25
  ;res@gsnHistogramBinIntervals=fspan(.6, 1.4,41)
  res@gsnHistogramBarWidthPercent=75
  plot = gsn_histogram(wks,(/array_lmr,array_lm/),res)
  delete([/wks,plot,res/])
end if
;Load in the Large Ensemble simulations. Need to take early part and late part, convert them back to absolute valeus and then stitch together
cesm1_lens_f1=addfile("data\CESM1-LargeEnsemble_all40.tas_global_avg_mon.1920-2017.nc","r")
lens_gmt_mon1=cesm1_lens_f1->tas_global_avg_mon;load in all global mean timeseries
lens_gmt_raw1=month_to_annual(lens_gmt_mon1(time|:,record|:),1);construction means of MC runs
lens_offset1_in=cesm1_lens_f1->tas_spatialmean_ann;load in spatial patterns
lat_wgts=NormCosWgtGlobe(lens_offset1_in&lat)
lens_offset1=wgt_areaave_Wrap(lens_offset1_in, lat_wgts,1.,1)
lens_gmt_raw1=lens_gmt_raw1+conform(lens_gmt_raw1, lens_offset1, 1)
cesm1_lens_f2=addfile("data\CESM1-LargeEnsemble_all40.tas_global_avg_mon.2018-2047.nc", "r")
lens_gmt_mon2=cesm1_lens_f2->tas_global_avg_mon;load in all global mean timeseries
lens_gmt_raw2=month_to_annual(lens_gmt_mon2(time|:,record|:),1);construction means of MC runs
lens_offset2_in=cesm1_lens_f2->tas_spatialmean_ann;load in spatial patterns
lens_offset2=wgt_areaave_Wrap(lens_offset2_in, lat_wgts,1.,1)
lens_gmt_raw2=lens_gmt_raw2+conform(lens_gmt_raw2, lens_offset2, 1)
;combine the two datasources to make coherent timeseries
lens_gmt_raw_combined=new((/98+83,40/), typeof(lens_gmt_raw2))
lens_gmt_raw_combined(:97,:)=lens_gmt_raw1
lens_gmt_raw_combined(98:,:)=lens_gmt_raw2
if ROUGHPLOTS then
  wks = gsn_open_wks("x11","Foo")
  plot = gsn_csm_xy(wks,ispan(1920,2100,1),lens_gmt_raw_combined(record|:,year|:),False)
  delete(wks)
end if
if FILTER then
  lens_notscaled=lens_gmt_raw_combined
  lens_gmt=lens_gmt_raw_combined
  ;apply band-pass filtering to Large Ensemble Runs
  lens_gmt=wgt_runave_n_Wrap(lens_gmt_raw_combined,wgts,1,0)
  lens_notscaled=lens_gmt_raw_combined-lens_gmt
  ;rebase w.r.t. 2013-2018 average
  lens_gmt=lens_gmt-conform(lens_gmt,dim_avg_n(lens_gmt(93:98,:),0),1)
  lens_notscaled=lens_notscaled-conform(lens_notscaled,dim_avg_n(lens_notscaled(93:98,:),0),1)
else
  lens_gmt=lens_gmt_raw_combined
  ;rebase w.r.t. 2013-2018 average
  lens_gmt=lens_gmt-conform(lens_gmt,dim_avg_n(lens_gmt(93:98,:),0),1)
end if
if ROUGHPLOTS then
  wks = gsn_open_wks("x11","Foo")
  plot = gsn_csm_xy(wks,ispan(2013,2047,1),lens_gmt(record|:,year|93:127),False)
  delete(wks)
end if
pred_20322047=dim_avg_n_Wrap(lens_gmt(112:127,:),0)
pred_20152018=dim_avg_n_Wrap(lens_gmt(93:98,:),0)
pred_20322047_mn=avg(pred_20322047)
pred_20322047_scaled=scaling_ratio*(pred_20322047-conform(pred_20322047,pred_20322047_mn,-1))^2+pred_20322047_mn
scaling_wTime=conform(lens_gmt(record|:,year|93:127),(pred_20322047_scaled-pred_20152018)/29.,0)*\
              conform(lens_gmt(record|:,year|93:127),ispan(-2,32,1),1)
total_pred=lens_gmt(record|:,year|93:127)
orig_pred=lens_gmt(record|:,year|93:127)
if FILTER then
  orig_pred=lens_gmt(record|:,year|93:127)+lens_notscaled(record|:,year|93:127)
  total_pred=lens_gmt(record|:,year|93:127)*scaling_wTime+lens_notscaled(record|:,year|93:127)
else
  orig_pred=lens_gmt(record|:,year|93:127)
  total_pred=lens_gmt(record|:,year|93:127)*scaling_wTime
end if
if PLOT_COMPONENTS then
  wks = gsn_open_wks("pdf","past2soon.Total_pred")
  res=True   
  res@tiMainString = " "
  res@tiXAxisString          = "Year"
  res@tiYAxisString          = "Global Mean Temp."
  res@xyDashPatterns         = (/0.,0.,0.,0.,0.,0./)
  res@xyLineThicknesses      = (/3.,3.,2.,2.,2.,2./)
  res@xyLineColors         = (/"Red","Blue","Red","Red","Blue","Blue"/)
  res@trXMaxF=2045
  res@trXMinF=2015
  array=orig_pred(0:5,:)
  if RUN_AVG then
    orig_pred=runave(orig_pred, 5, 0)
    total_pred=runave(total_pred, 5, 0)
  end if
  array(0,:)=dim_avg_n(orig_pred,0)
  array(1,:)=dim_avg_n(total_pred,0)
  array(2,:)=dim_avg_n(orig_pred,0)+2*dim_stddev_n(orig_pred,0)
  array(3,:)=dim_avg_n(orig_pred,0)-2*dim_stddev_n(orig_pred,0)
  array(4,:)=dim_avg_n(total_pred,0)+2*dim_stddev_n(total_pred,0)
  array(5,:)=dim_avg_n(total_pred,0)-2*dim_stddev_n(total_pred,0)
  plot = gsn_csm_xy(wks,ispan(2013,2047,1),array,res)
  delete(wks)
end if
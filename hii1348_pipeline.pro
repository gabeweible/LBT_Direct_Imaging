pro HII1348_pipeline, rho=rho, theta=theta, planet_x=planet_x, planet_y=planet_y, contrast=contrast, pre_inj_stuff=pre_inj_stuff, neg_inj=neg_inj, trial=trial, headless=headless, outpath=outpath, coadd=coadd, use_gauss=use_gauss, uncert=uncert

compile_opt IDL2; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
newline = string(10B)
start_time = systime(/JULIAN); Get the current time in a Julian date format (number of days since January 1, 4713 BC, with decimals)

;------------------------------[ Start User Input ]---------------------------------

; General/Combine Parameters
obj = 'HII1348'
raw_path = '/Users/gabeweible/OneDrive/research/HII1348/kevin/raw';'/Users/gabe/Data/HII1348/raw' ;'~/OneDrive/research/HII1348/HII1348/raw'
cube_start_frame = 0
output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_' + strcompress(coadd, /r) + '/';'/Users/gabe/reduction/testing_coadd_25/' ;'~/OneDrive/research/HII1349/testing_coadd_5/'
if keyword_set(outpath) then output_path=outpath

; Sky sub parameters
allowed_coadds = [1, 2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 200]
if total(allowed_coadds eq coadd) ne 1 then begin
  print, 'Error: invalid coadd value (must be a factor of 200)'
  stop
endif
fs_start = 400 / coadd & fs_end = 599 / coadd
print, 'Coadd:', coadd

; Centering parameters
do_block_right = 1
do_block_left = 1
do_block_bottom = 1
do_block_top = 1
do_cen_filter = 1
do_cen_median_sub = 0

; Clean parameters
; 0.993 is for KLIP, 0.995 is for ADI
corr_thresh = 0.994 ; Split the difference

; Planet injection parameters
;use_gauss can be set with a keyword
pxscale=0.0107 ;arcsec/pixel
; Passed values from keyword arguments (for contrast curve generation or astrometry/photometry)
if (keyword_set(rho) and keyword_set(theta)) or (keyword_set(planet_x) and keyword_set(planet_y)) then begin
   use_injection=1
endif else begin
   use_injeciton=0
endelse
if keyword_set(use_gauss) then use_gauss=[use_gauss]
if keyword_set(rho) then planet_theta=[theta]
if keyword_set(rho) then planet_r=[rho]
if keyword_set(planet_x) then planet_x=[planet_x]
if keyword_set(planet_y) then planet_y=[planet_y]
if keyword_set(contrast) then begin
	planet_contrast=[contrast]
	if keyword_set(rho) then n_planets=n_elements(planet_r)
	if keyword_set(planet_x) then n_planets=n_elements(planet_x)
endif


; Rotate/KLIP parameters (after testing/refinement)
silent=1; Don't print so much in adi.pro "Rotating by ..."
norm = 1; Normalize ADI frames to reduce residuals around the star
fs = 0; Run find_sources within ADI.pro
if keyword_set(rho) or keyword_set(planet_x) then use_injection=1 else use_injection = 0
do_destripe = 1
filter = 17.
bin = 3
bin_type = 'mean'
do_hyper = 0
if keyword_set(rho) or keyword_set(planet_x) then do_annmode=1 else do_annmode=0
combine_type = 'nwadi'
klip_fraction = 0
klip_start_frame = 75
klip_end_frame = 150
fill = 0
k_klip = 7
angsep= 1.
anglemax = 360.
nrings = 4.
n_ang = 2
annmode_inout = [83, 143]; Annulus inner and outer radii already setup for HII 1348 b; In px?

; Astrometry Parameters (Dewarp solns from https://scienceops.lbto.org/lbti/data-retrieval-reduction-publication/distortion-correction-and-astrometric-solution/)
Kx_dx = [[-1.13034544e+01,   1.45852226e-02,  -1.13372175e-05,   1.32602063e-09],
[ 1.03220943e+00,  -1.93058980e-05,   1.55798844e-08,  -3.86115281e-12],
[-2.57352199e-05,   2.70371257e-09,  -6.62650011e-12,   3.25017078e-15],
[ 8.02004325e-09,  -5.59037685e-13,   1.60256679e-15,  -8.18749145e-19]]

Ky_dx = [[-9.37023860e-01,   9.89973161e-01,  -4.32284634e-06,   7.08000564e-09],
[-1.59438718e-02,  -1.95099497e-05,   1.55801035e-09,   2.13743170e-13],
[ 9.34251811e-06,   1.26473736e-08,  -3.71968684e-12,   5.88384784e-16],
[ 3.41071678e-10,  -1.82060569e-12,   1.59690189e-15,  -3.14810895e-19]]

Kx_sx = [[-7.96016109e+00,   9.26584096e-03,  -7.93676069e-06,   5.13414639e-10],     [ 1.02925896e+00,  -1.59974177e-05,   9.57769272e-09,  -1.14409822e-12],
[-2.30169348e-05,  -3.45351550e-09,   1.89621010e-12,   6.72971750e-17],
[ 7.07041647e-09,   2.20511200e-12,  -1.66433082e-15,   2.55463711e-19]]

Ky_sx = [[-2.26409123e+00,   9.93175401e-01,  -6.67169688e-06,   8.08275391e-09],
[-1.38521740e-02,  -2.27910031e-05,   4.72367613e-09,  -1.64738716e-12],
[ 8.17060299e-06,   1.35240460e-08,  -5.52374318e-12,   2.14966954e-15],
[ 9.25982725e-10,  -2.50607186e-12,   2.58675626e-15,  -9.82463036e-19]]


;------------------------------[ End User Input ]---------------------------------

if keyword_set(planet_x) then begin
   planet_r = sqrt((250.-planet_x)^2 + (250.-planet_y)^2)
endif

if use_injection and not neg_inj then begin; Custom annulus inner and outer radii for artificial injections at arbitrary locations
   annmode=1
   ;Thicc rings
   annmode_inout=round([max([0.,planet_r/pxscale-25.+1.]),planet_r/pxscale+25.+2.])
   if annmode_inout[1] gt 16 then annmode_inout=round([max([0.,planet_r/pxscale-30.+1.]),planet_r/pxscale+30.+2.])
   ; Thinn rings
   ;annmode_inout=round([max([0.,planet_r/pxscale-10.+1.]),planet_r/pxscale+10.+2.])
   ;if annmode_inout[1] gt 16 then annmode_inout=round([max([0.,planet_r/pxscale-12.+1.]),planet_r/pxscale+12.+2.])
endif

wr=fix(250./float(nrings))
suffix=strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(ANGSEP,4))+'as_'+string(sigfig(anglemax,2))+'am_'+string(sigfig(nrings,2))+'rings_'+string(wr)+'wr_'+String(n_ang)+'nang_'+string(sigfig(filter,2))+'filter_'+string(sigfig(bin,2))+'bin'+string(corr_thresh)+'corrthresh'), /remove_all)
if use_injection then suffix=strcompress(suffix+'_inj',/rem)
if neg_inj then suffix=strcompress(suffix+'neg_inj',/rem)
if uncert then suffix=strcompress(suffix+'uncert',/rem)
;if neg_inj eq 1 then begin & annmode=0 & do_annmode=0 & endif

;if keyword_set(headless) then begin
;   if headless eq 1 then fs = 0
;   if headless eq 0 then fs = 1
;   if headless ne 1 and headless ne 0 then begin & print, 'Error: Invalid headless keyword value' & stop & endif
;   print, 'Headless:', string(headless), 'fs:', fs
;endif; headless if

;------------------------------[ Data Reduction ]---------------------------------

if pre_inj_stuff eq 1 then begin
   
   print, 'pre_inj_stuff:', pre_inj_stuff

   create_cube, obj, raw_path, cube_start_frame, coadd, output_path
   ;
   bad_pixels, output_path, obj
   ;
   sky_sub, obj, coadd, output_path, fs_start=fs_start, fs_end=fs_end
   
   dewarp, obj, output_path
   
   split, obj, output_path
   ;
   center, obj, output_path, do_block_right, do_block_left, do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub
   ;
   clean, obj, output_path, corr_thresh, do_cen_filter

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets, planet_contrast, pxscale, corr_thresh, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta, use_gauss=use_gauss, silent=silent
   if keyword_set(planet_x) then inject_planets, obj, output_path, n_planets, planet_contrast, pxscale, corr_thresh, do_cen_filter, planet_y=planet_y, planet_x=planet_x, use_gauss=use_gauss, silent=silent
   
   ;Change output folder manually in rotate.pro and klip.pro !!!!!!!!!!!!!!! (right now as long as it's testing_coadd_ whatever it's fine)
   ; I'm having trouble with find_sources here. (Everything is working now, but note that I might need to adjust the correction
   ; factor to get acurate values)
   ;klip, obj, output_path, use_injection, do_destripe, filter, bin, bin_type, do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame, klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout, suffix, corr_thresh, do_cen_filter, coadd, rho=rho, theta=theta, contrast=contrast, fs=fs, neg_inj=neg_inj
   
   adi, obj, output_path, use_injection, do_destripe, filter, suffix, corr_thresh, do_cen_filter, coadd, fs=fs, neg_inj=neg_inj,norm=norm, uncert=uncert, silent=silent
   
endif
if pre_inj_stuff eq 0 then begin
   
   print, 'pre_inj_stuff:', pre_inj_stuff

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets, planet_contrast, pxscale, corr_thresh, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta, use_gauss=use_gauss, silent=silent
   if keyword_set(planet_x) then begin
      inject_planets, obj, output_path, n_planets, planet_contrast, pxscale, corr_thresh, do_cen_filter, planet_y=planet_y, planet_x=planet_x, use_gauss=use_gauss
   endif

   ;klip, obj, output_path, use_injection, do_destripe, filter, bin, bin_type, do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame, klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout, suffix, corr_thresh, do_cen_filter, coadd, rho=rho, theta=theta, contrast=contrast, trial=trial, fs=fs, neg_inj=neg_inj

   adi, obj, output_path, use_injection, do_destripe, filter, suffix, corr_thresh, do_cen_filter, coadd, fs=fs, neg_inj=neg_inj,norm=norm, uncert=uncert, silent=silent
   
endif

;-----------------------------------[ El Fin ]--------------------------------------

print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'; 1440 minutes per day (Julian dates are measured in days)

end; That's all, folks!
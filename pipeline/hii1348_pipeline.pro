pro HII1348_pipeline, rho=rho, theta=theta, planet_x=planet_x, planet_y=planet_y,$
	contrast=contrast, pre_inj=pre_inj, neg_inj=neg_inj, trial=trial,$
   outpath=outpath, coadd=coadd, use_gauss=use_gauss, uncert=uncert, klip=klip,$
   fs=fs, two_soln_override=two_soln_override, extra=extra, nod=nod, cds=cds,$
   corr_thresh=corr_thresh, pre_clean=pre_clean, k_klip=k_klip, adi=adi, do_annmode=do_annmode,$
   szbin=szbin, angsep=angsep, bin_type=bin_type, filter=filter, pa_off=pa_off,$
	n_ang=n_ang,combine_type=combine_type,sky_type=sky_type,full_frame=full_frame

; hii1348_pipeline, pre_inj=1, neg_inj=0, coadd=25, uncert=0, klip=1 for normal use

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

;------------------------------[ Start User Input ]---------------------------------

; General/Combine Parameters
obj = 'HII1348'
band = 'Lp'; Std-L filter (3.70 microns)
annular_clean=0; bloc < 10 px and > 50 px (focus on cleaning based on 1st airy ring)
centered_clean=1
do_second_center=0
wl=3.7002 ; wavelength in microns
raw_path = '/Volumes/T7/hii1348_raw'
;'/Users/gabeweible/OneDrive/research/HII1348/kevin/raw'
; fwhm calculated with fullwid_halfmax(pupil_array read-in with readfits, /GAUSSIAN_FIT)
fwhm = 9.55415; Average of dith1 and dith2 SX and DX fwhm in x and in y (from pupil medians)
cube_start_frame = 0
coadd = 25
magnify = 1; Magnify
output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/HII1348/macbook_'+$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)
output_path += '/'

aperture='both'

;full_frame=0

; Passed values from keyword arguments (for contrast curve generation or astrometry/photometry)
if keyword_set(neg_inj) then begin
	if (neg_inj eq 1) then n_planets=1
endif else neg_inj = 0

if not keyword_set(use_gauss) then use_gauss=0; for cADI PSF model - don't use by default.
if not keyword_set(pre_clean) then pre_clean=0

; Planet injection parameters
; use_gauss can be set with a keyword
pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared)
pxscale_dx = 0.010700 ; +0.000042 or -0.000051 arcsec/pixel (from Steve and Jared)

if not keyword_set(pa_off) then pa_off=0.; let's not automatically set this.
;pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared)
pxscale_dx = 0.010610 ; +0.000042 or -0.000051 arcsec/pixel (from Steve and Jared)
truenorth_sx = -1.278-pa_off ; +0.131 or -0.225 deg (from Steve and Jared)
truenorth_dx = 1.001-pa_off ; + 0.254 or -0.237 deg (from Steve and Jared)

if keyword_set(outpath) then output_path=outpath

; Sky sub parameters
;allowed_coadds = [1, 2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 200]
;if total(allowed_coadds eq coadd) ne 1 then begin
;  print, 'Error: invalid coadd value (must be a factor of 200)'
;  stop
;endif
if not keyword_set(sky_type) then sky_type='median'
fs_start = 400 / coadd & fs_end = 599 / coadd
nnod = fix(200/coadd/2.)
;print, 'Coadd:', coadd

; Centering parameters
do_block_right = 1
do_block_left = 1
do_block_bottom = 1
do_block_top = 1
do_cen_filter = 1
do_cen_median_sub = 0

; CDS params
cds = 'midcds' ; second - first (all I have for HII 1348)

; Clean parameters
; 0.993 is for KLIP, 0.995 is for ADI
corr_thresh = 0.92;0.994 ; Split the difference

; Passed values from keyword arguments (for contrast curve generation or astrometry/photometry)
if neg_inj eq 1 then n_planets=1
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
if not keyword_set(nod) then nod='total'

; ADI/KLIP parameters (after testing/refinement)
silent = 0; Don't print so much in adi.pro "Rotating by ..."
normal = 1; Normalize ADI frames to reduce residuals around the star

;fs = 1; Run find_sources within ADI.pro; Set in a kwarg now...
if not keyword_set(fs) then fs = 0; Assume it's NOT fine?

if keyword_set(rho) or keyword_set(planet_x) then begin
	use_injection=1
	do_annmode=1
endif else begin
	use_injection=0
	;do_annmode=0
endelse

; Bad 
; bad pixel parameters
sigma_threshold=3.0 ; sigma_threshold*stdev outside of surrounding values are out!
boxhsize=8; 2*boxhsize x 2*boxhsize block of pixels for each comparisons params
stripe = 'center1024'

do_destripe = 1
filter = 30.;17.

bin = 2.
bin_type = 'median';'mean'
do_hyper = 0
combine_type = 'median';'nwadi'

; Start and end frame if only using KLIP on a fraction of dataset (start and end)
klip_fraction = 0; 1 to KLIP only a fraction, 0 to KLIP everything.
klip_start_frame = 75
klip_end_frame = 150

fill = 0
if not keyword_set(k_klip) then k_klip=7
if keyword_set(k_klip) then k_klip=k_klip
;k_klip = k_klip ; Number of KLIP components (eigenimages) to keep


; angsep is relevant for annuli - each annulus will have a unique
; parallactic-angle rotation constraint.
angsep= 0.8;1. ; # of (1.5x)FWHM  parallactic rotation required to use frames for KLIP-ADI
anglemax = 360.
nrings = 1.;4. ; Number of annuli for KLIP
n_ang = 2.;2. ; Number of segments in each annulus for KLIP

align_radius_thresh = 50

peak_clean=1	;reject frames where the peak flux is beneath a certain
;fraction compared to the average peak value
peak_thresh=0.90	;normalized units, so 1 is the average peak strength

stddev_clean=1
stddev_thresh=1.25; max. ratio of pixel stddev to median pixel stddev

; Annulus inner and outer radii already setup for HII 1348 b; In px?
annmode_inout_sx = [3,18];[83, 143]
; Annulus inner and outer radii already setup for HII 1348 b; In px?
annmode_inout_dx = [3,18];[83, 143]

; Astrometry Parameters (Dewarp solns from 
; https://scienceops.lbto.org/lbti/data-retrieval-reduction-publication/distortion-correction-and-astrometric-solution/)
Kx_dx = [[-1.13034544e+01, 1.45852226e-02, -1.13372175e-05, 1.32602063e-09],$
 	  		[1.03220943e+00, -1.93058980e-05, 1.55798844e-08, -3.86115281e-12],$
 	  		[-2.57352199e-05, 2.70371257e-09, -6.62650011e-12, 3.25017078e-15],$
 	  		[8.02004325e-09, -5.59037685e-13, 1.60256679e-15, -8.18749145e-19]]

Ky_dx = [[-9.37023860e-01, 9.89973161e-01, -4.32284634e-06, 7.08000564e-09],$
 	  		[-1.59438718e-02, -1.95099497e-05, 1.55801035e-09, 2.13743170e-13],$
 	  		[9.34251811e-06, 1.26473736e-08, -3.71968684e-12, 5.88384784e-16],$
 	  		[3.41071678e-10, -1.82060569e-12, 1.59690189e-15, -3.14810895e-19]]

Kx_sx = [[-7.96016109e+00, 9.26584096e-03, -7.93676069e-06, 5.13414639e-10],$
         [1.02925896e+00, -1.59974177e-05, 9.57769272e-09, -1.14409822e-12],$
         [-2.30169348e-05, -3.45351550e-09, 1.89621010e-12, 6.72971750e-17],$
         [7.07041647e-09, 2.20511200e-12, -1.66433082e-15, 2.55463711e-19]]

Ky_sx = [[-2.26409123e+00, 9.93175401e-01, -6.67169688e-06, 8.08275391e-09],$
         [-1.38521740e-02, -2.27910031e-05, 4.72367613e-09, -1.64738716e-12],$
         [8.17060299e-06, 1.35240460e-08, -5.52374318e-12, 2.14966954e-15],$
         [9.25982725e-10, -2.50607186e-12, 2.58675626e-15, -9.82463036e-19]]
         
half_cropped_sz = 250.;250. ; px (used to be 250, i.e., 500 px by 500 px cropped images)
annmode_ringw = 15;FWHM ; on either side of the object (so, 2x FWHM total)
wr=annmode_ringw

;------------------------------[ TWO SOLUTION OVERRIDE ]--------------------------
 
;if keyword_set(two_soln_override) && two_soln_override eq 1 then begin
;	pxscale_dx = pxscale_sx
 ;  Kx_dx = Kx_sx
  ; Ky_dx = Ky_sx
 ;  print, 'Two astrometric solution override active, using SX soln for both sides'
 ;  wait, 2; wait for two seconds to see the message
;endif; two_soln_overrride if

;------------------------------[ End User Input ]---------------------------------

if keyword_set(planet_x) then begin
   planet_r = sqrt((half_cropped_sz-planet_x)^2 + (half_cropped_sz-planet_y)^2)
endif

if use_injection and not neg_inj then begin; Custom annulus inner and outer radii for artificial injections at arbitrary locations
   annmode=1
   ;Thicc rings
   annmode_inout_sx = round([max([0.,planet_r/pxscale_sx-annmode_ringw+1.]),planet_r/pxscale_sx+annmode_ringw+2.])
   if annmode_inout_sx[1] gt 16 then BEGIN
   	annmode_inout_sx = round([max([0.,planet_r/pxscale_sx-annmode_ringw*1.2+1.]),planet_r/pxscale_sx+annmode_ringw*1.2+2.])
   endif
   
    annmode_inout_dx = round([max([0.,planet_r/pxscale_dx-annmode_ringw+1.]),planet_r/pxscale_dx+annmode_ringw+2.])
   if annmode_inout_dx[1] gt 16 then BEGIN 
   	annmode_inout_dx = round([max([0.,planet_r/pxscale_dx-annmode_ringw*1.2+1.]),planet_r/pxscale_dx+annmode_ringw*1.2+2.])
   endif
   ; Thinn rings
   ;annmode_inout=round([max([0.,planet_r/pxscale-10.+1.]),planet_r/pxscale+10.+2.])
   ;if annmode_inout[1] gt 16 then annmode_inout=round([max([0.,planet_r/pxscale-12.+1.]),planet_r/pxscale+12.+2.])
endif

if full_frame eq 1 then wr=fix((half_cropped_sz)/float(nrings)); width of ring for KLIP

if klip eq 1 then begin
	suffix=strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(angsep,4))+$
		'as_'+string(sigfig(anglemax,2))+'am_'+string(sigfig(nrings,2))+'rings_'+$
		string(wr)+'wr_'+String(n_ang)+'nang_'+string(sigfig(filter,3))+'filter_'+$
		string(sigfig(bin,2))+'bin'+string(corr_thresh)+'corrthresh'+'_inrad_'+string(annmode_inout_dx[0])), /remove_all)
endif
	
if use_injection then suffix=strcompress(suffix+'_inj',/rem)
if neg_inj then suffix=strcompress(suffix+'neg_inj',/rem)

if not keyword_set(uncert) then begin; not uncert (default)
	uncert=0
endif else begin; uncert
	suffix=strcompress(suffix+'uncert',/rem)
endelse

;------------------------------[ Pipeline ]---------------------------------

if pre_inj eq 1 then begin
   
   print, 'pre_inj:', pre_inj
   
   if pre_clean eq 1 then begin

		;create_cube, obj, raw_path, cube_start_frame, coadd, output_path, cds
		;
		;sky_sub, obj, coadd, output_path, fs_start=fs_start, fs_end=fs_end, fpn=fpn,$
		;	sky_type=sky_type, pre_bad_px=1, nnod=nnod
		
		;bad_pixels, output_path, obj, stripe, sigma_threshold=sigma_threshold,$
		;	boxhsize=boxhsize, pre_sky_sub=0, bad_px_arr=bad_px_arr
	
		;dewarp, output_path, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx, run=run,$
		 ;	skysub_first=1
	
		;split, obj, stripe, output_path, half_cropped_sz, aperture
		
		;center, obj, stripe, output_path, do_block_right, do_block_left,$
		;	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub,$
		;	half_cropped_sz, aperture, corr_thresh, good_framei=good_framei,$
		;	align_radius_thresh=align_radius_thresh
			
		center_old, obj, stripe, output_path, do_block_right, do_block_left,$
	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub, half_cropped_sz
			
	endif; pre_clean if
   ;
   clean, obj, stripe, output_path, corr_thresh, half_cropped_sz, aperture,$
   		annular_clean=annular_clean, centered_clean=centered_clean, do_second_center=do_second_center,$
   		peak_clean=peak_clean, peak_thresh=peak_thresh,stddev_clean=stddev_clean,$
   		stddev_thresh=stddev_thresh, fwhm=fwhm

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, corr_thresh, do_cen_filter,$
   	planet_r=rho, planet_theta=theta, use_gauss=use_gauss,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
   	
   if keyword_set(planet_x) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, corr_thresh, do_cen_filter,$
   	planet_y=planet_y, planet_x=planet_x, use_gauss=use_gauss, silent=silent,$
   	truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
   
   ; Change output folder manually in adi.pro and klip.pro !!!!!!!!!!!!!!! 
   ; (right now as long as it's macbook_<coadd> it's fine, and ssh is set
   ; appropriately)
   
   ; I'm having trouble with find_sources here. (Everything is working now, but
   ; note that I might need to adjust the correction factor to get acurate values)
   
   if adi eq 1 then begin
   		adi, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, corr_thresh,$
   			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj,normal=normal, uncert=uncert,$
   			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band, fwhm=fwhm,$
   			combine_type=combine_type
   	endif
   	
	if klip eq 1 then begin
		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, bin, bin_type,$
			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
			annmode_inout_sx, annmode_inout_dx, suffix, corr_thresh, do_cen_filter, coadd,$
			fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl
   endif; klip eq 1 if
   
endif
if pre_inj eq 0 then begin
   
   print, 'pre_inj:', pre_inj

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, corr_thresh, do_cen_filter,$
   	planet_r=rho, planet_theta=theta, use_gauss=use_gauss,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod,$
   	combine_type=combine_type
    
   if keyword_set(planet_x) then begin
   		print, planet_x, planet_y
      inject_planets, obj, output_path, n_planets, contrast, pxscale_sx,$
      pxscale_dx, corr_thresh, do_cen_filter, planet_y=planet_y, planet_x=planet_x,$
      use_gauss=use_gauss, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod,$
      silent=silent
   endif
	
	if adi eq 1 then begin
		print,'starting adi'
   		adi, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, corr_thresh,$
   			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
   			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band, fwhm=fwhm,$
   			combine_type=combine_type
   	endif
   	
	if klip eq 1 then begin
		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, bin, bin_type,$
			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
			annmode_inout_sx, annmode_inout_dx, suffix, corr_thresh, do_cen_filter, coadd,$
			trial=trial, fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx,$
			truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx,$
			magnify=magnify, fwhm=fwhm, wl=wl
	endif; klip eq 1 if
   
endif

;-----------------------------------[ El Fin ]--------------------------------------

; 1440 minutes per day (Julian dates are measured in days)
print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!

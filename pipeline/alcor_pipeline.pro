pro alcor_pipeline, rho=rho, theta=theta, planet_x=planet_x, planet_y=planet_y,$
	contrast=contrast, pre_inj=pre_inj, neg_inj=neg_inj, trial=trial,$
   outpath=outpath, coadd=coadd, use_gauss=use_gauss, uncert=uncert, klip=klip,$
   fs=fs, two_soln_override=two_soln_override, extra=extra, nod=nod, cds=cds,$
   corr_thresh=corr_thresh, pre_clean=pre_clean

; alcor_pipeline, pre_inj=1, neg_inj=0, coadd=25, uncert=0, klip=1 for normal use

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

;------------------------------[ Start User Input ]---------------------------------

half_cropped_sz = 190.;250. ; px (used to be 250, i.e., 500 px by 500 px cropped images)

; General/Combine Parameters
obj = 'Alcor'
band = 'M'
raw_path = '/Volumes/T7/alcor_raw'
stripe = 'second512'; not center 1024 pixel rows, but second set of 512 from
run = 1; sx/left aperture only (for dewarp.pro)
; the top of the full 2048x2048 frame
;'/Users/gabeweible/OneDrive/research/HII1348/kevin/raw'
cube_start_frame = 0
; fwhm calculated with fullwid_halfmax(pupil_array read-in with readfits, /GAUSSIAN_FIT)
fwhm = 13.0773; Average of dith1 and dith2 fwhm in x and in y (from pupil medians)
; CHANGE LATER????
magnify = 0; Magnify
output_path = '/Users/gabeweible/OneDrive/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)
output_path += '/'

if not keyword_set(pre_clean) then pre_clean=0

; Planet injection parameters; KEEP FOR NOW AND CHANGE LATER?????
; use_gauss can be set with a keyword
pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared)
pxscale_dx = 0.010700 ; +0.000042 or -0.000051 arcsec/pixel (from Steve and Jared)
truenorth_sx = -1.278 ; +0.131 or -0.225 deg (from Steve and Jared)
truenorth_dx = 1.001 ; + 0.254 or -0.237 deg (from Steve and Jared)

;'/Users/gabeweible/OneDrive/research/HII1348/macbook_' +$
;	+ strcompress(coadd, /r) + '/'
if keyword_set(outpath) then output_path=outpath

; Sky sub parameters
;allowed_coadds = [1, 2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 200]
;if total(allowed_coadds eq coadd) ne 1 then begin
;  print, 'Error: invalid coadd value (must be a factor of 200)'
;  stop
;endif

fpn = 1000; frames per nod
; inclusive start and end frames (IDL indices start at 0) for the first set
; of images at the second nodding position.
fs_start = fpn / coadd & fs_end = (2*fpn - 1) / coadd
;print, 'Coadd:', coadd

; Centering parameters; CHANGE THESE FOR SATURATED PSF????
; might be fine still? This may only matter for injections, really
do_block_right = 1
do_block_left = 1
do_block_bottom = 1
do_block_top = 1
do_cen_filter = 1
do_cen_median_sub = 0; subtract the median before centering images

; Clean parameters
; 0.993 is for KLIP, 0.995 is for ADI (HII 1348)
; 	CHANGE FOR ALCOR????
; With 'endcds', it looks like we might want a higher threshhold
; set in kwarg now;corr_thresh = 0.997 ; 

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
if not keyword_set(cds) then cds='endcds'; max exposure time.

; ADI/KLIP parameters (after testing/refinement)
if not keyword_set(use_gauss) then use_gauss = 0
silent = 0; Don't print so much in adi.pro "Rotating by ..."
normal = 0; Normalize ADI frames to reduce residuals around the star
; normal = 1 did not result in good star subtraction. I'm not sure yet if
; use_gauss should be 1 or 0 yet.

;fs = 1; Run find_sources within ADI.pro; Set in a kwarg now...
if not keyword_set(fs) then fs = 0; Assume it's NOT fine?

if keyword_set(rho) or keyword_set(planet_x) then begin
	use_injection=1
	do_annmode=1
endif else begin
	use_injection=0
	do_annmode=0
endelse

do_destripe = 1
; subtracts a smoothed version of the image
filter = 0; high-pass filter kernel size for ADI and KLIP; or 0 for no filter

; Can be re-optimized
bin = 3
bin_type = 'mean'; 'mean'
do_hyper = 0
combine_type = 'nwadi'; 'nwadi'

; Start and end frame if only using KLIP on a fraction of dataset (start and end)
klip_fraction = 0
klip_start_frame = 75
klip_end_frame = 150

fill = 0
; These KLIP parameters could be re-optimized
k_klip = 7; 5 ; 7; Number of KLIP components (eigenimages) to keep; 5? (in Kevin's code)
angsep= 1.;0.25; 1.; 0.25? (in Kevin's code)
anglemax = 360.
nrings = 4.;10;4. ; Number of annuli for KLIP; 10? (in Kevin's code)
n_ang = 2;1;2 ; Number of segments in each annulus for KLIP; 1? (in Kevin's code)

; CHANGE THESE VALUES FOR ALCOR B!!!!!!!!!
; Annulus inner and outer radii already setup for HII 1348 b; In px?
annmode_inout_sx = [83, 143]
; Annulus inner and outer radii already setup for HII 1348 b; In px?
annmode_inout_dx = [83, 143]


; Astrometry Parameters (Dewarp solns from 
; https://scienceops.lbto.org/lbti/data-retrieval-reduction-publication/distortion-correction-and-astrometric-solution/)

Kx_sx = [[-7.96016109e+00, 9.26584096e-03, -7.93676069e-06, 5.13414639e-10],$
         [1.02925896e+00, -1.59974177e-05, 9.57769272e-09, -1.14409822e-12],$
         [-2.30169348e-05, -3.45351550e-09, 1.89621010e-12, 6.72971750e-17],$
         [7.07041647e-09, 2.20511200e-12, -1.66433082e-15, 2.55463711e-19]]

Ky_sx = [[-2.26409123e+00, 9.93175401e-01, -6.67169688e-06, 8.08275391e-09],$
         [-1.38521740e-02, -2.27910031e-05, 4.72367613e-09, -1.64738716e-12],$
         [8.17060299e-06, 1.35240460e-08, -5.52374318e-12, 2.14966954e-15],$
         [9.25982725e-10, -2.50607186e-12, 2.58675626e-15, -9.82463036e-19]]
         
Kx_dx = [[-1.13034544e+01, 1.45852226e-02, -1.13372175e-05, 1.32602063e-09],$
 	  		[1.03220943e+00, -1.93058980e-05, 1.55798844e-08, -3.86115281e-12],$
 	  		[-2.57352199e-05, 2.70371257e-09, -6.62650011e-12, 3.25017078e-15],$
 	  		[8.02004325e-09, -5.59037685e-13, 1.60256679e-15, -8.18749145e-19]]

Ky_dx = [[-9.37023860e-01, 9.89973161e-01, -4.32284634e-06, 7.08000564e-09],$
 	  		[-1.59438718e-02, -1.95099497e-05, 1.55801035e-09, 2.13743170e-13],$
 	  		[9.34251811e-06, 1.26473736e-08, -3.71968684e-12, 5.88384784e-16],$
 	  		[3.41071678e-10, -1.82060569e-12, 1.59690189e-15, -3.14810895e-19]]

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
   annmode_inout_sx = round([max([0.,planet_r/pxscale_sx-25.+1.]),planet_r/pxscale_sx+25.+2.])
   if annmode_inout_sx[1] gt 16 then BEGIN
   	annmode_inout_sx = round([max([0.,planet_r/pxscale_sx-30.+1.]),planet_r/pxscale_sx+30.+2.])
   endif
   
    annmode_inout_dx = round([max([0.,planet_r/pxscale_dx-25.+1.]),planet_r/pxscale_dx+25.+2.])
   if annmode_inout_dx[1] gt 16 then BEGIN
   	annmode_inout_dx = round([max([0.,planet_r/pxscale_dx-30.+1.]),planet_r/pxscale_dx+30.+2.])
   endif
   ; Thinn rings
   ;annmode_inout=round([max([0.,planet_r/pxscale-10.+1.]),planet_r/pxscale+10.+2.])
   ;if annmode_inout[1] gt 16 then annmode_inout=round([max([0.,planet_r/pxscale-12.+1.]),planet_r/pxscale+12.+2.])
endif

wr=fix(half_cropped_sz/float(nrings)); width of ring for KLIP

suffix=strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(ANGSEP,4))+$
	'as_'+string(sigfig(anglemax,2))+'am_'+string(sigfig(nrings,2))+'rings_'+$
	string(wr)+'wr_'+String(n_ang)+'nang_'+string(sigfig(filter,2))+'filter_'+$
	string(sigfig(bin,2))+'bin'+string(corr_thresh)+'corrthresh'), /remove_all)
	
if use_injection then suffix=strcompress(suffix+'_inj',/rem)
if neg_inj then suffix=strcompress(suffix+'neg_inj',/rem)
if uncert then suffix=strcompress(suffix+'uncert',/rem)
;if neg_inj eq 1 then begin & annmode=0 & do_annmode=0 & endif

;------------------------------[ Pipeline ]---------------------------------

if pre_inj eq 1 then begin
   
   print, 'pre_inj:', pre_inj
   
   if pre_clean eq 1 then begin

		create_cube, obj, raw_path, cube_start_frame, coadd, output_path, cds
		;
		bad_pixels, output_path, obj, stripe
		;
		sky_sub, obj, coadd, output_path, fs_start=fs_start, fs_end=fs_end, fpn=fpn
	
		dewarp, output_path, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx, run=run
	
		split, obj, stripe, output_path, half_cropped_sz
		;
		center, obj, stripe, output_path, do_block_right, do_block_left,$
			do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub,$
			half_cropped_sz
			
	endif; pre_clean if
   ;
   clean, obj, stripe, output_path, corr_thresh, half_cropped_sz

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, corr_thresh, do_cen_filter,$
   	planet_r=planet_r, planet_theta=planet_theta, use_gauss=use_gauss,$
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
   
   adi, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, corr_thresh,$
   	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj,normal=normal, uncert=uncert,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band
   	
	if klip eq 1 then begin
		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, bin, bin_type,$
			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
			annmode_inout_sx, annmode_inout_dx, suffix, corr_thresh, do_cen_filter, coadd,$
			fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm
   endif; klip eq 1 if
   
endif
if pre_inj eq 0 then begin
   
   print, 'pre_inj:', pre_inj

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, corr_thresh, do_cen_filter,$
   	planet_r=planet_r, planet_theta=planet_theta, use_gauss=use_gauss,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
    
   if keyword_set(planet_x) then begin
      inject_planets, obj, output_path, n_planets, contrast, pxscale_sx,$
      pxscale_dx, corr_thresh, do_cen_filter, planet_y=planet_y, planet_x=planet_x,$
      use_gauss=use_gauss, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
   endif
	
   adi, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, corr_thresh,$
   	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band
   	
	if klip eq 1 then begin
		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, bin, bin_type,$
			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
			annmode_inout_sx, annmode_inout_dx, suffix, corr_thresh, do_cen_filter, coadd,$
			trial=trial, fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx,$
			truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx,$
			magnify=magnify, fwhm=fwhm
	endif; klip eq 1 if
   
endif

;-----------------------------------[ El Fin ]--------------------------------------

; 1440 minutes per day (Julian dates are measured in days)
print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!

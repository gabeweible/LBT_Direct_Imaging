pro pleiad3_pipeline, rho=rho, theta=theta, planet_x=planet_x, planet_y=planet_y,$
	contrast=contrast, pre_inj=pre_inj, neg_inj=neg_inj, trial=trial,$
   outpath=outpath, coadd=coadd, use_gauss=use_gauss, uncert=uncert, klip=klip,$
   fs=fs, two_soln_override=two_soln_override, extra=extra, nod=nod, cds=cds,$
   pre_clean=pre_clean, k_klip=k_klip, adi=adi, do_annmode=do_annmode,$
   szbin=szbin, angsep=angsep, bin_type=bin_type, filter=filter,$
	n_ang=n_ang,combine_type=combine_type,do_rdi=do_rdi,refcube=refcube,$
	annmode_inout=annmode_inout,basis=basis,ref_ang=ref_ang,do_new_binning=do_new_binning,$
	do_create_basis=do_create_basis, do_dewarp=do_dewarp, save_basis=save_basis, anglemax=anglemax,$; n_groups=n_groups,$
	normal=normal, rdi_extra=rdi_extra, mask_pt=mask_pt,$
	dxklip=dxklip, peak_clean=peak_clean, peak_thresh=peak_thresh, stddev_clean=stddev_clean,$
	stddev_thresh=stddev_thresh, set_coadd_type=set_coadd_type, $
	keep_number=keep_number, nod_filter=nod_filter, cube_indices=cube_indices
	
; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Juliasn date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)	
	
cenframewidth=40.
inner_rad_cen=7 & outer_rad_cen=20

; smooth before dewarping?
do_smooth = 1.; FWHM for Gaussian smoothing before and after dewarping.
; data.
; how to bin again before dewarping?
dewarp_bin_type='mean'
; doing sky subtraction after bad-px correction, now.
skysub_first = 0

; could also be set as 'dith1' or 'dith2' - for dewarping and splitting only.
if not keyword_set(nod_filter) then nod_filter = 'both'

; create master dark and master flat?
create_master_masks = 0
; create bad-pixel mast from master dark and master flat?
create_badpix_mask = 0

; when do darks start in the co-added cube?
dark_frame_start = 'None';

mask=0; radius for final ADI image masking before convolution
outmask=150

; new cube creation scheme
; not actually using Fixen et al. algo. b/c CDS frames only here.
fixsen=1

; for bad-pixel masks
do_destripe = 1; done in bad_pixels_fast.pro now...

; AUTOMATIC VALUES
;-----------------
if keyword_set(normal) then begin
	print, 'normalize ADI: ', normal
endif else begin
	normal = 0
endelse

; NEW 01/05/2026
if not keyword_set(filter) then filter = 0.0
if not keyword_set(keep_number) then keep_number=2400
if not keyword_set(do_rdi) then do_rdi = 0

; Normalize ADI frames to reduce residuals around the star
if not keyword_set(neg_inj) then neg_inj = 0
if not keyword_set(coadd) then coadd = 1
if not keyword_set(klip) then klip = 0
if not keyword_set (fs) then fs = 0
if not keyword_set(adi) then adi = 0
if not keyword_set(do_annmode) then do_annmode = 1
;if not keyword_set(n_groups) then n_groups = 1
if not keyword_set(szbin) then szbin = 1
if not keyword_set(do_create_basis) then do_create_basis=1
if not keyword_set(rdi_extra) then rdi_extra = ''
if not keyword_set(annmode_inout) then annmode_inout_dx=[0,350];annmode_inout_dx = [0, wr]
if keyword_set(annmode_inout) then annmode_inout_dx = annmode_inout
if not keyword_set(do_dewarp) then do_dewarp=1 ; assume that we should dewarp
if not keyword_set(do_new_binning) then do_new_binning = 0
if not keyword_set(nod) then nod = 'both'
if not keyword_set(cds) then cds = 'endcds'
if not keyword_set(use_gauss) then use_gauss=0; for cADI PSF model - don't use by default.
if not keyword_set(pre_clean) then pre_clean=0
if not keyword_set(save_basis) then save_basis=0
if not keyword_set(k_klip) then k_klip=4; how many KLIP components to retain?
if not keyword_set(angsep) then angsep=0.0; how many multiples of 1.5xfwhm for PSF frames?
if not keyword_set(n_ang) then n_ang=1; how many segments per annulus?
if not keyword_set(bin_type) then bin_type = 'mean'; 'mean' or 'median'
if not keyword_set(combine_type) then combine_type = 'median'; 'mean', 'median', or 'nwadi'
if keyword_set(neg_inj) then begin
	if (neg_inj eq 1) then n_planets=1
endif else neg_inj = 0
if keyword_set(use_gauss) then use_gauss=[use_gauss] else use_gauss=0
if keyword_set(rho) then planet_theta=[theta]
if keyword_set(rho) then planet_r=[rho]
if keyword_set(planet_x) then planet_x=[planet_x]
if keyword_set(planet_y) then planet_y=[planet_y]
if keyword_set(contrast) then begin
	planet_contrast=[contrast]
	if keyword_set(rho) then n_planets=n_elements(planet_r)
	if keyword_set(planet_x) then n_planets=n_elements(planet_x)
endif
if keyword_set(outpath) then output_path=outpath
if keyword_set(rho) or keyword_set(planet_x) then begin
	use_injection=1
	do_annmode=1
endif else begin
	use_injection=0
	if not keyword_set(do_annmode) then do_annmode = 0; default to no annuli
endelse

if ~keyword_set(dxklip) then dxklip=0
; normal = 1 did not result in good star subtraction. I'm not sure yet if
; use_gauss should be 1 or 0 yet.
; BUG FIXED TO NORMAL KWARG: now normalize to peak pixel value in FWHM/2 radius, not
; in the entire frame (which may have a bad pixel, for example.)
;------------------

; Annulus inner and outer radii in px
;annmode_inout_sx = [0, wr]
; Annulus inner and outer radii in px
wr = annmode_inout_dx[1]-annmode_inout_dx[0] ; px
; for cube creation
if not keyword_set(set_coadd_type) then begin
	coadd_type='median'; seems like it may be a better option than mean, at least for
	; even coadds (where really the mean of the two central values is taken, which
	; is effectively a pretty good outlier-resistant mean for, e.g., coadd=4 or coadd=6)
endif else begin
	coadd_type = set_coadd_type
endelse

;------------------------------[ Start User Input ]---------------------------------
; frame half-size in pixels. Trying as large as I can before hitting frame edges.
half_cropped_sz = 481; 

; General/Combine Parameters
obj = 'HIP17034'; OBJNAME in FITS header
band = 'L'
annular_clean=0; bloc < 10 px and > 50 px (focus on cleaning based on 1st airy ring)
centered_clean=1
wl = 3.672131; from SVO (at 77 K, effective, not central)
raw_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD3/raw'
darks_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD3/darks'

stripe = 'Full_Image' ; NEW: I'll be working with full 2048 x 2048 frames up until cropping.
aperture = 'both'

; choose a "run" for dewarp.pro if single-sided
if aperture eq 'right' then run = 2; DX/right aperture (only for dewarp.pro)
if aperture eq 'left' then run = 1; DX/right aperture (only for dewarp.pro)
if aperture eq 'both' then run = 5; for what used to be 'both'

cube_start_frame = 0
; fwhm calculated with fullwid_halfmax(pupil_array read-in with readfits, /GAUSSIAN_FIT)

;______________________________________
; OLD FROM TYC5709 - UPDATE FOR PLEIAD3
;______________________________________
fwhm = 10.6245 ; from VIP for TYC 5709 - to be updated for PLEIAD 3
lambda_over_d = 8.487; temporary - estimated for 10.61 mas/px nominal platescale

; bad pixel parameters
; sigma_threshold no longer relevant here - using master mask
; from the dark and flat + some sigma-clipping.
; I should define the sigma-clipping sigma level
; here, though (currently hard-coded to 4.0)
;sigma_threshold=50.0 ; sigma_threshold*stdev outside of surrounding values are out!
boxhsize=1; 2*boxhsize x 2*boxhsize block of pixels for each comparison
sigma_clip=3.5; for 3x3, 5x5, and 7x7 filters in sky_sub.pro

; magnification
magnify = 0; Magnify or not (not will keep native DX pxscale)

; update directory
output_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD3/gabe_macbook6_median_coadd/'

; use_gauss can be set with a keyword
; These are required rotations E of N, i.e., CCW.

;______________________________________
; OLD FROM TYC5709 - UPDATE FOR PLEIAD3
;______________________________________
pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared — OLD)
pxscale_dx = 0.010624 ; +0.000035 or -0.000038 arcsec/pixel (Jared 06/20/25)
truenorth_sx = -1.278 ; +0.131 or -0.225 deg (from Steve and Jared — OLD)
truenorth_dx = 0.874 ; + 0.176 or -0.292 deg (Jared 06/20/25)

; Sky sub parameters
fpn = 500; 'frames' per nod (really, ramps per nod position)

; There does not seem to be a lot of large-scale evolution between each nod
nnod = fix(fpn/coadd); nearest 1 nod position of frames (50) should be okay? (100 coadded by 20)

; inclusive start and end frames (IDL indices start at 0) for the first set
; of images at the second nodding position.
fs_start = fpn / coadd & fs_end = (2*fpn - 1) / coadd
;print, 'Coadd:', coadd

; Centering parameters;
; This may only matter for injections, really
do_block_right = 0
do_block_left = 0
do_block_bottom = 0
do_block_top = 0
do_cen_filter = 1
do_cen_median_sub = 0; subtract the median before centering images
silent = 0; Don't print so much in adi.pro "Rotating by ..."

align_radius_thresh = 20
align_inner_radius=7

if not keyword_set(peak_clean) then peak_clean=0	;reject frames where the peak flux is beneath a certain
;fraction compared to the average peak value
if not keyword_set(peak_thresh) then peak_thresh=0	;normalized units, so 1 is the average peak strength
	
if not keyword_set(stddev_clean) then stddev_clean=0
if not keyword_set(stddev_thresh) then stddev_thresh=0; max. ratio of pixel stddev to median pixel stddev

do_hyper = 0

fill = 0.; I think this just allows the rest of the image to follow along with the KLIP subtraction?
; These KLIP parameters could be re-optimized
nrings = 1;10;4. ; Number of annuli for KLIP; 10? (in Kevin's code)
full_frame=0; only klip in annuli.

;______________________________________
; OLD FROM TYC5709 - UPDATE FOR PLEIAD3
;______________________________________
; Astrometry Parameters (Dewarp solns from 
; https://scienceops.lbto.org/lbti/data-retrieval-reduction-publication/distortion-correction-and-astrometric-solution/)
; NEW from Jared ( added 06/20/2025 )
Kx_dx = [[-1.63626928e+01, 1.51722203e-02, -1.07497658e-05, 1.23226160e-09],$
		[1.02743011e+00, -1.11105659e-05, 1.17477374e-08, -3.68547617e-12],$
		[-2.35063091e-05, -4.97851470e-09, -3.48722610e-12, 3.55627029e-15],$
		[ 6.98282439e-09, 2.11881053e-12, 6.03681585e-16, -9.39558827e-19]]

Ky_dx = [[1.99522924e+01, 9.88214317e-01, -4.95231155e-06, 7.63350467e-09],$
		[-1.72489250e-02, -2.63781860e-05, 1.02979303e-08, -3.05000352e-12],$
		[5.62099668e-06, 2.60899303e-08, -2.00090215e-11, 6.47595533e-15],$
		[2.15105717e-09, -8.36555891e-12, 9.24058823e-15, -3.00298317e-18]]
 
 ; NEW from Jared ( added 06/20/2025 )
 Kx_sx = [[-1.13361367e+01, 1.76598270e-02, -1.40547004e-05, 1.95792550e-09],$
 			[1.03383062e+00, -3.17945681e-05, 2.58087494e-08, -4.92697861e-12],$
 			[-2.82524642e-05, 1.09594539e-08, -1.27950019e-11, 3.68000336e-15],$
 			[8.43745199e-09, -1.93507394e-12, 2.73218037e-15, -8.75155836e-19]]

Ky_sx = [[1.69544825e+01, 9.93121265e-01, -6.41091149e-06, 8.39074091e-09],$
		[-1.57576842e-02, -2.91980779e-05, 1.08894491e-08, -3.65243996e-12],$
		[6.75019653e-06, 2.32402580e-08, -1.60611833e-11, 5.55601798e-15],$
		[1.58471068e-09, -6.39813009e-12, 7.01919061e-15, -2.44373489e-18]]

annmode_ringw = wr ; on either side of the object (so, 2x FWHM total)
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

if klip eq 1 or do_rdi eq 1 then begin
	suffix=strcompress(reform('_'+string(k_klip)+'_k_'+string(sigfig(angsep,4))+$
		'_as_'+string(sigfig(anglemax,2))+'_am_'+string(sigfig(nrings,2))+'_rings_'+$
		string(wr)+'_wr_'+String(n_ang)+'_nang_'+string(sigfig(filter,3))+'_filter_'+$
		string(sigfig(szbin,2))+'_bin_'+string(keep_number)+'_keepnumber'), /remove_all)
endif
	
if use_injection then suffix=strcompress(suffix+'_inj',/rem)
if neg_inj then suffix=strcompress(suffix+'neg_inj',/rem)

if not keyword_set(uncert) then begin; not uncert (default)
	uncert=0
endif else begin; uncert
	suffix=strcompress(suffix+'uncert',/rem)
endelse

if not keyword_set(anglemax) then anglemax=34.

;------------------------------[ Pipeline ]---------------------------------

if pre_inj eq 1 then begin
   
   print, 'pre_inj:', pre_inj
   
   if pre_clean eq 1 then begin
        ; not actually using fixen algorithm up-the-ramp, just CDS in this case.
    	cosmic_ray_rejection_cube, obj, raw_path, cube_start_frame, coadd,$
    		output_path, read_header=1, coadd_type=coadd_type, legacy_mode='cds',$
    		debug=1, dark_frame_start=dark_frame_start, max_frames_per_group=fix(fpn-1) / coadd,$
    		min_r_squared=0.9, full_well=4095*0.98,skip_first_read=0,$
    		skip_last_read=0, fpn=fpn, verbose=1, del_nodframe1=1
		
		; create darks cube
		cosmic_ray_rejection_cube, obj, darks_path, cube_start_frame, coadd,$
				output_path, read_header=1, coadd_type=coadd_type, legacy_mode='cds',$
				debug=1, dark_frame_start=0,$
				max_frames_per_group=fix(2*fpn - 1) / coadd, dark_cube=1,$
				min_r_squared=0.9, full_well=4095*0.98,skip_first_read=0,$
				skip_last_read=0, fpn=fpn, verbose=1, del_nodframe1=1
		
		; not using bad_px_arr right now (find them all via computation)
		;bad_pixels_fast, output_path, obj, stripe,$
		;	boxhsize=boxhsize, pre_sky_sub=1, type='res_mean',$
		;	run=run, do_destripe=do_destripe,$
		;	fwhm=fwhm, just_destripe=0, do_second_round=0,$
		;	do_first_round=1, just_second_round=0, debug=1,$
		;	create_master_masks=create_master_masks,$
		;	create_badpix_mask=create_badpix_mask, coadd=coadd,$
		;	darks_filename='HIP17034_darks_cube.fits', vapp=0
		
		; NEW!! Sigma-clipping after bad-pixel-mask interpolation
		; this is frame-to-frame
		;sigma_clip_cubes, obj, output_path, cube_folder=output_path, $
         ;             sigma_clip=sigma_clip, debug=1, $
          ;            overwrite=1, frame_min=-10,$
           ;           cube_indices=[15,16,17,18,19,20,21]
      
      ; NEW!!! fill remaining NaNs with frame medians + MAD-estimated Gaussian
      ; noise
   ; fill_nan_cubes, obj, output_path, cube_folder=output_path, $
    ;                  debug=debug, overwrite=1, cube_indices=[7]
	;	
		; SKY SUBTRACTION IN 'pleiad3_pca_bkg_subtr_1cell.ipynb' !!!!
		
		; this isn't actually bad-px correction, I'm just destriping
		; the sky-subtracted frames, now.
	;	bad_pixels_fast, output_path, obj, stripe,$
	;		boxhsize=boxhsize, pre_sky_sub=0, type='res_mean',$
	;		run=run, do_destripe=do_destripe,$
	;		fwhm=fwhm, just_destripe=1, do_second_round=0,$
	;		do_first_round=1, just_second_round=0, debug=1,$
	;		create_master_masks=create_master_masks,$
	;		create_badpix_mask=create_badpix_mask, coadd=coadd,$
	;		darks_filename='HIP17034_darks_cube.fits', vapp=0,$
	;		destripe_skysub=1, pca_skysub=1
		
		;fill_nan_cubes, obj, output_path, cube_folder=output_path,$
		;		debug=debug, overwrite=0, input_suffix='_10comp_skysub_destriped_cube.fits',$
		;		output_suffix='_pca_skysub_filled_cube.fits',$
		;		input_prefix='test_pca_skysub_cube'
			
		;if do_dewarp eq 1 then begin
		;
		;	; new combined dewarping and splitting to reduce memory requirements.
		;	dewarp_split_combined, output_path, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx,$
		;		do_smooth=do_smooth, half_cropped_sz=half_cropped_sz, aperture=aperture, run=run,$
		;		fwhm=fwhm, hp_width=0, nod_filter=nod_filter, debug=1,$
		;		destripe_skysub=1, pca_skysub=1, filled=1
				
		;endif; else begin; dewarp if
		;	; 1 is for skysub_first
		;	split, obj, stripe, output_path, half_cropped_sz, aperture, skysub_first, do_dewarp,$
		;		fwhm=fwhm
		;endelse
		
    ; Centering moved to VIP pipeline.
			
	endif; pre_clean if
   ;
  ; clean, obj, stripe, output_path, keep_number, half_cropped_sz, aperture,$
   ;		annular_clean=annular_clean, centered_clean=centered_clean,$
   	;	peak_clean=peak_clean, peak_thresh=peak_thresh,stddev_clean=stddev_clean,$
   	;	stddev_thresh=stddev_thresh, fwhm=fwhm, do_dewarp=do_dewarp,$
   	;	do_block_right, do_block_left,$
	;	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub,$
	;	do_smooth=do_smooth

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, keep_number, do_cen_filter,$
   	planet_r=rho, planet_theta=theta, use_gauss=use_gauss,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
   	
   if keyword_set(planet_x) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, keep_number, do_cen_filter,$
   	planet_y=planet_y, planet_x=planet_x, use_gauss=use_gauss, silent=silent,$
   	truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
   
   ; Change output folder manually in adi.pro and klip.pro !!!!!!!!!!!!!!! 
   ; (right now as long as it's macbook_<coadd> it's fine, and ssh is set
   ; appropriately)
   
   ; I'm having trouble with find_sources here. (Everything is working now, but
   ; note that I might need to adjust the correction factor to get acurate values)
   
  ; 	if adi eq 1 then begin
   ;		adi, obj, half_cropped_sz, nod, output_path, use_injection, filter, keep_number,$
   	;		do_cen_filter, coadd, fs=fs, neg_inj=neg_inj,normal=normal, uncert=uncert,$
   	;		silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   	;		pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band, fwhm=fwhm,$
   	;		combine_type=combine_type, szbin=szbin, bin_type=bin_type,$
   	;		peak_thresh=peak_thresh, stddev_thresh=stddev_thresh,$
   	;		mask=mask, outmask=outmask, do_smooth=do_smooth
   	;endif
   	
	;if klip eq 1 then begin
	
	;		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, szbin, bin_type,$
	;			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
	;			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
	;			annmode_inout_sx, annmode_inout_dx, suffix, keep_number, do_cen_filter, coadd,$
	;			fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	;			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl,do_rdi=do_rdi,refcube=refcube,$
	;			basis=basis,ref_angles=ref_angles, do_new_binning=do_new_binning
		
   	;endif; klip eq 1 if
   
	;if do_rdi eq 1 then begin
	;	rdi, obj, half_cropped_sz, nod, output_path, use_injection, filter, szbin, bin_type,$
	;		do_annmode, combine_type, k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout_sx, annmode_inout_dx,$
	;		suffix, keep_number, do_cen_filter, coadd, trial=trial, fs=fs, neg_inj=neg_inj,$
	;		truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx,$
	;		pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl, refcube=refcube, ref_ang=ref_ang,$
	;		save_basis=save_basis, rdi_extra=rdi_extra,$
	;		mask_pt=mask_pt, do_new_binning=do_new_binning,$
	;		do_create_basis=do_create_basis, dxklip=dxklip,$
	;		normal=normal, do_smooth=do_smooth, keep_number=keep_number,$
	;		peak_thresh=peak_thresh, stddev_thresh=stddev_thresh
	;endif; rdi if
   
endif
if pre_inj eq 0 then begin
   
	print, 'pre_inj:', pre_inj

   ; two options, one for radius and angle and the other for x and y
   if keyword_set(rho) then inject_planets, obj, output_path, n_planets,$
   	contrast, pxscale_sx, pxscale_dx, keep_number, do_cen_filter,$
   	planet_r=rho, planet_theta=theta, use_gauss=use_gauss,$
   	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod,$
   	combine_type=combine_type
    
   if keyword_set(planet_x) then begin
   		print, planet_x, planet_y
      inject_planets, obj, output_path, n_planets, contrast, pxscale_sx,$
      pxscale_dx, keep_number, do_cen_filter, planet_y=planet_y, planet_x=planet_x,$
      use_gauss=use_gauss, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod,$
      silent=silent
   endif
	
	if adi eq 1 then begin
		print,'starting adi'
   		adi, obj, half_cropped_sz, nod, output_path, use_injection, filter, keep_number,$
   			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
   			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
   			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band, fwhm=fwhm,$
   			combine_type=combine_type, szbin=szbin, bin_type=bin_type,$
   			peak_thresh=peak_thresh, stddev_thresh=stddev_thresh,$
   			mask=mask, outmask=outmask, do_smooth=do_smooth
   	endif
   	
	if klip eq 1 then begin
		klip, obj, half_cropped_sz, nod, output_path, use_injection, do_destripe, filter, szbin, bin_type,$
			do_hyper, do_annmode, combine_type, klip_fraction, klip_start_frame,$
			klip_end_frame, fill, k_klip, angsep, anglemax, nrings, wr, n_ang,$
			annmode_inout_sx, annmode_inout_dx, suffix, keep_number, do_cen_filter, coadd,$
			trial=trial, fs=fs, neg_inj=neg_inj, truenorth_sx=truenorth_sx,$
			truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx,$
			magnify=magnify, fwhm=fwhm, wl=wl, do_rdi=do_rdi, refcube=refcube, ref_angles=ref_angles,basis=basis,$
			do_new_binning=do_new_binning
	endif; klip eq 1 if
	
	if do_rdi eq 1 then begin
		rdi, obj, half_cropped_sz, nod, output_path, use_injection, filter, szbin, bin_type,$
		do_annmode, combine_type, k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout_sx, annmode_inout_dx,$
		suffix, keep_number, do_cen_filter, coadd, trial=trial, fs=fs, neg_inj=neg_inj,$
		truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx,$
		pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl, refcube=refcube, ref_ang=ref_ang,$
		save_basis=save_basis, rdi_extra=rdi_extra,$
		mask_pt=mask_pt, do_create_basis=do_create_basis, do_new_binning=do_new_binning,$
		dxklip=dxklip, normal=normal, do_smooth=do_smooth, keep_number=keep_number,$
		peak_thresh=peak_thresh, stddev_thresh=stddev_thresh
	endif; rdi if
   
endif

;-----------------------------------[ El Fin ]--------------------------------------

; 1440 minutes per day (Julian dates are measured in days)
print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!

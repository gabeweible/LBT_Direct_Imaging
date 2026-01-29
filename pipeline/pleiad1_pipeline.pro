pro pleiad1_pipeline, outpath=outpath, coadd=coadd, extra=extra, nod=nod,$
	do_dewarp=do_dewarp, set_coadd_type=set_coadd_type, nod_filter=nod_filter,$
	cube_indices=cube_indices
	
; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
; Get the current time in a Juliasn date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

;--------------------
; BEGIN USER OPTIONS:
;--------------------

; CUBE creation parameters
if not keyword_set(set_coadd_type) then begin
	coadd_type='median'; seems like it may be a better option than mean, at least for
	; even coadds (where really the mean of the two central values is taken, which
	; is effectively a pretty good outlier-resistant mean for, e.g., coadd=4 or coadd=6)
endif else begin
	coadd_type = set_coadd_type
endelse
cube_start_frame = 0
; when do darks start in the co-added cube? 'None' if they are in a separate folder.
dark_frame_start = 'None';

; BADPIX parameters
; create master dark and master flat?
create_master_masks = 0
; create bad-pixel mast from master dark and master flat?
create_badpix_mask = 0
; removing column and row offsets.
do_destripe = 1; done in bad_pixels_fast.pro now...
boxhsize=1; 2*boxhsize x 2*boxhsize block of pixels for each comparison
sigma_clip=3.5; for 3x3, 5x5, and 7x7 filters in sky_sub.pro

; SKYSUB parameters
fpn = 1000; 'frames' per nod (really, ramps per nod position)
; doing sky subtraction after bad-px correction, now (0)
skysub_first = 0

; DEWARP/split parameters
do_smooth = 1.; FWHM for Gaussian smoothing before and after dewarping.
dewarp_bin_type='mean' ; how to bin again before dewarping, if relevant?

; frame half-size in pixels for initial cropping.
half_cropped_sz = 295

; General/Combine Parameters
obj = 'HIP17900'; OBJNAME in FITS header
band = 'L'
wl = 3.672131; from SVO (at 77 K, effective, not central)
raw_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD1/raw'
darks_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD1/darks'
stripe = 'Full_Image' ; NEW: I'll be working with full 2048 x 2048 frames up until cropping.
aperture = 'both'

;______________________________________
; OLD FROM TYC5709 - UPDATE FOR PLEIAD3
;______________________________________
fwhm = 10.6245 ; from VIP for TYC 5709 - to be updated for PLEIAD 3
; "In practice, D = 8.36 and Bmax = 22.65 m, because the adaptive secondaries in each LBT are purposefully undersized to control the thermal background." - Patru et al. 2017
lambda_over_d = 8.487; temporary - estimated for 10.61 mas/px nominal platescale

; update directory (coadd = 12 discards 4 frames per nod position, including the first at each position.)
; loss is less than 1 minute of integration time, so 99.6% efficiency. Even coadd is preferred so that
; the two middle values are actually averaged in the median.
output_path = '/Users/gweible/OneDrive - University of Arizona/research/PLEIAD_LBTI_data/PLEIAD1/gabe_macbook1/'

;______________________________________
; OLD FROM TYC5709 - UPDATE FOR PLEIAD3
;______________________________________
pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared — OLD)
pxscale_dx = 0.010624 ; +0.000035 or -0.000038 arcsec/pixel (Jared 06/20/25)
; These are required rotations E of N, i.e., CCW.
truenorth_sx = -1.278 ; +0.131 or -0.225 deg (from Steve and Jared — OLD)
truenorth_dx = 0.874 ; + 0.176 or -0.292 deg (Jared 06/20/25)
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
        
; AUTOMATIC VALUES (if not set to something otherwise)
;-----------------------------------------------------
if not keyword_set(coadd) then coadd = 1
if not keyword_set(do_dewarp) then do_dewarp=1 ; assume that we should dewarp
if not keyword_set(nod) then nod = 'both'
if keyword_set(outpath) then output_path=outpath
; choose a "run" for dewarp.pro if single-sided
if aperture eq 'right' then run = 2; DX/right aperture (only for dewarp.pro)
if aperture eq 'left' then run = 1; DX/right aperture (only for dewarp.pro)
if aperture eq 'both' then run = 5; for what used to be 'both'
if not keyword_set(nod_filter) then nod_filter = 'both' ; could also be set as 'dith1' or 'dith2'
;-----------------------------------------------------

;------------------------------[ Pipeline ]---------------------------------

; not actually using fixen algorithm up-the-ramp, just CDS in this case.
;cosmic_ray_rejection_cube, obj, raw_path, cube_start_frame, coadd,$
 ;   output_path, read_header=1, coadd_type=coadd_type, legacy_mode='cds',$
  ;    debug=1, dark_frame_start=dark_frame_start, max_frames_per_group=fix(fpn-1) / coadd,$
   ; min_r_squared=0.9, full_well=4095*0.98,skip_first_read=0,$
    ;skip_last_read=0, fpn=fpn, verbose=1, del_nodframe1=1, nod_counter_add=0

; create darks cube
;cosmic_ray_rejection_cube, obj, darks_path, cube_start_frame, coadd,$
 ;       output_path, read_header=1, coadd_type=coadd_type, legacy_mode='cds',$
  ;      debug=1, dark_frame_start=0,$
   ;     max_frames_per_group=fix(2*fpn - 1) / coadd, dark_cube=1,$
    ;    min_r_squared=0.9, full_well=4095*0.98,skip_first_read=0,$
     ;   skip_last_read=0, fpn=fpn, verbose=1, del_nodframe1=1

; not using bad_px_arr right now (find them all via computation)
;bad_pixels_fast, output_path, obj, stripe,$
;	boxhsize=boxhsize, pre_sky_sub=1, type='res_mean',$
;	run=run, do_destripe=do_destripe,$
;	fwhm=fwhm, just_destripe=0, debug=0,$
;	create_master_masks=create_master_masks,$
;	create_badpix_mask=create_badpix_mask, coadd=coadd,$
;	darks_filename='HIP17900_darks_cube.fits', vapp=0, post_pca_crop=0,$
;	start_nod=0

; NEW!! Sigma-clipping after bad-pixel-mask interpolation
; this is frame-to-frame
;sigma_clip_cubes, obj, output_path, cube_folder=output_path, $
 ;             sigma_clip=sigma_clip, debug=1, $
  ;            overwrite=1, frame_min=-40,$ ; ADJUST FRAME_MIN TO MATCH bad_pixels_fast (post destripe)!!!
   ;           cube_indices=[7,8,9]

; NEW!!! fill remaining NaNs with frame medians + MAD-estimated Gaussian
; noise
;fill_nan_cubes, obj, output_path, cube_folder=output_path,$
 ;   debug=debug, overwrite=1;, cube_indices=[18,19,20,21,22,23,24,25]
;	
; SKY SUBTRACTION IN 'pleiad1_pca_bkg_subtr_1cell.ipynb' !!!!

; this isn't actually bad-px correction, I'm just destriping
; the sky-subtracted frames, now.
;bad_pixels_fast, output_path, obj, stripe,$
 ;   boxhsize=boxhsize, pre_sky_sub=0, type='res_mean',$
  ;  run=run, do_destripe=do_destripe,$
   ; fwhm=fwhm, just_destripe=0, debug=0,$
;    create_master_masks=create_master_masks,$
 ;   create_badpix_mask=create_badpix_mask, coadd=coadd,$
  ;  darks_filename='HIP17900_darks_cube.fits', vapp=0,$
   ; destripe_skysub=1, pca_skysub=1, post_pca_crop=1, start_nod=0,$
    ;skip_badpix_mask=1

fill_nan_cubes, obj, output_path, cube_folder=output_path,$
		debug=debug, overwrite=1, input_suffix='_30comp_corrected_cube.fits',$
		output_suffix='_pca_skysub_filled_cube.fits',$
		input_prefix='test_pca_skysub_cube'
    
;if do_dewarp eq 1 then begin

	; new combined dewarping and splitting to reduce memory requirements.
;	dewarp_split_combined, output_path, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx,$
;		do_smooth=do_smooth, half_cropped_sz=half_cropped_sz, aperture=aperture, run=run,$
;		fwhm=fwhm, hp_width=0, nod_filter=nod_filter, debug=1,$
;		destripe_skysub=1, pca_skysub=1, filled=1
        
;endif; else begin; dewarp if
;	; 1 is for skysub_first
;	split, obj, stripe, output_path, half_cropped_sz, aperture, skysub_first, do_dewarp,$
;		fwhm=fwhm
;endelse

; Centering, cleaning, PSF subtr. moved to VIP pipeline.

;-----------------------------------[ El Fin ]--------------------------------------

; 1440 minutes per day (Julian dates are measured in days)
print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!

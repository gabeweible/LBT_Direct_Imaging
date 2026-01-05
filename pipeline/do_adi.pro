FUNCTION DO_ADI, obj_name, half_cropped_sz, cube_folder, use_injection, filter, keep_number,$
	do_cen_filter, coadd, adi_cube_nframes, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs, fwhm=fwhm,$
	combine_type=combine_type, szbin=szbin, bin_type=bin_type,$
	peak_thresh=peak_thresh, stddev_thresh=stddev_thresh,$
	adi_suffix=adi_suffix, do_smooth=do_smooth
	
compile_opt idl2
newline = string(10B)

sz = 2*half_cropped_sz
xs = sz & ys = sz
xhs = xs / 2. & yhs = ys / 2.

if do_smooth gt 0 then lp_PSF = psf_Gaussian(NPIX=11, FWHM=[do_smooth, do_smooth], /normalize, /double)

; Pre-compute coordinate grid once
xx_grid = rebin(indgen(fix(xs)), fix(xs), fix(ys))
yy_grid = rebin(reform(indgen(fix(ys)), 1, fix(ys)), fix(xs), fix(ys))

if not keyword_set(szbin) then szbin=1
if not keyword_set(bin_type) then bin_type='median'

;Do this for runs eq 1 and runs eq 2
if runs lt 3 then begin
	output_folder = cube_folder + 'processed_left/'
	truenorth = truenorth_sx
endif else begin
	truenorth = truenorth_dx
	output_folder = cube_folder + 'processed_right/'
endelse

; Do this for runs eq 1 and runs eq 3
if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'

if use_injection then begin
	print, 'Reading post-injection cube...'
	obj_cube = readfits_fast(strcompress(output_folder + dither_folder + obj_name +$
   		'_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+'_stddev_thresh_'+$
   		string(stddev_thresh) + '_cube_skysub_cen_clean_inj.fits', /r))
   		
   inj_string = '_inj'
endif else begin
	print, 'Reading post-clean cube...'
   	obj_cube = readfits_fast(strcompress(output_folder + dither_folder + obj_name +$
   		'_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+'_stddev_thresh_'+$
   		string(stddev_thresh) + '_cube_skysub_cen_clean.fits', /r))
   	
   inj_string = ''
endelse


; restore parallactic angles
restore, filename = strcompress(output_folder + dither_folder + obj_name + '_' + string(keep_number) +$
	'_peak_thresh_'+string(peak_thresh)+'_stddev_thresh_'+string(stddev_thresh) +$
	'_parang_clean.sav', /r)

; Efficient binning
processed_obj_cube = obj_cube

if keyword_set(szbin) then begin

	if szbin ge 2 then begin
		; Optimized binning using array reshaping
		print, 'Binning science frames...'
		n_frames = (size(processed_obj_cube, /dim))[2]
		n_bins = floor(n_frames / szbin)
		
		; Reshape and reduce reference cube
		reshaped_cube = reform(processed_obj_cube[*,*,0:n_bins*szbin-1], $
							   (size(processed_obj_cube, /dim))[0], $
							   (size(processed_obj_cube, /dim))[1], $
							   szbin, n_bins)
		
		; Choose binning method
		if bin_type eq 'median' then begin
			processed_obj_cube = median(reshaped_cube, dimension=3, /double, /even)
		endif
		if bin_type eq 'mean' then begin
			processed_obj_cube = mean(reshaped_cube, dimension=3, /nan, /double)
		endif
		if bin_type eq 'res_mean' then begin
			RESISTANT_Mean, reshaped_cube, 3.5, processed_obj_cube, dim=3
		endif
		
		; Average angles for binned frames
		processed_obj_angles = total(reform(angles[0:n_bins*szbin-1], szbin, n_bins), 1) / szbin
	endif
	
	delvar, reshaped_cube; free up some memory
endif

n_frames = (size(processed_obj_cube))[3]
adi_cube_nframes = n_frames

	if filter ge 1 then begin; just high-pass filtering
		hp_sz = round(5*filter)
		if odd(hp_sz) eq 0 then hp_sz += 1
		hp_PSF = psf_Gaussian(NPIX=hp_sz, FWHM=[filter, filter], /double, /normalize)
		
		print, 'High-pass filtering science frames...'
		for iii = 0, (size(processed_obj_cube, /dim))[2]-1 do begin
			print, 'Filtering frame ', iii
			frame_iii = processed_obj_cube[*,*,iii]
			
			; Replace NaNs with zeros
    		nan_indices = where(finite(frame_iii) ne 1, nan_count)
    		if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
			
			; high-pass/unsharp mask
			processed_obj_cube[*,*,iii] -= convolve(frame_iii, hp_PSF)
			
			; try a median instead?
			;processed_obj_cube[*,*,iii] -= filter_image(frame_iii,$
			;	smooth=filter,$ ;smooth=filter,$
         	;	/ALL_PIXELS)
		endfor
	endif
;endelse

;if keyword_set(normal) then begin
;	; grab pixel values within a radius of 
;   		center_pixels = []
;		for dx=-fwhm*(1.5),fwhm*(1.5) do begin
;		for dy=-fwhm*(1.5),fwhm*(1.5) do begin
;			if sqrt( dx^2. + dy^2. ) le fwhm*1.5 then center_pixels = [[center_pixels],$
;				medframe[xhs+dx, yhs+dy]]
;		endfor
;		endfor
;		
;		med_peak_mean = mean(center_pixels, /double, /nan)
;endif

;Initialize our cube for Angular Differential Imaging
adi_cube = processed_obj_cube
n_processed_frames = (size(processed_obj_cube))[3]

if normal eq 1 then begin; only normalize to the smoothed peak	
	for ii=0, n_processed_frames - 1 do begin
	   
		; grab frame for PSF subtraction
	   	adi_frame_ii = adi_cube[*,*,ii]
	   
		; try a more robust normalization?
		if ii mod 100 eq 0 then print, 'normalzing adi_frame' + string(ii) + '...'
		; grab pixel values within a radius of 0.5 FWHM
		center_pixels = []
		for dx=-fwhm*3.0,fwhm*3.0 do begin
			for dy=-fwhm*3.0,fwhm*3.0 do begin
				if sqrt( dx^2. + dy^2. ) le 3.0*fwhm then center_pixels = [[center_pixels],$
					adi_frame_ii[xhs+dx, yhs+dy]]
			endfor
		endfor
			
		ii_peak_mean = mean(center_pixels, /double, /nan)
		
		; scale a unique medframe normalization to frame_ii mean flux in IWA
		; 100. just rescales to make the numbers nicer.
		scaled_frame_ii = adi_frame_ii / (ii_peak_mean / 100.)
		adi_cube[*,*,ii] = scaled_frame_ii
	endfor
	;medarr, scaled_medframe_cube, medframe; new medframe from scaling one for each 
	; ADI image
endif; normalization if

medarr, adi_cube, medframe; to be subtracted from the image (cADI)

; PSF SUBTRACTION
for ii=0, n_processed_frames - 1 do begin

	frame_ii = adi_cube[*,*, ii]
	
   	; Subtract the median frame (star PSF) from the new adi_cube frame
   	frame_ii -= medframe
   	
   	median_ii = median(frame_ii, /double, /even)
   	; Rotate the (PSF-subtracted) adi_cube frame to truenorth
  	; (negatives because these are E of N, a.k.a. CCW, ROT does CW)
  	if do_smooth gt 0 then begin
  		nan_indices = where(finite(frame_ii) ne 1, nan_count)
		if nan_count gt 0 then frame_ii[nan_indices] = median_ii
  		frame_ii = convolve(temporary(frame_ii), lp_PSF)
	endif
  	
  	;print, 'Rotating adi_cube about', half_cropped_sz, half_cropped_sz
	adi_cube[*,*,ii] = rot(frame_ii, -angles[ii] - truenorth,$
   		1.0, half_cropped_sz, half_cropped_sz, cubic=-1.0, missing=median_ii, /pivot)
   
    ; Rotate obj_cube frame to truenorth
   	; (negatives because these are E of N, a.k.a. CCW, ROT does CW)
   	;print, 'Rotating processed_obj_cube about', half_cropped_sz, half_cropped_sz
   	processed_obj_cube[*,*,ii] = rot(temporary(processed_obj_cube[*,*,ii]),$
   		-angles[ii] - truenorth, 1.0, half_cropped_sz, half_cropped_sz, cubic=-1.0,$
   		missing=median(processed_obj_cube[*,*,ii], /double, /even), /pivot)
   		
endfor; rotate if

if do_smooth gt 0 then begin; just low-pass filtering

	print, 'Low-pass filtering science frames...'
	for iii = 0, (size(processed_obj_cube, /dim))[2]-1 do begin
		if iii mod 100 eq 0 then print, 'Filtering: ', iii, '/', (size(processed_obj_cube, /dim))[2]-1, ' in do_adi'
		frame_iii = processed_obj_cube[*,*,iii]
		adi_frame_iii = adi_cube[*,*,iii]
		median_iii = median(adi_frame_iii)
		
		; Replace NaNs with zeros
		nan_indices = where(finite(frame_iii) ne 1, nan_count)
		if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
		
		nan_indices = where(finite(adi_frame_iii) ne 1, nan_count)
		if nan_count gt 0 then adi_frame_iii[nan_indices] = median_iii
		
		; low-pass
		processed_obj_cube[*,*,iii] = convolve(frame_iii, lp_PSF)
		adi_cube[*,*,iii] = convolve(adi_frame_iii, lp_PSF)
	endfor
endif

; Write the de-rotated obj_cube (No ADI)
;writefits, strcompress(adi_suffix + '_cube_skysub_cen_filt_derot.fits', /rem),$
;	processed_obj_cube
	
; something weird is messing up the centering just between these two steps...
medarr, processed_obj_cube, derot_medframe

if do_smooth gt 0 then begin
	nan_indices = where(finite(derot_medframe) ne 1, nan_count)
	if nan_count gt 0 then derot_medframe[nan_indices] = median(derot_medframe, /double, /even)
	
	derot_medframe = convolve(temporary(derot_medframe), lp_PSF)
endif

writefits, strcompress(adi_suffix + '_median_derot.fits', /rem), derot_medframe

if combine_type eq 'median' then medarr, adi_cube, adiframe
if combine_type eq 'mean' then adiframe = mean(adi_cube, dim=3, /double, /nan)
if combine_type eq 'nwadi' then adiframe = nw_ang_comb(adi_cube, angles+truenorth,$
	do_smooth=do_smooth)

return, adiframe

END; do_adi function end
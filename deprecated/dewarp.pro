pro dewarp, dir, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx, run=run,$
	skysub_first=skysub_first, do_smooth=do_smooth,$ ; sigma_clip=sigma_clip,$
	fwhm=fwhm, half_cropped_sz=half_cropped_sz, second_coadd=second_coadd,$
	bin_type=bin_type, dewarp_filter=dewarp_filter, hp_width=hp_width


; PURPOSE: dewarps LBTI LMIRCAM images using coefficients Kx_dx, Ky_dx, Kx_sx,
; and Ky_sx (update, if needed) for two-mirror, non-interferometric observations,
; resulting in two new cubes being written, one that was dewarped with the sx 
; mirror coefficients, and the other with the dx mirror coefficients, which then
; will need to be differentiated in split.pro.
;
; ARGUMENTS:
; dir: string directory to find the sky-subtraced cube to read in and for 
; writing our dewarped cubes to
;
; obj: string target object name (probably just the host star if directly imaging
; an exoplanet or exoplanets
; (or brown dwarfs! :))
;
; RETURNS: None, but will write two dewarped datacubes to dir (one for the sx
; mirror and one for the dx mirror)
;
; AUTHOR: Gabe Weible, Jan 2023, Undisclosed location

; run=1 for SX only, run=2 for DX only

COMPILE_OPT IDL2; Strictarr, 32 bit integers by default
newline = string(10B)

print, 'Dewarping cube...'

if not keyword_set(second_coadd) then second_coadd = 0
if not keyword_set(skysub_first) then skysub_first = 0
if not keyword_set(bin_type) then bin_type = 'res_mean'

; dir =  '/Users/gabeweible/OneDrive/research/HII1348/macbook_25/'
; obj = 'HII1348'

; Setup
; Read in our sky-subtracted cube to dewarp
if skysub_first eq 0 then begin
	print, 'reading skysub cube...'
	path = dir + obj + '_cube_skysub.fits'
	print, 'cube at:', path 
	obj_cube = readfits_fast(path)
	
	;writefits, '~/Desktop/skysub_cube_frame_0.fits', obj_cube[*,*,0]
endif

;if skysub_first eq 1 then print, 'reading_bad_px_cube...' & obj_cube = readfits_fast(dir + obj + '_cube_no_bad_pixels.fits')

print, size(obj_cube)
frames = (size(obj_cube))[3]; Total # of images in each datacube
xdim = (size(obj_cube))[1]
ydim = (size(obj_cube))[2]

; read-in old angles (UPDATE TO SKY_SUB PARANG AFTER RUNNING SKY_SUB AGAIN...)
restore, strcompress(dir+obj+'_parang_sky_sub.sav',/rem)
if second_coadd ge 2 then begin

	new_flags = list(); new flag list to fill when binning
	coadd_frame = [[]]; will be filled with each coadd
	coadd_angle = 0.; will be added to with each coadd
	binned_obj_cube = list()
	new_angles = list()
	
	szbin = second_coadd
	; Optimized binning using array reshaping
	print, 'Binning science frames...'
	
	k = 1; frame counter
	for ii = 0, frames-1 do begin
		print, 'File index', ii, '/', frames
		
		; grab the frame and flag and parallactic angle for index ii
		frame = obj_cube[*,*,ii]
		flag = flags[ii]
		angle = angles[ii]
		
		; Add to our current coadd frame as a sum, and the same with the current coadd angle
		coadd_frame = [ [[coadd_frame]], [[frame]] ]; add the frame to the small group
		coadd_angle = temporary(coadd_angle) + angle
		if k mod second_coadd eq 0 then begin
			
			new_flags.Add, flag
			coadd_angle = temporary(coadd_angle) * (1. / second_coadd)
			new_angles.Add, coadd_angle
			coadd_angle = 0.
			
			;Add the mean frame to the cube and reset the frame, and add the mean angle and reset the angle
			; Choose binning method
			if bin_type eq 'median' then begin
				coadd_frame = median(temporary(coadd_frame), dimension=3, /even, /double)
			endif
			
			if bin_type eq 'mean' then begin
				coadd_frame = mean(temporary(coadd_frame), dimension=3, /nan, /double)
			endif
			
			if bin_type eq 'res_mean' then begin
				resistant_mean, temporary(coadd_frame), 5.0, coadd_frame, dim=3, /double
			endif
			
			; add to the binned cube
			binned_obj_cube.Add, [[coadd_frame]]
		
			; reset for coadding
			coadd_frame = [[]]
			
		endif
		
		; keep track of binning and increment the counter
	   	print, 'Number of Frames in Cube: ', n_elements(binned_obj_cube)
	   	k += 1
	   	if ii eq 0 then writefits, '~/Desktop/dewarp_bin_test.fits', coadd_frame
		
	endfor ; ii frame loop for
	
	; replace and then delete the duplicates.
	obj_cube = binned_obj_cube.toArray(/TRANSPOSE, /NO_COPY)
	angles = new_angles.toArray(/NO_COPY)
	flags = new_flags.toArray(/NO_COPY)
	
	delvar, binned_obj_cube, angles, flags
	
endif

; save new angles
save, filename=strcompress(dir+obj+'_parang_dewarp.sav',/rem), angles, flags
print, 'Save file written!'

; Pre-compute normalized Gaussian for low-pass filtering
if do_smooth gt 0 then begin
	lp_npix = fix(11*do_smooth)
	if ~ODD(lp_npix) then lp_npix += 1
	lp_PSF = psf_Gaussian(NPIX=lp_npix, FWHM=[do_smooth, do_smooth], /normalize)
endif

if hp_width gt 0 then begin
	hp_npix = fix(11*hp_width)
	if ~ODD(hp_npix) then hp_npix += 1
	hp_PSF = psf_Gaussian(NPIX=hp_npix, FWHM=[hp_width, hp_width], /normalize)
endif

; second sigma-clipping after the second binning.
if do_smooth eq 1 then begin
	for ii = 0,(size(obj_cube))[3]-1 do begin
		print, 'low-pass filtering frame', ii
	
		frame_ii = obj_cube[*,*,ii]
	
		nan_indices = where(finite(frame_ii) ne 1, nan_count)
		if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /nan)
		frame_ii_smooth = convolve(frame_ii, lp_PSF)
		
		; High-pass filter
    	if hp_width gt 0 then frame_ii_smooth -= convolve(frame_ii_smooth, hp_PSF)
    	; low-pass again just to be safe.
    	frame_ii_smooth = convolve(frame_ii_smooth, lp_PSF)
    	
		
		; put the filtered image back into the cube
		obj_cube[*,*,ii] = frame_ii_smooth
		; debugging
		if ii eq 0 then writefits, '~/Desktop/smooth_clip_test.fits', obj_cube[*,*,ii]
	endfor
endif

if not keyword_set(run) then run = 'both'

if run eq 'both' then begin
	; 2 runs because we haven't split the data yet, but do have two mirrors with 
	; separate distortion coefficients
	for run = 1,2 do begin
		print, 'Starting run: ' + string(run)
		; different distortion coefficient matrices for each mirror
		if run eq 1 then begin; left mirror (top PSF)
			mirror = 'SX'
			Kx = Kx_sx
			Ky = Ky_sx
		endif
	
		if run eq 2 then begin; right mirror (bottom PSF)
			mirror = 'DX'
			Kx = Kx_dx
			Ky = Ky_dx
		endif

		print,"Everything's loaded in - let's dewarp!"
	
		; Initialize lists for our cube and angles
		new_cube = list()
		warped_cube = list()
	
		for ii=0, frames-1 do begin
		  print, 'Working on frame: ' + string(ii + 1) + ' / ' + string(frames)
		  print, 'side: ' + mirror + newline
	  
		  ;Grab a frame from our cube to work on
		  warped = obj_cube[*,*,ii]
		  padded = replicate(!VALUES.F_NAN, 2048, 2048); full frame
		  ; insert warped center stripe into full-frame padding
		  ; 512 px (up to index 511) are NaN, 513th px (index 512) is the first with
		  ; our actual data and not NaNs
		  if stripe eq 'center1024' then begin
		  		padded[0, 512] = warped; location of bottom-left corner of observations
		
		  endif else if stripe eq 'second512' then begin
		  		padded[0, 1024] = warped; second stripe of 512 down from the top?
		  endif
	  
		  ; Test before dewarping
		  warped_cube.add, [[padded]]
	  		
	  		;shifting_frame = filter_image(padded, SMOOTH=2., /ALL_PIXELS)
		  dewarped = POLY_2D(temporary(padded), Kx, Ky, 2, cubic=-1.0); dewarp w/coefficient matrices
		  
		  ; low-pass filter. Should help cubic convolution.
	  		if do_smooth eq 1 then begin
	  			nan_indices = where(finite(dewarped) ne 1, nan_count)
				if nan_count gt 0 then dewarped[nan_indices] = 0.
				dewarped = convolve(temporary(dewarped), lp_PSF)
			endif
	  
		  ; Add to lists
		  new_cube.Add, [[dewarped]]
		endfor; ii for

	;--------------------------------------------------------------------------------------------------------

		print, 'Converting new cube ' + string(run) + ' to array...'
		new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
	
		; Test
		warped_cube = warped_cube.toArray(/TRANSPOSE, /NO_COPY)
	
		print, 'New cube converted to array! Writing cube...' + newline
		writefits, dir + obj + '_dewarped_' + mirror + '.fits', new_cube
	
		; Test
		writefits, dir + obj + '_warped' + '.fits', warped_cube
	
		print, mirror + ' cube written!'

	endfor; run for (1 and 2 for dual-sided observations)
	
; run=1 for SX (left) aperture, run=2 for DX (right) aperture
endif else if (run eq 1) or (run eq 2) then begin; for single-sided observations

	; frames needs re-defined after binning.
	print, size(obj_cube)
	frames = (size(obj_cube))[3]; Total # of images in each datacube
	
	if run eq 1 then begin; left mirror (top PSF)
		mirror = 'SX'
		Kx = Kx_sx
		Ky = Ky_sx
	endif

	if run eq 2 then begin; right mirror (bottom PSF)
		mirror = 'DX'
		Kx = Kx_dx
		Ky = Ky_dx
	endif

	print,"Everything's loaded in - let's dewarp!"

	; Initialize lists for our cube and angles
	new_cube = list()
	;warped_cube = list()

	for ii=0, frames-1 do begin
	  	print, 'De-warping working on frame: ' + string(ii + 1) + ' / ' + string(frames)
	  	print, 'side: ' + mirror + newline
  
	  	padded = replicate(!VALUES.F_NAN, 2048, 2048); full frame
	  	; insert warped center stripe into full-frame padding
	  	; 512 px (up to index 511) are NaN, 513th px (index 512) is the first with
	  	; our actual data and not NaNs
	  	if stripe eq 'center1024' then begin
		  	padded[0, 512] = obj_cube[*,*,ii]; location of bottom-left corner of observations
	
	  	endif else if stripe eq 'second512' then begin
			padded[0, 1024] = obj_cube[*,*,ii]; second stripe of 512 down from the top?
	  	endif
	  	
	  	if ii eq 0 then writefits, '~/Desktop/padded_test.fits', padded
  
	  	; Test before dewarping
	  	;warped_cube.add, [[padded]]
  
  		;shifting_frame = filter_image(padded, SMOOTH=2., /ALL_PIXELS)
	  	dewarped = POLY_2D(padded, Kx, Ky, 2, cubic=-0.5,$
	  		missing=!values.f_nan); dewarp w/coefficient matrices
	  
	  	; low-pass filter. Should help cubic convolution.
	  	if do_smooth eq 1 then begin
	  		nan_indices = where(finite(dewarped) ne 1, nan_count)
			if nan_count gt 0 then dewarped[nan_indices] = 0.
			dewarped = convolve(temporary(dewarped), lp_PSF)
		endif
  
	  ; Add to lists
	  	new_cube.Add, [[dewarped]]
	  	
	  	if ii eq 0 then writefits, '~/Desktop/dewarped_test.fits', dewarped
	endfor; ii for
	
	delvar, obj_cube; free up some memory

;--------------------------------------------------------------------------------------------------------

	print, 'Converting new cube ' + string(run) + ' to array...'
	new_cube = (temporary(new_cube)).toArray(/TRANSPOSE, /NO_COPY)

	; Test
	;warped_cube = warped_cube.toArray(/TRANSPOSE, /NO_COPY)

	print, 'New cube converted to array! Writing cube...' + newline
	writefits, dir + obj + '_dewarped_' + mirror + '.fits', new_cube

	; Test
	;writefits, dir + obj + '_warped' + '.fits', warped_cube

	print, mirror + ' cube written!'
endif; single-sided if (run eq 1) or (run eq 2)

end

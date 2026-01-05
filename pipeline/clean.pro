pro clean, obj_name, stripe, cube_folder, keep_number, half_cropped_sz, aperture,$
	annular_clean=annular_clean,centered_clean=centered_clean, peak_clean=peak_clean, peak_thresh=peak_thresh,$
	stddev_clean=stddev_clean, stddev_thresh=stddev_thresh, fwhm=fwhm,$
	do_dewarp=do_dewarp, do_block_right, do_block_left,$
	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub, do_smooth=do_smooth
	
newline = string(10B)
compile_opt idl2

; high-pass unsharp mask size
new_filter=27.0

; 7x7 should be plenty for fwhm =1.5
;if do_smooth gt 0 then lp_PSF = psf_Gaussian(NPIX=11, FWHM=[do_smooth, do_smooth], /normalize, /double)

if not keyword_set(annular_clean) then annular_clean = 0;
if not keyword_set(centered_clean) then centered_clean=1
if not keyword_set(stddev_clean) then stddev_clean=0
if not keyword_set(peak_clean) then peak_clean=0
if not keyword_set(stddev_thresh) then stddev_thresh=0
if not keyword_set(peak_thresh) then peak_thresh=0

sz = 2*half_cropped_sz + 1

;HII 1348 (SX AND DX)
if (stripe eq 'center1024') and aperture eq ('both') then for runs=1,4 do begin

	; Do this for runs eq 1 and runs eq 3
	if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'

	;Do this for runs eq 1 and runs eq 2
	if runs lt 3 then output_folder = cube_folder + '/processed_left/' else output_folder = cube_folder + '/processed_right/'

	print, 'Finding bad frames...'
	cube = readfits_fast(output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', hdr)
	
	restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'

	;192,195

	bad_pixels = []
	goods = findgen((size(cube))[3])

	;initialize a new cube to clean
	cube2 = cube
	
	;remove saturated center
	;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.
	
	if annular_clean eq 1 then begin
		
		for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) le 10. then cube2[xx,yy,*]=0.
		for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) ge 50. then cube2[xx,yy,*]=0.
	endif
	
	; PSF core + 1st Airy ring
	if centered_clean eq 1 then begin
		for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) ge 35. then cube2[xx,yy,*]=0.
	endif

	;compute pupil median
	lmed = median(cube2, dim=3, /double, /even)

	;First frame special case
	maxcor = max(crosscorr(cube2[*,*,0], lmed, /cut))
	print, 'Frame ', 0, ' has max correlation value of ', maxcor, ' with pupil median.'
	;First frame vs. the rest
	corrs = maxcor
	
	; crop everything once
	peak_cropped = cube2[half_cropped_sz-fwhm*3.0:half_cropped_sz+fwhm*3.0,half_cropped_sz-fwhm*3.0:half_cropped_sz+fwhm*3.0,*]
	
	; compute first fwhm x fwhm resistant mean and standard deviation (resistant mean, not resistant standard deviation)
	frame0 = peak_cropped[*,*, 0]
	resistant_mean, frame0, 3.5, peaks, /double
	frame_stddevs = stddev(frame0, /double, /nan)
	
	; stdev cleaning
	;frame_stddevs = stddev(cube2[*,*,0], /double, /nan)

	for jj=1, (size(cube))[3] - 1 do begin
	
		frame_jj = cube2[*,*,jj]
		
		maxcor = max(crosscorr(frame_jj, lmed, /cut))
		print, newline, 'Frame ', jj, ' has max correlation value of ', maxcor, ' with pupil median.'
		;First frame vs. the rest
		corrs = [corrs, maxcor]
		
		; resistant mean, but a non-resistant standard deviation
		resistant_mean, frame_jj, 3.5, peaks_jj, /double
		frame_stddevs_jj = stddev(frame_jj, /double, /nan)
		
		peaks = [peaks, peaks_jj]
		frame_stddevs = [frame_stddevs, frame_stddevs_jj]
		
		if peak_clean eq 1 then print, 'Peak flux / median peak flux = ',peaks_jj/median(peaks)
		if stddev_clean eq 1 then print, 'Frame stddev / median frame stddev', frame_stddevs_jj/median(frame_stddevs)
	endfor; max correlation for
	
	delvar, cube2 ; no longer needed.

	; select by max. cross correlation
	bad_pixels = where(corrs lt keep_number)
	
	; select by ratio of peak flux to median peak flux
	if peak_clean eq 1 then bad_pixels=[bad_pixels,where(peaks/median(peaks, /double, /even) lt peak_thresh, /null)]
	
	; select by ratio of frade stddev to mean frame stddev
	if stddev_clean eq 1 then bad_pixels=[bad_pixels, where( frame_stddevs/median(frame_stddevs, /double, /even) gt stddev_thresh, /null )]
	
	bad_pixels = bad_pixels[rem_dup(temporary(bad_pixels))]

	print, bad_pixels
	print, 'Found ', n_elements(bad_pixels), ' bad frames.'
	print, 'Average correlation = ', mean(corrs, /double, /nan)
	print, 'Median correlation = ', median(corrs, /double, /even)
	print, 'Standard deviation = ', stddev(corrs, /double, /nan)
	;hak
	print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'

	if n_elements(bad_pixels) gt 1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', left_bad_pixels_cube
	endif else if n_elements(bad_pixels) eq 1 and bad_pixels ne -1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', left_bad_pixels_cube
	endif

	cube = (temporary(cube))[*,*,goods]; Reassign our cube to only have our good pixels
	medarr, cube, medframe; Get a median frame for our new bad-pixel-free cube

	;Write our clean cube, our bad pixels, and our median frame
	writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_cube_skysub_cen_clean.fits', cube, hdr
	writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_cube_skysub_cen_bad_pixels.fits', left_bad_pixels_cube, hdr
	writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) + '_pupil.fits', medframe
	save, filename = output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_parang_clean.sav', angles

	print, 'Run: ' + string(runs) + ' complete
endfor	; runs for

; ALCOR (DX ONLY)
if stripe eq 'second512' then begin
	if aperture eq 'left' then for runs=1,2 do begin
		; Do this for runs eq 1 and runs eq 3
		if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
	
		output_folder = cube_folder + '/processed_left/'
	
		print, 'Finding bad frames...'
		cube = readfits_fast(output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', hdr)
		restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
	
		bad_pixels = []
		goods = findgen((size(cube))[3])
	
		;initialize a new cube to clean
		cube2 = cube
		;remove saturated center
		;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.
	
		if annular_clean eq 1 then begin
			for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) le 10. then cube2[xx,yy,*] = 0.
			for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) ge 50. then cube2[xx,yy,*] = 0.
		endif
		
		; PSF core + 1st Airy ring
		if centered_clean eq 1 then begin
			for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) ge 50. then cube2[xx,yy,*] = 0.
		endif
	
		;compute pupil median
		lmed = median(cube2, dim=3, /double, /even)
	
		;writefits, output_folder + dither_folder + obj_name + '_test-med_' + '_' + string(keep_number) +  '_corrthresh.fits', lmed
	
		;First frame special case
		maxcor = max(crosscorr(cube2[*,*,0], lmed, /cut))
		print, 'Frame ', 0, ' has max correlation value of ', maxcor, ' with pupil median.'
		;First frame vs. the rest
		corrs = maxcor
		
		; crop everything once
		peak_cropped = cube2[half_cropped_sz-fwhm*3.:half_cropped_sz+fwhm*3.,half_cropped_sz-fwhm*3.:half_cropped_sz+fwhm*3.,*]
		
		; compute first fwhm x fwhm resistant mean and standard deviation (resistant mean, not resistant standard deviation)
		frame0 = peak_cropped[*,*, 0]
		resistant_mean, frame0, 3.5, peaks, /double
		frame_stddevs = stddev(frame0, /double, /nan)
		
		; stdev cleaning
		;frame_stddevs = stddev(cube2[*,*,0], /double, /nan)
	
		for jj=1, (size(cube))[3] - 1 do begin
	
			frame_jj = cube2[*,*,jj]
			
			maxcor = max(crosscorr(frame_jj, lmed, /cut))
			if jj mod 100 eq 0 then print, newline, 'Frame ', jj, ' has max correlation value of ', maxcor, ' with pupil median.'
			;First frame vs. the rest
			corrs = [corrs, maxcor]
			
			; resistant mean, but a non-resistant standard deviation
			resistant_mean, frame_jj, 3.5, peaks_jj, /double
			frame_stddevs_jj = stddev(frame_jj, /double, /nan)
			
			peaks = [peaks, peaks_jj]
			frame_stddevs = [frame_stddevs, frame_stddevs_jj]
			
			if peak_clean eq 1 then print, 'Peak flux / median peak flux = ',peaks_jj/median(peaks)
			if stddev_clean eq 1 then print, 'Frame stddev / median frame stddev', frame_stddevs_jj/median(frame_stddevs)
		endfor; max correlation for
	
		delvar, cube2; no longer needed.
	
		bad_pixels = where(corrs lt keep_number)
		if peak_clean eq 1 then bad_pixels=[bad_pixels,where(peaks/median(peaks, /double, /even) lt peak_thresh, /null)]
		
		; select by ratio of frade stddev to mean frame stddev
		if stddev_clean eq 1 then bad_pixels=[bad_pixels ,where( frame_stddevs/median(frame_stddevs, /double, /even) gt stddev_thresh, /null )]
	
		bad_pixels = bad_pixels[rem_dup(temporary(bad_pixels))]
		print, bad_pixels
		print, 'Found ', n_elements(bad_pixels), ' bad frames.'
		print, 'Average correlation = ', mean(corrs, /double, /nan)
		print, 'Median correlation = ', median(corrs, /double, /even)
		print, 'Standard deviation = ', stddev(corrs, /double, /nan)
		;hak
		print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'
	
		if n_elements(bad_pixels) gt 1 then begin
			remove, bad_pixels, goods, angles
			left_bad_pixels_cube = cube[*,*,bad_pixels]
			;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', left_bad_pixels_cube
		endif else if n_elements(bad_pixels) eq 1 and bad_pixels ne -1 then begin
			remove, bad_pixels, goods, angles
			left_bad_pixels_cube = cube[*,*,bad_pixels]
			;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', left_bad_pixels_cube
		endif
	
		cube = (temporary(cube))[*,*,goods]; Reassign our cube to only have our good pixels
		medarr, cube, medframe; Get a median frame for our new bad-pixel-free cube
	
		;Write our clean cube, our bad pixels, and our median frame
		writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_cube_skysub_cen_clean.fits', cube, hdr
		writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_cube_skysub_cen_bad_pixels.fits', left_bad_pixels_cube, hdr
		writefits, output_folder + dither_folder + obj_name + '_' + string(keep_number) + '_pupil.fits', medframe
		save, filename = output_folder + dither_folder + obj_name + '_' + string(keep_number) +  '_parang_clean.sav', angles
	
		print, 'Run: ' + string(runs) + ' complete
	endfor; SX only
	
	if aperture eq 'right' then begin
		dith1_cube = readfits_fast(cube_folder + '/processed_right/' + 'dith1/' + obj_name + '_cube_skysub_cen.fits', hdr)
		dith2_cube = readfits_fast(cube_folder + '/processed_right/' + 'dith2/' + obj_name + '_cube_skysub_cen.fits', hdr)
	endif
	
	if aperture eq 'right' then for runs=3,4 do begin
		; Do this for runs eq 1 and runs eq 3
		if runs mod 2 then begin
			dither_folder = 'dith1/'
			not_dither_folder = 'dith2/'
			
			cube = dith1_cube[0:sz-1,0:sz-1,*]
			;not_cube = dith2_cube
		endif else begin
			dither_folder = 'dith2/'
			not_dither_folder = 'dith1/'
			
			cube = dith2_cube[0:sz-1,0:sz-1,*]
			;not_cube = dith1_cube
		endelse
		sz -= 1; cropped so there are an even # of pixels in each image.
	
		output_folder = cube_folder + '/processed_right/'
		
		print, 'Finding bad frames...'
		
		restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
		if runs eq 3 then begin
			angles = dx_dith1_angles_array
		endif else begin
			angles = dx_dith2_angles_array
		endelse
	
		bad_pixels = []
		n_cube_frames = (size(cube))[3]
		;n_not_cube_frames = (size(not_cube))[3]
		goods = findgen(n_cube_frames)
	
		;initialize a new cube to clean
		;cube2 = dblarr(2*half_cropped_sz, 2*half_cropped_sz, n_cube_frames + n_not_cube_frames, /nozero)
		;cube2[*,*,0:n_cube_frames-1] = cube
		;cube2[*,*,n_cube_frames:n_cube_frames+n_not_cube_frames-1] = not_cube
		;if ~((size(cube2))[3] eq ((size(cube))[3] + (size(not_cube))[3])) then stop, (size(cube2))[3], (size(cube))[3] + (size(not_cube))[3]
		cube2 = cube
		
		;remove saturated center
		;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.
	
		if annular_clean eq 1 then begin
			for ii = 0,(size(cube2))[3]-1 do begin
				; print periodic progress tracking
				if ii mod 100 eq 0 then begin
					print, 'Annular cropping on frame ', ii, '/', (size(cube2))[3]-1
				endif
				
				med_ii = median(cube2[*,*,ii], /double, /even)
				for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) le 10. then cube2[xx,yy,ii] = med_ii
				for xx=0,sz-1 do for yy=0,sz-1 do if sqrt( (( xx - half_cropped_sz)^2.) + (( yy - half_cropped_sz)^2.) ) ge 40. then cube2[xx,yy,ii] = med_ii
			endfor
		endif
		
		; this is all inside of the IWA, more or less.
		if centered_clean eq 1 then begin
			xx = (findgen(sz)# replicate(1, sz)) - half_cropped_sz
			yy = (replicate(1, sz)# findgen(sz)) - half_cropped_sz
			mask = sqrt(xx^2 + yy^2) ge 40.
			
			; Apply the mask to each frame
			for ii = 0,(size(cube2))[3]-1 do begin
				; print periodic progress tracking
				if ii mod 100 eq 0 then begin
					print, 'Center-cropping on frame ', ii, '/', (size(cube2))[3]-1
				endif
				
				med_ii = median(cube2[*,*,ii], /double, /even)
				indices = where(mask eq 1)
				if indices[0] ne -1 then begin
					frame = cube2[0:sz-1, 0:sz-1, ii]
					frame[indices] = med_ii
					cube2[0:sz-1, 0:sz-1, ii] = frame
				endif
			endfor
		endif
	
		; Compute pupil median
		lmed = median(cube2, dim=3, /double, /even)
		
		; Create mask for peak_cropped once outside the loop
		peak_cropped_sz = 2*fix(fwhm) + 1
		half_peak_cropped_sz = peak_cropped_sz / 2.
		xx = (findgen(peak_cropped_sz)# replicate(1, peak_cropped_sz)) - half_peak_cropped_sz
		yy = (replicate(1, peak_cropped_sz)# findgen(peak_cropped_sz)) - half_peak_cropped_sz
		peak_mask = sqrt(xx^2 + yy^2) ge fix(fwhm)
		
		; Pre-allocate arrays
		num_frames = n_cube_frames
		corrs = fltarr(num_frames)
		peaks = fltarr(num_frames)
		frame_stddevs = fltarr(num_frames)
		
		; Extract the cropped region all at once
		x_range = [half_cropped_sz-fix(fwhm), half_cropped_sz+fix(fwhm)]
		y_range = [half_cropped_sz-fix(fwhm), half_cropped_sz+fix(fwhm)]
		peak_cropped = cube2[x_range[0]:x_range[1], y_range[0]:y_range[1], 0:n_cube_frames-1]
		
		; Apply circular mask to peak_cropped all at once
		for i=0, num_frames-1 do begin
			mask_indices = where(peak_mask eq 1)
			if mask_indices[0] ne -1 then begin
				temp = peak_cropped[*, *, i]
				temp[mask_indices] = !values.f_nan
				peak_cropped[*, *, i] = temp
			endif
		endfor
		
		; Process all frames with minimal repeated code
		for jj=0, num_frames-1 do begin
			frame_jj = cube2[*,*,jj]
			frame_jj_peak = peak_cropped[*,*,jj]
			
			corrs[jj] = max(crosscorr(frame_jj, lmed, /cut))
			
			; Calculate statistics
			peaks[jj] = mean(frame_jj_peak, /double, /nan)
			frame_stddevs[jj] = stddev(frame_jj_peak, /double, /nan)
			
			; Print progress and stats
			if jj mod 100 eq 0 then begin
				print, 'Frame ', jj, ' has max correlation value of ', corrs[jj], ' with pupil median.'
				if jj gt 0 then begin
					if peak_clean eq 1 then print, 'Peak flux / median peak flux = ', peaks[jj]/median(peaks[0:jj])
					if stddev_clean eq 1 then print, 'Frame stddev', frame_stddevs[jj]
				endif
			endif
		endfor
		
		delvar, cube2; no longer needed.
		
		; select by ratio of peak flux to median peak flux
		if peak_clean eq 1 then begin
			total_bad_peaks = 0
			; below threshold
			peaks_too_low = where(peaks/median(peaks, /double, /even) lt peak_thresh, /null)
			total_bad_peaks = temporary(total_bad_peaks) + (size(peaks_too_low))[1]
			if n_elements(bad_pixels) eq 0 then bad_pixels = peaks_too_low else bad_pixels=[bad_pixels, peaks_too_low]
			
			; above threshold
			;peaks_too_high = where(peaks/median(peaks, /double, /even) gt 2.-peak_thresh, /null)
			;bad_pixels=[bad_pixels, peaks_too_high]
			;total_bad_peaks += (size(peaks_too_high))[1]
		endif
		;print, 'bad_pixels line 395: ', bad_pixels
		
		; select by ratio of frade stddev to mean frame stddev
		if stddev_clean eq 1 then begin
		
			; get an outlier-resistant mean and standard deviation of the standard deviations
			resistant_mean, frame_stddevs, 3.5, mean_std, std_std, /double, goodvec=goodvec
			std_std *= sqrt(n_elements(goodvec) - 1)
			;mean_std = mean(frame_stddevs, /double, /nan)
			;std_std = stddev(frame_stddevs, /double, /nan)
			print, newline, 'mean, std of frame_stddevs: ', mean_std, std_std
			
			;stddev above threshold
			;bad_stddev_is = where(frame_stddevs/median(frame_stddevs, /double, /even) gt stddev_thresh, /null)
			bad_stddev_is = where((frame_stddevs - mean_std) gt (stddev_thresh*std_std)); trim 3-sigma outliers in pixel-to-pixel standard deviation of peak
			;print, 'bad_stddev_is:', bad_stddev_is
			; none found.
			if bad_stddev_is eq [-1] then bad_stddev_is = !NULL
			
			total_bad_stddevs = (size(bad_stddev_is))[1]
			if n_elements(bad_pixels) eq 0 then bad_pixels = bad_stddev_is else bad_pixels=[bad_pixels, bad_stddev_is]
			bad_pixels = bad_pixels[rem_dup(bad_pixels)]
		endif
		;print, 'bad_pixels line 412: ', bad_pixels
		
		; reject a certain minimum number of frames
		if keyword_set(keep_number) then begin
		
			; See how many frames we should reject in total
			total_number = (size(cube))[3]
			reject_number = total_number - keep_number
			
			current_reject_number = 0
			; remove duplicates and count current # of rejects
			if n_elements(bad_pixels) gt 0 then bad_pixels = bad_pixels[rem_dup(bad_pixels)]
			rejects_pre_corr = n_elements(bad_pixels)
			
			current_reject_number = rejects_pre_corr
			print, current_reject_number, ' rejected before cross-corr comparison'
			
			; use cross-correlation for the rest.
			; initialize based on current difference in rejection from
			; target.
			left_to_reject = reject_number - current_reject_number
			while current_reject_number lt reject_number do begin
				low_corr_is = N_LOWEST_INDICES(corrs, left_to_reject)
				
				; add lowest crosscorss to bad_pixels and remove duplicates
				if n_elements(bad_pixels) eq 0 then bad_pixels = low_corr_is else bad_pixels=[bad_pixels, low_corr_is]
				;print, 'bad_pixels line 438: ', bad_pixels
				bad_pixels = bad_pixels[rem_dup(bad_pixels)]
				
				; update with non-redundant # of rejects.
				current_reject_number = n_elements(bad_pixels)
				
				; increase number to reject based on current difference from
				; the target.
				left_to_reject += reject_number - current_reject_number
			endwhile
			
			total_bad_corrs = current_reject_number - rejects_pre_corr
		endif
		
		; double-check on rejects
		;print, 'bad_pixels line 453: ', bad_pixels
		if n_elements(bad_pixels) gt 0 then bad_pixels = bad_pixels[rem_dup(bad_pixels)]
	
		print, newline, 'Found ', n_elements(bad_pixels), ' bad frames.'
		print, 'Average correlation = ', mean(corrs, /double, /nan)
		print, 'Median correlation = ', median(corrs, /double, /even)
		print, 'Standard deviation = ', stddev(corrs, /double, /nan)
		print, newline, total_bad_corrs, 'frames rejected for crosscor'
		if peak_clean eq 1 then print, total_bad_peaks, 'frames rejected for peak flux'
		if stddev_clean eq 1 then print, total_bad_stddevs, 'frames rejected for stddev', newline
		;hak
		print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'
	
		if n_elements(bad_pixels) gt 1 then begin
			;print, 'bad_pixels: ', bad_pixels
			remove, bad_pixels, goods, angles
			right_bad_pixels_cube = cube[*,*,bad_pixels]
			;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', right_bad_pixels_cube
		endif else if n_elements(bad_pixels) eq 1 and bad_pixels ne -1 then begin
			remove, bad_pixels, goods, angles
			right_bad_pixels_cube = cube[*,*,bad_pixels]
			;writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + '_' + string(keep_number) +  '_corrthresh.fits', right_bad_pixels_cube
		endif
	
		cube = (temporary(cube))[*,*,goods]; Reassign our cube to only have our good pixels
		
		; third shift - now that the images are reasonably co-registered and common-centered
		; to the frame center, we can use a new 'med' to further co-register
		; remove NaNs (again)
		for jj=0,(size(cube))[3]-1 do begin
			frame_jj = cube[*,*, jj]
			frame_jj[~finite(frame_jj)] = median(frame_jj, /double, /even)
			cube[*,*,jj] = frame_jj
		endfor
		;initialize a new cube to modify with the same data
		obj_cube2 = cube
		print, 'cube copy created! Blocking out what needs blocked...'
	   
		; set to zero if we're blocking a side out
		if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *, *] /= 1000000000000000000.
		if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] /= 1000000000000000000.
		if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] /= 1000000000000000000.
		if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *] /= 1000000000000000000.
		print, 'Requested blocking complete! Changing some stuff in the second cube...'
	
		obj_cube2[where(obj_cube2 lt 0)] /= 1000000000000000000. ;reduce negative values so CC doesn't get confused on subtracted PSF.
	
		n_frames = (size(obj_cube2))[3]
		; Pre-compute constants
		half_sizex = half_cropped_sz
		half_sizey = half_cropped_sz
		
		; Pre-allocate shift arrays
		dx_array = fltarr(n_frames)
		dy_array = fltarr(n_frames)
		
		; ---- FIRST OPTIMIZATION: More efficient filtering of frames ----
		if do_cen_filter then begin
			print, 'starting cen_filter for third shift...'
			
			; Since filter_image only accepts 2D arrays, we need to process frames individually
			; But we can still optimize by reducing temporary array creation
			for ii=0, n_frames - 1 do begin
				; periodic progress tracking
				if ii mod 100 eq 0 then begin
					print, 'Cen_filter on frame:' + string(ii)
					print, 'Run:' + string(runs)
				endif
				
				; Extract frame only once to avoid repeated memory allocation
				current_frame = obj_cube2[*,*,ii]
				
				; Filter in steps to avoid creating intermediate temporary arrays
				filtered_frame = filter_image(current_frame, smooth=new_filter, /ALL_PIXELS)
				current_frame -= filtered_frame
				
				; Apply second filter
				filtered_frame = filter_image(current_frame, smooth=3., /ALL_PIXELS)
				
				; Store back to cube
				obj_cube2[*,*,ii] = filtered_frame
			endfor
		endif; do_cen_filter if
		
	;	; ---- SECOND OPTIMIZATION: Calculate reference frame once ----
	;	ref_frame = median(obj_cube2, dim=3, /even, /double)
	;	
	;	; ---- THIRD OPTIMIZATION: Calculate all shifts first, then apply them ----
	;	; Calculate all shifts first
	;	for ii=0, n_frames-1 do begin
	;		frame_ii = obj_cube2[*,*,ii]
	;		
	;		corr = crosscorr(ref_frame, frame_ii, pmax, /cut)
	;		if pmax[0] ge 0. then dx = half_sizex-pmax[0] else dx = -half_sizex+abs(pmax[0])
	;		if pmax[1] ge 0. then dy = half_sizey-pmax[1] else dy = -half_sizey+abs(pmax[1])
	;		
	;		dx_array[ii] = -dx
	;		dy_array[ii] = -dy
	;		
	;		if ii mod 100 eq 0 then begin
	;			print, 'Calculated shift for frame', ii, ' / ', n_frames - 1, ' by: ', String(dx_array[ii]), String(dy_array[ii])
	;		endif
	;	endfor
		
		; Apply all shifts
		;for ii=0, n_frames-1 do begin
		;	cube[*,*,ii] = shift_sub(temporary(cube[*,*,ii]), dx_array[ii], dy_array[ii], cubic=-1.0)
		;	
		;	; Apply smoothing if needed
		;	if do_smooth eq 1 then begin
		;		frame_ii = cube[*,*,ii]
		;		nan_indices = where(finite(frame_ii) ne 1, nan_count)
		;		if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
		;		cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		;	endif
		;	
		;	if ii mod 100 eq 0 then begin
		;		print, 'Applied third shift to frame', ii, ' / ', n_frames - 1
		;		print, 'Run: ' + string(runs)
		;	endif
		;endfor
		;print, 'Third frame shifting complete! Starting the fourth shift...'
		
		; ---- FOURTH OPTIMIZATION: Calculate final common shift once ----
	;	med = median(cube, dim=3, /even, /double)
	;	peak = mpfit2dpeak(med, A, /tilt)
	;	xx = A[4] & yy = A[5]
	;	common_dx = half_cropped_sz-xx
	;	common_dy = half_cropped_sz-yy
		
	;	print, 'Common shift for all frames: ', string(common_dx), string(common_dy)
		
		; Apply the common shift to all frames
	;	for ii=0, n_frames-1 do begin
	;		cube[*,*,ii] = shift_sub(temporary(cube[*,*,ii]), common_dx, common_dy, cubic=-1.0)
			
			; Apply smoothing if needed
	;		if do_smooth eq 1 then begin
	;			frame_ii = cube[*,*,ii]
	;			nan_indices = where(finite(frame_ii) ne 1, nan_count)
	;			if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
	;			cube[*,*,ii] = convolve(frame_ii, lp_PSF)
	;		endif
			
	;		if ii mod 100 eq 10 then begin
	;			print, 'Fourth shift applied to frame: ' + string(ii)
	;			print, 'Run: ' + string(runs)
	;		endif
	;	endfor
	;	print, 'Fourth shifting complete!'
		   
		medarr, cube, medframe; Get a median frame
		
		;Write our clean cube, our bad pixels, and our median frame
		writefits, strcompress(output_folder + dither_folder + obj_name + '_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
			'_stddev_thresh_'+string(stddev_thresh) + '_cube_skysub_cen_clean.fits', /r), cube, hdr
			
		writefits, strcompress(output_folder + dither_folder + obj_name + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
			'_stddev_thresh_'+string(stddev_thresh)+  '_cube_skysub_cen_bad_pixels.fits', /r), right_bad_pixels_cube, hdr
			
		writefits, strcompress(output_folder + dither_folder + obj_name + '_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
			'_stddev_thresh_'+string(stddev_thresh) + '_pupil.fits', /r), medframe
			
		save, filename = strcompress(output_folder + dither_folder + obj_name + '_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
			'_stddev_thresh_'+string(stddev_thresh) +  '_parang_clean.sav', /r), angles
	
		;cgcleanup; Found online, closes any graphics or widget windows from our plots
		print, 'Run: ' + string(runs) + ' complete
	endfor; DX only
	
endif; second stripe if

print, 'Finished!'

end

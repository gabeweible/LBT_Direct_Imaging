FUNCTION bin_refcube, refcube, bin, ref_ang, filter, lp_filter, output_path,$
	obj, short_suffix, processed_ref_angles, runs=runs, bin_type=bin_type,$
	half_cropped_sz=half_cropped_sz, normal=normal, fwhm=fwhm
	
compile_opt IDL2, logical_predicate
newline = string(10B)

if ~keyword_set(bin_type) then bin_type='median'

sz = 2*half_cropped_sz
xs = sz & ys = sz
xhs = xs / 2. & yhs = ys / 2.

if not keyword_set(runs) then begin

	print, 'Binning references...'
	; Optimized binning using array reshaping
	n_frames = (size(refcube, /dim))[2]
	n_bins = floor(n_frames / bin)
	
	; Reshape and reduce reference cube
	reshaped_cube = reform(refcube[*,*,0:n_bins*bin-1], $
						   (size(refcube, /dim))[0], $
						   (size(refcube, /dim))[1], $
						   bin, n_bins)
	
	; Choose binning method
	if bin_type eq 'median' then $
		processed_refcube = median(reshaped_cube, dimension=3, /even, /double) $
	else $
		processed_refcube = mean(reshaped_cube, dimension=3, /nan, /double)
	
	; Average angles for binned frames
	processed_ref_angles = total(reform(ref_ang[0:n_bins*bin-1], bin, n_bins), 1) / bin
	
	; Normalize all of the frames so that they each have a max value of 1 (should help with residuals around the star)
	if keyword_set(normal) then begin
		print, 'Normalizing frames...'
	   if normal eq 1 then for iv=0, n_frames - 1 do begin; only normalize to the smoothed peak		
			; try a more robust normalization?
			
			frame_iv = processed_refcube[*,*,iv]
			
			; grab pixel values within a radius of fwhm/2.
			center_pixels = []
			for dx=-fwhm/2.,fwhm/2. do begin
			for dy=-fwhm/2.,fwhm/2. do begin
				if sqrt( dx^2. + dy^2. ) lt fwhm/2. then center_pixels = [[center_pixels],$
					frame_iv[xhs+dx, yhs+dy]]
			endfor
			endfor
			
			processed_refcube[*,*,iv] /= (median(center_pixels, /even, /double) / 100.)
				
	   endfor; double-check on normalization
	endif
	
	; by turning bad pixels to 0, then filtering, I effectively replace them with
	; nearby values.
	
	; High-pass filtering AND/OR low-pass filtering
	if filter gt 0 and lp_filter gt 0 then begin
		
		lp_PSF = psf_Gaussian(NPIX=xs, FWHM=[lp_filter,lp_filter], /normalize)
		hp_PSF = psf_Gaussian(NPIX=xs, FWHM=[filter, filter], /normalize)
		
		print, 'Filtering science frames...'
		
		for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
			; Replace NaNs with zeros
			nan_indices = where(finite(frame_iv) ne 1, nan_count)
			if nan_count gt 0 then frame_iv[nan_indices] = 0.
		
			frame_iii = processed_refcube[*,*,iii]
			; high-pass
			processed_refcube[*,*,iii] -= convolve(frame_iii, hp_PSF)
			
			frame_iii = processed_refcube[*,*,iii]
			; low-pass
			processed_refcube[*,*,iii] = convolve(frame_iii, lp_PSF)
		endfor
		
	endif else begin; only one filtering or the other.
	
		if filter gt 0 then begin; just high-pass filtering
		
			hp_PSF = psf_Gaussian(NPIX=xs, FWHM=[filter, filter], /normalize)
			print, 'Filtering science frames...'
			
			for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
				frame_iii = processed_refcube[*,*,iii]
				
				; Replace NaNs with zeros
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = 0.
				
				; high-pass
				processed_refcube[*,*,iii] -= convolve(frame_iii, hp_PSF)
			endfor
			
		endif
		
		if lp_filter gt 0 then begin; just low-pass filtering
		
			lp_PSF = psf_Gaussian(NPIX=xs, FWHM=[lp_filter,lp_filter], /normalize)
			print, 'Filtering science frames...'
			
			for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
				frame_iii = processed_refcube[*,*,iii]
				
				; Replace NaNs with zeros
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = 0.
				
				; low-pass
				processed_refcube[*,*,iii] = convolve(frame_iii, lp_PSF)
			endfor
			
		endif
		
	endelse
	
	; Save preprocessed reference cube
	print, 'deleting old refs and saving new refs...'
	refs_filename = output_path+obj+short_suffix+'_processed_refs.sav'
	file_delete, refs_filename, /allow_nonexistent
	save, processed_refcube, processed_ref_angles, $
		  filename=refs_filename
endif else begin

	print, 'Binning references...'
	; Optimized binning using array reshaping
	n_frames = (size(refcube, /dim))[2]
	n_bins = floor(n_frames / bin)
	
	; Reshape and reduce reference cube
	reshaped_cube = reform(refcube[*,*,0:n_bins*bin-1], $
						   (size(refcube, /dim))[0], $
						   (size(refcube, /dim))[1], $
						   bin, n_bins)
	
	; Choose binning method
	if bin_type eq 'median' then $
		processed_refcube = median(reshaped_cube, dimension=3, /even, /double) $
	else $
		processed_refcube = mean(reshaped_cube, dimension=3, /nan, /double)
	
	; Average angles for binned frames
	processed_ref_angles = total(reform(ref_ang[0:n_bins*bin-1], bin, n_bins), 1) / bin
	
	; Normalize all of the frames so that they each have a max value of 1 (should help with residuals around the star)
	if keyword_set(normal) then begin
		print, 'Normalizing frames...'
	   if normal eq 1 then for iv=0, n_frames - 1 do begin; only normalize to the smoothed peak		
			; try a more robust normalization?
			
			frame_iv = processed_refcube[*,*,iv]
			
			; grab pixel values within a radius of fwhm/2.
			center_pixels = []
			for dx=-fwhm/2.,fwhm/2. do begin
			for dy=-fwhm/2.,fwhm/2. do begin
				if sqrt( dx^2. + dy^2. ) lt fwhm/2. then center_pixels = [[center_pixels],$
					frame_iv[xhs+dx, yhs+dy]]
			endfor
			endfor
			
			processed_refcube[*,*,iv] /= (median(center_pixels, /even, /double) / 100.)
				
	   endfor; double-check on normalization
	endif
	
	; by turning bad pixels to 0, then filtering, I effectively replace them with
	; nearby values.
	
	; High-pass filtering AND/OR low-pass filtering
	if filter gt 0 and lp_filter gt 0 then begin
		
		lp_PSF = psf_Gaussian(NPIX=xs, FWHM=[lp_filter,lp_filter], /normalize)
		hp_PSF = psf_Gaussian(NPIX=xs, FWHM=[filter, filter], /normalize)
		
		print, 'Filtering science frames...'
		
		for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
			; Replace NaNs with zeros
			nan_indices = where(finite(frame_iv) ne 1, nan_count)
			if nan_count gt 0 then frame_iv[nan_indices] = 0.
		
			frame_iii = processed_refcube[*,*,iii]
			; high-pass
			processed_refcube[*,*,iii] -= convolve(frame_iii, hp_PSF)
			
			frame_iii = processed_refcube[*,*,iii]
			; low-pass
			processed_refcube[*,*,iii] = convolve(frame_iii, lp_PSF)
		endfor
		
	endif else begin; only one filtering or the other.
	
		if filter gt 0 then begin; just high-pass filtering
		
			hp_PSF = psf_Gaussian(NPIX=xs, FWHM=[filter, filter], /normalize)
			print, 'Filtering science frames...'
			
			for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
				frame_iii = processed_refcube[*,*,iii]
				
				; Replace NaNs with zeros
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = 0.
				
				; high-pass
				processed_refcube[*,*,iii] -= convolve(frame_iii, hp_PSF)
			endfor
			
		endif
		
		if lp_filter gt 0 then begin; just low-pass filtering
		
			lp_PSF = psf_Gaussian(NPIX=xs, FWHM=[lp_filter,lp_filter], /normalize)
			print, 'Filtering science frames...'
			
			for iii = 0, (size(processed_refcube, /dim))[2]-1 do begin
				frame_iii = processed_refcube[*,*,iii]
				
				; Replace NaNs with zeros
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = 0.
				
				; low-pass
				processed_refcube[*,*,iii] = convolve(frame_iii, lp_PSF)
			endfor
			
		endif
		
	endelse
	
	; Save preprocessed reference cube
	print, 'deleting old refs and saving new refs...'
	refs_filename = strcompress(output_path+obj+short_suffix+$
		'_processed_refs_runs_'+string(runs)+'.sav', /r)
	file_delete, refs_filename, /allow_nonexistent
	save, processed_refcube, processed_ref_angles, $
		  filename=refs_filename
endelse

return, processed_refcube
	  
end
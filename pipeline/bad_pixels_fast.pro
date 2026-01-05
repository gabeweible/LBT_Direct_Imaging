pro bad_pixels_fast, output_folder, obj_name, stripe,$; sigma_threshold=sigma_threshold, $
    boxhsize=boxhsize, pre_sky_sub=pre_sky_sub, bad_px_arr=bad_px_arr, $; n_groups=n_groups,$
	type=type, run=run, do_destripe=do_destripe,$
    fwhm=fwhm, just_destripe=just_destripe, do_second_round=do_second_round,$
    do_first_round=do_first_round, just_second_round=just_second_round, debug=debug,$
    create_master_masks=create_master_masks,$
    create_badpix_mask=create_badpix_mask, use_flat_corrected=use_flat_corrected,$
    darks_filename=darks_filename, coadd=coadd, vapp=vapp, destripe_skysub=destripe_skysub,$
    pca_skysub=pca_skysub
    
    compile_opt idl2, logical_predicate
    newline = string(10B)
    
   ; if not keyword_set(pca_skysub) then pca_skysub=0
    
    if vapp eq 1 then begin
    
    	boxehs = 60.; host star exclusion radius
    	
		g_boxehs = boxehs/5.; ghost exclusion radius
		
		; OWA bright ring masking radii (the whole big chunk...)
		; make a little thinner for better stats?
		owa_boxehs_in = 216;175+10.
		owa_boxehs_out = 277;278-10.
		
		; 2nd ring (distant - faint)
		;owa_boxehs_in2 = 424
		;owa_boxehs_out2 = 450
	endif
    
    if obj_name eq 'Alcor' then begin
		; mask hot column segments that couldn't easily be fixed
		hc1_x1 = 385.-1.
		hc1_y1 = 250.-1.
		hc1_y2 = 278.-1.
		; start and end (coadded) indices to mask for this hot column
		hc1_f1 = 0 * (20/coadd)
		hc1_f2 = 1545 * (20/coadd)
		
		; mask hot column segments that couldn't easily be fixed
		hc2_x1 = 1281.-1.
		hc2_y1 = 210.-1.
		hc2_y2 = 235.-1.
		; start and end (coadded) indices to mask for this hot column
		hc2_f1 = 1550 * (20/coadd)
		hc2_f2 = 2849 * (20/coadd)
		
		; these have been multiplied x 1000 to match counts/s instead of
		; counter/ms
		; NEW: updated based on the cubes pre-destriping.
		frame_max = 6550
		frame_min = 2000
		
		; min and max values after destriping for additional masking.
		; Peak values can vary quite a lot, so frame_max_post_destripe
		; has had to be increased a couple of times... hot pixels
		; are not really the concern, however. There are many more which are dead.
		; NEW: UPDATED BASED ON POST-DESTRIPE RESULTS
		frame_max_post_destripe = 4000
		frame_min_post_destripe = -130
		
		hot_sigma = 2.0
		flat_thresh = [0.48, 1.15]
		
		; Configuration parameters
		max_sources = 2 ; Maximum number of sources to search for (pos and negative)
		min_separation = 100.0  ; px
		min_peak_significance = 3.0  ;Minimum significance above noise for a valid peak
		
		max_cluster_size = 1000
	endif
	
	; start by trying the same values as for Alcor,
	; except a smaller PSF exclusion size. (half the size as before)
	; since there is no super bright IWA w/the vAPP - just a first Airy
	; ring that looks pretty normal.
	if obj_name eq 'tyc5709' then begin
	
		boxehs = 80.; host star exclusion radius. Should exclude the disk, as well
		; as the optical ghost, ideally.
		
		hot_sigma = 2.0
		flat_thresh = [0.73, 1.14]
		
		; these have been multiplied x 1000 to match counts/s instead of
		; counter/ms
		; UPDATE AFTER CHECKING ALL NEW CUBES...
		frame_max = 1400
		; shrink a bit to accomodate more frames not like the example that
		; set 0.88 as a good minimum.
		frame_min = 665
		
		; min and max values after destriping for additional masking.
		; Peak values can vary quite a lot, so frame_max_post_destripe
		; has had to be increased a couple of times... hot pixels
		; are not really the concern, however. There are many more which are dead.
		frame_max_post_destripe = 575
		frame_min_post_destripe = -8.25
		
		; Configuration parameters
		max_sources = 4  ; Maximum number of sources to search for (pos and negative)
		min_separation = 100.0  ; px
		min_peak_significance = 3.0  ;Minimum significance above noise for a valid peak
		
		max_cluster_size = 1000
	endif
    
    new_filter = 27.
    nantr = 10.0; max stdev for final bad-pix ident post-normal correction - only for destriping purposes

    ; Set default parameters more efficiently
    pre_sky_sub = keyword_set(pre_sky_sub) ? pre_sky_sub : 0
    ;sigma_threshold = keyword_set(sigma_threshold) ? sigma_threshold : 5.0
    boxhsize = keyword_set(boxhsize) ? boxhsize : 1; 2 => 5x5 filter
    
    if not keyword_set(type) then type = 'mean'
    ;if not keyword_set(darks_filename) then darks_filename = output_folder + obj_name + '_darks.fits'
    
    if just_destripe ne 1 then begin
    
    	if do_first_round eq 1 then begin
    
    	print, 'Reading in cubes for bad-px correction'
		
		; Find all NOD cube files
		print, output_folder + obj_name + '_NOD_?_nod??_grp00_cube.fits'
		nod_files = file_search(output_folder + obj_name + '_NOD_?_nod??_grp00_cube.fits', count=n_nod_files)
		if n_nod_files eq 0 then begin
			print, 'ERROR: No NOD cube files found!'
			return
		endif
		
		print, 'Found ', n_nod_files, ' NOD cube files'
		
		; Sort files by nod number to ensure proper chronological order
		nod_numbers = intarr(n_nod_files)
		for i = 0, n_nod_files-1 do begin
		    nod_numbers[i] = extract_nod_number(nod_files[i])
		endfor
		
		; Sort the files by nod number
		sort_indices = sort(nod_numbers)
		nod_files = nod_files[sort_indices]
		
		print, 'Files sorted by nod number:'
		for i = 0, n_nod_files-1 do begin
		    print, '  ', file_basename(nod_files[i])
		endfor
		
		; Separate NOD_A and NOD_B files
		nod_a_files = []
		nod_b_files = []
		
		for i = 0, n_nod_files-1 do begin
			if strpos(nod_files[i], '_NOD_A_') ne -1 then begin
				nod_a_files = [nod_a_files, nod_files[i]]
			endif else if strpos(nod_files[i], '_NOD_B_') ne -1 then begin
				nod_b_files = [nod_b_files, nod_files[i]]
			endif
		endfor
		
		n_nod_a = n_elements(nod_a_files)
		n_nod_b = n_elements(nod_b_files)
		
		print, 'Found ', n_nod_a, ' NOD_A files and ', n_nod_b, ' NOD_B files'
		
		if (n_nod_a eq 0) or (n_nod_b eq 0) then begin
		  print, 'ERROR: Missing one or both nodding positions'
		  return
		endif
		
		; Read first NOD file to get dimensions
		first_cube = readfits_fast(nod_a_files[0])
		cube_size = size(first_cube, /dimensions)
		x_dim = cube_size[0]
		y_dim = cube_size[1]
		frames_per_nod = cube_size[2]
		delvar, first_cube  ; Free memory immediately
		
		print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'
		
		; Replaced the older flat-field creation section with this outlier-resistant version
		if create_master_masks eq 1 then begin
		  print, 'Generating master dark...'
		  darks = readfits_fast(darks_filename)
		  resistant_mean, darks, 3.5, master_dark, dim=3, /double
		  delvar, darks  ; Free memory
		
		  print, 'Saving master dark to FITS file...'
		  writefits, output_folder + obj_name + '_master_dark.fits', master_dark
		  print, 'Master dark FITS saved.'
		
		  print, 'Creating outlier-resistant master flat fields for each nod position...'
		  
		  ; --- Create NOD_A master flat (outlier-resistant) ---
		  print, 'Creating outlier-resistant master flat for NOD_A...'
		  
		  ; First, create individual nod flats for each NOD_A file
		  nod_a_individual_flats = []
		  
		  for i = 0, n_nod_a-1 do begin
			print, 'Processing NOD_A file ', i+1, '/', n_nod_a, ' for individual flat'
			nod_cube = readfits_fast(nod_a_files[i])
			
			; Subtract master dark from each frame
			for kk = 0, frames_per_nod-1 do begin
			  nod_cube[*,*,kk] -= master_dark
			endfor
			
			; Create array to hold normalized frames for this nod
			normalized_frames = dblarr(x_dim, y_dim, frames_per_nod)
			valid_frame_count = 0
			
			; Normalize each frame by its median
			for kk = 0, frames_per_nod-1 do begin
			  frame = nod_cube[*,*,kk]
			  frame_median = median(frame, /even, /double)
			  ; don't normalize if median is < 0 (would flip PSF sign)
			  if frame_median gt 0 then begin
				normalized_frames[*,*,valid_frame_count] = frame / frame_median
				valid_frame_count += 1
			  endif
			endfor
			
			; Trim the array to only valid frames
			if valid_frame_count gt 0 then begin
			  normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
			  
			  ; Create outlier-resistant mean for this individual nod
			  individual_nod_flat = dblarr(x_dim, y_dim)
			  resistant_mean, normalized_frames, 3.5, individual_nod_flat, dim=3, /double
			  
			  ; Normalize this individual flat to its median
			  individual_median = median(individual_nod_flat, /even, /double)
			  ; don't normalize if median is < 0 (would flip PSF sign)
			  if individual_median gt 0 then individual_nod_flat /= individual_median
			  
			  ; Add to collection
			  if n_elements(nod_a_individual_flats) eq 0 then begin
				nod_a_individual_flats = dblarr(x_dim, y_dim, n_nod_a)
			  endif
			  nod_a_individual_flats[*,*,i] = individual_nod_flat
			endif
			
			delvar, nod_cube, normalized_frames, individual_nod_flat  ; Free memory
		  endfor
		  
		  ; Now create final NOD_A master flat from individual nod flats
		  print, 'Combining individual NOD_A flats with outlier resistance...'
		  nod_a_flat = dblarr(x_dim, y_dim)
		  resistant_mean, nod_a_individual_flats, 3.5, nod_a_flat, dim=3, /double
		  
		  ; Final normalization to median
		  nod_a_median = median(nod_a_flat, /even, /double)
		  ; don't normalize if median is < 0 (would flip PSF sign)
		  if nod_a_median gt 0 then nod_a_flat /= nod_a_median
		  
		  delvar, nod_a_individual_flats  ; Free memory
		  
		  ; --- Create NOD_B master flat (outlier-resistant) ---
		  print, 'Creating outlier-resistant master flat for NOD_B...'
		  
		  ; First, create individual nod flats for each NOD_B file
		  nod_b_individual_flats = []
		  
		  for i = 0, n_nod_b-1 do begin
			print, 'Processing NOD_B file ', i+1, '/', n_nod_b, ' for individual flat'
			nod_cube = readfits_fast(nod_b_files[i])
			
			; Subtract master dark from each frame
			for kk = 0, frames_per_nod-1 do begin
			  nod_cube[*,*,kk] -= master_dark
			endfor
			
			; Create array to hold normalized frames for this nod
			normalized_frames = dblarr(x_dim, y_dim, frames_per_nod)
			valid_frame_count = 0
			
			; Normalize each frame by its median
			for kk = 0, frames_per_nod-1 do begin
			  frame = nod_cube[*,*,kk]
			  frame_median = median(frame, /even, /double)
			  ; don't normalize if median is < 0 (would flip PSF sign)
			  if frame_median gt 0 then begin
				normalized_frames[*,*,valid_frame_count] = frame / frame_median
				valid_frame_count += 1
			  endif
			endfor
			
			; Trim the array to only valid frames
			if valid_frame_count gt 0 then begin
			  normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
			  
			  ; Create outlier-resistant mean for this individual nod
			  individual_nod_flat = dblarr(x_dim, y_dim)
			  resistant_mean, normalized_frames, 3.5, individual_nod_flat, dim=3, /double
			  
			  ; Normalize this individual flat to its median
			  individual_median = median(individual_nod_flat, /even, /double)
			  ; don't normalize if median is < 0 (would flip PSF sign)
			  if individual_median gt 0 then individual_nod_flat /= individual_median
			  
			  ; Add to collection
			  if n_elements(nod_b_individual_flats) eq 0 then begin
				nod_b_individual_flats = dblarr(x_dim, y_dim, n_nod_b)
			  endif
			  nod_b_individual_flats[*,*,i] = individual_nod_flat
			endif
			
			delvar, nod_cube, normalized_frames, individual_nod_flat  ; Free memory
		  endfor
		  
		  ; Now create final NOD_B master flat from individual nod flats
		  print, 'Combining individual NOD_B flats with outlier resistance...'
		  nod_b_flat = dblarr(x_dim, y_dim)
		  resistant_mean, nod_b_individual_flats, 3.5, nod_b_flat, dim=3, /double
		  
		  ; Final normalization to median
		  nod_b_median = median(nod_b_flat, /even, /double)
		  ; don't normalize if median is < 0 (would flip PSF sign)
		  if nod_b_median gt 0 then nod_b_flat /= nod_b_median
		  
		  delvar, nod_b_individual_flats  ; Free memory
		
		  ; Save the outlier-resistant master flats
		  writefits, output_folder + obj_name + '_nod_a_master_flat_resistant.fits', nod_a_flat
		  writefits, output_folder + obj_name + '_nod_b_master_flat_resistant.fits', nod_b_flat
		
		  ; Stitch the two signal-less halves together:
		  ; (left half of b, right half of a)
		  ; mid-way point between the two nods is ~pixel index 897
		  ; This is true for Alcor data - may need to be adjusted for TYC5709?
		  super_flat = nod_a_flat ; make a copy that has the right dimensions
		  if obj_name eq 'Alcor' then super_flat[0:897, *] = nod_b_flat[0:897, *]
		  if obj_name eq 'tyc5709' then super_flat[0:970, *] = nod_b_flat[0:970, *]
		  delvar, nod_a_flat, nod_b_flat; don't need these anymore...
		  writefits, output_folder + obj_name + '_super_master_flat_resistant_raw.fits', super_flat
		
		  ; Final normalization of the stitched flat by its median
		  print, 'Final normalization of stitched master flat...'
		  super_median = median(super_flat, /even, /double)
		  ; don't normalize if median is < 0 (would flip PSF sign)
		  if super_median gt 0 then super_flat /= super_median
		
		  writefits, output_folder + obj_name + '_super_master_flat_resistant_normalized.fits', super_flat
		  print, 'Outlier-resistant master flat correction complete.'
		  
		  ; Create master bad-pixel mask
		  if create_badpix_mask eq 1 then begin
			print, 'Creating master bad-pixel mask from master dark and outlier-resistant master flat'
			
			; THESE VALUES WERE SET FOR ALCOR - MAY NEED ADJUSTED FOR TYC5709
			gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask,$
					   hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
					   
			writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
		  endif
		  
		endif else begin; read-in pre-generated masks
		
		  print, 'Reading in pre-generated master files...'
		  ; Try to read the resistant version first, fall back to original if needed
		  resistant_flat_file = output_folder + obj_name + '_super_master_flat_resistant_normalized.fits'
		  original_flat_file = output_folder + obj_name + '_super_master_flat_normalized.fits'
		  
		  if file_test(resistant_flat_file) then begin
			print, 'Reading outlier-resistant master flat...'
			super_flat = readfits_fast(resistant_flat_file)
		  endif else begin
			print, 'Outlier-resistant flat not found, reading original flat...'
			super_flat = readfits_fast(original_flat_file)
		  endelse
		  
		  master_dark = readfits_fast(output_folder + obj_name + '_master_dark.fits')
		  
		  if create_badpix_mask eq 1 then begin
			print, 'Creating master bad-pixel mask from master dark and master flat'
			
			; THESE VALUES WERE SET FOR ALCOR - MAY NEED ADJUSTED FOR TYC5709
			gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask,$
							   hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
							   
			writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
		  endif else begin
			print, 'Reading-in existing bad pixel mask...'
            master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask.fits')
		  endelse
		endelse
		
		; ===========================================================================
		; MEMORY EFFICIENT PROCESSING: Process each NOD file separately
		; ===========================================================================
		print, 'Processing NOD files individually for bad pixel correction and flat fielding...'
		
		; Initialize frame counting for hot column masking
		total_frame_count = 0L
		
		; Process each NOD file separately (already sorted)
		for i = 0, n_nod_files-1 do begin
			print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
			
			; Read the current NOD cube
			nod_cube = readfits_fast(nod_files[i])
			nod_frames = (size(nod_cube, /dimensions))[2]
			
			; Apply bad pixel correction
			n_bad_pix = n_elements(where(master_badpix_mask eq 0))
			if n_bad_pix gt 0 then begin
				print, '  Applying bad pixel correction to ', nod_frames, ' frames'
				print, '  Total bad pixels: ', n_bad_pix, ' (', $
					   string((100.* n_bad_pix) / (x_dim * y_dim), format='(F5.2)'), '% of each frame)'
				
				; NEW 06/20/25: add Gaussian noise to interpolated pixels
                ; Process each frame in this NOD cube
                for kk = 0, nod_frames-1 do begin
                    if (kk mod 100) eq 0 then print, '    Processing frame ', kk+1, '/', nod_frames
                    
                    frame = nod_cube[*,*,kk]
                    
                    ; Subtract master dark
                    frame -= master_dark
                    
                    ; Apply flat field correction
                    frame /= super_flat
                    
                    if obj_name eq 'Alcor' then begin
                        ; Mask hot column segments based on global frame index
                        global_frame_idx = total_frame_count + kk
                        if (global_frame_idx ge hc1_f1) and (global_frame_idx le hc1_f2) then begin
                            master_badpix_mask[hc1_x1, hc1_y1:hc1_y2] = 0
                        endif
                        if (global_frame_idx ge hc2_f1) and (global_frame_idx le hc2_f2) then begin
                            master_badpix_mask[hc2_x1, hc2_y1:hc2_y2] = 0
                        endif
                        
                        ; Extra bad column in the middle
                        master_badpix_mask[931, *] = 0
                    endif
                    
                    ; Mask edges
                    master_badpix_mask[0:4, *] = 0; left 5 columns
                    master_badpix_mask[x_dim-4:x_dim-1, *] = 0; right 5 columns
                    master_badpix_mask[*, y_dim-3:y_dim-1] = 0; top 3 rows
                    master_badpix_mask[*, 0] = 0; bottom row
                    
                    ; I need to find that artifact in the rd2 data
                    ; which appears most prevalent in the first half
                    ; of the cleaned observations.
                    if obj_name eq 'tyc5709' then begin
                        ; bad readout channel in center
                        master_badpix_mask[896:959, *] = 0
                        
                        ; crosstalk column segment in bottom right
                        ; (enlarged by 10 px on either side...)
                        master_badpix_mask[1408, 228-10:244+10] = 0
                        
                        ; bad center rows near cold-stop artifact
                        master_badpix_mask[*, 434:536] = 0
                        
                        ; bad column to the very right
                        master_badpix_mask[1984, *] = 0
                    endif
                    
                    ; overwrite bad-pixel mask with hard-coded regions.
                    if i eq 0 then writefits, output_folder+obj_name+'_master_badpix_mask.fits',$
                        master_badpix_mask
                    
                    ; Store original frame before masking for noise estimation
                    original_frame = frame
                    
                    ; Mask master bad-pixels as NaNs
                    frame[where(master_badpix_mask eq 0)] = !values.f_nan
                    
                    ; turn pure zeros in NaN for interpolation, also.
                    frame[where(frame eq 0)] = !values.f_nan
                    
                    ; mask crazy values for interpolation:
                    ; get rid of crazy values before bad-px interpolation
                    frame[where(frame gt frame_max)] = !values.f_nan
                    frame[where(frame lt frame_min)] = !values.f_nan
                    
                    ; Store which pixels were originally NaN (before interpolation)
                    original_nan_mask = ~finite(frame)
                    
                    ; NEW:: Efficiently fix only small NaN clusters
                    efficient_nan_correction, frame, fixed_frame, $
                        max_cluster_size=max_cluster_size, npix=24, /weight, /silent
                        
                    ; get rid of crazy values after bad-px interpolation
                    fixed_frame[where(fixed_frame gt frame_max)] = !values.f_nan
                    fixed_frame[where(fixed_frame lt frame_min)] = !values.f_nan
                    
                    ; Identify pixels that were interpolated (were NaN, now finite)
                    interpolated_mask = (original_nan_mask and finite(fixed_frame))
                    
                    ; Robust noise estimation using median and MAD
                    ; Use only finite values for noise estimation
                    good_pixels = where(finite(original_frame) and master_badpix_mask eq 1, n_good)
                    if n_good gt 100 then begin  ; Need sufficient pixels for robust estimate
                        good_values = original_frame[good_pixels]
                        
                        ; Calculate robust noise estimate: 1.48 * MAD
                        median_val = median(good_values, /even)
                        mad_val = median(abs(good_values - median_val), /even)
                        noise_sigma = 1.4826 * mad_val
                        
                        ; Generate Gaussian noise for interpolated pixels
                        interpolated_indices = where(interpolated_mask, n_interpolated)
                        if n_interpolated gt 0 then begin
                            ; Generate random Gaussian noise with estimated sigma
                            gaussian_noise = randomn(seed, n_interpolated) * noise_sigma
                            
                            ; Add noise to interpolated pixels
                            fixed_frame[interpolated_indices] += gaussian_noise
                        endif
                        
                        if kk eq 0 and i eq 0 then begin
                            print, '    Noise sigma estimate: ', noise_sigma
                            print, '    Number of interpolated pixels: ', n_interpolated
                        endif
                    endif else begin
                        if kk eq 0 and i eq 0 then print, '    Warning: Insufficient good pixels for noise estimation'
                    endelse
                    
                    if debug eq 1 and kk eq 2 and i eq 0 then writefits,$
                        strcompress('~/Desktop/nod_'+string(i)+'_post_px_fix.fits', /r),$
                        fixed_frame
                    
                    nod_cube[*,*,kk] = fixed_frame
                endfor
				
				; Update total frame count
				total_frame_count += nod_frames
				
			endif else begin
				print, '  No bad pixels, applying dark/flat correction only'
				
				; Still need to apply dark and flat corrections
				for kk = 0, nod_frames-1 do begin
					frame = nod_cube[*,*,kk]
					frame -= master_dark
					frame /= super_flat
					
					; Handle potential division by zero
					bad_values = where(~finite(frame))
					if bad_values[0] ne -1 then frame[bad_values] = 0.0
					
					nod_cube[*,*,kk] = frame
				endfor
				
				total_frame_count += nod_frames
			endelse
			
			; Apply destriping if requested
			if (do_destripe eq 1) and (run eq 2) and (obj_name ne 'tyc5709') then begin
				print, '  Applying destriping to NOD cube...'
				
				; Pre-calculate arrays for masking operations
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				; Process each frame in the cube
				for kk = 0, nod_frames-1 do begin
					if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_ii_med = median(frame_ii, /double, /even)
					
					;===========================================
					; Find positive PSF peak
					;===========================================
					frame_ii_cen = frame_ii
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					mx = MAX(frame_ii_cen, location)
					ind = ARRAY_INDICES(frame_ii_cen, location)
					xx = ind[0] & yy = ind[1]
					
					;===========================================
					; Create mask using vectorized operations
					;===========================================
					rad_pos = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
					ghost_rad_pos = sqrt((x_arr-(xx+11))^2 + (y_arr-(yy+58.5))^2)
					
					if pre_sky_sub eq 0 then begin
						; Find negative PSF peak
						frame_ii_cen = -1.*frame_ii
						frame_ii_cen = frame_ii_cen > (-1.*frame_ii_med)
						frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
						
						neg_mx = MAX(frame_ii_cen, location)
						neg_ind = ARRAY_INDICES(frame_ii_cen, location)
						xx_neg = neg_ind[0] & yy_neg = neg_ind[1]
						
						rad_neg = sqrt((x_arr-xx_neg)^2 + (y_arr-yy_neg)^2)
						ghost_rad_neg = sqrt((x_arr-(xx_neg+11))^2 + (y_arr-(yy_neg+58.5))^2)
						
						mask = (rad_pos le boxehs) OR (rad_neg le boxehs) OR $
							   (ghost_rad_pos le g_boxehs) OR (ghost_rad_neg le g_boxehs) OR $
							   ((rad_pos ge owa_boxehs_in) AND (rad_pos le owa_boxehs_out)) OR $
							   ((rad_neg ge owa_boxehs_in) AND (rad_neg le owa_boxehs_out))
					endif else begin
						mask = (rad_pos le boxehs) OR $
							   (ghost_rad_pos le g_boxehs) OR $
							   ((rad_pos ge owa_boxehs_in) AND (rad_pos le owa_boxehs_out))
					endelse
					
					; Create masked frame for stripe calculation
					frame_for_stripes = frame_ii
					frame_for_stripes[where(mask)] = !values.f_nan
					
					; Additional masking based on pixel values
					frame_ii_median = median(frame_for_stripes, /double, /even)
					frame_ii_stddev = stddev(frame_for_stripes, /nan, /double)
					
					value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
								(frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; write the first mask as debugging test...
					if debug eq 1 and kk eq 2 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
						frame_for_stripes
					
					; Apply destriping (vertical then horizontal)
					destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					frame_ii = destriped_frame
					frame_for_stripes = frame_ii
					
					frame_for_stripes[where(mask)] = !values.f_nan
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                    
                    nod_cube[*,*,kk] = destriped_frame
				endfor
			endif
			
			; NEW!! We need two PSF masks for the TYC5709 data, and they are not
			; the vAPP masks.
			; Apply destriping if requested
			if (do_destripe eq 1) and obj_name eq 'tyc5709' then begin
				print, '  Applying destriping to tyc5709 NOD cube...'
				
				; Pre-calculate arrays for masking operations
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				; Process each frame in the cube
				for kk = 0, nod_frames-1 do begin
					if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_ii_med = median(frame_ii, /double, /even)
					
					;===========================================
					; Prepare frame for peak detection
					;===========================================
					frame_ii_cen = frame_ii
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					; Calculate noise statistics for peak significance testing
					noise_stats = frame_ii_cen[where(finite(frame_ii_cen))]
					noise_mean = mean(noise_stats, /double)
					noise_stddev = stddev(noise_stats, /double)
					min_peak_value = noise_mean + min_peak_significance * noise_stddev
					
					; write the first processed frame as a debugging test...
					if debug eq 1 and kk eq 2 and i eq 0 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_cen_test.fits', /r),$
						frame_ii_cen
					
					; Enhanced source detection for both positive and negative sources					
					;===========================================
					; Find all significant PSF peaks (both positive and negative)
					;===========================================
					peak_x = []  ; Array to store x coordinates of all peaks
					peak_y = []  ; Array to store y coordinates of all peaks
					peak_type = []  ; Array to store peak type: 1=positive, -1=negative
					
					; First find positive sources
					temp_frame_pos = frame_ii_cen  ; Working copy for positive peak finding
					for source_num = 0, max_sources-1 do begin
						; Find the maximum in the current temp_frame
						mx = MAX(temp_frame_pos, location, /nan)
						
						; Check if this peak is significant enough
						if ~finite(mx) or mx lt min_peak_value then break
						
						; Get the indices of the peak
						ind = ARRAY_INDICES(temp_frame_pos, location)
						xx = ind[0] & yy = ind[1]
						
						; Check minimum separation from all previously found peaks
						valid_peak = 1
						if n_elements(peak_x) gt 0 then begin
							for prev_peak = 0, n_elements(peak_x)-1 do begin
								separation = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
								if separation lt min_separation then begin
									valid_peak = 0
									break
								endif
							endfor
						endif
						
						if valid_peak then begin
							; Add this positive peak to our list
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, 1]  ; 1 for positive
							
							if debug eq 1 and kk eq 2 then begin
								print, '    Found positive peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', mx
							endif
						endif
						
						; Mask out this peak region to find the next one
						rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
						mask_peak = where(rad_temp le boxehs)
						if mask_peak[0] ne -1 then temp_frame_pos[mask_peak] = !values.f_nan
					endfor
					
					; Now find negative sources (dips below background)
					; Create inverted frame for negative peak detection
					frame_ii_neg = -1.0 * (frame_ii - frame_ii_med)  ; Invert and center around median
					frame_ii_neg = frame_ii_neg > 0  ; Keep only negative deviations (now positive)
					frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
					
					; Calculate minimum significance threshold for negative peaks
					; (same logic but for negative deviations)
					min_neg_peak_value = min_peak_significance * noise_stddev
					
					temp_frame_neg = frame_ii_neg  ; Working copy for negative peak finding
					for source_num = 0, max_sources-1 do begin
						; Find the maximum in the inverted frame (= most negative in original)
						mx = MAX(temp_frame_neg, location, /nan)
						
						; Check if this negative peak is significant enough
						if ~finite(mx) or mx lt min_neg_peak_value then break
						
						; Get the indices of the peak
						ind = ARRAY_INDICES(temp_frame_neg, location)
						xx = ind[0] & yy = ind[1]
						
						; Check minimum separation from all previously found peaks (pos and neg)
						valid_peak = 1
						if n_elements(peak_x) gt 0 then begin
							for prev_peak = 0, n_elements(peak_x)-1 do begin
								separation = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
								if separation lt min_separation then begin
									valid_peak = 0
									break
								endif
							endfor
						endif
						
						if valid_peak then begin
							; Add this negative peak to our list
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, -1]  ; -1 for negative
							
							if debug eq 1 and kk eq 2 then begin
								; Show the actual negative value in original frame
								actual_neg_value = frame_ii[xx, yy]
								print, '    Found negative peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', actual_neg_value
							endif
						endif
						
						; Mask out this peak region to find the next one
						rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
						mask_peak = where(rad_temp le boxehs)
						if mask_peak[0] ne -1 then temp_frame_neg[mask_peak] = !values.f_nan
					endfor
					
					n_sources = n_elements(peak_x)
					n_pos_sources = n_elements(where(peak_type eq 1))
					n_neg_sources = n_elements(where(peak_type eq -1))
					
					if debug eq 1 and kk eq 2 then begin
						print, '    Total sources found in frame ', kk+1, ': ', n_sources
						print, '    Positive sources: ', n_pos_sources, ', Negative sources: ', n_neg_sources
					endif
					
					;===========================================
					; Create combined mask for all PSF peaks (positive and negative)
					;===========================================
					if n_sources gt 0 then begin
						; Initialize mask as all false
						mask = bytarr(x_dim, y_dim)
						
						; Add each source to the mask (same masking radius for pos and neg)
						for source_num = 0, n_sources-1 do begin
							rad_pos = sqrt((x_arr - peak_x[source_num])^2 + (y_arr - peak_y[source_num])^2)
							source_mask = rad_pos le boxehs
							mask = mask OR source_mask
							
							if debug eq 1 and kk eq 2 then begin
								source_type_str = peak_type[source_num] eq 1 ? 'positive' : 'negative'
								print, '    Masking ', source_type_str, ' source at (', peak_x[source_num], ',', peak_y[source_num], ')'
							endif
						endfor
					endif else begin
						; No sources found, create empty mask
						mask = bytarr(x_dim, y_dim)
						if debug eq 1 then print, '    Warning: No significant sources found in frame ', kk+1
					endelse
					
					; Optional: Create separate masks if you need different treatment
					if n_pos_sources gt 0 then begin
						pos_indices = where(peak_type eq 1)
						mask_pos = bytarr(x_dim, y_dim)
						for j = 0, n_elements(pos_indices)-1 do begin
							idx = pos_indices[j]
							rad_pos = sqrt((x_arr - peak_x[idx])^2 + (y_arr - peak_y[idx])^2)
							source_mask = rad_pos le boxehs
							mask_pos = mask_pos OR source_mask
						endfor
					endif
					
					if n_neg_sources gt 0 then begin
						neg_indices = where(peak_type eq -1)
						mask_neg = bytarr(x_dim, y_dim)
						for j = 0, n_elements(neg_indices)-1 do begin
							idx = neg_indices[j]
							rad_pos = sqrt((x_arr - peak_x[idx])^2 + (y_arr - peak_y[idx])^2)
							source_mask = rad_pos le boxehs
							mask_neg = mask_neg OR source_mask
						endfor
					endif
					
					; Create masked frame for stripe calculation
					frame_for_stripes = frame_ii
					masked_pixels = where(mask)
					if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
					
					; Additional masking based on pixel values
					frame_ii_median = median(frame_for_stripes, /double, /even)
					frame_ii_stddev = stddev(frame_for_stripes, /nan, /double)
					
					value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
								(frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; write the first mask as debugging test...
					if debug eq 1 and kk eq 2 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
						frame_for_stripes
					
					; Apply destriping (vertical then horizontal)
					destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					frame_ii = destriped_frame
					frame_for_stripes = frame_ii
					
					if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					nod_cube[*,*,kk] = destriped_frame
				
				endfor
			endif
			
			; Save the processed NOD cube
			; FIXED: Properly extract base name from current file being processed
			current_file_basename = file_basename(nod_files[i])  ; Gets full filename without path
			print, 'Debug: Processing file: ', current_file_basename

			; Remove '_cube.fits' suffix more robustly
			suffix_pos = strpos(current_file_basename, '_cube.fits')
			if suffix_pos ne -1 then begin
    				current_base_name = strmid(current_file_basename, 0, suffix_pos)
			endif else begin
    				; Fallback: remove .fits extension if _cube.fits not found
    				suffix_pos = strpos(current_file_basename, '.fits')
    				if suffix_pos ne -1 then begin
        				current_base_name = strmid(current_file_basename, 0, suffix_pos)
    				endif else begin
        				current_base_name = current_file_basename  ; No extension found
    				endelse
			endelse

			print, 'Debug: Base name after suffix removal: ', current_base_name

			; Construct unique output filename for this specific NOD file
			output_filename = output_folder + current_base_name + '_corrected_cube.fits'
			print, 'Debug: Full output path: ', output_filename
			print, '  Saving corrected cube: ', file_basename(output_filename)

			; Verify the output filename is unique for this iteration
			if i gt 0 then begin
    				print, '  (This should be different from previous files)'
			endif

			writefits, output_filename, nod_cube	
			
			; Free memory
			delvar, obj_cube
		endfor
		
		print, 'All NOD cubes processed individually and saved.'
		print, 'Total frames processed: ', total_frame_count
		
		; Create a summary file listing all the corrected cubes
		summary_file = output_folder + obj_name + '_corrected_cubes_list.txt'
		openw, lun, summary_file, /get_lun
		printf, lun, '# List of corrected cube files for object: ' + obj_name
		printf, lun, '# Generated on: ' + systime()
		printf, lun, '# Total NOD files processed: ' + string(n_nod_files)
		printf, lun, '# Total frames processed: ' + string(total_frame_count)
		printf, lun, '#'
		printf, lun, '# Format: cube_file'
		
		for i = 0, n_nod_files-1 do begin
			base_name = file_basename(nod_files[i], '_cube.fits')
			cube_file = base_name + '_corrected_cube.fits'
			printf, lun, cube_file
		endfor
		
		free_lun, lun
		print, 'Summary file created: ', file_basename(summary_file)
		
		endif; do_first_round if
		
	endif else begin; just destriping
	
		print, 'Just destriping mode - processing individual NOD cubes...'
		
		print, 'Reading-in existing bad pixel mask...'
        master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask.fits')
		
		; Find all corrected cube files (or pre-destripe files)
		if keyword_set(use_flat_corrected) then begin
			nod_files = file_search(output_folder + obj_name + '_NOD_?_nod??_grp??_corrected_cube.fits', count=n_nod_files)
		endif else if keyword_set(destripe_skysub) and (pca_skysub eq 0) then begin
			nod_files = file_search(output_folder + obj_name + '_NOD_?_nod??_grp??_skysub_cube.fits', count=n_nod_files)
		endif else if keyword_set(destripe_skysub) and (pca_skysub eq 1) then begin
		    print, 'searching: ', output_folder + 'test_pca_skysub_cube_nod??_10comp.fits'
		    nod_files = file_search(output_folder + 'test_pca_skysub_cube_nod??_10comp.fits', count=n_nod_files)
		endif else begin
			nod_files = file_search(output_folder + obj_name + '_NOD_?_nod??_no_bad_px_pre_destripe_cube.fits', count=n_nod_files)
		endelse
		
		; Sort files by nod number to ensure proper chronological order
		nod_numbers = intarr(n_nod_files)
		for i = 0, n_nod_files-1 do begin
		    nod_numbers[i] = extract_nod_number(nod_files[i])
		endfor
		
		; Sort the files by nod number
		sort_indices = sort(nod_numbers)
		nod_files = nod_files[sort_indices]
		
		print, 'Files sorted by nod number:'
		for i = 0, n_nod_files-1 do begin
		    print, '  ', file_basename(nod_files[i])
		endfor
		
		; Separate NOD_A and NOD_B files
		nod_a_files = []
		nod_b_files = []
		
		; Automatic detection via filenames
		if pca_skysub eq 0 then begin
            for i = 0, n_nod_files-1 do begin
                if strpos(nod_files[i], '_NOD_A_') ne -1 then begin
                    nod_a_files = [nod_a_files, nod_files[i]]
                endif else if strpos(nod_files[i], '_NOD_B_') ne -1 then begin
                    nod_b_files = [nod_b_files, nod_files[i]]
                endif
            endfor
        ; more specific detection via nod index being even or odd.
        endif else begin
            ; loop over nod indices and manually determine whether we are on
            ; NOD_A or NOD_B 
            for i = 0, n_nod_files-1 do begin
                ; even nod indices = NOD_A
		        if (i mod 2) eq 0 then begin
		            nod_a_files = [nod_a_files, nod_files[i]]
		        endif else begin; odd nod indices = NOD_B
		            nod_b_files = [nod_b_files, nod_files[i]]
		        endelse
		    endfor
        endelse
		
		n_nod_a = n_elements(nod_a_files)
		n_nod_b = n_elements(nod_b_files)
		
		print, 'Found ', n_nod_a, ' NOD_A files and ', n_nod_b, ' NOD_B files'
		
		if (n_nod_a eq 0) or (n_nod_b eq 0) then begin
		  print, 'ERROR: Missing one or both nodding positions'
		  return
		endif
		
		; Read first NOD file to get dimensions
		first_cube = readfits_fast(nod_a_files[0])
		cube_size = size(first_cube, /dimensions)
		x_dim = cube_size[0]
		y_dim = cube_size[1]
		frames_per_nod = cube_size[2]
		delvar, first_cube  ; Free memory immediately
		
		print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'
		
		; NEW!! We need two PSF masks for the TYC5709 data, and they are not
		; the vAPP masks.
		; Apply destriping if requested
		if (do_destripe eq 1) and obj_name eq 'tyc5709' then begin
		
			; Initialize frame counting for hot column masking
			total_frame_count = 0L
		
			; Process each NOD file separately (already sorted)
			for i = 0, n_nod_files-1 do begin
				print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
			
				; Read the current NOD cube
				nod_cube = readfits_fast(nod_files[i])
				nod_frames = (size(nod_cube, /dimensions))[2]
				print, '  Applying destriping to tyc5709 NOD cube...'
				
				; Pre-calculate arrays for masking operations
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				; Process each frame in the cube
				; Process each frame in the cube
				for kk = 0, nod_frames-1 do begin
					if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_ii_med = median(frame_ii, /double, /even)
					
					;===========================================
					; Prepare frame for peak detection
					;===========================================
					frame_ii_cen = frame_ii
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					; Calculate noise statistics for peak significance testing
					noise_stats = frame_ii_cen[where(finite(frame_ii_cen))]
					noise_mean = mean(noise_stats, /double)
					noise_stddev = stddev(noise_stats, /double)
					min_peak_value = noise_mean + min_peak_significance * noise_stddev
					
					; write the first processed frame as a debugging test...
					if debug eq 1 and kk eq 2 and i eq 0 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_cen_test.fits', /r),$
						frame_ii_cen
					
					; Enhanced source detection for both positive and negative sources					
					;===========================================
					; Find all significant PSF peaks (both positive and negative)
					;===========================================
					peak_x = []  ; Array to store x coordinates of all peaks
					peak_y = []  ; Array to store y coordinates of all peaks
					peak_type = []  ; Array to store peak type: 1=positive, -1=negative
					
					; First find positive sources
					temp_frame_pos = frame_ii_cen  ; Working copy for positive peak finding
					for source_num = 0, max_sources-1 do begin
						; Find the maximum in the current temp_frame
						mx = MAX(temp_frame_pos, location, /nan)
						
						; Check if this peak is significant enough
						if ~finite(mx) or mx lt min_peak_value then break
						
						; Get the indices of the peak
						ind = ARRAY_INDICES(temp_frame_pos, location)
						xx = ind[0] & yy = ind[1]
						
						; Check minimum separation from all previously found peaks
						valid_peak = 1
						if n_elements(peak_x) gt 0 then begin
							for prev_peak = 0, n_elements(peak_x)-1 do begin
								separation = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
								if separation lt min_separation then begin
									valid_peak = 0
									break
								endif
							endfor
						endif
						
						if valid_peak then begin
							; Add this positive peak to our list
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, 1]  ; 1 for positive
							
							if debug eq 1 and kk eq 2 then begin
								print, '    Found positive peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', mx
							endif
						endif
						
						; Mask out this peak region to find the next one
						rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
						mask_peak = where(rad_temp le boxehs)
						if mask_peak[0] ne -1 then temp_frame_pos[mask_peak] = !values.f_nan
					endfor
					
					; Now find negative sources (dips below background)
					; Create inverted frame for negative peak detection
					frame_ii_neg = -1.0 * (frame_ii - frame_ii_med)  ; Invert and center around median
					frame_ii_neg = frame_ii_neg > 0  ; Keep only negative deviations (now positive)
					frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
					
					; Calculate minimum significance threshold for negative peaks
					; (same logic but for negative deviations)
					min_neg_peak_value = min_peak_significance * noise_stddev
					
					temp_frame_neg = frame_ii_neg  ; Working copy for negative peak finding
					for source_num = 0, max_sources-1 do begin
						; Find the maximum in the inverted frame (= most negative in original)
						mx = MAX(temp_frame_neg, location, /nan)
						
						; Check if this negative peak is significant enough
						if ~finite(mx) or mx lt min_neg_peak_value then break
						
						; Get the indices of the peak
						ind = ARRAY_INDICES(temp_frame_neg, location)
						xx = ind[0] & yy = ind[1]
						
						; Check minimum separation from all previously found peaks (pos and neg)
						valid_peak = 1
						if n_elements(peak_x) gt 0 then begin
							for prev_peak = 0, n_elements(peak_x)-1 do begin
								separation = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
								if separation lt min_separation then begin
									valid_peak = 0
									break
								endif
							endfor
						endif
						
						if valid_peak then begin
							; Add this negative peak to our list
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, -1]  ; -1 for negative
							
							if debug eq 1 and kk eq 2 then begin
								; Show the actual negative value in original frame
								actual_neg_value = frame_ii[xx, yy]
								print, '    Found negative peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', actual_neg_value
							endif
						endif
						
						; Mask out this peak region to find the next one
						rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
						mask_peak = where(rad_temp le boxehs)
						if mask_peak[0] ne -1 then temp_frame_neg[mask_peak] = !values.f_nan
					endfor
					
					n_sources = n_elements(peak_x)
					n_pos_sources = n_elements(where(peak_type eq 1))
					n_neg_sources = n_elements(where(peak_type eq -1))
					
					if debug eq 1 and kk eq 2 then begin
						print, '    Total sources found in frame ', kk+1, ': ', n_sources
						print, '    Positive sources: ', n_pos_sources, ', Negative sources: ', n_neg_sources
					endif
					
					;===========================================
					; Create combined mask for all PSF peaks (positive and negative)
					;===========================================
					if n_sources gt 0 then begin
						; Initialize mask as all false
						mask = bytarr(x_dim, y_dim)
						
						; Add each source to the mask (same masking radius for pos and neg)
						for source_num = 0, n_sources-1 do begin
							rad_pos = sqrt((x_arr - peak_x[source_num])^2 + (y_arr - peak_y[source_num])^2)
							source_mask = rad_pos le boxehs
							mask = mask OR source_mask
							
							if debug eq 1 and kk eq 2 then begin
								source_type_str = peak_type[source_num] eq 1 ? 'positive' : 'negative'
								print, '    Masking ', source_type_str, ' source at (', peak_x[source_num], ',', peak_y[source_num], ')'
							endif
						endfor
					endif else begin
						; No sources found, create empty mask
						mask = bytarr(x_dim, y_dim)
						if debug eq 1 then print, '    Warning: No significant sources found in frame ', kk+1
					endelse
					
					; Optional: Create separate masks if you need different treatment
					if n_pos_sources gt 0 then begin
						pos_indices = where(peak_type eq 1)
						mask_pos = bytarr(x_dim, y_dim)
						for j = 0, n_elements(pos_indices)-1 do begin
							idx = pos_indices[j]
							rad_pos = sqrt((x_arr - peak_x[idx])^2 + (y_arr - peak_y[idx])^2)
							source_mask = rad_pos le boxehs
							mask_pos = mask_pos OR source_mask
						endfor
					endif
					
					if n_neg_sources gt 0 then begin
						neg_indices = where(peak_type eq -1)
						mask_neg = bytarr(x_dim, y_dim)
						for j = 0, n_elements(neg_indices)-1 do begin
							idx = neg_indices[j]
							rad_pos = sqrt((x_arr - peak_x[idx])^2 + (y_arr - peak_y[idx])^2)
							source_mask = rad_pos le boxehs
							mask_neg = mask_neg OR source_mask
						endfor
					endif
					
					; Create masked frame for stripe calculation
					frame_for_stripes = frame_ii
					masked_pixels = where(mask)
					if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
					
					; Additional masking based on pixel values
					frame_ii_median = median(frame_for_stripes, /double, /even)
					frame_ii_stddev = stddev(frame_for_stripes, /nan, /double)
					
					value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
								(frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; mask bad pixels before looking for stripes.
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					; write the first mask as debugging test...
					if debug eq 1 and kk eq 2 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
						frame_for_stripes
					
					; Apply destriping (vertical then horizontal)
					destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					frame_ii = destriped_frame
					frame_for_stripes = frame_ii
					
					if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; mask bad pixels before looking for stripes.
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					nod_cube[*,*,kk] = destriped_frame
				endfor
			
				; Save the processed NOD cube
				; FIXED: Properly extract base name from current file being processed
				current_file_basename = file_basename(nod_files[i])  ; Gets full filename without path
				print, 'Debug: Processing file: ', current_file_basename
		
				; Remove '_cube.fits' suffix more robustly
				suffix_pos = strpos(current_file_basename, '_cube.fits')
				if suffix_pos ne -1 then begin
						current_base_name = strmid(current_file_basename, 0, suffix_pos)
				endif else begin
						; Fallback: remove .fits extension if _cube.fits not found
						suffix_pos = strpos(current_file_basename, '.fits')
						if suffix_pos ne -1 then begin
							current_base_name = strmid(current_file_basename, 0, suffix_pos)
						endif else begin
							current_base_name = current_file_basename  ; No extension found
						endelse
				endelse
		
				print, 'Debug: Base name after suffix removal: ', current_base_name
		
				; Construct unique output filename for this specific NOD file
				if keyword_set(destripe_skysub) then begin
					output_filename = output_folder + current_base_name + '_skysub_destriped_cube.fits'
				endif else begin
					output_filename = output_folder + current_base_name + '_corrected_cube.fits'
				endelse
				
				print, 'Debug: Full output path: ', output_filename
				print, '  Saving corrected cube: ', file_basename(output_filename)
		
				; Verify the output filename is unique for this iteration
				if i gt 0 then begin
						print, '  (This should be different from previous files)'
				endif
		
				writefits, output_filename, nod_cube	
				
				; Free memory
				delvar, obj_cube
			endfor; nod_cubes for
		endif else if (do_destripe eq 1) and obj_name eq 'Alcor' then begin;tyc5709 destriping if
	        ; Initialize frame counting for hot column masking
			total_frame_count = 0L
		
			; Process each NOD file separately (already sorted)
			for i = 0, n_nod_files-1 do begin
				print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
			
				; Read the current NOD cube
				nod_cube = readfits_fast(nod_files[i])
				nod_frames = (size(nod_cube, /dimensions))[2]
				print, '  Applying destriping to Alcor NOD cube...'
				
				; Pre-calculate arrays for masking operations
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				; Process each frame in the cube
				; Process each frame in the cube
				for kk = 0, nod_frames-1 do begin
					if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_ii_med = median(frame_ii, /double, /even)
					
					;===========================================
					; Find positive PSF peak
					;===========================================
					frame_ii_cen = frame_ii
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					mx = MAX(frame_ii_cen, location)
					ind = ARRAY_INDICES(frame_ii_cen, location)
					xx = ind[0] & yy = ind[1]
					
					;===========================================
					; Create mask using vectorized operations
					;===========================================
					rad_pos = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
					ghost_rad_pos = sqrt((x_arr-(xx+11))^2 + (y_arr-(yy+58.5))^2)
					
					if pre_sky_sub eq 0 then begin
						; Find negative PSF peak
						frame_ii_cen = -1.*frame_ii
						frame_ii_cen = frame_ii_cen > (-1.*frame_ii_med)
						frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
						
						neg_mx = MAX(frame_ii_cen, location)
						neg_ind = ARRAY_INDICES(frame_ii_cen, location)
						xx_neg = neg_ind[0] & yy_neg = neg_ind[1]
						
						rad_neg = sqrt((x_arr-xx_neg)^2 + (y_arr-yy_neg)^2)
						ghost_rad_neg = sqrt((x_arr-(xx_neg+11))^2 + (y_arr-(yy_neg+58.5))^2)
						
						mask = (rad_pos le boxehs) OR (rad_neg le boxehs) OR $
							   (ghost_rad_pos le g_boxehs) OR (ghost_rad_neg le g_boxehs) OR $
							   ((rad_pos ge owa_boxehs_in) AND (rad_pos le owa_boxehs_out)) OR $
							   ((rad_neg ge owa_boxehs_in) AND (rad_neg le owa_boxehs_out))
					endif else begin
						mask = (rad_pos le boxehs) OR $
							   (ghost_rad_pos le g_boxehs) OR $
							   ((rad_pos ge owa_boxehs_in) AND (rad_pos le owa_boxehs_out))
					endelse
					
					; Create masked frame for stripe calculation
					frame_for_stripes = frame_ii
					frame_for_stripes[where(mask)] = !values.f_nan
					
					; Additional masking based on pixel values
					frame_ii_median = median(frame_for_stripes, /double, /even)
					frame_ii_stddev = stddev(frame_for_stripes, /nan, /double)
					
					value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
								(frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; mask bad pixels before looking for stripes.
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					; write the first mask as debugging test...
					if debug eq 1 and kk eq 2 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
						frame_for_stripes
					
					; Apply destriping (vertical then horizontal)
					destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					frame_ii = destriped_frame
					frame_for_stripes = frame_ii
					
					frame_for_stripes[where(mask)] = !values.f_nan
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					; mask bad pixels before looking for stripes.
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
											   /no_fit, /nodisp)
											   
					; get rid of crazy values after bad-px interpolation
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                    
                    nod_cube[*,*,kk] = destriped_frame
				endfor
			
				; Save the processed NOD cube
				; FIXED: Properly extract base name from current file being processed
				current_file_basename = file_basename(nod_files[i])  ; Gets full filename without path
				print, 'Debug: Processing file: ', current_file_basename
		
				; Remove '_cube.fits' suffix more robustly
				suffix_pos = strpos(current_file_basename, '_cube.fits')
				if suffix_pos ne -1 then begin
						current_base_name = strmid(current_file_basename, 0, suffix_pos)
				endif else begin
						; Fallback: remove .fits extension if _cube.fits not found
						suffix_pos = strpos(current_file_basename, '.fits')
						if suffix_pos ne -1 then begin
							current_base_name = strmid(current_file_basename, 0, suffix_pos)
						endif else begin
							current_base_name = current_file_basename  ; No extension found
						endelse
				endelse
		
				print, 'Debug: Base name after suffix removal: ', current_base_name
		
				; Construct unique output filename for this specific NOD file
				if keyword_set(destripe_skysub) then begin
					output_filename = output_folder + current_base_name + '_skysub_destriped_cube.fits'
				endif else begin
					output_filename = output_folder + current_base_name + '_corrected_cube.fits'
				endelse
				
				print, 'Debug: Full output path: ', output_filename
				print, '  Saving corrected cube: ', file_basename(output_filename)
		
				; Verify the output filename is unique for this iteration
				if i gt 0 then begin
						print, '  (This should be different from previous files)'
				endif
		
				writefits, output_filename, nod_cube	
				
				; Free memory
				delvar, obj_cube
			endfor; nod_cubes for
	    endif
	endelse; just destriping else

end

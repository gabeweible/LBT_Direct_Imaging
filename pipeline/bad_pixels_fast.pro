pro bad_pixels_fast, output_folder, obj_name, stripe,$
    boxhsize=boxhsize, pre_sky_sub=pre_sky_sub, bad_px_arr=bad_px_arr,$
	type=type, run=run, do_destripe=do_destripe,$
    fwhm=fwhm, just_destripe=just_destripe, do_second_round=do_second_round,$
    do_first_round=do_first_round, just_second_round=just_second_round, debug=debug,$
    create_master_masks=create_master_masks,$
    create_badpix_mask=create_badpix_mask, use_flat_corrected=use_flat_corrected,$
    darks_filename=darks_filename, coadd=coadd, vapp=vapp, destripe_skysub=destripe_skysub,$
    pca_skysub=pca_skysub, post_pca_crop=post_pca_crop
    
    compile_opt idl2, logical_predicate
    newline = string(10B)
    
    if vapp eq 1 then begin
    	boxehs = 60.
		g_boxehs = boxehs/5.
		owa_boxehs_in = 216
		owa_boxehs_out = 277
	endif
    
    if obj_name eq 'Alcor' then begin
		hc1_x1 = 385.-1.
		hc1_y1 = 250.-1.
		hc1_y2 = 278.-1.
		hc1_f1 = 0 * (20/coadd)
		hc1_f2 = 1545 * (20/coadd)
		
		hc2_x1 = 1281.-1.
		hc2_y1 = 210.-1.
		hc2_y2 = 235.-1.
		hc2_f1 = 1550 * (20/coadd)
		hc2_f2 = 2849 * (20/coadd)
		
		frame_max = 6550
		frame_min = 2000
		frame_max_post_destripe = 4000
		frame_min_post_destripe = -130
		hot_sigma = 2.0
		flat_thresh = [0.48, 1.15]
		max_sources = 2
		min_separation = 100.0
		min_peak_significance = 3.0
		max_cluster_size = 500000
	endif
	
	if obj_name eq 'tyc5709' then begin
		boxehs = 80.
		hot_sigma = 2.0
		flat_thresh = [0.73, 1.14]
		frame_max = 1400
		frame_min = 665
		frame_max_post_destripe = 575
		frame_min_post_destripe = -8.25
		max_sources = 4
		min_separation = 100.0
		min_peak_significance = 3.0
		max_cluster_size = 500000
	endif
	
	if obj_name eq 'HIP17034' then begin
		boxehs = 80.
		hot_sigma = 1.75
		flat_thresh = [0.73, 1.22]
		frame_max = 3100
		frame_min = 665
		frame_max_post_destripe = 2350
		frame_min_post_destripe = -14.
		max_sources = 8
		min_separation = 100.0
		min_peak_significance = 1.5
		max_cluster_size = 500000
	endif
	
	if obj_name eq 'HIP17900' then begin
		boxehs = 90.
		; BELOW MAY NEED TWEAKED!!
		hot_sigma = 1.75
		flat_thresh = [0.73, 1.16]
		frame_max = 6600
		frame_min = 690
		
		; 650 and 1000 are round min and max counts/s from the bkg?
		frame_max_post_destripe = 6100
		frame_min_post_destripe = -14
		
		max_sources = 30
		min_separation = 75.0
		min_peak_significance = 1.25
		max_cluster_size = 500000
	endif
    
    new_filter = 27.
    nantr = 20.0

    pre_sky_sub = keyword_set(pre_sky_sub) ? pre_sky_sub : 0
    boxhsize = keyword_set(boxhsize) ? boxhsize : 1
    
    if not keyword_set(type) then type = 'mean'
    
    if just_destripe ne 1 then begin
    
    	if do_first_round eq 1 then begin
    
            print, 'Reading in cubes for bad-px correction'
            
            print, output_folder + obj_name + '_NOD_?_nod??_grp00_cube.fits'
            nod_files = file_search(output_folder + obj_name + '_NOD_?_nod??_grp00_cube.fits', count=n_nod_files)
            if n_nod_files eq 0 then begin
                print, 'ERROR: No NOD cube files found!'
                return
            endif
            
            print, 'Found ', n_nod_files, ' NOD cube files'
            
            nod_numbers = intarr(n_nod_files)
            for i = 0, n_nod_files-1 do begin
                nod_numbers[i] = extract_nod_number(nod_files[i])
            endfor
            
            sort_indices = sort(nod_numbers)
            nod_files = nod_files[sort_indices]
            
            print, 'Files sorted by nod number:'
            for i = 0, n_nod_files-1 do begin
                print, '  ', file_basename(nod_files[i])
            endfor
            
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
            
            first_cube = readfits_fast(nod_a_files[0])
            cube_size = size(first_cube, /dimensions)
            x_dim = cube_size[0]
            y_dim = cube_size[1]
            frames_per_nod = cube_size[2]
            delvar, first_cube
            
            print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'
            
            if create_master_masks eq 1 then begin
              print, 'Generating master dark...'
              darks = readfits_fast(darks_filename)
              resistant_mean, double(darks), 3.5, master_dark, dim=3, /double
              delvar, darks
              master_dark = float(master_dark)
            
              print, 'Saving master dark to FITS file...'
              writefits, output_folder + obj_name + '_master_dark.fits', master_dark
              print, 'Master dark FITS saved.'
            
              print, 'Creating outlier-resistant master flat fields for each nod position...'
              
              print, 'Creating outlier-resistant master flat for NOD_A...'
              nod_a_individual_flats = []
              
              for i = 0, n_nod_a-1 do begin
                print, 'Processing NOD_A file ', i+1, '/', n_nod_a, ' for individual flat'
                nod_cube = readfits_fast(nod_a_files[i])
                
                for kk = 0, frames_per_nod-1 do begin
                  nod_cube[*,*,kk] -= master_dark
                endfor
                
                normalized_frames = fltarr(x_dim, y_dim, frames_per_nod)
                valid_frame_count = 0
                
                for kk = 0, frames_per_nod-1 do begin
                  frame = nod_cube[*,*,kk]
                  frame_median = float(median(double(frame), /even, /double))
                  if frame_median gt 0 then begin
                    normalized_frames[*,*,valid_frame_count] = frame / frame_median
                    valid_frame_count += 1
                  endif
                endfor
                
                if valid_frame_count gt 0 then begin
                  normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
                  
                  individual_nod_flat = dblarr(x_dim, y_dim)
                  resistant_mean, double(normalized_frames), 3.5, individual_nod_flat, dim=3, /double
                  
                  individual_median = float(median(individual_nod_flat, /even, /double))
                  individual_nod_flat = float(individual_nod_flat)
                  
                  if individual_median gt 0 then individual_nod_flat /= individual_median
                  
                  if n_elements(nod_a_individual_flats) eq 0 then begin
                    nod_a_individual_flats = fltarr(x_dim, y_dim, n_nod_a)
                  endif
                  nod_a_individual_flats[*,*,i] = individual_nod_flat
                endif
                
                delvar, nod_cube, normalized_frames, individual_nod_flat
              endfor
              
              print, 'Combining individual NOD_A flats with outlier resistance...'
              nod_a_flat = dblarr(x_dim, y_dim)
              resistant_mean, double(nod_a_individual_flats), 3.5, nod_a_flat, dim=3, /double
              
              nod_a_median = float(median(nod_a_flat, /even, /double))
              nod_a_flat = float(nod_a_flat)
              
              if nod_a_median gt 0 then nod_a_flat /= nod_a_median
              
              delvar, nod_a_individual_flats
              
              print, 'Creating outlier-resistant master flat for NOD_B...'
              nod_b_individual_flats = []
              
              for i = 0, n_nod_b-1 do begin
                print, 'Processing NOD_B file ', i+1, '/', n_nod_b, ' for individual flat'
                nod_cube = readfits_fast(nod_b_files[i])
                
                for kk = 0, frames_per_nod-1 do begin
                  nod_cube[*,*,kk] -= master_dark
                endfor
                
                normalized_frames = fltarr(x_dim, y_dim, frames_per_nod)
                valid_frame_count = 0
                
                for kk = 0, frames_per_nod-1 do begin
                  frame = nod_cube[*,*,kk]
                  frame_median = float(median(double(frame), /even, /double))
                  if frame_median gt 0 then begin
                    normalized_frames[*,*,valid_frame_count] = frame / frame_median
                    valid_frame_count += 1
                  endif
                endfor
                
                if valid_frame_count gt 0 then begin
                  normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
                  
                  individual_nod_flat = dblarr(x_dim, y_dim)
                  resistant_mean, double(normalized_frames), 3.5, individual_nod_flat, dim=3, /double
                  
                  individual_median = float(median(individual_nod_flat, /even, /double))
                  individual_nod_flat = float(individual_nod_flat)
                  
                  if individual_median gt 0 then individual_nod_flat /= individual_median
                  
                  if n_elements(nod_b_individual_flats) eq 0 then begin
                    nod_b_individual_flats = fltarr(x_dim, y_dim, n_nod_b)
                  endif
                  nod_b_individual_flats[*,*,i] = individual_nod_flat
                endif
                
                delvar, nod_cube, normalized_frames, individual_nod_flat
              endfor
              
              print, 'Combining individual NOD_B flats with outlier resistance...'
              nod_b_flat = dblarr(x_dim, y_dim)
              resistant_mean, double(nod_b_individual_flats), 3.5, nod_b_flat, dim=3, /double
              
              nod_b_median = float(median(nod_b_flat, /even, /double))
              nod_b_flat = float(nod_b_flat)
              
              if nod_b_median gt 0 then nod_b_flat /= nod_b_median
              
              delvar, nod_b_individual_flats
            
              writefits, output_folder + obj_name + '_nod_a_master_flat_resistant.fits', nod_a_flat
              writefits, output_folder + obj_name + '_nod_b_master_flat_resistant.fits', nod_b_flat
            
              super_flat = nod_a_flat
              if obj_name eq 'Alcor' then super_flat[0:897, *] = nod_b_flat[0:897, *]
              if obj_name eq 'tyc5709' then super_flat[0:970, *] = nod_b_flat[0:970, *]
              ; Full_Image splitting
              if (obj_name eq 'HIP17034') or (obj_name eq 'HIP17900') then super_flat[0:1023, *] = nod_b_flat[0:1023, *]
              delvar, nod_a_flat, nod_b_flat
              writefits, output_folder + obj_name + '_super_master_flat_resistant_raw.fits', super_flat
            
              print, 'Final normalization of stitched master flat...'
              super_median = float(median(double(super_flat), /even, /double))
              if super_median gt 0 then super_flat /= super_median
            
              writefits, output_folder + obj_name + '_super_master_flat_resistant_normalized.fits', super_flat
              print, 'Outlier-resistant master flat correction complete.'
              
              if create_badpix_mask eq 1 then begin
                print, 'Creating master bad-pixel mask from master dark and outlier-resistant master flat'
                
                gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask,$
                           hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
                           
                writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
              endif
              
            endif else begin
            
              print, 'Reading in pre-generated master files...'
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
                
                gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask,$
                                   hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
                                   
                writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
              endif else begin
                print, 'Reading-in existing bad pixel mask...'
                master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask.fits')
              endelse
            endelse
            
            ; ===========================================================================
            ; NEW: BUILD COMPLETE BAD PIXEL MASK BEFORE PROCESSING ANY CUBES
            ; ===========================================================================
            print, '=========================================='
            print, 'Building complete bad pixel mask...'
            print, '=========================================='
            
            ; Apply all hard-coded mask regions
            if obj_name eq 'Alcor' then begin
                master_badpix_mask[hc1_x1, hc1_y1:hc1_y2] = 0
                master_badpix_mask[hc2_x1, hc2_y1:hc2_y2] = 0
                master_badpix_mask[931, *] = 0
            endif
            
            ; Mask edges by 4 px, always.
            master_badpix_mask[0:4, *] = 0
            master_badpix_mask[x_dim-4:x_dim-1, *] = 0
            master_badpix_mask[*, y_dim-4:y_dim-1] = 0
            master_badpix_mask[*, 0:4] = 0
            
            ;bad readout channel, which is always bad.
            master_badpix_mask[896:959, *] = 0
            
            ; known bad column
            master_badpix_mask[1984, *] = 0
            
            if obj_name eq 'tyc5709' then begin
                master_badpix_mask[1408, 228-10:244+10] = 0
                master_badpix_mask[*, 434:536] = 0
            endif
            
            if obj_name eq 'HIP17034' then begin
                master_badpix_mask[*, 870:1100] = 0
            endif
            
            if obj_name eq 'HIP17900' then begin
                ; cross-talk column issues to mask
                master_badpix_mask[1408, 1500:1570] = 0
                master_badpix_mask[1344:1345, 470:540] = 0
                master_badpix_mask[447, 440:510] = 0
                master_badpix_mask[575, 1525:1595] = 0
            endif
            
            print, 'Saving complete bad pixel mask...'
            writefits, output_folder+obj_name+'_master_badpix_mask_complete.fits', master_badpix_mask
            
            n_bad_pix = n_elements(where(master_badpix_mask eq 0))
            print, 'Complete bad pixel mask created with ', n_bad_pix, ' bad pixels (', $
                   string((100.* n_bad_pix) / (x_dim * y_dim), format='(F5.2)'), '% of each frame)'
            print, '=========================================='
            
            ; ===========================================================================
            ; MEMORY EFFICIENT PROCESSING: Process each NOD file separately
            ; ===========================================================================
            print, 'Processing NOD files individually for bad pixel correction and flat fielding...'
            
            total_frame_count = 0L
            seed = 12345L
            for i = 0, n_nod_files-1 do begin
                print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
                
                nod_cube = readfits_fast(nod_files[i])
                nod_frames = (size(nod_cube, /dimensions))[2]
                
                working_mask = master_badpix_mask
                
                print, '  Applying bad pixel correction to ', nod_frames, ' frames'
                print, '  Using complete bad pixel mask uniformly across all cubes'
                
                for kk = 0, nod_frames-1 do begin
                    if (kk mod 41) eq 0 then print, '    Processing frame ', kk+1, '/', nod_frames
                    
                    frame = nod_cube[*,*,kk]
                    frame -= master_dark
                    frame /= super_flat
                    
                    if obj_name eq 'Alcor' then begin
                        global_frame_idx = total_frame_count + kk
                        frame_mask = working_mask
                        
                        if (global_frame_idx lt hc1_f1) or (global_frame_idx gt hc1_f2) then begin
                            frame_mask[hc1_x1, hc1_y1:hc1_y2] = 1
                        endif
                        if (global_frame_idx lt hc2_f1) or (global_frame_idx gt hc2_f2) then begin
                            frame_mask[hc2_x1, hc2_y1:hc2_y2] = 1
                        endif
                    endif else begin
                        frame_mask = working_mask
                    endelse
                    
                    original_frame = frame
                    if (i eq 0) and (kk eq 0) then writefits, '~/Desktop/original_frame_test.fits', original_frame
                    
                    frame[where(frame_mask eq 0)] = !values.f_nan
                    frame[where(frame eq 0)] = !values.f_nan
                    frame[where(frame gt frame_max)] = !values.f_nan
                    frame[where(frame lt frame_min)] = !values.f_nan
                    
                    original_nan_mask = ~finite(frame)
                    
                    efficient_nan_correction, frame, fixed_frame, $
                        max_cluster_size=max_cluster_size, npix=24, /weight, /silent
                        
                    fixed_frame[where(fixed_frame gt frame_max)] = !values.f_nan
                    fixed_frame[where(fixed_frame lt frame_min)] = !values.f_nan
                    
                    interpolated_mask = (original_nan_mask and finite(fixed_frame))
                    
                    good_pixels = where(finite(original_frame) and frame_mask eq 1, n_good)
                    if n_good gt 100 then begin
                        good_values = original_frame[good_pixels]
                        
                        median_val = median(double(good_values), /even, /double)
                        mad_val = float(median(abs(double(good_values) - median_val), /even, /double))
                        median_val = float(median_val)
                        
                        ; double the noise to down-weight intepolations!
                        noise_sigma = 2 * 1.4826 * mad_val
                        
                        interpolated_indices = where(interpolated_mask, n_interpolated)
                        if n_interpolated gt 0 then begin
                            gaussian_noise = randomn(seed, n_interpolated) * noise_sigma
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
                
                total_frame_count += nod_frames
                
                if (do_destripe eq 1) and (run eq 2) and (obj_name ne 'tyc5709') then begin
                    print, '  Applying destriping to NOD cube...'
                    
                    x_arr = indgen(x_dim) # replicate(1, y_dim)
                    y_arr = replicate(1, x_dim) # indgen(y_dim)
                    
                    for kk = 0, nod_frames-1 do begin
                        if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
                        frame_ii = nod_cube[*,*,kk]
                        frame_ii_med = float(median(double(frame_ii), /double, /even))
                        
                        frame_ii_cen = frame_ii
                        frame_ii_cen = frame_ii_cen > frame_ii_med
                        frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
                        
                        mx = MAX(frame_ii_cen, location)
                        ind = ARRAY_INDICES(frame_ii_cen, location)
                        xx = ind[0] & yy = ind[1]
                        
                        rad_pos = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                        ghost_rad_pos = sqrt((x_arr-(xx+11))^2 + (y_arr-(yy+58.5))^2)
                        
                        if pre_sky_sub eq 0 then begin
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
                        
                        frame_for_stripes = frame_ii
                        frame_for_stripes[where(mask)] = !values.f_nan
                        
                        frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                        frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
                        
                        value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
                                    (frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
                        
                        nanmask = where(value_mask)
                        if nanmask[0] ne -1 then begin
                            frame_for_stripes[nanmask] = !values.f_nan
                            frame_ii[nanmask] = !values.f_nan
                        endif
                        
                        if debug eq 1 and kk eq 2 then writefits,$
                            strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
                            frame_for_stripes
                        
                        destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
                                                   /no_fit, /nodisp)
                                                   
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
                                                   
                        destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                        destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                        
                        nod_cube[*,*,kk] = destriped_frame
                    endfor
                endif
                
                if (do_destripe eq 1) and ((obj_name eq 'tyc5709') or (obj_name eq 'HIP17034') or (obj_name eq 'HIP17900')) then begin
                    print, '  Applying destriping to NOD cube...'
                    
                    x_arr = indgen(x_dim) # replicate(1, y_dim)
                    y_arr = replicate(1, x_dim) # indgen(y_dim)
                    
                    for kk = 0, nod_frames-1 do begin
                        if (kk mod 41) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
                        frame_ii = nod_cube[*,*,kk]
                        frame_ii_med = float(median(double(frame_ii), /double, /even))
                        
                        frame_ii_cen = frame_ii
                        frame_ii_cen = frame_ii_cen > frame_ii_med
                        frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
                        
                        noise_stats = frame_ii_cen[where(finite(frame_ii_cen))]
                        noise_mean = float(mean(double(noise_stats), /double))
                        noise_stddev = float(stddev(double(noise_stats), /double))
                        min_peak_value = noise_mean + min_peak_significance * noise_stddev
                        
                        if debug eq 1 and kk eq 2 and i eq 0 then writefits,$
                            strcompress('~/Desktop/nod_'+string(i)+'_cen_test.fits', /r),$
                            frame_ii_cen
                        
                        peak_x = []
                        peak_y = []
                        peak_type = []
                        
                        temp_frame_pos = frame_ii_cen
                        for source_num = 0, max_sources-1 do begin
                            mx = MAX(temp_frame_pos, location, /nan)
                            
                            if ~finite(mx) or mx lt min_peak_value then break
                            
                            ind = ARRAY_INDICES(temp_frame_pos, location)
                            xx = ind[0] & yy = ind[1]
                            
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
                                peak_x = [peak_x, xx]
                                peak_y = [peak_y, yy]
                                peak_type = [peak_type, 1]
                                
                                if debug eq 1 and kk eq 2 then begin
                                    print, '    Found positive peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', mx
                                endif
                            endif
                            
                            rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                            mask_peak = where(rad_temp le boxehs)
                            if mask_peak[0] ne -1 then temp_frame_pos[mask_peak] = !values.f_nan
                        endfor
                        
                        frame_ii_neg = -1.0 * (frame_ii - frame_ii_med)
                        frame_ii_neg = frame_ii_neg > 0
                        frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
                        
                        min_neg_peak_value = min_peak_significance * noise_stddev
                        
                        temp_frame_neg = frame_ii_neg
                        for source_num = 0, max_sources-1 do begin
                            mx = MAX(temp_frame_neg, location, /nan)
                            
                            if ~finite(mx) or mx lt min_neg_peak_value then break
                            
                            ind = ARRAY_INDICES(temp_frame_neg, location)
                            xx = ind[0] & yy = ind[1]
                            
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
                                peak_x = [peak_x, xx]
                                peak_y = [peak_y, yy]
                                peak_type = [peak_type, -1]
                                
                                if debug eq 1 and kk eq 2 then begin
                                    actual_neg_value = frame_ii[xx, yy]
                                    print, '    Found negative peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', actual_neg_value
                                endif
                            endif
                            
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
                        
                        if n_sources gt 0 then begin
                            mask = bytarr(x_dim, y_dim)
                            
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
                            mask = bytarr(x_dim, y_dim)
                            if debug eq 1 then print, '    Warning: No significant sources found in frame ', kk+1
                        endelse
                        
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
                        
                        frame_for_stripes = frame_ii
                        masked_pixels = where(mask)
                        if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
                        
                        frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                        frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
                        
                        low_nantr_thresh = frame_ii_median - nantr*frame_ii_stddev
                        high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                        if debug eq 1 and kk mod 41 then begin
                            print, 'low_nantr_thresh: ', low_nantr_thresh
                            print, 'high_nantr_thresh: ', high_nantr_thresh
                        endif
                        
                        value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                                    (frame_for_stripes gt high_nantr_thresh)
                        
                        nanmask = where(value_mask)
                        if nanmask[0] ne -1 then begin
                            frame_for_stripes[nanmask] = !values.f_nan
                            frame_ii[nanmask] = !values.f_nan
                        endif
                            
                        if debug eq 1 and kk eq 2 then writefits,$
                            strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
                            frame_for_stripes
                            
                        if i eq 0 and kk eq 2 then begin
                            print, '  Checking stripe statistics for nod00, frame 2:'
                            print, '  frame_for_stripes finite pixels: ', n_elements(where(finite(frame_for_stripes)))
                            print, '  frame_ii finite pixels: ', n_elements(where(finite(frame_ii)))
                            print, '  frame_for_stripes median: ', median(frame_for_stripes, /double)
                            print, '  frame_ii median: ', median(frame_ii, /double)
                        endif
                        
                        if i eq 2 and kk eq 2 then begin
                            print, '  Checking stripe statistics for nod02, frame 2:'
                            print, '  frame_for_stripes finite pixels: ', n_elements(where(finite(frame_for_stripes)))
                            print, '  frame_ii finite pixels: ', n_elements(where(finite(frame_ii)))
                            print, '  frame_for_stripes median: ', median(frame_for_stripes, /double)
                            print, '  frame_ii median: ', median(frame_ii, /double)
                        endif
                        
                        ; Create quadrant-specific frames for vertical destriping
                        frame_for_stripes_top = frame_for_stripes
                        frame_for_stripes_bottom = frame_for_stripes
                        
                        ; Mask out the bottom half for top quadrant processing
                        frame_for_stripes_top[*, 0:1023] = !values.f_nan
                        
                        ; Mask out the top half for bottom quadrant processing  
                        frame_for_stripes_bottom[*, 1024:2047] = !values.f_nan
                        
                        if debug eq 1 and kk eq 2 then BEGIN
                            writefits,$
                                strcompress('~/Desktop/nod_'+string(i)+'_masked_test_top.fits', /r),$
                                frame_for_stripes_top
                            
                            writefits,$
                                strcompress('~/Desktop/nod_'+string(i)+'_masked_test_bottom.fits', /r),$
                                frame_for_stripes_bottom
                        endif
                        
                        ; Destripe top quadrant (rows 1024-2047)
                        destriped_frame_top = destripe(frame_for_stripes_top, frame_ii, 90., fraction=0.02, $
                                                       /no_fit, /nodisp)
                                                       
                        destriped_frame_top[where(destriped_frame_top gt frame_max_post_destripe)] = !values.f_nan
                        destriped_frame_top[where(destriped_frame_top lt frame_min_post_destripe)] = !values.f_nan
                        
                        ; Destripe bottom quadrant (rows 0-1023)
                        destriped_frame_bottom = destripe(frame_for_stripes_bottom, frame_ii, 90., fraction=0.02, $
                                                          /no_fit, /nodisp)
                                                          
                        destriped_frame_bottom[where(destriped_frame_bottom gt frame_max_post_destripe)] = !values.f_nan
                        destriped_frame_bottom[where(destriped_frame_bottom lt frame_min_post_destripe)] = !values.f_nan
                        
                        ; Combine the two destriped quadrants
                        destriped_frame = frame_ii  ; Start with original frame
                        destriped_frame[*, 0:1023] = destriped_frame_bottom[*, 0:1023]
                        destriped_frame[*, 1024:2047] = destriped_frame_top[*, 1024:2047]
                        
                        if i eq 0 and kk eq 2 then begin
                            print, '  After vertical destripe (90deg) with quadrant split:'
                            print, '  destriped_frame finite pixels: ', n_elements(where(finite(destriped_frame)))
                            print, '  destriped_frame median: ', median(destriped_frame, /double)
                            print, '  Number of NaNs added: ', n_elements(where(~finite(destriped_frame))) - n_elements(where(~finite(frame_ii)))
                        endif
                        
                        ; Continue with the rest of the processing
                        frame_ii = destriped_frame
                        frame_for_stripes = frame_ii
                        
                        frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                        frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
                        
                        low_nantr_thresh = frame_ii_median - nantr*frame_ii_stddev
                        high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                        if debug eq 1 and kk mod 41 then begin
                            print, 'low_nantr_thresh: ', low_nantr_thresh
                            print, 'high_nantr_thresh: ', high_nantr_thresh
                        endif
                        
                        value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                                    (frame_for_stripes gt high_nantr_thresh)
                        
                        nanmask = where(value_mask)
                        if nanmask[0] ne -1 then begin
                            frame_for_stripes[nanmask] = !values.f_nan
                            frame_ii[nanmask] = !values.f_nan
                        endif
                        
                        if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
                        
                        ; Horizontal destriping (0 degree) - this remains unchanged
                        destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
                                                   /no_fit, /nodisp)
                                                   
                        destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                        destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                        
                        nod_cube[*,*,kk] = destriped_frame
                    
                    endfor
                endif
                
                current_file_basename = file_basename(nod_files[i])
                print, 'Debug: Processing file: ', current_file_basename
    
                suffix_pos = strpos(current_file_basename, '_cube.fits')
                if suffix_pos ne -1 then begin
                        current_base_name = strmid(current_file_basename, 0, suffix_pos)
                endif else begin
                        suffix_pos = strpos(current_file_basename, '.fits')
                        if suffix_pos ne -1 then begin
                            current_base_name = strmid(current_file_basename, 0, suffix_pos)
                        endif else begin
                            current_base_name = current_file_basename
                        endelse
                endelse
    
                print, 'Debug: Base name after suffix removal: ', current_base_name
    
                output_filename = output_folder + current_base_name + '_corrected_cube.fits'
                print, 'Debug: Full output path: ', output_filename
                print, '  Saving corrected cube: ', file_basename(output_filename)
    
                if i gt 0 then begin
                        print, '  (This should be different from previous files)'
                endif
    
                writefits, output_filename, nod_cube	
                
                delvar, nod_cube
            endfor
            
            print, 'All NOD cubes processed individually and saved.'
            print, 'Total frames processed: ', total_frame_count
            
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
		
		endif

	endif else begin ; begin for just_destripe = 1 (e.g., post–PCA skysub)
	
		print, 'Just destriping mode - processing individual NOD cubes...'
		
		print, 'Reading-in existing bad pixel mask...'
        master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask.fits')
		
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
		
		nod_numbers = intarr(n_nod_files)
		for i = 0, n_nod_files-1 do begin
		    nod_numbers[i] = extract_nod_number(nod_files[i])
		endfor
		
		sort_indices = sort(nod_numbers)
		nod_files = nod_files[sort_indices]
		
		print, 'Files sorted by nod number:'
		for i = 0, n_nod_files-1 do begin
		    print, '  ', file_basename(nod_files[i])
		endfor
		
		nod_a_files = []
		nod_b_files = []
		
		if pca_skysub eq 0 then begin
            for i = 0, n_nod_files-1 do begin
                if strpos(nod_files[i], '_NOD_A_') ne -1 then begin
                    nod_a_files = [nod_a_files, nod_files[i]]
                endif else if strpos(nod_files[i], '_NOD_B_') ne -1 then begin
                    nod_b_files = [nod_b_files, nod_files[i]]
                endif
            endfor
        endif else begin
            for i = 0, n_nod_files-1 do begin
		        if (i mod 2) eq 0 then begin
		            nod_a_files = [nod_a_files, nod_files[i]]
		        endif else begin
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
		
		first_cube = readfits_fast(nod_a_files[0])
		cube_size = size(first_cube, /dimensions)
		x_dim = cube_size[0]
		y_dim = cube_size[1]
		frames_per_nod = cube_size[2]
		delvar, first_cube
		
		print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'
		
		; Below is for L-band direct imaging (no vAPP, like for WISPIT 2 and HIP 17034)
		if (do_destripe eq 1) and ((obj_name eq 'tyc5709') or (obj_name eq 'HIP17034') or (obj_name eq 'HIP17900')) then begin
		
			total_frame_count = 0L
		
			for i = 0, n_nod_files-1 do begin
				print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
			
				nod_cube = readfits_fast(nod_files[i])
				nod_frames = (size(nod_cube, /dimensions))[2]
				print, 'Applying destriping to tyc5709 or PLEIAD NOD cube...'
				
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				for kk = 0, nod_frames-1 do begin
					if (kk mod 41) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_for_stripes = frame_ii
					
					if post_pca_crop eq 1 then begin
					    if ODD(i) then begin; odd nod index <=> NOD_B (right)
					        frame_for_stripes[0:1023,*] = !values.f_nan
					    endif else begin; even nod index <=> NOD_A (left)
					        frame_for_stripes[1023:2047,*] = !values.f_nan
					    endelse
					endif
					
					frame_ii_med = float(median(double(frame_for_stripes), /double, /even))
					
					frame_ii_cen = frame_for_stripes
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					noise_stats = frame_ii_cen[where(finite(frame_ii_cen))]
					noise_mean = float(mean(double(noise_stats), /double))
					noise_stddev = float(stddev(double(noise_stats), /double))
					min_peak_value = noise_mean + min_peak_significance * noise_stddev
					
					if debug eq 1 and kk eq 2 and i eq 0 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_cen_test.fits', /r),$
						frame_ii_cen
					
					peak_x = []
					peak_y = []
					peak_type = []
					
					temp_frame_pos = frame_ii_cen
					for source_num = 0, max_sources-1 do begin
						mx = MAX(temp_frame_pos, location, /nan)
						
						if ~finite(mx) or mx lt min_peak_value then break
						
						ind = ARRAY_INDICES(temp_frame_pos, location)
						xx = ind[0] & yy = ind[1]
						
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
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, 1]
							
							if debug eq 1 and kk eq 2 then begin
								print, '    Found positive peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', mx
							endif
						endif
						
						rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
						mask_peak = where(rad_temp le boxehs)
						if mask_peak[0] ne -1 then temp_frame_pos[mask_peak] = !values.f_nan
					endfor
					
					frame_ii_neg = -1.0 * (frame_for_stripes - frame_ii_med)
					frame_ii_neg = frame_ii_neg > 0
					frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
					
					min_neg_peak_value = min_peak_significance * noise_stddev
					
					temp_frame_neg = frame_ii_neg
					for source_num = 0, max_sources-1 do begin
						mx = MAX(temp_frame_neg, location, /nan)
						
						if ~finite(mx) or mx lt min_neg_peak_value then break
						
						ind = ARRAY_INDICES(temp_frame_neg, location)
						xx = ind[0] & yy = ind[1]
						
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
							peak_x = [peak_x, xx]
							peak_y = [peak_y, yy]
							peak_type = [peak_type, -1]
							
							if debug eq 1 and kk eq 2 then begin
								actual_neg_value = frame_ii[xx, yy]
								print, '    Found negative peak ', n_elements(peak_x), ' at (', xx, ',', yy, ') with value ', actual_neg_value
							endif
						endif
						
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
					
					if n_sources gt 0 then begin
						mask = bytarr(x_dim, y_dim)
						
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
						mask = bytarr(x_dim, y_dim)
						if debug eq 1 then print, '    Warning: No significant sources found in frame ', kk+1
					endelse
					
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
					
					masked_pixels = where(mask)
					if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan
					
					frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
					frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
					
					low_nantr_thresh = frame_ii_median - nantr*frame_ii_stddev
					high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
					if debug eq 1 and kk mod 41 then begin
                        print, 'low_nantr_thresh: ', low_nantr_thresh
                        print, 'high_nantr_thresh: ', high_nantr_thresh
                    endif
					
					value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
								(frame_for_stripes gt high_nantr_thresh)
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
                        frame_ii[nanmask] = !values.f_nan
					endif
					
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
                    ; Create quadrant-specific frames for vertical destriping
                    frame_for_stripes_top = frame_for_stripes
                    frame_for_stripes_bottom = frame_for_stripes
                    
                    ; Mask out the bottom half for top quadrant processing
                    frame_for_stripes_top[*, 0:1023] = !values.f_nan
                    
                    ; Mask out the top half for bottom quadrant processing  
                    frame_for_stripes_bottom[*, 1024:2047] = !values.f_nan
                    
                    if debug eq 1 and kk eq 2 then BEGIN
                        writefits,$
						    strcompress('~/Desktop/nod_'+string(i)+'_masked_test_top.fits', /r),$
						    frame_for_stripes_top
						
					    writefits,$
						    strcompress('~/Desktop/nod_'+string(i)+'_masked_test_bottom.fits', /r),$
						    frame_for_stripes_bottom
                    endif
                    
                    ; Destripe top quadrant (rows 0-1023)
                    destriped_frame_top = destripe(frame_for_stripes_top, frame_ii, 90., fraction=0.02, $
                                                   /no_fit, /nodisp)
                                                   
                    destriped_frame_top[where(destriped_frame_top gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame_top[where(destriped_frame_top lt frame_min_post_destripe)] = !values.f_nan
                    
                    ; Destripe bottom quadrant (rows 1024-2047)
                    destriped_frame_bottom = destripe(frame_for_stripes_bottom, frame_ii, 90., fraction=0.02, $
                                                      /no_fit, /nodisp)
                                                      
                    destriped_frame_bottom[where(destriped_frame_bottom gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame_bottom[where(destriped_frame_bottom lt frame_min_post_destripe)] = !values.f_nan
                    
                    ; Combine the two destriped quadrants
                    destriped_frame = frame_ii  ; Start with original frame
                    destriped_frame[*, 0:1023] = destriped_frame_bottom[*, 0:1023]
                    destriped_frame[*, 1024:2047] = destriped_frame_top[*, 1024:2047]
                    
                    ; Continue with the rest of the processing
                    frame_ii = destriped_frame
                    frame_for_stripes = frame_ii
                    
                    frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
					frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
					
					low_nantr_thresh = frame_ii_median - nantr*frame_ii_stddev
					high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
					if debug eq 1 and kk mod 41 then begin
                        print, 'low_nantr_thresh: ', low_nantr_thresh
                        print, 'high_nantr_thresh: ', high_nantr_thresh
                    endif
					
					value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
								(frame_for_stripes gt high_nantr_thresh)
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
                    if masked_pixels[0] ne -1 then frame_for_stripes[masked_pixels] = !values.f_nan

                    frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
                    
                    ; Horizontal destriping (0 degree) - this remains unchanged
                    destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
                                               /no_fit, /nodisp)
                                               
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                    
                    nod_cube[*,*,kk] = destriped_frame
				endfor
			
				current_file_basename = file_basename(nod_files[i])
				print, 'Debug: Processing file: ', current_file_basename
		
				suffix_pos = strpos(current_file_basename, '_cube.fits')
				if suffix_pos ne -1 then begin
						current_base_name = strmid(current_file_basename, 0, suffix_pos)
				endif else begin
						suffix_pos = strpos(current_file_basename, '.fits')
						if suffix_pos ne -1 then begin
							current_base_name = strmid(current_file_basename, 0, suffix_pos)
						endif else begin
							current_base_name = current_file_basename
						endelse
				endelse
		
				print, 'Debug: Base name after suffix removal: ', current_base_name
		
				if keyword_set(destripe_skysub) then begin
					output_filename = output_folder + current_base_name + '_skysub_destriped_cube.fits'
				endif else begin
					output_filename = output_folder + current_base_name + '_corrected_cube.fits'
				endelse
				
				print, 'Debug: Full output path: ', output_filename
				print, '  Saving corrected cube: ', file_basename(output_filename)
		
				if i gt 0 then begin
						print, '  (This should be different from previous files)'
				endif
		
				writefits, output_filename, nod_cube	
				
				delvar, nod_cube
			endfor
			
			; Below is for just_destripe=1 and M-band vAPP data (like for Alcor)
		endif else if (do_destripe eq 1) and obj_name eq 'Alcor' then begin
	        total_frame_count = 0L
		
			for i = 0, n_nod_files-1 do begin
				print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])
			
				nod_cube = readfits_fast(nod_files[i])
				nod_frames = (size(nod_cube, /dimensions))[2]
				print, '  Applying destriping to Alcor NOD cube...'
				
				x_arr = indgen(x_dim) # replicate(1, y_dim)
				y_arr = replicate(1, x_dim) # indgen(y_dim)
				
				for kk = 0, nod_frames-1 do begin
					if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
					frame_ii = nod_cube[*,*,kk]
					frame_ii_med = median(frame_ii, /double, /even)
					
					frame_ii_cen = frame_ii
					frame_ii_cen = frame_ii_cen > frame_ii_med
					frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)
					
					mx = MAX(frame_ii_cen, location)
					ind = ARRAY_INDICES(frame_ii_cen, location)
					xx = ind[0] & yy = ind[1]
					
					rad_pos = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
					ghost_rad_pos = sqrt((x_arr-(xx+11))^2 + (y_arr-(yy+58.5))^2)
					
					if pre_sky_sub eq 0 then begin
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
					
					frame_for_stripes = frame_ii
					frame_for_stripes[where(mask)] = !values.f_nan
					
					frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
					frame_ii_stddev = float(stddev(double(frame_for_stripes), /nan, /double))
					
					value_mask = (frame_for_stripes lt (frame_ii_median - nantr*frame_ii_stddev)) OR $
								(frame_for_stripes gt (frame_ii_median + nantr*frame_ii_stddev))
					
					nanmask = where(value_mask)
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					if debug eq 1 and kk eq 2 then writefits,$
						strcompress('~/Desktop/nod_'+string(i)+'_masked_test.fits', /r),$
						frame_for_stripes
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 90., fraction=0.02, $
											   /no_fit, /nodisp)
											   
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
					
					frame_ii = destriped_frame
					frame_for_stripes = frame_ii
					
					frame_for_stripes[where(mask)] = !values.f_nan
					if nanmask[0] ne -1 then begin
					    frame_for_stripes[nanmask] = !values.f_nan
					    frame_ii[nanmask] = !values.f_nan
					endif
					
					frame_for_stripes[where(master_badpix_mask eq 0)] = !values.f_nan
					
					destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, $
											   /no_fit, /nodisp)
											   
                    destriped_frame[where(destriped_frame gt frame_max_post_destripe)] = !values.f_nan
                    destriped_frame[where(destriped_frame lt frame_min_post_destripe)] = !values.f_nan
                    
                    nod_cube[*,*,kk] = destriped_frame
				endfor
			
				current_file_basename = file_basename(nod_files[i])
				print, 'Debug: Processing file: ', current_file_basename
		
				suffix_pos = strpos(current_file_basename, '_cube.fits')
				if suffix_pos ne -1 then begin
						current_base_name = strmid(current_file_basename, 0, suffix_pos)
				endif else begin
						suffix_pos = strpos(current_file_basename, '.fits')
						if suffix_pos ne -1 then begin
							current_base_name = strmid(current_file_basename, 0, suffix_pos)
						endif else begin
							current_base_name = current_file_basename
						endelse
				endelse
		
				print, 'Debug: Base name after suffix removal: ', current_base_name
		
				if keyword_set(destripe_skysub) then begin
					output_filename = output_folder + current_base_name + '_skysub_destriped_cube.fits'
				endif else begin
					output_filename = output_folder + current_base_name + '_corrected_cube.fits'
				endelse
				
				print, 'Debug: Full output path: ', output_filename
				print, '  Saving corrected cube: ', file_basename(output_filename)
		
				if i gt 0 then begin
						print, '  (This should be different from previous files)'
				endif
		
				writefits, output_filename, nod_cube	
				
				delvar, nod_cube
			endfor
	    endif
	endelse

end
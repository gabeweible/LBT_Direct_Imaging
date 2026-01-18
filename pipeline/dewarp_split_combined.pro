;+ 
; NAME:
;   dewarp_split_combined
; PURPOSE:
;   Dewarps and crops sky-subtracted nod-based cubes produced by sky_sub.pro
;   into left/right aperture output cubes for each NOD position.
; 
; INPUTS:
;   dir      - Directory containing sky-subtracted nod cubes and parang files
;   obj      - Target object name
;   stripe   - 'center1024' or 'second512'
;   Kx_sx, Ky_sx - Polynomial distortion coefficients for SX mirror
;   Kx_dx, Ky_dx - Polynomial distortion coefficients for DX mirror
;
; KEYWORDS:
;   run              - 1 (SX), 2 (DX), or 5 (both) [default: 5]
;   do_smooth        - Apply Gaussian smoothing [default: 0]
;   fwhm             - FWHM of smoothing kernel [default: 13.0]
;   half_cropped_sz  - Half size of cropped region [default: 190]
;   aperture         - 'left', 'right', or 'both' [default: 'both']
;   hp_width         - FWHM of optional high-pass filter [default: 0]
;   nod_filter       - Process only specific nod position: 'dith1', 'dith2', or 'both' [default: 'both']
;
; AUTHOR: GPT (2024), adapted from user-provided IDL pipeline
;-
pro dewarp_split_combined, dir, obj, stripe, Kx_sx, Ky_sx, Kx_dx, Ky_dx, $
    run=run, do_smooth=do_smooth, fwhm=fwhm, half_cropped_sz=half_cropped_sz, $
    aperture=aperture, hp_width=hp_width, nod_filter=nod_filter, debug=debug,$
    destripe_skysub=destripe_skysub, pca_skysub=pca_skysub, filled=filled

  COMPILE_OPT IDL2
  newline = STRING(10B)

  if ~keyword_set(run) then run = 5
  if ~keyword_set(do_smooth) then do_smooth = 1
  if ~keyword_set(half_cropped_sz) then half_cropped_sz = 190
  if ~keyword_set(aperture) then aperture = 'both'
  if ~keyword_set(hp_width) then hp_width = 0
  if ~keyword_set(nod_filter) then nod_filter = 'both'
  if ~keyword_set(filled) then filled=0

  ; Validate nod_filter parameter
  if nod_filter ne 'both' and nod_filter ne 'dith1' and nod_filter ne 'dith2' then begin
    message, 'Invalid nod_filter: use "dith1", "dith2", or "both"'
  endif

  ; Define PSF positions and pre-compute cropping indices
  crop_size = 2*half_cropped_sz+1
  
  
	if stripe eq 'center1024' and obj eq 'hii1348' then begin
		SX_PSF_Y = 701
		DX_PSF_Y = 259
		SX_PSF_X = 547
		DX_PSF_X = 1475
		LEFT_PSF_X = SX_PSF_X - 200
		RIGHT_PSF_X = SX_PSF_X + 200
		
		; Pre-compute SX cropping indices
		sx_left_x0 = LEFT_PSF_X - half_cropped_sz
		sx_left_x1 = LEFT_PSF_X + half_cropped_sz
		sx_right_x0 = RIGHT_PSF_X - half_cropped_sz
		sx_right_x1 = RIGHT_PSF_X + half_cropped_sz
		sx_y0 = SX_PSF_Y - half_cropped_sz
		sx_y1 = SX_PSF_Y + half_cropped_sz
    
    ; TYC 5709 - START WITH HII 1348 GUESSES
	endif else if stripe eq 'center1024' and obj eq 'tyc5709' then begin
		print, 'TYC 5709 Center1024 correction enabled.', newline
		
		; in full-frmae dewarped coordinates.
		; "1" is for the "left" nod. DX is top, SX is bottom (opposite of HII 1348)
		DX_PSF_X1 = 476-2
		DX_PSF_Y1 = 1331-18-3
		
		SX_PSF_X1 = 550+5
		SX_PSF_Y1 = 746-15
		
		DX_PSF_X2 = 1390+4
		DX_PSF_Y2 = 1294-8
		
		SX_PSF_X2 = 1488-19+5
		SX_PSF_Y2 = 754-5-6+6
		
		; Pre-compute cropping indices
		; bottom-left
		sx_left_x0 = SX_PSF_X1- half_cropped_sz
		sx_left_x1 = SX_PSF_X1 + half_cropped_sz
		sx_left_y0 = SX_PSF_Y1 - half_cropped_sz
		sx_left_y1 = SX_PSF_Y1 + half_cropped_sz
		
		; top-left
		dx_left_x0 = DX_PSF_X1- half_cropped_sz
		dx_left_x1 = DX_PSF_X1 + half_cropped_sz
		dx_left_y0 = DX_PSF_Y1 - half_cropped_sz
		dx_left_y1 = DX_PSF_Y1 + half_cropped_sz
		
		; bottom-right
		sx_right_x0 = SX_PSF_X2- half_cropped_sz
		sx_right_x1 = SX_PSF_X2 + half_cropped_sz
		sx_right_y0 = SX_PSF_Y2 - half_cropped_sz
		sx_right_y1 = SX_PSF_Y2 + half_cropped_sz
		
		; top-right
		dx_right_x0 = DX_PSF_X2- half_cropped_sz
		dx_right_x1 = DX_PSF_X2 + half_cropped_sz
		dx_right_y0 = DX_PSF_Y2 - half_cropped_sz
		dx_right_y1 = DX_PSF_Y2 + half_cropped_sz
    
	  endif else if stripe eq 'second512' then begin
		; ALCOR
		DX_PSF_X1 = 449
		DX_PSF_X2 = 1372
		DX_PSF_Y1 = 1279
		DX_PSF_Y2 = 1258
		
		; top-left
		dx_left_x0 = DX_PSF_X1- half_cropped_sz
		dx_left_x1 = DX_PSF_X1 + half_cropped_sz
		dx_left_y0 = DX_PSF_Y1 - half_cropped_sz
		dx_left_y1 = DX_PSF_Y1 + half_cropped_sz
			
		; top-right
		dx_right_x0 = DX_PSF_X2- half_cropped_sz
		dx_right_x1 = DX_PSF_X2 + half_cropped_sz
		dx_right_y0 = DX_PSF_Y2 - half_cropped_sz
		dx_right_y1 = DX_PSF_Y2 + half_cropped_sz
		
	  ; start w/ WISPIT 2 values, but now we are full-frame.
	  endif else if stripe eq 'Full_Image' and (obj eq 'HIP17034') then begin
		print, 'HIP 17034 Full_Image correction enabled.', newline
		
		; in full-frmae dewarped coordinates.
		; "1" is for the "left" nod = NOD_A = dith1. DX is top, SX is bottom
		; (opposite of HII 1348)
		; below have been updated once for NOD_A (dith1) and nod00 only.
		; '2' is for the "right" nod = NOD_B = dith2
		DX_PSF_X1 = 481
		DX_PSF_Y1 = 1560
		
		SX_PSF_X1 = 502
		SX_PSF_Y1 = 557
		
		DX_PSF_X2 = 1394 + 20
		DX_PSF_Y2 = 1286 + 245
		
		SX_PSF_X2 = 1474 - 41
		SX_PSF_Y2 = 749 - 165
		
		; Pre-compute cropping indices
		; bottom-left
		sx_left_x0 = SX_PSF_X1- half_cropped_sz
		sx_left_x1 = SX_PSF_X1 + half_cropped_sz
		sx_left_y0 = SX_PSF_Y1 - half_cropped_sz
		sx_left_y1 = SX_PSF_Y1 + half_cropped_sz
		
		; top-left
		dx_left_x0 = DX_PSF_X1- half_cropped_sz
		dx_left_x1 = DX_PSF_X1 + half_cropped_sz
		dx_left_y0 = DX_PSF_Y1 - half_cropped_sz
		dx_left_y1 = DX_PSF_Y1 + half_cropped_sz
		
		; bottom-right
		sx_right_x0 = SX_PSF_X2- half_cropped_sz
		sx_right_x1 = SX_PSF_X2 + half_cropped_sz
		sx_right_y0 = SX_PSF_Y2 - half_cropped_sz
		sx_right_y1 = SX_PSF_Y2 + half_cropped_sz
		
		; top-right
		dx_right_x0 = DX_PSF_X2- half_cropped_sz
		dx_right_x1 = DX_PSF_X2 + half_cropped_sz
		dx_right_y0 = DX_PSF_Y2 - half_cropped_sz
		dx_right_y1 = DX_PSF_Y2 + half_cropped_sz
			
	  endif else if stripe eq 'Full_Image' and (obj eq 'HIP17900') then begin
		print, 'HIP 17900 Full_Image correction enabled.', newline
		
		; in full-frmae dewarped coordinates.
		; "1" is for the "left" nod = NOD_A = dith1. DX is top, SX is bottom
		; (opposite of HII 1348)
		; below have been updated once for NOD_A (dith1) and nod00 only.
		; '2' is for the "right" nod = NOD_B = dith2
		DX_PSF_X1 = 481+50+7
		DX_PSF_Y1 = 1560-5-2
		
		SX_PSF_X1 = 502-20
		SX_PSF_Y1 = 557-76-3
		
		DX_PSF_X2 = 1394 + 20-49+98+3
		DX_PSF_Y2 = 1286 + 245-10+5
		
		SX_PSF_X2 = 1474 - 41+22-45-3
		SX_PSF_Y2 = 749 - 165-85+2
		
		; Pre-compute cropping indices
		; bottom-left
		sx_left_x0 = SX_PSF_X1- half_cropped_sz
		sx_left_x1 = SX_PSF_X1 + half_cropped_sz
		sx_left_y0 = SX_PSF_Y1 - half_cropped_sz
		sx_left_y1 = SX_PSF_Y1 + half_cropped_sz
		
		; top-left
		dx_left_x0 = DX_PSF_X1- half_cropped_sz
		dx_left_x1 = DX_PSF_X1 + half_cropped_sz
		dx_left_y0 = DX_PSF_Y1 - half_cropped_sz
		dx_left_y1 = DX_PSF_Y1 + half_cropped_sz
		
		; bottom-right
		sx_right_x0 = SX_PSF_X2- half_cropped_sz
		sx_right_x1 = SX_PSF_X2 + half_cropped_sz
		sx_right_y0 = SX_PSF_Y2 - half_cropped_sz
		sx_right_y1 = SX_PSF_Y2 + half_cropped_sz
		
		; top-right
		dx_right_x0 = DX_PSF_X2- half_cropped_sz
		dx_right_x1 = DX_PSF_X2 + half_cropped_sz
		dx_right_y0 = DX_PSF_Y2 - half_cropped_sz
		dx_right_y1 = DX_PSF_Y2 + half_cropped_sz
			
	  endif else begin
		message, 'Invalid stripe.'
	  endelse

  ; Pre-compute smoothing kernels once
  use_lp_smooth = (do_smooth gt 0)
  use_hp_smooth = (hp_width gt 0)

  ; Get input cubes and sort them by nod sequence number
	if not keyword_set(pca_skysub) then begin
		if keyword_set(destripe_skysub) and (filled eq 0) then begin
			cube_files = file_search(dir + obj + '_NOD_?_nod??_grp??_skysub_skysub_destriped_cube.fits', count=n_cubes)
		endif else begin
			cube_files = file_search(dir + obj + '_NOD_?_nod??_grp??_skysub_cube.fits', count=n_cubes)
		endelse
	endif else begin ; PCA skysub cubes from Jupyter Notebook
	
		; destriping after PCA sky-subtraction
		if keyword_set(destripe_skysub) and (filled eq 0) then begin
			cube_files = file_search(dir + 'test_pca_skysub_cube_nod??_10comp_skysub_destriped_cube.fits', count=n_cubes)
		endif else if keyword_set(destripe_skysub) and (filled eq 1) then begin
			search_dir = dir + 'test_pca_skysub_cube_nod??_pca_skysub_filled_cube.fits'
			print, 'search_dir: ', search_dir
			cube_files = file_search(search_dir, count=n_cubes)
		endif else begin
			cube_files = file_search(dir + 'test_pca_skysub_cube_nod??_10comp.fits', count=n_cubes)
		endelse
	
	endelse
  	
  if n_cubes eq 0 then begin
    print, 'No sky-subtracted nod cubes found.'
    return
  endif

  ; Sort files by nod sequence number (nod00, nod01, nod02, etc.)
  ; hopefully this works for the PCA skysub cubes too...?
  nod_numbers = intarr(n_cubes)
  for i=0, n_cubes-1 do begin
    base = file_basename(cube_files[i])
    ; Extract nod number from filename (e.g., "nod05" -> 5)
    nod_pos = strpos(base, '_nod')
    if nod_pos ne -1 then begin
      nod_str = strmid(base, nod_pos+4, 2)
      nod_numbers[i] = fix(nod_str)
    endif else begin
      nod_numbers[i] = i  ; fallback
    endelse
  endfor
  
  ; Sort files by nod sequence number
  sort_idx = sort(nod_numbers)
  cube_files = cube_files[sort_idx]
  nod_numbers = nod_numbers[sort_idx]

  print, 'Found ', n_cubes, ' cubes to process (sorted by nod sequence)'
  
  ; Filter cubes based on nod_filter parameter
  if nod_filter ne 'both' then begin
		valid_indices = []
		for i=0, n_cubes-1 do begin
		  ; Determine NOD position based on sequence number
		  if (nod_numbers[i] mod 2) eq 0 then begin
			current_nod = 'dith1'
		  endif else begin
			current_nod = 'dith2'
		  endelse
		  
		  if current_nod eq nod_filter then begin
			valid_indices = [valid_indices, i]
		  endif
		endfor
		
		if n_elements(valid_indices) eq 0 then begin
		  print, 'No cubes match the specified nod_filter: ', nod_filter
		  return
		endif
		
		; Update arrays to only include filtered cubes
		cube_files = cube_files[valid_indices]
		nod_numbers = nod_numbers[valid_indices]
		n_cubes = n_elements(cube_files)
		
		print, 'Filtered to ', n_cubes, ' cubes matching nod_filter: ', nod_filter
  endif
  
  for i=0, n_cubes-1 do begin
    print, '  ', file_basename(cube_files[i]), ' (nod', string(nod_numbers[i], format='(I02)'), ')'
  endfor
  
  ; Pre-allocate padded array to avoid repeated allocation
  pad = fltarr(2048, 2048)
  
  ; Determine what needs to be processed
  print, 'Run: ', run
  process_sx = (run eq 1 or run eq 5)
  process_dx = (run eq 2 or run eq 5)
  print, 'Process_sx: ', process_sx
  print, 'Process_dx: ', process_dx

  for c=0, n_cubes-1 do begin
    path = cube_files[c]
    base = file_basename(path)
    
    ; FIXED: Determine NOD position based on sequence number
    ; Even nod numbers (00, 02, 04...) go to dith1 = NOD_A
    ; Odd nod numbers (01, 03, 05...) go to dith2 = NOD_B
    if (nod_numbers[c] mod 2) eq 0 then begin
      nod = 'dith1'
    endif else begin
      nod = 'dith2'
    endelse
    
    print, 'Processing cube ', c+1, '/', n_cubes, ': ', base, ' -> ', nod

    ; Read entire cube at once for better I/O performance
    ;t0 = systime(1)
    cube_data = readfits(path, header, /silent)
    ;print, 'Cube read time: ', systime(1)-t0, ' seconds'
    
    sz = size(cube_data)
    xdim = sz[1]
    ydim = sz[2] 
    frames = sz[3]

    ; Initialize output arrays only for what we need
    if process_sx then begin
      cube_sx = fltarr(crop_size, crop_size, frames)
    endif
    if process_dx then begin
      cube_dx = fltarr(crop_size, crop_size, frames)
    endif

    ; Process frames
    ;t0 = systime(1)
    for i=0, frames-1 do begin
      if (i mod 50) eq 0 then print, 'Frame ', i, '/', frames
      
      ; Extract frame and calculate median once
      data = cube_data[*,*,i]
      med = median(data, /even)
      
      ; Apply smoothing operations efficiently
      ; (do this before padding!)
      if use_lp_smooth then begin
        ; Handle NaNs efficiently
        bad_mask = ~finite(data)
        if total(bad_mask) gt 0 then data[where(bad_mask)] = med
        
        data = filter_image(data, FWHM=do_smooth, PSF=lp_PSF, /all_pixels)
        
        if use_hp_smooth then begin
          hp_filtered = filter_image(data, FWHM=hp_width, PSF=hp_PSF, /all_pixels)
          data -= hp_filtered
          data = filter_image(data, FWHM=do_smooth, PSF=lp_PSF, /all_pixels)
        endif
      endif
      
      ; Efficient padding operation
      pad[*] = med  ; Fill entire array with median
      if stripe eq 'center1024' then begin
        pad[0:xdim-1, 512:512+ydim-1] = data
      endif else if stripe eq 'second512' then begin
        pad[0:xdim-1, 1024:1024+ydim-1] = data
      endif else if stripe eq 'Full_Image' then begin
      	pad[0:xdim-1, 0:ydim-1] = data
      endif

      ; Process SX mirror data
      if process_sx then begin
        dew_sx = poly_2d(pad, Kx_sx, Ky_sx, 2, cubic=-1.0, missing=med)
        if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_sx_test.fits', dew_sx
        
        ; SX DITH 2
        if nod eq 'dith2' then begin
        
			if ((stripe eq 'center1024') and (obj eq 'tyc5709')) or ((stripe eq 'Full_Image') and ((obj eq 'HIP17034') or (obj eq 'HIP17900'))) then begin
				cropped = dew_sx[sx_right_x0:sx_right_x1, sx_right_y0:sx_right_y1]
			endif else begin
				cropped = dew_sx[sx_left_x0:sx_left_x1, sx_y0:sx_y1]
			endelse
        
        	; smooth again after dewarping, if set
        	if use_lp_smooth then begin
        		
        		bad_mask = ~finite(cropped)
        		if total(bad_mask) gt 0 then cropped[where(bad_mask)] = med
        		
				cube_sx[*,*,i] = filter_image(cropped, FWHM=do_smooth, PSF=lp_PSF,$
												/all_pixels)
			endif else begin
				cube_sx[*,*,i] = cropped
			endelse
			
						if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_sx2_cropped_test.fits', cube_sx[*,*,i]
			
        endif
        
        ; SX DITH 1
        if nod eq 'dith1' then begin
        	;print, 'working on SX dith 1'
        
			if ((stripe eq 'center1024') and (obj eq 'tyc5709')) or ((stripe eq 'Full_Image') and ((obj eq 'HIP17034') or (obj eq 'HIP17900'))) then begin
				;print, 'TUC 5709 cropping SX dith1 ...'
				cropped = dew_sx[sx_left_x0:sx_left_x1, sx_left_y0:sx_left_y1]
			endif else begin
				cropped = dew_sx[sx_right_x0:sx_right_x1, sx_y0:sx_y1]
			endelse
		
        	; smooth again after dewarping, if set
        	if use_lp_smooth then begin
        		
        		bad_mask = ~finite(cropped)
        		if total(bad_mask) gt 0 then cropped[where(bad_mask)] = med
        		
				cube_sx[*,*,i] = filter_image(cropped, FWHM=do_smooth, PSF=lp_PSF,$
												/all_pixels)
			endif else begin
				cube_sx[*,*,i] = cropped
			endelse
			
			if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_sx1_cropped_test.fits', cube_sx[*,*,i]
			
        endif
        
      endif

      ; Process DX mirror data
      if process_dx then begin
        dew_dx = poly_2d(pad, Kx_dx, Ky_dx, 2, cubic=-1.0, missing=med)
        if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_dx_test.fits', dew_dx
        
        ; SX DITH 2
        if nod eq 'dith2' then begin
        
			if ((stripe eq 'center1024') and (obj eq 'tyc5709')) or ((stripe eq 'Full_Image') and ((obj eq 'HIP17034') or (obj eq 'HIP17900'))) then begin
				cropped = dew_dx[dx_right_x0:dx_right_x1, dx_right_y0:dx_right_y1]
			endif else begin; Alcor
				cropped = dew_dx[dx_right_x0:dx_right_x1, dx_right_y0:dx_right_y1]
			endelse
        
        	; smooth again after dewarping, if set
        	if use_lp_smooth then begin
        		
        		bad_mask = ~finite(cropped)
        		if total(bad_mask) gt 0 then cropped[where(bad_mask)] = med
        		
				cube_dx[*,*,i] = filter_image(cropped, FWHM=do_smooth, PSF=lp_PSF,$
												/all_pixels)
			endif else begin
				cube_dx[*,*,i] = cropped
			endelse
			
			if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_dx2_cropped_test.fits', cube_dx[*,*,i]
			
        endif
        
        ; SX DITH 1
        if nod eq 'dith1' then begin
        	;print, 'working on SX dith 1'
        
			if ((stripe eq 'center1024') and (obj eq 'tyc5709')) or ((stripe eq 'Full_Image') and ((obj eq 'HIP17034') or (obj eq 'HIP17900'))) then begin
				;print, 'TUC 5709 cropping SX dith1 ...'
				cropped = dew_dx[dx_left_x0:dx_left_x1, dx_left_y0:dx_left_y1]
			endif else begin; Alcor
				cropped = dew_dx[dx_left_x0:dx_left_x1, dx_left_y0:dx_left_y1]
			endelse
		
        	; smooth again after dewarping, if set
        	if use_lp_smooth then begin
        		
        		bad_mask = ~finite(cropped)
        		if total(bad_mask) gt 0 then cropped[where(bad_mask)] = med
        		
				cube_dx[*,*,i] = filter_image(cropped, FWHM=do_smooth, PSF=lp_PSF,$
												/all_pixels)
			endif else begin
				cube_dx[*,*,i] = cropped
			endelse
			
			if debug eq 1 and i eq 0 then writefits, '~/Desktop/dewarped_dx1_cropped_test.fits', cube_dx[*,*,i]
			
        endif
        
      endif
      
    endfor
    
    ;print, 'Frame processing time: ', systime(1)-t0, ' seconds'

    ; Efficient file writing with timing - uses original cube index for numbering
    ;t0 = systime(1)
    original_cube_index = nod_numbers[c]  ; Use the original nod number for file naming
    
    if process_sx then begin
        outfile = dir + 'processed_left/' + nod + '/' + obj + '_nod' + string(original_cube_index, format='(I02)') + '_cube_sx.fits'
        print, 'Writing SX to: ', outfile
        writefits, outfile, cube_sx, /silent
    endif
    
    if process_dx then begin
        outfile = dir + 'processed_right/' + nod + '/' + obj + '_nod' + string(original_cube_index, format='(I02)') + '_cube_dx.fits'
        print, 'Writing DX to: ', outfile
        writefits, outfile, cube_dx, /silent
    endif
    
    ;print, 'File writing time: ', systime(1)-t0, ' seconds'
    print, 'Completed cube: ', base, ' in folder: ', nod
    print, '----------------------------------------'
  endfor

  print, 'All cubes processed successfully.'
end
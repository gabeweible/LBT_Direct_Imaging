; Efficient bad pixel correction for frames where bad pixels are already NaN
; Only applies fixpix to small clusters/isolated NaN pixels
; Large NaN regions remain as NaN

pro efficient_nan_correction, frame, fixed_frame, $
    max_cluster_size=max_cluster_size, npix=npix, weight=weight, silent=silent
    
    compile_opt idl2

    ; Set defaults
    if not keyword_set(max_cluster_size) then max_cluster_size = 50
    if not keyword_set(npix) then npix = 24
    if not keyword_set(weight) then weight = 0
    if not keyword_set(silent) then silent = 0
    
    ; Get image dimensions
    sz = size(frame)
    x_dim = sz[1]
    y_dim = sz[2]
    
    ; Create mask of NaN pixels (0 = NaN/bad, 1 = good)
    nan_mask = finite(frame)
    
    ; Find connected components of NaN pixels
    nan_regions = label_region(1 - nan_mask, /all_neighbors, /ulong)
    
    ; Get statistics on each NaN region
    h = histogram(nan_regions, min=1, reverse_indices=ri)
    n_regions = n_elements(h)
    
    if not keyword_set(silent) then $
        print, 'Found ', strtrim(n_regions,2), ' NaN regions'
    
    ; Create a working copy of the frame
    temp_frame = frame
    
    ; Process each region
    small_nan_count = 0
    large_nan_count = 0
    
    for i = 0, n_regions-1 do begin
        if ri[i+1] gt ri[i] then begin  ; check if region exists
            region_size = h[i]
            region_indices = ri[ri[i]:ri[i+1]-1]
            
            if region_size gt max_cluster_size then begin
                ; Large region - keep as NaN, but mark for exclusion from fixpix
                large_nan_count += region_size
                ; Set these to a special value that fixpix will ignore
                ; We'll restore them to NaN after fixpix
                temp_frame[region_indices] = -999999.0
            endif else begin
                ; Small region - leave as NaN for fixpix to handle
                small_nan_count += region_size
            endelse
        endif
    endfor
    
    if not keyword_set(silent) then begin
        print, 'Small NaN clusters/isolated pixels: ', strtrim(small_nan_count,2)
        print, 'Large NaN regions (kept as NaN): ', strtrim(large_nan_count,2)
    endif
    
    ; Apply fixpix only to the remaining NaN pixels (small clusters)
    if small_nan_count gt 0 then begin
        if not keyword_set(silent) then $
            print, 'Applying fixpix to small NaN regions...'
        
        dummy = replicate(!values.f_nan, x_dim, y_dim)
        fixpix, temp_frame, dummy, fixed_frame, npix=npix, $
                weight=weight, /nan, silent=silent
    endif else begin
        fixed_frame = temp_frame
        if not keyword_set(silent) then $
            print, 'No small NaN regions found for fixpix correction'
    endelse
    
    ; Restore large regions back to NaN
    if large_nan_count gt 0 then begin
        wlarge = where(fixed_frame eq -999999.0, count)
        if count gt 0 then fixed_frame[wlarge] = !values.f_nan
    endif
    
end

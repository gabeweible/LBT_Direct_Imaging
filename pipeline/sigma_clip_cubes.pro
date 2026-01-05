pro sigma_clip_cubes, obj_name, output_folder, cube_folder=cube_folder, $
                      sigma_clip=sigma_clip, input_suffix=input_suffix, $
                      output_suffix=output_suffix, debug=debug, $
                      overwrite=overwrite, cube_indices=cube_indices,$
                      frame_min=frame_min

; Procedure to apply sigma clipping to all NOD cube files
; This should be run before sky_sub.pro to preprocess all cubes
;
; INPUTS:
;   obj_name - Object name (string)
;   output_folder - Output directory path
;
; KEYWORDS:
;   cube_folder - Input cube directory (default: same as output_folder)
;   sigma_clip - Sigma threshold for clipping (default: 5.0)
;   input_suffix - Input file suffix pattern (default: '_corrected_cube.fits')
;   output_suffix - Output file suffix (default: '_sigma_clipped_cube.fits')
;   debug - Enable debug output (default: 0)
;   overwrite - Overwrite existing files (default: 0)
;   cube_indices - Array of cube indices to process (0-based). If not set, processes all cubes.

compile_opt idl2
newline = string(10B)

; mask values less than this threshold before local-median replacement
; (for very stubborn pixels...)
if not keyword_set(frame_min) then frame_min = -8.5

; Set defaults
if not keyword_set(cube_folder) then cube_folder = output_folder
if not keyword_set(sigma_clip) then sigma_clip = 3.5
if not keyword_set(input_suffix) then input_suffix = '_corrected_cube.fits'
if not keyword_set(output_suffix) then output_suffix = '_sigma_clipped_cube.fits'
if not keyword_set(debug) then debug = 0
if not keyword_set(overwrite) then overwrite = 0

print, 'Starting sigma clipping of NOD cubes...'
print, 'Sigma clipping threshold: ', sigma_clip
print, 'Input suffix: ', input_suffix
print, 'Output suffix: ', output_suffix
print, newline

; Find all input cube files
search_pattern = cube_folder + obj_name + '_NOD_?_nod??*' + input_suffix
cube_files = file_search(search_pattern, count=n_cube_files)

if n_cube_files eq 0 then begin
    print, 'ERROR: No cube files found with pattern: ', search_pattern
    return
endif

print, 'Found ', n_cube_files, ' total cube files'

; Sort files by nod number to ensure proper order
nod_numbers = intarr(n_cube_files)
for i = 0, n_cube_files-1 do begin
    nod_numbers[i] = extract_nod_number(cube_files[i])
endfor

sort_indices = sort(nod_numbers)
cube_files = cube_files[sort_indices]

; Handle cube_indices selection
if keyword_set(cube_indices) then begin
    ; Validate cube_indices
    cube_indices = long(cube_indices)  ; Ensure long integers
    valid_indices = where(cube_indices ge 0 and cube_indices lt n_cube_files, n_valid)
    
    if n_valid eq 0 then begin
        print, 'ERROR: No valid cube indices specified. Valid range is 0 to ', n_cube_files-1
        return
    endif
    
    if n_valid lt n_elements(cube_indices) then begin
        invalid_indices = where(cube_indices lt 0 or cube_indices ge n_cube_files, n_invalid)
        print, 'WARNING: ', n_invalid, ' invalid cube indices ignored (valid range: 0-', n_cube_files-1, ')'
        cube_indices = cube_indices[valid_indices]
    endif
    
    ; Select only the specified cubes
    selected_cube_files = cube_files[cube_indices]
    selected_cube_indices = cube_indices
    n_process_cubes = n_elements(selected_cube_files)
    
    print, 'Processing ', n_process_cubes, ' selected cubes (indices: ', $
           strjoin(string(cube_indices, format='(I0)'), ', '), ')'
endif else begin
    ; Process all cubes
    selected_cube_files = cube_files
    selected_cube_indices = lindgen(n_cube_files)
    n_process_cubes = n_cube_files
    
    print, 'Processing all ', n_process_cubes, ' cube files'
endelse

print, 'Files to be processed:'
for i = 0, n_process_cubes-1 do begin
    cube_idx = selected_cube_indices[i]
    print, '  [', cube_idx, '] ', file_basename(selected_cube_files[i])
endfor
print, newline

; Process each selected cube file
for i = 0, n_process_cubes-1 do begin
    input_file = selected_cube_files[i]
    cube_idx = selected_cube_indices[i]
    
    ; Generate output filename
    base_name = file_basename(input_file, input_suffix)
    output_file = output_folder + base_name + output_suffix
    
    ; Check if output already exists
    if file_test(output_file) and (overwrite eq 0) then begin
        print, 'Skipping (already exists): ', file_basename(output_file)
        continue
    endif
    
    print, 'Processing cube [', cube_idx, '] ', i+1, '/', n_process_cubes, ': ', file_basename(input_file)
    
    ; Load the cube
    print, '  Loading cube...'
    obj_cube = readfits_fast(input_file)
    
    ; Get cube dimensions
    cube_size = size(obj_cube, /dimensions)
    x_dim = cube_size[0]
    y_dim = cube_size[1]
    total_frames = cube_size[2]
    
    print, '  Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', total_frames
    print, '  Applying sigma clipping to all frames...'
    
    ; Track statistics
    total_pixels_processed = 0L
    total_pixels_clipped = 0L
    
    ; Process each frame
    for frame_idx = 0, total_frames-1 do begin
        ; Extract frame
        frame_ii = obj_cube[*,*,frame_idx]
        original_frame = frame_ii  ; Keep copy for statistics
        
        ; replace low pixel values with local median (see next step...)
        frame_ii[where(frame_ii lt frame_min)] = !values.f_nan
        
        ; Handle NaNs first - replace with median-filtered values
        nan_mask = where(finite(frame_ii) ne 1, n_nans)
        if n_nans gt 0 then begin
            frame_ii_med = filter_image(frame_ii, median=5., /all)
            frame_ii[nan_mask] = frame_ii_med[nan_mask]
            
            ; *** NEW: Add Gaussian noise to replaced pixels ***
            ; Robust noise estimation using median and MAD on original frame
            ; Use only finite values from original frame for noise estimation
            good_pixels = where(finite(original_frame), n_good)
            if n_good gt 100 then begin  ; Need sufficient pixels for robust estimate
                good_values = original_frame[good_pixels]
                
                ; Calculate robust noise estimate: 1.4826 * MAD
                median_val = median(good_values, /even)
                mad_val = median(abs(good_values - median_val), /even)
                noise_sigma = 1.4826 * mad_val
                
                ; Generate Gaussian noise for replaced pixels
                if n_nans gt 0 then begin
                    ; Generate random Gaussian noise with estimated sigma
                    gaussian_noise = randomn(seed, n_nans) * noise_sigma
                    
                    ; Add noise to the replaced pixels
                    frame_ii[nan_mask] += gaussian_noise
                endif
                
                ; Debug output for noise estimation (first frame only)
                if frame_idx eq 0 and debug eq 1 then begin
                    print, '    Noise sigma estimate: ', noise_sigma
                    print, '    Number of pixels replaced with median+noise: ', n_nans
                endif
            endif else begin
                if frame_idx eq 0 and debug eq 1 then print, '    Warning: Insufficient good pixels for noise estimation'
            endelse
        endif
        
        ; Apply progressive sigma filtering
        ; First pass: 3x3 filter
        frame_ii_clipped = sigma_filter(frame_ii, 3, N_sigma=sigma_clip, /all, /iter)
        n_clipped_3x3 = total(finite(frame_ii) and ~finite(frame_ii_clipped))
        frame_ii = frame_ii_clipped
        
        ; Second pass: 5x5 filter  
        frame_ii_clipped = sigma_filter(frame_ii, 5, N_sigma=sigma_clip, /all, /iter)
        n_clipped_5x5 = total(finite(frame_ii) and ~finite(frame_ii_clipped))
        frame_ii = frame_ii_clipped
        
        ; Third pass: 7x7 filter
        frame_ii_clipped = sigma_filter(frame_ii, 7, N_sigma=sigma_clip, /all, /iter)
        n_clipped_7x7 = total(finite(frame_ii) and ~finite(frame_ii_clipped))
        
        ; Store the final clipped frame
        obj_cube[*,*,frame_idx] = frame_ii_clipped
        
        ; Update statistics
        frame_pixels = x_dim * y_dim
        frame_clipped = n_clipped_3x3 + n_clipped_5x5 + n_clipped_7x7
        total_pixels_processed += frame_pixels
        total_pixels_clipped += frame_clipped
        
        ; Debug output for first frame
        if (frame_idx eq 0) and (debug eq 1) then begin
            debug_file = output_folder + base_name + '_sigma_clip_debug_frame0.fits'
            writefits, debug_file, frame_ii_clipped
            print, '    Debug: First frame saved to ', file_basename(debug_file)
            print, '    Debug: NaNs replaced: ', n_nans
            print, '    Debug: Clipped in 3x3: ', n_clipped_3x3
            print, '    Debug: Clipped in 5x5: ', n_clipped_5x5  
            print, '    Debug: Clipped in 7x7: ', n_clipped_7x7
        endif
        
        ; Progress reporting
        if (frame_idx mod 100) eq 0 then begin
            percent_complete = 100.0 * frame_idx / (total_frames-1)
            print, '    Progress: ', string(percent_complete, format='(F5.1)'), '% (frame index', $
                   frame_idx, '/', total_frames-1, ')'
        endif
    endfor
    
    ; Report statistics for this cube
    clip_percentage = 100.0 * total_pixels_clipped / total_pixels_processed
    print, '  Statistics:'
    print, '    Total pixels processed: ', total_pixels_processed
    print, '    Total pixels clipped: ', total_pixels_clipped, $
           ' (', string(clip_percentage, format='(F5.2)'), '%)'
    
    ; Save the sigma-clipped cube
    print, '  Saving sigma-clipped cube: ', file_basename(output_file)
    writefits, output_file, obj_cube
    
    ; Free memory
    delvar, obj_cube
    
    print, '  Cube processing complete.', newline
endfor

print, 'Sigma clipping complete for selected cubes!'
print, 'Processed ', n_process_cubes, ' cube files'

; Create a summary file
summary_file = output_folder + obj_name + '_sigma_clipped_cubes_list.txt'
openw, lun, summary_file, /get_lun
printf, lun, '# List of sigma-clipped cube files for object: ' + obj_name
printf, lun, '# Generated on: ' + systime()
printf, lun, '# Sigma clipping threshold: ' + string(sigma_clip)
if keyword_set(cube_indices) then begin
    printf, lun, '# Selected cube indices: ' + strjoin(string(cube_indices, format='(I0)'), ', ')
    printf, lun, '# Cubes processed: ' + string(n_process_cubes) + ' of ' + string(n_cube_files) + ' total'
endif else begin
    printf, lun, '# Total cubes processed: ' + string(n_process_cubes)
endelse
printf, lun, '#'
printf, lun, '# Format: cube_file'

for i = 0, n_process_cubes-1 do begin
    base_name = file_basename(selected_cube_files[i], input_suffix)
    output_name = base_name + output_suffix
    printf, lun, output_name
endfor

free_lun, lun
print, 'Summary file created: ', file_basename(summary_file)

end
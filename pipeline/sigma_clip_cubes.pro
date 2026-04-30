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
;   cube_indices - Array of nod NUMBERS to process (i.e. the XX in 'nodXX' as
;                  it appears in the filename, NOT positional indices into the
;                  file list). If not set, all cubes are processed. Requesting
;                  a nod number whose file does not exist produces a warning
;                  and that nod is skipped.

compile_opt idl2
newline = string(10B)

; mask values less than this threshold before local-median replacement
; (for very stubborn pixels...)
if not keyword_set(frame_min) then frame_min = -40

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

; Extract nod numbers from each filename, then sort files & nod_numbers together
nod_numbers = intarr(n_cube_files)
for i = 0, n_cube_files-1 do begin
    nod_numbers[i] = extract_nod_number(cube_files[i])
endfor

sort_indices = sort(nod_numbers)
cube_files   = cube_files[sort_indices]
nod_numbers  = nod_numbers[sort_indices]   ; keep aligned with cube_files

; Handle cube_indices selection
; cube_indices refers to ACTUAL nod numbers (e.g. the 02 in 'nod02'),
; not to positions in the file list. This guarantees correct selection
; even when some nod cubes are missing.
if keyword_set(cube_indices) then begin
    cube_indices = long(cube_indices)
    n_requested  = n_elements(cube_indices)

    ; For each requested nod number, find its position in cube_files (-1 if missing)
    matched_file_positions = lonarr(n_requested) - 1L
    for j = 0, n_requested-1 do begin
        match_pos = where(nod_numbers eq cube_indices[j], n_match)
        if n_match gt 0 then matched_file_positions[j] = match_pos[0]
    endfor

    ; Split into found vs. missing
    found = where(matched_file_positions ge 0, n_found, $
                  complement=missing, ncomplement=n_missing)

    if n_missing gt 0 then begin
        missing_nods = cube_indices[missing]
        print, 'WARNING: ', n_missing, $
               ' requested nod number(s) have no matching file and will be skipped: ', $
               strjoin(string(missing_nods, format='(I02)'), ', ')
        print, '         Available nod numbers in folder: ', $
               strjoin(string(nod_numbers,  format='(I02)'), ', ')
    endif

    if n_found eq 0 then begin
        print, 'ERROR: None of the requested nod numbers were found among existing files.'
        return
    endif

    selected_cube_files   = cube_files[matched_file_positions[found]]
    selected_cube_indices = cube_indices[found]   ; actual nod numbers
    n_process_cubes       = n_found

    print, 'Processing ', n_process_cubes, ' selected cubes (nod numbers: ', $
           strjoin(string(selected_cube_indices, format='(I02)'), ', '), ')'
endif else begin
    ; Process all cubes
    selected_cube_files   = cube_files
    selected_cube_indices = nod_numbers          ; actual nod numbers, not positions
    n_process_cubes       = n_cube_files

    print, 'Processing all ', n_process_cubes, ' cube files'
endelse

print, 'Files to be processed:'
for i = 0, n_process_cubes-1 do begin
    cube_idx = selected_cube_indices[i]
    print, '  [nod ', string(cube_idx, format='(I02)'), '] ', $
           file_basename(selected_cube_files[i])
endfor
print, newline

; Process each selected cube file
for i = 0, n_process_cubes-1 do begin
    input_file = selected_cube_files[i]
    cube_idx   = selected_cube_indices[i]   ; this is now the real nod number

    ; Generate output filename
    base_name   = file_basename(input_file, input_suffix)
    output_file = output_folder + base_name + output_suffix

    ; Check if output already exists
    if file_test(output_file) and (overwrite eq 0) then begin
        print, 'Skipping (already exists): ', file_basename(output_file)
        continue
    endif

    print, 'Processing cube [nod ', string(cube_idx, format='(I02)'), '] ', $
           i+1, '/', n_process_cubes, ': ', file_basename(input_file)

    ; Load the cube
    print, '  Loading cube...'
    obj_cube = readfits_fast(input_file)

    ; Get cube dimensions
    cube_size    = size(obj_cube, /dimensions)
    x_dim        = cube_size[0]
    y_dim        = cube_size[1]
    total_frames = cube_size[2]

    print, '  Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', total_frames
    print, '  Applying sigma clipping to all frames...'

    ; Track statistics
    total_pixels_processed = 0L
    total_pixels_clipped   = 0L

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
                median_val = median(double(good_values), /even, /double)
                mad_val = float(median(abs(double(good_values) - median_val), /even, /double))
                median_val = float(median_val)

                ; double the noise so interpolated pixels are further downweighted!!
                noise_sigma = 1.4826 * mad_val

                ; Generate Gaussian noise for replaced pixels
                seed = 12345L
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
        frame_ii = frame_ii_clipped

        ; Second pass: 5x5 filter
        frame_ii_clipped = sigma_filter(frame_ii, 5, N_sigma=sigma_clip, /all, /iter)
        frame_ii = frame_ii_clipped

        ; Third pass: 7x7 filter
        frame_ii_clipped = sigma_filter(frame_ii, 7, N_sigma=sigma_clip, /all, /iter)

        ; Store the final clipped frame
        obj_cube[*,*,frame_idx] = frame_ii_clipped

        ; Update statistics
        frame_pixels = x_dim * y_dim
        total_pixels_processed += frame_pixels

        ; Progress reporting
        if (frame_idx mod 41) eq 0 then begin
            percent_complete = 100.0 * frame_idx / (total_frames-1)
            print, '    Progress: ', string(percent_complete, format='(F5.1)'), '% (frame index', $
                   frame_idx, '/', total_frames-1, ')'
        endif
    endfor

    ; Report statistics for this cube
    print, '  Statistics:'
    print, '    Total pixels processed: ', total_pixels_processed

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
    printf, lun, '# Requested nod numbers: ' + $
            strjoin(string(cube_indices,         format='(I02)'), ', ')
    printf, lun, '# Processed nod numbers: ' + $
            strjoin(string(selected_cube_indices, format='(I02)'), ', ')
    printf, lun, '# Cubes processed: ' + string(n_process_cubes) + $
                 ' of ' + string(n_cube_files) + ' total available'
endif else begin
    printf, lun, '# Total cubes processed: ' + string(n_process_cubes)
endelse
printf, lun, '#'
printf, lun, '# Format: cube_file'

for i = 0, n_process_cubes-1 do begin
    base_name   = file_basename(selected_cube_files[i], input_suffix)
    output_name = base_name + output_suffix
    printf, lun, output_name
endfor

free_lun, lun
print, 'Summary file created: ', file_basename(summary_file)

end
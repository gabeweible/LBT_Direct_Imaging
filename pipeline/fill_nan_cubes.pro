pro fill_nan_cubes, obj_name, output_folder, cube_folder=cube_folder, $
                      input_suffix=input_suffix, output_suffix=output_suffix, $
                      debug=debug, overwrite=overwrite, cube_indices=cube_indices,$
                      input_prefix=input_prefix

; Procedure to replace NaN values in sigma-clipped cube files with frame median + Gaussian noise
; This should be run after sigma_clip_cubes.pro to ensure all pixels have finite values
;
; INPUTS:
;   obj_name - Object name (string)
;   output_folder - Output directory path
;
; KEYWORDS:
;   cube_folder - Input cube directory (default: same as output_folder)
;   input_suffix - Input file suffix pattern (default: '_sigma_clipped_cube.fits')
;   output_suffix - Output file suffix (default: '_filled_cube.fits')
;   debug - Enable debug output (default: 0)
;   overwrite - Overwrite existing files (default: 0)
;   cube_indices - Array of nod NUMBERS to process (i.e. the XX in 'nodXX' as
;                  it appears in the filename, NOT positional indices into the
;                  file list). If not set, all cubes are processed. Requesting
;                  a nod number whose file does not exist produces a warning
;                  and that nod is skipped.

compile_opt idl2, logical_predicate, strictarrsubs
newline = string(10B)

; Set defaults
if ~ keyword_set(cube_folder) then cube_folder = output_folder
if ~keyword_set(input_suffix) then input_suffix = '_sigma_clipped_cube.fits'
if ~keyword_set(output_suffix) then output_suffix = '_filled_cube.fits'
if ~keyword_set(debug) then debug = 0
if ~keyword_set(overwrite) then overwrite = 0
if ~keyword_set(input_prefix) then input_prefix = obj_name + '_NOD_?'

print, 'Starting NaN filling of sigma-clipped cubes...'
print, 'Input suffix: ', input_suffix
print, 'Output suffix: ', output_suffix
print, newline

; Find input cube files. select_nod_cubes handles globbing, nod-number
; extraction, sorting, and (optional) cube_indices filtering. Selection by
; cube_indices matches against the actual nod number in each filename, not
; against a positional index, so missing nods don't shift the selection.
search_pattern = cube_folder + input_prefix + '_nod??*' + input_suffix
selected_cube_files = select_nod_cubes(search_pattern, $
                                       cube_indices=cube_indices, $
                                       nod_numbers=selected_cube_indices, $
                                       count=n_process_cubes, $
                                       /verbose)

if n_process_cubes eq 0 then begin
    print, 'ERROR: No cube files matched (pattern + cube_indices filter).'
    return
endif

print, 'Files to be processed:'
for i = 0, n_process_cubes-1 do begin
    cube_idx = selected_cube_indices[i]
    print, '  [nod ', string(cube_idx, format='(I02)'), '] ', $
           file_basename(selected_cube_files[i])
endfor
print, newline

; Process each selected cube file. seed is intentionally left undefined:
; randomn() will initialize it from the system clock on its first call and
; advance the same sequence across all cubes/frames, giving non-reproducible
; but uncorrelated noise. (Previously seed was reset to 12345L at the start
; of each cube, which made the same noise pattern repeat from cube to cube.)
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
    obj_cube = double(readfits_fast(input_file))

    ; Get cube dimensions
    cube_size    = size(obj_cube, /dimensions)
    x_dim        = cube_size[0]
    y_dim        = cube_size[1]
    total_frames = cube_size[2]

    print, '  Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', total_frames
    print, '  Filling NaN values with frame median + noise...'

    ; Track statistics
    total_pixels_processed = 0L
    total_nans_filled      = 0L

    ; Process each frame
    for frame_idx = 0, total_frames-1 do begin
        ; Extract frame
        frame_ii = obj_cube[*,*,frame_idx]

        ; Find NaN pixels
        nan_mask = where(~finite(frame_ii), n_nans)

        if n_nans gt 0 then begin
            ; Find finite pixels for statistics
            finite_mask = where(finite(frame_ii), n_finite)

            if n_finite gt 10 then begin  ; Need sufficient pixels for robust statistics
                finite_values = frame_ii[finite_mask]

                ; Calculate frame median
                frame_median = median(double(finite_values), /even, /double)

                ; Calculate robust noise estimate using MAD
                mad_val     = median(abs(double(finite_values) - frame_median), /even, /double)
                noise_sigma = 1.4826 * mad_val

                ; Generate Gaussian noise for NaN pixels
                gaussian_noise = randomn(seed, n_nans, /double) * noise_sigma

                ; Replace NaN values with median + noise
                frame_ii[nan_mask] = frame_median + gaussian_noise

                ; Store the filled frame back to cube
                obj_cube[*,*,frame_idx] = frame_ii

                ; Update statistics
                total_nans_filled += n_nans

                ; Debug output for first frame
                if (frame_idx eq 0) and (debug eq 1) then begin
                    debug_file = output_folder + base_name + '_nan_fill_debug_frame0.fits'
                    writefits, debug_file, frame_ii
                    print, '    Debug: First frame saved to ', file_basename(debug_file)
                    print, '    Debug: Frame median: ', frame_median
                    print, '    Debug: Noise sigma (MAD-based): ', noise_sigma
                    print, '    Debug: NaNs filled: ', n_nans, ' of ', x_dim*y_dim, ' pixels'

                    ; Additional statistics for debug
                    filled_values = frame_median + gaussian_noise
                    print, '    Debug: Filled values range: [', min(filled_values), ', ', max(filled_values), ']'
                endif

            endif else begin
                ; Not enough finite pixels - fill with zeros and warn
                print, '    Warning: Frame ', frame_idx, ' has insufficient finite pixels (', n_finite, ')'
                print, '             Filling NaNs with zeros instead of median+noise'
                frame_ii[nan_mask] = 0.0
                obj_cube[*,*,frame_idx] = frame_ii
                total_nans_filled += n_nans
            endelse

        endif ; if n_nans gt 0

        ; Update total pixels processed
        total_pixels_processed += (x_dim * y_dim)

        ; Progress reporting
        if (frame_idx mod 100) eq 0 then begin
            percent_complete = 100.0 * frame_idx / (total_frames-1)
            print, '    Progress: ', string(percent_complete, format='(F5.1)'), '% (frame ', $
                   frame_idx, '/', total_frames-1, ')'
        endif
    endfor

    ; Final validation - check for any remaining NaNs
    remaining_nans = total(~finite(obj_cube))
    if remaining_nans gt 0 then begin
        print, '  WARNING: ', remaining_nans, ' NaN values remain in cube after processing!'
    endif else begin
        print, '  SUCCESS: All NaN values have been replaced with finite values'
    endelse

    ; Report statistics for this cube
    fill_percentage = 100.0 * total_nans_filled / total_pixels_processed
    print, '  Statistics:'
    print, '    Total pixels processed: ', total_pixels_processed
    print, '    Total NaNs filled: ', total_nans_filled, $
           ' (', string(fill_percentage, format='(F5.2)'), '%)'

    ; Save the filled cube
    print, '  Saving filled cube: ', file_basename(output_file)
    writefits, output_file, float(obj_cube)

    ; Free memory
    delvar, obj_cube

    print, '  Cube processing complete.', newline
endfor

print, 'NaN filling complete for selected cubes!'
print, 'Processed ', n_process_cubes, ' cube files'

; Create a summary file
summary_file = output_folder + obj_name + '_filled_cubes_list.txt'
openw, lun, summary_file, /get_lun
printf, lun, '# List of NaN-filled cube files for object: ' + obj_name
printf, lun, '# Generated on: ' + systime()
printf, lun, '# Input files were sigma-clipped cubes'
if keyword_set(cube_indices) then begin
    printf, lun, '# Requested nod numbers: ' + $
            strjoin(string(cube_indices,         format='(I02)'), ', ')
    printf, lun, '# Processed nod numbers: ' + $
            strjoin(string(selected_cube_indices, format='(I02)'), ', ')
    printf, lun, '# Cubes processed: ' + string(n_process_cubes)
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
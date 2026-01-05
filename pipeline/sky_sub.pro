pro sky_sub, obj_name, coadd, output_folder, cube_folder=cube_folder,$
	fs_start=fs_start, fs_end=fs_end, fpn=fpn,$
	pre_bad_px=pre_bad_px, nnod=nnod, debug=debug,$
	start_cube=start_cube, cube_indices=cube_indices

; Optimized version - sky subtraction only with mean combination
compile_opt idl2
newline = string(10B)

; Set defaults
if not keyword_set(pre_bad_px) then pre_bad_px = 0
if not keyword_set(nnod) then nnod = 50
if keyword_set(cube_folder) then cube_folder = cube_folder else cube_folder = output_folder
if not keyword_set(debug) then debug = 0
if not keyword_set(start_cube) then start_cube = 0

print, 'Starting optimized sky subtraction with mean combination only...'
print, 'Processing individual NOD cubes sequentially...'

; Print processing parameters
print, newline
print, 'Using mean sky combination'
print, nnod, ' nearest frames for each sky model'
print, fpn, ' frames per nod before co-addition'
print, 'coadd = ', coadd
print, 'means ', fpn/coadd, ' co-added frames per nod'
print, newline

; Find all corrected cube files
if pre_bad_px eq 0 then begin
	cube_files = file_search(cube_folder + obj_name + '_NOD_?_nod??_grp??_sigma_clipped_cube.fits', count=n_cube_files)
endif else begin
	cube_files = file_search(cube_folder + obj_name + '_NOD_?_nod??_cube.fits', count=n_cube_files)
endelse

if n_cube_files eq 0 then begin
	print, 'ERROR: No NOD cube files found!'
	return
endif

print, 'Found ', n_cube_files, ' NOD cube files'

; Sort files by nod number
nod_numbers = intarr(n_cube_files)
for i = 0, n_cube_files-1 do begin
    nod_numbers[i] = extract_nod_number(cube_files[i])
endfor

sort_indices = sort(nod_numbers)
cube_files = cube_files[sort_indices]
nod_numbers = nod_numbers[sort_indices]

print, 'Files sorted by nod number'

; Separate NOD_A and NOD_B files with their indices
nod_a_files = []
nod_b_files = []
nod_a_indices = []
nod_b_indices = []

for i = 0, n_cube_files-1 do begin
	if strpos(cube_files[i], '_NOD_A_') ne -1 then begin
		nod_a_files = [nod_a_files, cube_files[i]]
		nod_a_indices = [nod_a_indices, i]
	endif else if strpos(cube_files[i], '_NOD_B_') ne -1 then begin
		nod_b_files = [nod_b_files, cube_files[i]]
		nod_b_indices = [nod_b_indices, i]
	endif
endfor

n_nod_a = n_elements(nod_a_files)
n_nod_b = n_elements(nod_b_files)

print, 'Found ', n_nod_a, ' NOD_A files and ', n_nod_b, ' NOD_B files'

if (n_nod_a eq 0) or (n_nod_b eq 0) then begin
  print, 'ERROR: Missing one or both nodding positions'
  return
endif

; Read first cube to get dimensions
first_cube = readfits_fast(cube_files[0])
cube_size = size(first_cube, /dimensions)
; After reading the first cube
print, 'First input cube dimensions: ', cube_size

x_dim = cube_size[0]
y_dim = cube_size[1]
delvar, first_cube

;print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' y per cube'

; Determine which cubes to process
if keyword_set(cube_indices) then begin
	; Validate cube_indices
	valid_indices = where(cube_indices ge 0 and cube_indices lt n_cube_files, n_valid)
	if n_valid eq 0 then begin
		print, 'ERROR: No valid cube indices specified!'
		print, 'Valid range is 0 to ', n_cube_files-1
		return
	endif
	if n_valid lt n_elements(cube_indices) then begin
		print, 'WARNING: Some cube indices are out of range and will be ignored'
		print, 'Valid range is 0 to ', n_cube_files-1
	endif
	process_indices = cube_indices[valid_indices]
	print, 'Processing only specified cube indices: ', process_indices
endif else begin
	; Process all cubes from start_cube to end
	process_indices = indgen(n_cube_files - start_cube) + start_cube
	print, 'Processing cubes from index ', start_cube, ' to ', n_cube_files-1
endelse

n_process = n_elements(process_indices)
print, 'Total cubes to process: ', n_process

; Process each specified cube
for proc_idx = 0, n_process-1 do begin
	cube_idx = process_indices[proc_idx]
	print, newline, 'Processing cube ', proc_idx+1, '/', n_process, ' (index ', cube_idx, '): ', file_basename(cube_files[cube_idx])
	
	; Determine if current cube is NOD_A or NOD_B
	current_is_nod_a = strpos(cube_files[cube_idx], '_NOD_A_') ne -1
	
	; Load current cube
	print, '  Loading current cube...'
	current_cube = readfits_fast(cube_files[cube_idx])
	
	frames_per_cube = (size(current_cube))[3]
	
	; Find which cubes to use for sky computation
	sky_cube_indices = find_sky_cubes(cube_idx, current_is_nod_a, nod_a_indices, nod_b_indices, $
	                                  n_cube_files, max_sky_cubes=2)
	
	print, '  Using cubes for sky: ', sky_cube_indices
	
	; Load sky cubes (up to 2 additional cubes)
	n_sky_cubes = n_elements(sky_cube_indices)
	sky_cubes = ptrarr(n_sky_cubes)
	
	for sky_idx = 0, n_sky_cubes-1 do begin
		print, '    Loading sky cube: ', file_basename(cube_files[sky_cube_indices[sky_idx]])
		sky_cubes[sky_idx] = ptr_new(readfits_fast(cube_files[sky_cube_indices[sky_idx]]))
	endfor
	
	; Initialize sky-subtracted cube
	sky_sub_cube = dblarr(x_dim, y_dim, frames_per_cube)
	
	; Process each frame in the current cube
	print, '  Processing frames for sky subtraction...'
	for frame_idx = 0, frames_per_cube-1 do begin
		; Initialize arrays for sky frame collection
		sky_frame_sum = dblarr(x_dim, y_dim)
		sky_frame_count = 0L
		
		; Collect sky frames from all sky cubes
		for sky_cube_idx = 0, n_sky_cubes-1 do begin
			sky_cube_data = *sky_cubes[sky_cube_idx]
			n_sky_frames = (size(sky_cube_data, /dimensions))[2]
			
			; Calculate distances from current frame to all frames in this sky cube
			frame_distances = abs(indgen(n_sky_frames) - frame_idx)
			sorted_frame_indices = sort(frame_distances)
			
			; Take up to nnod nearest frames from this sky cube
			n_frames_to_use = (n_sky_frames lt nnod) ? n_sky_frames : nnod
			
			; Add selected frames to running sum (vectorized)
			for sel_idx = 0, n_frames_to_use-1 do begin
				sky_frame_sum += sky_cube_data[*,*,sorted_frame_indices[sel_idx]]
				sky_frame_count++
			endfor
		endfor
		
		; Compute mean sky and subtract
		if sky_frame_count gt 0 then begin
			sky = sky_frame_sum / sky_frame_count
			sky_sub_cube[*,*,frame_idx] = current_cube[*,*,frame_idx] - sky
		endif else begin
			print, 'WARNING: No sky frames available for frame ', frame_idx
			sky_sub_cube[*,*,frame_idx] = current_cube[*,*,frame_idx]
		endelse
		
		; Debug output for first frame of first processed cube
		if proc_idx eq 0 and frame_idx eq 0 and debug eq 1 then begin
			writefits, output_folder + 'skysub_cube_' + strtrim(cube_idx,2) + '_frame_0.fits', sky_sub_cube[*,*,frame_idx]
			writefits, output_folder + 'skysub_first_sky.fits', sky
		endif
		
		; Report progress periodically
		if (frame_idx mod 20) eq 0 then begin
			print, '    Processed frame index', frame_idx, ' of ', frames_per_cube-1
		endif
	endfor
	
	; Free memory used by sky cubes
	for sky_idx = 0, n_sky_cubes-1 do begin
		ptr_free, sky_cubes[sky_idx]
	endfor
	
	; Set median of all frames to zero (vectorized where possible)
	print, '  Setting median of all frames to zero...'
	for frame_idx = 0, frames_per_cube-1 do begin
		frame_med = median(sky_sub_cube[*,*,frame_idx], /double, /even)
		sky_sub_cube[*,*,frame_idx] -= frame_med
	endfor
	
	; Before saving each cube
	print, 'Output cube dimensions: ', size(sky_sub_cube, /dimensions)
	
	; Save sky-subtracted cube
	base_name = file_basename(cube_files[cube_idx])
	if pre_bad_px eq 0 then begin
		output_name = str_replace(base_name, '_sigma_clipped_cube.fits', '_skysub_cube.fits')
	endif else begin
		output_name = str_replace(base_name, '_cube.fits', '_skysub_cube.fits')
	endelse
	
	output_file = output_folder + output_name
	print, '  Saving sky-subtracted cube: ', file_basename(output_file)
	writefits, output_file, sky_sub_cube
	
	; Free memory
	delvar, current_cube, sky_sub_cube
	
	print, '  Cube processing complete.'
endfor

print, 'Sky subtraction processing complete!'
print, 'Processed ', n_process, ' individual NOD cubes.'

end
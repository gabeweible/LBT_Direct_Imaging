function find_sky_cubes, current_cube_idx, current_is_nod_a, nod_a_indices, nod_b_indices, $
                        n_total_cubes, max_sky_cubes=max_sky_cubes

; Function to find nearest opposite-nod cubes for sky computation
; 
; INPUTS:
;   current_cube_idx - Index of the current cube being processed
;   current_is_nod_a - Boolean: 1 if current cube is NOD_A, 0 if NOD_B
;   nod_a_indices - Array of indices for all NOD_A cubes
;   nod_b_indices - Array of indices for all NOD_B cubes
;   n_total_cubes - Total number of cubes
;
; KEYWORDS:
;   max_sky_cubes - Maximum number of sky cubes to return (default: 2)
;
; RETURNS:
;   Array of cube indices to use for sky computation

compile_opt idl2

; Set default for maximum sky cubes
if not keyword_set(max_sky_cubes) then max_sky_cubes = 2

; Determine which cubes to use for sky based on nodding position
if current_is_nod_a then begin
    ; Use NOD_B cubes for sky
    target_indices = nod_b_indices
endif else begin
    ; Use NOD_A cubes for sky
    target_indices = nod_a_indices
endelse

; Check if we have any target cubes
if n_elements(target_indices) eq 0 then begin
    print, 'ERROR: No target cubes available for sky computation'
    return, -1
endif

; Calculate distances to all target cubes
distances = abs(target_indices - current_cube_idx)
sorted_indices = sort(distances)

; Determine how many sky cubes to use based on position
; Strategy:
; - First cube: use 1 nearest cube (the one immediately following)
; - Last cube: use 1 nearest cube (the one immediately preceding)  
; - Middle cubes: use up to max_sky_cubes nearest cubes
; - But never exceed the number of available target cubes

n_available = n_elements(target_indices)

if (current_cube_idx eq 0) or (current_cube_idx eq n_total_cubes-1) then begin
    ; First or last cube - use only 1 nearest cube
    n_cubes_to_use = 1
endif else begin
    ; Middle cube - use up to max_sky_cubes nearest cubes
    n_cubes_to_use = (n_available lt max_sky_cubes) ? n_available : max_sky_cubes
endelse

; Select the nearest cubes
sky_cube_indices = target_indices[sorted_indices[0:n_cubes_to_use-1]]

; Sort the selected indices to maintain chronological order
sky_cube_indices = sky_cube_indices[sort(sky_cube_indices)]

return, sky_cube_indices

end
pro process_ramp_self_calibrating, ramp_data, fitted_image, cr_count, $
    full_well, debug=debug, obj=obj, timearr=timearr, dark_cube=dark_cube, $
    min_r_squared=min_r_squared, skip_first_read=skip_first_read, $
    sigma_threshold=sigma_threshold, max_iterations=max_iterations,$
    skip_last_read=skip_last_read

compile_opt idl2

if not keyword_set(debug) then debug=0
if not keyword_set(dark_cube) then dark_cube=0
if not keyword_set(skip_first_read) then skip_first_read=0
if not keyword_set(skip_las_read) then skip_last_read=0
if not keyword_set(min_r_squared) then min_r_squared=0.9
if not keyword_set(sigma_threshold) then sigma_threshold=2.5
if not keyword_set(max_iterations) then max_iterations=0

; Get dimensions
sz = size(ramp_data)
nx = sz[1]
ny = sz[2]
N = sz[3]

; Check if we have enough frames
if N lt 3 then begin
    print, 'Error: Need at least 3 frames in ramp data'
    return
endif

if skip_first_read and N lt 4 then begin
    print, 'Error: Need at least 4 frames in ramp data when skipping first read'
    return
endif

use_ramp_data = ramp_data

; Initialize output arrays
fitted_image = fltarr(nx, ny)
cr_count = bytarr(nx, ny)

; Main fitting loop with sigma-clipping
for test_y = 0, ny-1 do begin
    for test_x = 0, nx-1 do begin
        test_ramp = reform(use_ramp_data[test_x, test_y, *])
        
        ; Determine starting index based on skip_first_read flag
        if skip_first_read then begin
            use_ramp = test_ramp[1:*]
            use_times = timearr[1:*]
        endif else if skip_last_read then begin
            use_ramp = test_ramp[0:(size(timearr))[1]-1]
            use_times = timearr[0:(size(timearr))[1]-1]
        endif else begin
            use_ramp = test_ramp
            use_times = timearr
        endelse
        
        ; Exclude saturated reads
        good_reads = where(use_ramp lt full_well, n_good)
        
        if n_good ge 3 then begin
            good_ramp = use_ramp[good_reads]
            good_times = use_times[good_reads]

            ; Final fit with clipped points
            if n_elements(good_ramp) ge 3 then begin
                final_ramp = good_ramp
                final_times = good_times
                
                ; Linear least-squares fit
                sum_x = total(final_times)
                sum_y = total(final_ramp)
                sum_xy = total(final_times * final_ramp)
                sum_x2 = total(final_times^2.)
                
                denominator = n_good * sum_x2 - sum_x^2.
                if denominator ne 0.0 then begin
                    slope_value = (n_good * sum_xy - sum_x * sum_y) / denominator
                    intercept_value = (sum_y - slope_value * sum_x) / n_good
                    
                    ; Calculate R-squared
                    fit_values = slope_value * final_times + intercept_value
                    ss_res = total((final_ramp - fit_values)^2)
                    ss_tot = total((final_ramp - mean(final_ramp))^2)
                    
                    if ss_tot gt 0 then begin
                        r_squared = 1.0 - ss_res/ss_tot
                    endif else begin
                        r_squared = 0.0
                    endelse
                    
                    ; Check fit quality
                    if ((slope_value lt 0.0) or (r_squared lt min_r_squared)) and dark_cube ne 1 then begin
                        fitted_image[test_x, test_y] = !values.f_nan
                    endif else begin
                        fitted_image[test_x, test_y] = slope_value
                    endelse
                endif else begin
                    if dark_cube ne 1 then fitted_image[test_x, test_y] = !values.f_nan
                endelse
            endif else begin
                if dark_cube ne 1 then fitted_image[test_x, test_y] = !values.f_nan
            endelse
            
        endif else if n_good eq 2 then begin
            good_ramp = use_ramp[good_reads]
            good_times = use_times[good_reads]
            slope_value = (good_ramp[1] - good_ramp[0]) / (good_times[1] - good_times[0])
            
            if (slope_value lt 0.0) and dark_cube ne 1 then begin
                fitted_image[test_x, test_y] = !values.f_nan
            endif else begin
                fitted_image[test_x, test_y] = slope_value
            endelse
        endif else begin
            if dark_cube ne 1 then fitted_image[test_x, test_y] = !values.f_nan
        endelse
    endfor
endfor

; Debug output
if debug eq 1 then begin
    writefits, '~/Desktop/fitted_image_self_cal.fits', fitted_image
endif

cr_count = fitted_image
cr_count[where(finite(fitted_image))] = 0
cr_count[where(~finite(fitted_image))] = 1

end
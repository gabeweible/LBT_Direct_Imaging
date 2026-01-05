pro process_ramp_with_cr_rejection, ramp_data, fitted_image, cr_count, $
    full_well, debug=debug, obj=obj,$
    timearr=timearr, dark_cube=dark_cube, min_r_squared=min_r_squared, $
    skip_first_read=skip_first_read
    
; Process a ramp read with cosmic ray rejection
;
; INPUTS:
;   ramp_data - 3D array of ramp reads [x, y, frame]
;   full_well - Full well value in counts (should be scaled to ~end of linearity,
;               maybe 80% of full or so)
;   debug - Debug flag (optional)
;   obj - object name (string)
;   timearr - array of float times for the reads. For LMIRcam, these should be
;             evenly spaced after the first (min. exposure) read.
;   dark_cube - binary flag for creating a cube of dark frames, or not.
;   skip_first_read - binary flag to exclude the first read from fitting (optional)
;
; OUTPUTS:
;   fitted_image - Fitted image after CR rejection
;   cr_count - Image showing the number of cosmic rays detected per pixel

compile_opt idl2

if not keyword_set(debug) then debug=0
if not keyword_set(dark_cube) then dark_cube=0
if not keyword_set(skip_first_read) then skip_first_read=0
; good default min_r_squared for Alcor. TYC 5709 might like 0.8 better...
if not keyword_set(min_r_squared) then min_r_squared=0.97

; Get dimensions
sz = size(ramp_data)
nx = sz[1]
ny = sz[2]
N = sz[3]  ; Number of reads in the ramp

; Check if we have enough frames
if N lt 3 then begin
    print, 'Error: Need at least 3 frames in ramp data'
    return
endif

; If skipping first read, check we still have enough frames
if skip_first_read and N lt 4 then begin
    print, 'Error: Need at least 4 frames in ramp data when skipping first read'
    return
endif

; Initialize output arrays
fitted_image = fltarr(nx, ny)
cr_count = bytarr(nx, ny)

; NEW:: Replaced fixsen logic with a simple slope fit, plus some flagging (negative slopes, poor fits, etc.)
; Create a simple slope test image for using linear fit to all (non-saturated) reads
fitted_image = fltarr(nx, ny)
for test_y = 0, ny-1 do begin
    for test_x = 0, nx-1 do begin
        test_ramp = reform(ramp_data[test_x, test_y, *])
        
        ; Determine starting index based on skip_first_read flag
        if skip_first_read then begin
            start_idx = 1
            use_ramp = test_ramp[1:*]
            use_times = timearr[1:*]
        endif else begin
            start_idx = 0
            use_ramp = test_ramp
            use_times = timearr
        endelse
        
        ; Exclude saturated reads (at or above full_well)
        good_reads = where(use_ramp lt full_well, n_good)
        
        if n_good ge 3 then begin  ; Need at least 3 points to assess fit quality
            ; Use only non-saturated reads for the fit
            good_ramp = use_ramp[good_reads]
            good_times = use_times[good_reads]
            
            ; Perform linear least-squares fit using non-saturated reads
            ; y = mx + b, where y is the ramp data and x is the time
            sum_x = total(good_times)
            sum_y = total(good_ramp)
            sum_xy = total(good_times * good_ramp)
            sum_x2 = total(good_times^2.)
            
            ; Calculate slope and intercept using least-squares formula
            denominator = n_good * sum_x2 - sum_x^2.
            if denominator ne 0.0 then begin
                slope_value = (n_good * sum_xy - sum_x * sum_y) / denominator
                intercept_value = (sum_y - slope_value * sum_x) / n_good
                
                ; Calculate R-squared (coefficient of determination)
                fit_values = slope_value * good_times + intercept_value
                ss_res = total((good_ramp - fit_values)^2)  ; Sum of squares of residuals
                ss_tot = total((good_ramp - mean(good_ramp))^2)  ; Total sum of squares
                
                if ss_tot gt 0 then begin
                    r_squared = 1.0 - ss_res/ss_tot
                endif else begin
                    r_squared = 0.0  ; All values are identical
                endelse
                
                ; Check fit quality and slope sign
                if ((slope_value lt 0.0) or (r_squared lt min_r_squared)) and dark_cube ne 1 then begin
                    fitted_image[test_x, test_y] = !values.f_nan
                endif else begin
                    fitted_image[test_x, test_y] = slope_value
                endelse
            endif else begin
                ; Denominator is zero - poor fit
                if dark_cube ne 1 then fitted_image[test_x, test_y] = !values.f_nan
            endelse
        endif else if n_good eq 2 then begin
            ; Only 2 points - use simple slope but can't assess fit quality
            good_ramp = use_ramp[good_reads]
            good_times = use_times[good_reads]
            slope_value = (good_ramp[1] - good_ramp[0]) / (good_times[1] - good_times[0])
            
            ; Set negative slopes to NaN
            if (slope_value lt 0.0) and dark_cube ne 1 then begin
                fitted_image[test_x, test_y] = !values.f_nan
            endif else begin
                fitted_image[test_x, test_y] = slope_value
            endelse
        endif else begin
            ; Not enough good reads for a fit - set to NaN
            if dark_cube ne 1 then fitted_image[test_x, test_y] = !values.f_nan
        endelse
    endfor
endfor
;if debug eq 1 then writefits, '~/Desktop/fitted_image_test.fits', fitted_image

cr_count = fitted_image
cr_count[where(finite(fitted_image))] = 0
cr_count[where(~finite(fitted_image))] = 1

end
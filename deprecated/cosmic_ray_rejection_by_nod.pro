pro cosmic_ray_rejection_cube, object_name, raw_fits_path, start_frame, coadd, output_path, $
    read_header=read_header, coadd_type=coadd_type, $
    full_well=full_well, read_variance=read_variance, bits_to_keep=bits_to_keep, $
    sigma_cr=sigma_cr, gain=gain, legacy_mode=legacy_mode, cds=cds, end_frame=end_frame,$
    debug=debug
; Based on create_cube.pro but implementing the CR rejection algorithm from
; "Cosmic-Ray Rejection and Readout Efficiency for Large-Area Arrays" (D. J. FIXSEN et al. 2000)
;
; INPUTS (Additional to create_cube.pro):
;   full_well - Full well value for detector saturation (default: 50000)
;   read_variance - Read noise variance in counts^2 (default: 100)
;   bits_to_keep - Bits of noise to keep (default: 1.0)
;   sigma_cr - Sigma for cosmic ray cutoff (default: 5.0)
;   gain - Detector gain in photons/count (default: 1.0)
;   legacy_mode - Set to use legacy 'cds', 'midcds', etc. options instead of 
;                 the Fixsen algorithm (not recommended)
;
compile_opt idl2; 32-bit Integers and only square brackets for array indexing
newline = string(10B)

if not keyword_set(debug) then debug=0

; Set defaults for CR rejection parameters
; taken from Leisenring 2012 for LMIRCam
if n_elements(full_well) eq 0 then full_well = 32767; assumes +- 16-bit integers
if n_elements(read_variance) eq 0 then read_variance = 9.30^2. ; assumes CDS, not Fowler sampling
if n_elements(bits_to_keep) eq 0 then bits_to_keep = 1.0
if n_elements(sigma_cr) eq 0 then sigma_cr = 3.5
if n_elements(gain) eq 0 then gain = 8.9; LMIRCam fast mode electrons/ADU
if n_elements(legacy_mode) eq 0 then legacy_mode = 0

; Constants from the paper
E = 0.0024788  ; Smallest S/N ratio bin
K = 8          ; Number of signal (S) cuts
M = 64         ; Max number of reads

; Initialize the weight table for CR rejection
init_cr_weights, bits_to_keep, gain, read_variance, sigma_cr, debug=debug

starttime=systime(/JULIAN)
; Search for the raw images in the specified data path and give some output to the user
print, 'Searching for FITS files in', raw_fits_path, '...'
files = FILE_SEARCH(raw_fits_path, '*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'
if n_elements(end_frame) eq 0 then end_frame = filecount-1

; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
print, 'Initializing...'
obj_cube = list()
cr_counts = list()  ; Track number of detected cosmic rays

; default to reading this in (angles, nod/dither position flag, etc.)
if not keyword_set(read_header) then read_header=1
if not keyword_set(coadd_type) then coadd_type='res_mean'

if read_header eq 1 then begin
    angles = list()
    dits = list()
    flags = list()
endif; read_header if

; Initialize an empty coadd frame and our angle for the coadd as 0
coadd_frame = [[]]
cr_count_frame = [[]]
if read_header eq 1 then coadd_angle = 0.

; k counts what image we're on, which is used to know when to add the coadded frame to the cube
k = 1
; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'
for ii = start_frame, end_frame do begin
    print, 'File index', ii, '/', filecount-1
    print, 'Reading in:', files[ii], newline
    if read_header eq 1 then begin
        raw_frame = readfits(files[ii], head)
        angle = fxpar(head, 'LBT_PARA')
        dit = fxpar(head,'ITIME')
        flag = fxpar(head,'FLAG')
    endif else begin
        raw_frame = readfits(files[ii])
    endelse
    
    ; Check if we have multiple reads up the ramp in the FITS file
    sz = size(raw_frame)
    
    if sz[0] ge 3 then begin
        ; We have ramp data with multiple reads
        print, 'Processing frame with up-the-ramp data...'
        
        if legacy_mode then begin
            ; Use legacy CDS options if requested
            print, 'Using legacy CDS mode: ', cds
            if cds eq 'midcds' then begin
                frame = reform(raw_frame[*,*,1] - raw_frame[*,*,0])  ; second - first CDS
                cr_map = bytarr(sz[1], sz[2])  ; No CR detection in legacy mode
            endif else if cds eq 'endcds' then begin
                frame = reform(raw_frame[*,*,2] - raw_frame[*,*,0])  ; third - first CDS
                cr_map = bytarr(sz[1], sz[2])
            endif else if cds eq 'last' then begin
                frame = reform(raw_frame[*,*,2])  ; just the third frame
                cr_map = bytarr(sz[1], sz[2])
            endif else begin
                ; Default to endcds
                frame = reform(raw_frame[*,*,2] - raw_frame[*,*,0])
                cr_map = bytarr(sz[1], sz[2])
            endelse
        endif else begin
            ; Apply the Fixsen CR rejection algorithm to the ramp data
            print, 'Using Fixsen CR rejection algorithm...'
            frame = dblarr(sz[1], sz[2])
            cr_map = bytarr(sz[1], sz[2])
            
            ; Process ramp reads with CR rejection
            process_ramp_with_cr_rejection, raw_frame, frame, cr_map, $
                full_well, bits_to_keep, read_variance, gain, sigma_cr, E=E,$
                K=K, debug=debug
        endelse
    endif else begin
        ; No ramp data, just use the frame as is
        print, 'No ramp data found, using frame as-is'
        frame = float(raw_frame)
        cr_map = bytarr(sz[1], sz[2])
    endelse
    
    if debug eq 1 then writefits, '~/Desktop/Fixsen_test.fits', frame
     
    ; Add to our current coadd frame as a sum, and the same with the current coadd angle
    coadd_frame = [ [[coadd_frame]], [[frame]] ]  ; add the frame to the small group
    cr_count_frame = [ [[cr_count_frame]], [[cr_map]] ]  ; add the CR map to the small group
    
    if read_header eq 1 then coadd_angle = temporary(coadd_angle) + angle
    
    if k mod coadd eq 0 then begin
        if read_header eq 1 then begin
            dits.Add, dit
            flags.Add, flag
            coadd_angle = temporary(coadd_angle) * (1. / coadd)
            angles.Add, coadd_angle
            coadd_angle = 0.
        endif; header-reading if
        
        ; Add the mean frame to the cube and reset the frame, and add the mean angle and reset the angle
        ; Choose binning method
        if coadd_type eq 'median' then begin
            coadd_frame = median(temporary(coadd_frame), dimension=3, /even, /double)
            cr_count_frame = median(temporary(cr_count_frame), dimension=3, /even, /double)
        endif
        
        if coadd_type eq 'mean' then begin
            coadd_frame = mean(temporary(coadd_frame), dimension=3, /nan, /double)
            cr_count_frame = mean(temporary(cr_count_frame), dimension=3, /nan, /double)
        endif
        
        if coadd_type eq 'res_mean' then begin
            resistant_mean, temporary(coadd_frame), 3.0, res_mean_coadd_frame, dim=3, /double
            resistant_mean, temporary(cr_count_frame), 3.0, res_mean_cr_frame, dim=3, /double
            coadd_frame = res_mean_coadd_frame
            cr_count_frame = res_mean_cr_frame
        endif
        
        obj_cube.Add, coadd_frame
        cr_counts.Add, cr_count_frame
    
        ; reset for coadding
        coadd_frame = [[]]
        cr_count_frame = [[]]
        
    endif
    print, 'Number of Frames in Cube: ', n_elements(obj_cube)
    k += 1
endfor; ii = start_frame for loop

; Writefits takes an array as input, so we'll need to convert our list over to an array
print, 'Converting to array...'
obj_cube = (temporary(obj_cube)).toArray(/TRANSPOSE, /NO_COPY)
cr_counts = (temporary(cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
print, 'Writing cube FITS...'

; Write the cube
writefits, strcompress(output_path+object_name+'_cube.fits',/rem), obj_cube
writefits, strcompress(output_path+object_name+'_cr_count.fits',/rem), cr_counts
print, 'Cube FITS written! Writing save file...'

if read_header eq 1 then begin
    angles = (temporary(angles)).toArray(/NO_COPY)
    flags = (temporary(flags)).toArray(/NO_COPY)
    save, filename=strcompress(output_path+object_name+'_parang.sav',/rem), angles, flags
    print, 'Save file written!'
endif; read_header if
print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end

; Implementation of the Fixsen et al. (2000) CR rejection algorithm
pro init_cr_weights, bits_to_keep, gain, read_variance, sigma_cr, debug=debug
; Initialize the weight tables for CR rejection as per Fixsen et al. 2000
common cr_weights, WT, a, f, v, q
if not keyword_set(debug) then debug=0

; Constants from the paper
E = 0.0024788  ; Smallest S/N ratio bin
K = 8          ; Number of signal (S) cuts
M = 64         ; Max number of reads

; Initialize the WT array
WT = dblarr(K, M, M/2+1)

; Calculate parameters
; For LMIRCam data with 3 samples (reset + 2 reads)
N = 3  ; Number of reads

; Initialize variables used for calculations
Z = dblarr(M)
Y = dblarr(M)

; For each S/N ratio bin
for little_k = 0, K-1 do begin
    little_s = E * exp(little_k)  ; S/N ratio
    WT[little_k, 0, 0] = 0.0
    
    ; Reset variables for this iteration
    Z[*] = 0.0
    Y[*] = 0.0
    little_x = 0.0
    little_y = 0.0
    
    ; For all ramp lengths
    for little_n = 1, N-1 do begin
        W = reform(WT[little_k, little_n, *])
        
        ; New x, new last Y
        little_x = 1.0/(2.0 + little_s - little_x)
        Y[little_n] = little_x
        little_y++
        
        ; Update row and sum
        little_z = 0.0
        for i = 1, little_n-1 do begin
            Y[i] *= little_x
            Z[i] += little_y * Y[i]
            little_z += Z[i]
        endfor
        
        ; Sum weight
        little_y = little_x
        for i = 1, little_n-1 do little_y += Y[i]
        
        ; Differences - FIXED: using correct variable naming
        for i = 0, (little_n+1)/2-1 do W[i] = (Z[i] - Z[i+1]) * gain * (N-1)
        
        ; Final weight
        Z[little_n] = little_y
        W[(little_n+1)/2] = little_z + little_y
    endfor
    
    ; Renormalize last row
    little_x = 1.0/(little_z + little_y)
    for i = 0, N/2-1 do W[i] *= little_x
endfor

; Set the shared variables for use in the CR rejection routine
a = bits_to_keep * sqrt(gain)                       ; Renormalize Bits
f = 3.0 * read_variance / (N * gain * (N+1) * (N+2))  ; Output offset
q = sigma_cr / gain                                ; Renormalize cutoff
v = 2.0 * q * read_variance * (N*N + 1) / (gain * N * (N-1))  ; Renormalize variance

if debug eq 1 then begin
	print, '--- Weight Table Initialization ---'
	print, 'a value:', a
	print, 'f value:', f
	print, 'v value:', v
	print, 'q value:', q
	print, 'a * sqrt(f) =', a * sqrt(f)  ; Check if this equals 2.07954
	print, '--- End Weight Table Info ---'
endif

end

pro process_ramp_with_cr_rejection, ramp_data, fitted_image, cr_count, $
    full_well, bits_to_keep, read_variance, gain, sigma_cr, E=E, K=K, debug=debug
; Process a ramp read with cosmic ray rejection
; Implements the algorithm from Fixsen et al. 2000
; Inputs:
;   ramp_data - 3D array of ramp reads [x, y, frame]
;   full_well - Full well value in counts
;   bits_to_keep, read_variance, gain, sigma_cr - Algorithm parameters
; Outputs:
;   fitted_image - Fitted image after CR rejection
;   cr_count - Image showing the number of cosmic rays detected per pixel

if not keyword_set(debug) then debug=0

common cr_weights, WT, a, f, v, q

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

; Initialize output arrays
fitted_image = dblarr(nx, ny)
cr_count = bytarr(nx, ny)

; Print CR rejection parameters
if debug eq 1 then begin

	; Create a simple slope test image for debugging
	simple_slope = dblarr(nx, ny)
	for test_y = 0, ny-1 do begin
		for test_x = 0, nx-1 do begin
			test_ramp = reform(ramp_data[test_x, test_y, *])
			simple_slope[test_x, test_y] = (test_ramp[N-1] - test_ramp[0]) / (N-1)
		endfor
	endfor
	writefits, 'simple_slope_test.fits', simple_slope
	
	print, '--- CR Rejection Parameters ---'
	print, 'full_well:', full_well
	print, 'read_variance:', read_variance
	print, 'gain:', gain
	print, 'sigma_cr:', sigma_cr
	print, 'bits_to_keep:', bits_to_keep
	print, 'N (reads):', N
	print, '--- End CR Rejection Params ---'

endif

; Constants for the algorithm
EVe = E * read_variance * exp(1)  ; Smallest S/N ratio bin × readout variance × e
little_n = N - 1  ; Number of differences between reads

; Process each pixel
for little_y = 0, ny-1 do begin
    for little_x = 0, nx-1 do begin
        ; Extract the ramp for this pixel
        pixel_ramp = reform(ramp_data[little_x, little_y, *])
        
        ; Check if the pixel is saturated in any frame
        saturated = 0
        for i = 0, N-1 do begin
            if pixel_ramp[i] gt full_well then saturated = 1
        endfor
        
        if saturated then begin
            ; Mark as bad pixel and continue to next
            cr_count[little_x, little_y] = N
            fitted_image[little_x, little_y] = 0
            continue
        endif
        
        ; Convert to float array for processing
        Z = dblarr(N)  ; Array to store data points
        C = dblarr(N)  ; Array to store indices of CRs
        
        ; Convert data to float array
        for i = 0, N-1 do Z[i] = float(pixel_ramp[i])
        
        ; Calculate initial signal estimate from first and last points
        signal = (Z[N-1] - Z[0]) / (N-1)
        
        ; Look for cosmic rays
        little_c = -1  ; Count of cosmic rays, start at -1 for 0-indexed
        
        ; For 3-sample data (LMIRCam), we have a simplified case
        if N eq 3 then begin
            ; Calculate differences between frames
            diff1 = Z[1] - Z[0]  ; First difference
            diff2 = Z[2] - Z[1]  ; Second difference
            
            ; If these differences vary significantly, we have a CR
            thresh = sigma_cr * sqrt(2.0 * read_variance)
            if abs(diff2 - diff1) gt thresh then begin
                ; Found a cosmic ray
                little_c = 0  ; 1 CR detected
                
                ; Identify which segment has the CR
                if diff2 gt diff1 + thresh then begin
                    ; CR in second segment
                    C[0] = 1  ; Index where CR happens
                    ; Use only first segment for signal
                    signal = diff1
                endif else begin
                    ; CR in first segment
                    C[0] = 0  ; Index where CR happens
                    ; Use only second segment for signal
                    signal = diff2
                endelse
            endif else begin
                ; No CR detected, use optimal weighting of both segments
                ; For 3-sample case with no CR, optimal estimate is mean of differences
                signal = (diff1 + diff2) / 2.0
            endelse
            
            ; DIRECT APPROACH FOR 3-SAMPLE CASE
            if little_c lt 0 then begin
                ; No CR detected, use simple slope
                fitted_image[little_x, little_y] = signal
            endif else begin
                ; Use the clean segment
                if C[0] eq 0 then begin
                    fitted_image[little_x, little_y] = diff2
                endif else begin
                    fitted_image[little_x, little_y] = diff1
                endelse
            endelse
            
            ; Record number of cosmic rays
            cr_count[little_x, little_y] = little_c + 1  ; Add 1 because c is 0-indexed
            continue  ; Skip the complex weighting for 3-sample case
            
        endif else begin
            ; General case implementation for more than 3 samples
            ; This follows the paper's algorithm more closely
            
            ; Find worst point
            worst_X = 0
            x_val = v  ; Threshold for CR detection
            
            for i = 1, N-1 do begin
                dev_squared = (signal - (Z[i] - Z[i-1]))^2
                if dev_squared gt x_val then begin
                    worst_X = i
                    x_val = dev_squared
                endif
            endfor
            
            ; Check if it's a cosmic ray
            while x_val gt signal*q + v and worst_X gt 0 do begin
                ; Found a cosmic ray
                little_c++
                
                ; Adjust signal estimate
                signal += (signal + Z[worst_X] - Z[worst_X-1]) / (little_n - little_c - 1)
                
                ; Record CR position
                for jj = little_c, 1, -1 do begin
                    if C[jj-1] gt worst_X then begin
                        C[jj] = C[jj-1]
                    endif else begin
                        break
                    endelse
                endfor
                C[little_c] = worst_X  ; FIXED: using little_c instead of little_k
                
                ; Look for more CRs
                x_val = v
                worst_X = 0
                
                ; Check each segment between CRs
                for jj = 0, little_c do begin
                    start_Y = (jj eq 0) ? 0 : C[jj-1]  ; FIXED: using jj instead of k
                    end_Y = (jj eq little_c) ? N-1 : C[jj]
                    
                    for loop_Y = start_Y+1, end_Y do begin
                        dev_squared = (signal - (Z[loop_Y] - Z[loop_Y-1]))^2
                        if dev_squared gt x_val then begin
                            worst_X = loop_Y
                            x_val = dev_squared
                        endif
                    endfor
                endfor
            endwhile
        endelse
        
        ; Record number of cosmic rays
        cr_count[little_x, little_y] = little_c + 1  ; Add 1 because c is 0-indexed
        
        ; Find appropriate S/N bin for weighting
        little_k2 = 0
        x_val = EVe
        while signal gt x_val and little_k2 lt K-1 do begin
            little_k2++
            x_val *= exp(1)
        endwhile
        
        ; Apply weighting to get final signal
        final_signal = f  ; Start with offset
        
        if little_c ge 0 then begin
            ; With CRs, apply segment-by-segment weighting
            weight_sum = 0.0
            
            ; For each clean segment, apply appropriate weights
            for i = 0, little_c do begin
                start_idx = (i eq 0) ? 0 : C[i-1]
                end_idx = (i eq little_c) ? N-1 : C[i]
                segment_len = end_idx - start_idx
                
                if segment_len gt 0 then begin
                    ; Get weights for this segment length
                    W = reform(WT[little_k2, segment_len, *])
                    
                    ; Apply weights to differences
                    for j = 0, segment_len/2 do begin
                        if j lt segment_len/2 then begin
                            diff_weight = W[j]
                            final_signal += diff_weight * (Z[start_idx+j] - Z[end_idx-j])
                        endif else begin
                            weight_sum += W[j]
                        endelse
                    endfor
                endif
            endfor
            
            ; Normalize by weight sum
            if weight_sum gt 0 then final_signal /= weight_sum
        endif else begin
            ; No CRs, use optimal weighting for full ramp
            W = reform(WT[little_k2, little_n, *])
            weight_sum = 0.0
            
            ; Apply weights to differences (FIXED: removed negative weight check)
            for i = 0, little_n/2 do begin
                if i lt little_n/2 then begin
                    diff_weight = W[i]
                    final_signal += diff_weight * (Z[i] - Z[little_n-i])
                endif else begin
                    weight_sum += W[i]
                endelse
            endfor
            
            ; Normalize by weight sum
            if weight_sum gt 0 then final_signal /= weight_sum
        endelse
        
        ; For debugging, print some information about a sample pixel
        if debug and little_x eq nx/2 and little_y eq ny/2 then begin
            print, '--- Sample Pixel Info ---'
            print, '(x,y):', little_x, little_y
            print, 'Ramp values:', Z
            print, 'Simple slope:', (Z[N-1] - Z[0]) / (N-1)
            print, 'CRs detected:', little_c + 1
            print, 'Final signal:', final_signal
            print, 'a * sqrt(final_signal):', a * sqrt(final_signal)
            print, '--- End Sample Pixel ---'
        endif
        
        ; Convert to output format
        if final_signal lt 0.0 then begin
            fitted_image[little_x, little_y] = 0.0
        endif else begin
            fitted_image[little_x, little_y] = a * sqrt(final_signal)
        endelse
    endfor
endfor

; For debugging, check how many pixels have the constant value
if debug then begin
    const_count = 0
    for y = 0, ny-1 do begin
        for x = 0, nx-1 do begin
            if abs(fitted_image[x, y] - 2.07954) lt 0.001 then const_count++
        endfor
    endfor
    print, 'Pixels with value ~2.07954:', const_count, ' of ', nx*ny, ' (', float(const_count)/(nx*ny)*100.0, '%)'
end

end
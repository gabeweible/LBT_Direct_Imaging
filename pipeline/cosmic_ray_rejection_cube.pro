pro cosmic_ray_rejection_cube, object_name, raw_fits_path, start_frame, coadd, output_path, $
    read_header=read_header, coadd_type=coadd_type, $
    full_well=full_well, legacy_mode=legacy_mode, cds=cds, end_frame=end_frame,$
    debug=debug, dark_frame_start=dark_frame_start, max_frames_per_group=max_frames_per_group,$
    resume_from_nod=resume_from_nod, dark_cube=dark_cube,$
    min_r_squared=min_r_squared, skip_first_read=skip_first_read,$
    skip_last_read=skip_last_read
; Based on create_cube.pro but implementing the CR rejection algorithm 

; INPUTS (Additional to create_cube.pro):
;   full_well - Full well value for detector saturation
;   legacy_mode - Set to use legacy 'cds', 'midcds', etc. options (disabled currently)
;   dark_frame_start - Frame index at which dark frames start (frames >= this index are treated as darks)
;   max_frames_per_group - Maximum frames per group within each nod (default: 100)
;   resume_from_nod - Nod number to resume processing from (0-based indexing)
;
compile_opt idl2; 32-bit Integers and only square brackets for array indexing
newline = string(10B)


if not keyword_set(resume_from_nod) then resume_from_nod = 0
if not keyword_set(debug) then debug=0
if not keyword_set(max_frames_per_group) then max_frames_per_group=1000
if not keyword_set(skip_first_read) then skip_first_read=0

; Set defaults for CR rejection parameters
if n_elements(full_well) eq 0 then full_well = 4095
; good default for Alcor. TYC 5709 might like 0.8 better...
if not keyword_set(min_r_squared) then min_r_squared=0.9

starttime=systime(/JULIAN)

print, 'Maximum frames per group: ', max_frames_per_group
if resume_from_nod gt 0 then print, 'Resuming from nod: ', resume_from_nod

; Search for the raw images in the specified data path and give some output to the user
print, 'Searching for FITS files in', raw_fits_path, '...'
files = FILE_SEARCH(raw_fits_path, '*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'
if n_elements(end_frame) eq 0 then end_frame = filecount-1

; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
print, 'Initializing...'
prev_flag = ''
nod_counter = 0
group_counter = 0
frames_in_current_group = 0

; NEW: Variables to track resume functionality
resume_mode = (resume_from_nod gt 0)
skip_processing = resume_mode  ; Start by skipping if resuming
current_nod_for_resume = 0

nod_cube = list()
nod_cr_counts = list()
nod_angles = list()
nod_dits = list()

; Initialize dark frame storage
if n_elements(dark_frame_start) gt 0 and dark_frame_start ne 'None' then begin
    dark_cube_list = list()
    dark_cr_counts = list()
    dark_angles = list()
    dark_dits = list()
endif

; default to reading this in (angles, nod/dither position flag, etc.)
if not keyword_set(read_header) then read_header=0
if not keyword_set(coadd_type) then coadd_type='res_mean'

if read_header eq 1 then begin
    angles = list()
    dits = list()
    flags = list()
endif; read_header if

; Initialize variables
; Initialize an empty coadd frame and our angle for the coadd as 0
coadd_frame = list()
cr_count_frame = list()
if read_header eq 1 then coadd_angle = 0.

; k counts what image we're on, which is used to know when to add the coadded frame to the cube
k = 1
; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'

for ii = start_frame, end_frame do begin
    ; Reduce print statements for speed
    print, 'File index', ii, '/', filecount-1
    
    ; Check if this is a dark frame
    is_dark_frame = 0
    if n_elements(dark_frame_start) gt 0 and dark_frame_start ne 'None' then begin
        if ii ge dark_frame_start then begin
            is_dark_frame = 1
            print, 'Processing as dark frame (index >= ', dark_frame_start, ')'
        endif
    endif
    
    raw_frame = readfits(files[ii], head)
    angle = fxpar(head, 'LBT_PARA')
    ; Extract timing information from header
    dit = float(fxpar(head, 'EXPTIME'))
    ft = float(fxpar(head, 'FRAME'))
    sz = size(raw_frame)
    ngroup = sz[3]  ; Number of reads
    
    ; Calculate time array for each read
    exptimearr = ft + float(indgen(ngroup)) * (dit - ft) / (ngroup - 1)
    exptimearr /= 1000. ; ms to s
    print, 'exptimearr: ', exptimearr, 's'
    
    flag = fxpar(head,'FLAG')
    
    if ii eq 0 then prev_flag = flag; set first flag
    
    ; NEW: Handle resume logic - detect nod changes to track current nod number
    if resume_mode and prev_flag ne '' and flag ne prev_flag then begin
        current_nod_for_resume += 1
        print, 'Resume mode: Detected nod change, current_nod_for_resume now = ', current_nod_for_resume
        
        ; Check if we've reached the target nod to resume from
        if current_nod_for_resume ge resume_from_nod then begin
            skip_processing = 0
            nod_counter = current_nod_for_resume  ; Sync the actual nod counter
            print, 'Resume mode: Reached target nod ', resume_from_nod, ', starting processing...'
        endif
    endif
    
    ; NEW: Skip processing if we haven't reached the resume point yet
    if skip_processing then begin
       ; if debug then print, 'Resume mode: Skipping frame ', ii, ' (current nod ', current_nod_for_resume, ' < target ', resume_from_nod, ')'
        prev_flag = flag
        continue
    endif
    
    ; FULL PROCESSING MODE: Continue with original logic (modified for splitting)
    ; Check if we have multiple reads up the ramp in the FITS file
    sz = size(raw_frame)
    
    if sz[0] ge 3 then begin
        ; We have ramp data with multiple reads
        print, 'Processing frame with up-the-ramp data...'
        
        frame = fltarr(sz[1], sz[2])
        cr_map = bytarr(sz[1], sz[2])
        
        ; Process ramp reads with CR rejection
        process_ramp_self_calibrating, raw_frame, frame, cr_map, $
            full_well, debug=debug, obj=object_name, timearr=exptimearr,$
            dark_cube=dark_cube, min_r_squared=min_r_squared,$
            skip_first_read=skip_first_read, skip_last_read=skip_last_read
            
       ; endelse
    endif else begin
        ; No ramp data, just use the frame as is
        print, 'No ramp data found, using frame as-is'
        frame = double(raw_frame)
        cr_map = bytarr(sz[1], sz[2])
    endelse

    ; Handle dark frames separately
    if is_dark_frame then begin
        ; Add dark frame directly to dark cube (no coadding for darks)
        dark_cube_list.Add, frame
        dark_cr_counts.Add, cr_map
        if read_header eq 1 then begin
            dark_angles.Add, angle
            dark_dits.Add, dit
        endif
        print, 'Added frame to dark cube. Total dark frames: ', dark_cube_list.Count()
    endif else begin
        ; Process science frames with coadding as before
        ; Add to our current coadd frame as a sum, and the same with the current coadd angle
        coadd_frame.Add, frame  ; add the frame to the small group
        cr_count_frame.Add, cr_map  ; add the CR map to the small group
        
        if read_header eq 1 then coadd_angle = temporary(coadd_angle) + angle
        
        if k mod coadd eq 0 then begin
            print, 'coadding...'
            dits.Add, dit
            flags.Add, flag
            coadd_angle = temporary(coadd_angle) * (1. / coadd)
            angles.Add, coadd_angle
            averaged_angle = coadd_angle  ; Store it before reset
            coadd_angle = 0.
            
            ; Add the mean frame to the cube and reset the frame, and add the mean angle and reset the angle
            ; Choose binning method
            ; When ready to coadd (k mod coadd eq 0):
            if coadd_frame.Count() ge 2 then begin
                ; Convert list to array for processing
                temp_cube = coadd_frame.toArray(/TRANSPOSE, /NO_COPY)
                temp_cr_cube = cr_count_frame.toArray(/TRANSPOSE, /NO_COPY)
                
                ; Apply your coadding method
                if coadd_type eq 'median' then begin
                    coadd_frame = median(temp_cube, dimension=3, /even, /double)
                    cr_count_frame = median(temp_cr_cube, dimension=3, /even, /double)
                endif else if coadd_type eq 'mean' then begin
                    coadd_frame = mean(temp_cube, dimension=3, /nan, /double)
                    cr_count_frame = mean(temp_cr_cube, dimension=3, /nan, /double)
                endif else if coadd_type eq 'res_mean' then begin
                    resistant_mean, temp_cube, 3.5, res_mean_coadd_frame, dim=3, /double
                    resistant_mean, temp_cr_cube, 3.5, res_mean_cr_frame, dim=3, /double
                    coadd_frame = res_mean_coadd_frame
                    cr_count_frame = res_mean_cr_frame
                endif
            endif else begin
                ; Single frame case
                coadd_frame = coadd_frame[0]
                cr_count_frame = cr_count_frame[0]
            endelse
            
            ; Detect nod switch
            if prev_flag ne '' and flag ne prev_flag then begin
                print, 'Nod switch detected...'
                print, 'prev_flag: ', prev_flag
                print, 'flag: ', flag
                ; Write cube from previous group
                if nod_cube.Count() gt 0 then begin
                    outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
                
                    nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
                    nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
                    angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
                    dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)
                
                    writefits, outbase+'_cube.fits', nod_cube_arr
                    writefits, outbase+'_cr_count.fits', nod_cr_arr
                    save, filename=outbase+'_parang.sav', angles_arr, dits_arr
                
                    print, 'Saved nod cube for ', prev_flag, ' as nod ', nod_counter, ' group ', group_counter
                endif
            
                ; Clear for new nod
                nod_cube = list()
                nod_cr_counts = list()
                nod_angles = list()
                nod_dits = list()
                nod_counter += 1
                group_counter = 0
                frames_in_current_group = 0
            endif else if frames_in_current_group ge max_frames_per_group then begin
                ; Save current group and start new group within same nod
                if nod_cube.Count() gt 0 then begin
                    outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
					
                    nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
                    nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
                    angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
                    dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)
                
                    writefits, outbase+'_cube.fits', nod_cube_arr
                    writefits, outbase+'_cr_count.fits', nod_cr_arr
                    save, filename=outbase+'_parang.sav', angles_arr, dits_arr
                
                    print, 'Saved nod cube for ', prev_flag, ' as nod ', nod_counter-1, ' group ', group_counter
                endif
                
                ; Start new group within same nod
                nod_cube = list()
                nod_cr_counts = list()
                nod_angles = list()
                nod_dits = list()
                group_counter += 1
                frames_in_current_group = 0
            endif

            ; Append current coadd to current group
            nod_cube.Add, coadd_frame
            nod_cr_counts.Add, cr_count_frame
            nod_angles.Add, averaged_angle  ; Use the stored averaged angle
            nod_dits.Add, dit
            frames_in_current_group += 1
            
            prev_flag = flag
            
            ; Reset for next coadd group
            coadd_frame = list()
            cr_count_frame = list()
            coadd_frame_initialized = 0
            
        endif
        print, 'Number of Frames in Current Group: ', nod_cube.Count(), ', Group: ', group_counter
    endelse
    
    ; Only increment k for science frames (for coadding logic)
    if not is_dark_frame then k += 1

endfor; ii = start_frame for loop

;write any remaining science frames (only in full processing mode)
if nod_cube.Count() gt 0 then begin
    outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)

    nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
    nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
    angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
    dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)

    writefits, outbase+'_cube.fits', nod_cube_arr
    writefits, outbase+'_cr_count.fits', nod_cr_arr
    save, filename=outbase+'_parang.sav', angles_arr, dits_arr

    print, 'Final nod cube saved for ', prev_flag, ' as nod ', nod_counter, ' group ', group_counter
endif

; Write dark frames if any were collected
if dark_cube_list.Count() gt 0 then begin
    print, 'Writing dark frame cube with ', dark_cube_list.Count(), ' frames...'
    dark_outbase = output_path + strcompress(object_name + '_darks', /r)
    
    dark_cube_arr = (temporary(dark_cube_list)).toArray(/TRANSPOSE, /NO_COPY)
    dark_cr_arr = (temporary(dark_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
    
    writefits, dark_outbase+'.fits', dark_cube_arr
    writefits, dark_outbase+'_cr_count.fits', dark_cr_arr
    
    if read_header eq 1 then begin
        dark_angles_arr = (temporary(dark_angles)).toArray(/NO_COPY)
        dark_dits_arr = (temporary(dark_dits)).toArray(/NO_COPY)
        save, filename=dark_outbase+'_parang.sav', dark_angles_arr, dark_dits_arr
    endif
    
    print, 'Dark frame cube saved as: ', dark_outbase, '.fits'
endif

if dark_cube_list.Count() eq 0 then begin
    print, 'No dark frames found to write.'
endif

print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end
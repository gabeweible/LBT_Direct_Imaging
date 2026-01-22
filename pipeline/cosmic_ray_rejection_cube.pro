pro cosmic_ray_rejection_cube, object_name, raw_fits_path, start_frame, coadd, output_path, $
    read_header=read_header, coadd_type=coadd_type, $
    full_well=full_well, legacy_mode=legacy_mode, cds=cds, end_frame=end_frame,$
    debug=debug, dark_frame_start=dark_frame_start, max_frames_per_group=max_frames_per_group,$
    resume_from_nod=resume_from_nod, dark_cube=dark_cube,$
    min_r_squared=min_r_squared, skip_first_read=skip_first_read,$
    skip_last_read=skip_last_read, fpn=fpn, verbose=verbose,$
    del_nodframe1=del_nodframe1, nod_counter_add=nod_counter_add
; Based on create_cube.pro but implementing the CR rejection algorithm 
; OPTIMIZED for CDS mode with coadd=1 and large frame counts

; INPUTS (Additional to create_cube.pro):
;   full_well - Full well value for detector saturation
;   legacy_mode - Set to 'cds' to use correlated double sampling mode
;   dark_frame_start - Frame index at which dark frames start (frames >= this index are treated as darks)
;   max_frames_per_group - Maximum frames per group within each nod (default: 100)
;   resume_from_nod - Nod number to resume processing from (0-based indexing)
;   del_nodframe1 - If set, exclude first frame of each nod from coadding and save separately
;
compile_opt idl2; 32-bit Integers and only square brackets for array indexing
newline = string(10B)

if not keyword_set(resume_from_nod) then resume_from_nod = 0
if not keyword_set(debug) then debug=0
if not keyword_set(max_frames_per_group) then max_frames_per_group=fix(fpn - 1)/coadd
if not keyword_set(skip_first_read) then skip_first_read=0
if not keyword_set(legacy_mode) then legacy_mode=''
if not keyword_set(del_nodframe1) then del_nodframe1 = 0

; Set defaults for CR rejection parameters
; largest DN possible for a frame, not necessarily physical saturation level
if n_elements(full_well) eq 0 then full_well = 4095*0.98

; good default?
if not keyword_set(min_r_squared) then min_r_squared=0.9

; Optimization flag for CDS mode
is_cds_mode = (legacy_mode eq 'cds')

starttime=systime(/JULIAN)

print, 'Maximum frames per group: ', max_frames_per_group
if resume_from_nod gt 0 then print, 'Resuming from nod: ', resume_from_nod
if is_cds_mode then print, 'Using CDS (Correlated Double Sampling) mode - OPTIMIZED'
if del_nodframe1 eq 1 then print, 'del_nodframe1 enabled: First frame of each nod will be saved separately'

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

; Track first frame of each nod for del_nodframe1
if del_nodframe1 eq 1 then begin
    first_nod_frames = list()
    first_nod_cr_counts = list()
    first_nod_angles = list()
    first_nod_dits = list()
    first_nod_flags = list()
    is_first_frame_of_nod = 1  ; Flag to track if this is the first frame of current nod
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
k = 0

; Flag to track if we've updated the header for CDS mode (optimization)
cds_header_updated = 0

; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'
    
for ii = start_frame, end_frame do begin
    ; Reduce print statements for speed
    if (ii mod 100) eq 0 or ii eq start_frame or ii eq end_frame then begin
        print, 'File index', ii, '/', filecount-1
    endif
    
    ; Check if this is a dark frame
    is_dark_frame = 0
    if n_elements(dark_frame_start) gt 0 and dark_frame_start ne 'None' then begin
        if ii ge dark_frame_start then begin
            is_dark_frame = 1
        endif
    endif
    
    raw_frame = readfits(files[ii], head)
    angle = fxpar(head, 'LBT_PARA')
    ; Extract timing information from header
    dit = float(fxpar(head, 'EXPTIME'))
    ft = float(fxpar(head, 'FRAME'))
    sz = size(raw_frame)
    ngroup = sz[3]  ; Number of reads
    
    ; Only calculate time array if NOT in CDS mode (optimization)
    if not is_cds_mode then begin
        exptimearr = ft + float(indgen(ngroup)) * (dit - ft) / (ngroup - 1)
        exptimearr /= 1000. ; ms to s
    endif
    
    flag = fxpar(head,'FLAG')
    
    if ii eq 0 then prev_flag = flag; set first flag
    
    ; CRITICAL: Detect nod switch IMMEDIATELY after reading flag, before any processing
    ; This must happen before resume logic and before frame processing
    nod_switch_detected = 0
    if prev_flag ne '' and flag ne prev_flag and not skip_processing then begin
        nod_switch_detected = 1
        print, 'Nod switch detected at file index ', ii
        print, 'prev_flag: ', prev_flag
        print, 'flag: ', flag
    endif
    
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
        prev_flag = flag
        continue
    endif
    
    ; FULL PROCESSING MODE: Continue with original logic (modified for splitting)
    ; Check if we have multiple reads up the ramp in the FITS file
    sz = size(raw_frame)
    
    if sz[0] ge 3 then begin
        ; We have ramp data with multiple reads
        
        ; OPTIMIZED CDS MODE: Correlated Double Sampling
        if is_cds_mode then begin
            ; Direct CDS calculation in one operation: (last - first) * scale_factor
            ; This is 3x faster than separate operations
            frame = (raw_frame[*,*,ngroup-1] - raw_frame[*,*,0]) * (1000. / (dit - ft))
            
            ; Skip CR map allocation entirely for CDS - saves ~16 MB per 2048x2048 frame
            ; cr_map is not needed for CDS mode
            
            ; Only update header once at very first frame or after nod switch (not every frame!)
            if (not cds_header_updated) or (nod_switch_detected and k eq 0) then begin
                sxaddpar, head, 'CDSMODE', 'T', 'Correlated Double Sampling applied'
                sxaddpar, head, 'CDSINTEG', (dit - ft) / 1000., 'CDS integration time (seconds)'
                sxaddpar, head, 'BUNIT', 'counts/s', 'Data units after CDS'
                cds_header_updated = 1
            endif
            
        endif else begin
            ; Standard ramp fitting with CR rejection
            frame = fltarr(sz[1], sz[2])
            cr_map = bytarr(sz[1], sz[2])
            
            ; Process ramp reads with CR rejection
            process_ramp_self_calibrating, raw_frame, frame, cr_map, $
                full_well, debug=debug, obj=object_name, timearr=exptimearr,$
                dark_cube=dark_cube, min_r_squared=min_r_squared,$
                skip_first_read=skip_first_read, skip_last_read=skip_last_read
        endelse
            
    endif else begin
        ; No ramp data, just use the frame as is
        print, 'No ramp data found, using frame as-is'
        frame = raw_frame
        if not is_cds_mode then cr_map = bytarr(sz[1], sz[2])
    endelse

    ; Handle dark frames separately
    if is_dark_frame then begin
        ; Add dark frame directly to dark cube (no coadding for darks)
        dark_cube_list.Add, frame
        if not is_cds_mode then dark_cr_counts.Add, cr_map
        if read_header eq 1 then begin
            dark_angles.Add, angle
            dark_dits.Add, dit
        endif
        if (dark_cube_list.Count() mod 100) eq 0 then begin
            print, 'Added frame to dark cube. Total dark frames: ', dark_cube_list.Count()
        endif
    endif else begin
        ; HANDLE NOD SWITCH BEFORE PROCESSING NEW NOD'S FRAMES
        if nod_switch_detected then begin
            ; Discard any partial coadd at nod boundary
            if coadd_frame.Count() gt 0 then begin
                print, 'Discarding ', coadd_frame.Count(), ' frames from partial coadd at nod boundary'
            endif
            
            ; Write cube from previous nod
            if nod_cube.Count() gt 0 then begin
                if keyword_set(nod_counter_add) then BEGIN
                    outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter+nod_counter_add, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
                endif else BEGIN
                    outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
                endelse
                
                nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
                angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
                dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)
            
                ; Write with header preservation
                writefits, outbase+'_cube.fits', nod_cube_arr, head
                
                ; Only write CR counts if not in CDS mode
                if not is_cds_mode then begin
                    nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
                    writefits, outbase+'_cr_count.fits', nod_cr_arr
                endif
                
                save, filename=outbase+'_parang.sav', angles_arr, dits_arr
            
                print, 'Saved nod cube for ', prev_flag, ' as nod ', nod_counter, ' group ', group_counter
            endif else if nod_cube.Count() gt 0 then begin
                print, 'Discarding ', nod_cube.Count(), ' frames at nod boundary (too few for a complete group)'
            endif
        
            ; Clear for new nod
            nod_cube = list()
            nod_cr_counts = list()
            nod_angles = list()
            nod_dits = list()
            nod_counter += 1
            group_counter = 0
            frames_in_current_group = 0
            
            ; Reset coadd buffers (discard partial coadd)
            coadd_frame = list()
            cr_count_frame = list()
            coadd_angle = 0.
            
            ; Reset k to maintain consistent coadding alignment for each nod
            k = 0
            
            ; Update prev_flag
            prev_flag = flag
            
            ; NEW: Set flag for next nod's first frame
            if del_nodframe1 eq 1 then is_first_frame_of_nod = 1
            
            ; Reset CDS header flag for new nod
            cds_header_updated = 0
        endif
        
        ; Process science frames with coadding as before
        
        ; NEW: Handle first frame of nod if del_nodframe1 is set
        if del_nodframe1 eq 1 and is_first_frame_of_nod eq 1 then begin
            ; Save this frame to the special first-frame cube
            first_nod_frames.Add, frame
            if not is_cds_mode then first_nod_cr_counts.Add, cr_map
            if read_header eq 1 then begin
                first_nod_angles.Add, angle
                first_nod_dits.Add, dit
                first_nod_flags.Add, flag
            endif
            if (first_nod_frames.Count() mod 100) eq 0 then begin
                print, 'Saved first frame of nod to separate cube. Total first frames: ', first_nod_frames.Count()
            endif
            
            ; Mark that we've processed the first frame of this nod
            is_first_frame_of_nod = 0
            
            ; IMPORTANT: Increment k here so the coadding counter stays in sync
            ; This ensures consistent coadding behavior across nods
            k += 1
            
            ; Skip adding to coadd_frame - continue to next iteration
            prev_flag = flag
            continue
        endif
        
        ; Add to our current coadd frame as a sum, and the same with the current coadd angle
        coadd_frame.Add, frame  ; add the frame to the small group
        if not is_cds_mode then cr_count_frame.Add, cr_map  ; add the CR map only if not CDS
        
        if read_header eq 1 then coadd_angle = temporary(coadd_angle) + angle
        
        if k mod coadd eq 0 then begin
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
                
                ; Apply your coadding method
                if coadd_type eq 'median' then begin
                    coadd_frame = median(temp_cube, dimension=3, /even)
                    if not is_cds_mode then begin
                        temp_cr_cube = cr_count_frame.toArray(/TRANSPOSE, /NO_COPY)
                        cr_count_frame = median(temp_cr_cube, dimension=3, /even)
                    endif
                endif else if coadd_type eq 'mean' then begin
                    coadd_frame = mean(temp_cube, dimension=3, /nan)
                    if not is_cds_mode then begin
                        temp_cr_cube = cr_count_frame.toArray(/TRANSPOSE, /NO_COPY)
                        cr_count_frame = mean(temp_cr_cube, dimension=3, /nan)
                    endif
                endif else if coadd_type eq 'res_mean' then begin
                    resistant_mean, temp_cube, 3.5, res_mean_coadd_frame, dim=3
                    coadd_frame = res_mean_coadd_frame
                    if not is_cds_mode then begin
                        temp_cr_cube = cr_count_frame.toArray(/TRANSPOSE, /NO_COPY)
                        resistant_mean, temp_cr_cube, 3.5, res_mean_cr_frame, dim=3
                        cr_count_frame = res_mean_cr_frame
                    endif
                endif
            endif else begin
                ; Single frame case (coadd=1)
                coadd_frame = coadd_frame[0]
                if not is_cds_mode then cr_count_frame = cr_count_frame[0]
            endelse
            
            ; Append current coadd to current group
            nod_cube.Add, coadd_frame
            if not is_cds_mode then nod_cr_counts.Add, cr_count_frame
            nod_angles.Add, averaged_angle  ; Use the stored averaged angle
            nod_dits.Add, dit
            frames_in_current_group += 1
            
            ; Check if we need to split into a new group (within same nod) AFTER adding
            if frames_in_current_group ge max_frames_per_group then begin
                ; Save current group and start new group within same nod
                if nod_cube.Count() gt 0 then begin
                    if keyword_set(nod_counter_add) then BEGIN
                        outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter+nod_counter_add, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
                    endif else BEGIN
                        outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
                    endelse
					
                    nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
                    angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
                    dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)
                
                    ; Write with header preservation
                    writefits, outbase+'_cube.fits', nod_cube_arr, head
                    
                    ; Only write CR counts if not in CDS mode
                    if not is_cds_mode then begin
                        nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
                        writefits, outbase+'_cr_count.fits', nod_cr_arr
                    endif
                    
                    save, filename=outbase+'_parang.sav', angles_arr, dits_arr
                
                    print, 'Saved nod cube for ', prev_flag, ' as nod ', nod_counter, ' group ', group_counter
                endif
                
                ; Start new group within same nod
                nod_cube = list()
                nod_cr_counts = list()
                nod_angles = list()
                nod_dits = list()
                group_counter += 1
                frames_in_current_group = 0
            endif
            
            ; Reset for next coadd group
            coadd_frame = list()
            cr_count_frame = list()
            coadd_frame_initialized = 0
            
        endif
        if verbose eq 1 then print, 'Number of Frames in Current Group: ', nod_cube.Count(), ', Group: ', group_counter
        
        ; Update prev_flag for next iteration (only if not already updated by nod switch)
        if not nod_switch_detected then prev_flag = flag
        
        k += 1 ; Only increment k for science frames (for coadding logic)
    endelse

endfor; ii = start_frame for loop

;write any remaining science frames (only in full processing mode)
if nod_cube.Count() gt 0 then begin
    if keyword_set(nod_counter_add) then BEGIN
        outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter+nod_counter_add, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
    endif else BEGIN
        outbase = output_path + strcompress(object_name + '_' + prev_flag + '_nod' + string(nod_counter, format='(I02)') + '_grp' + string(group_counter, format='(I02)'), /r)
    endelse

    nod_cube_arr = (temporary(nod_cube)).toArray(/TRANSPOSE, /NO_COPY)
    angles_arr = (temporary(nod_angles)).toArray(/NO_COPY)
    dits_arr = (temporary(nod_dits)).toArray(/NO_COPY)

    ; Write with header preservation
    writefits, outbase+'_cube.fits', nod_cube_arr, head
    
    ; Only write CR counts if not in CDS mode
    if not is_cds_mode then begin
        nod_cr_arr = (temporary(nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
        writefits, outbase+'_cr_count.fits', nod_cr_arr
    endif
    
    save, filename=outbase+'_parang.sav', angles_arr, dits_arr

    print, 'Final nod cube saved for ', prev_flag, ' as nod ', nod_counter, ' group ', group_counter
endif

; Write dark frames if any were collected
if dark_cube_list.Count() gt 0 then begin
    print, 'Writing dark frame cube with ', dark_cube_list.Count(), ' frames...'
    dark_outbase = output_path + strcompress(object_name + '_darks', /r)
    
    dark_cube_arr = (temporary(dark_cube_list)).toArray(/TRANSPOSE, /NO_COPY)
    
    ; Write with header preservation
    writefits, dark_outbase+'.fits', dark_cube_arr, head
    
    ; Only write CR counts if not in CDS mode
    if not is_cds_mode then begin
        dark_cr_arr = (temporary(dark_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
        writefits, dark_outbase+'_cr_count.fits', dark_cr_arr
    endif
    
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

; Write first nod frames if del_nodframe1 was set and frames were collected
if del_nodframe1 eq 1 and first_nod_frames.Count() gt 0 then begin
    print, 'Writing first nod frames cube with ', first_nod_frames.Count(), ' frames...'
    first_frames_outbase = output_path + strcompress(object_name + '_first_nod_frames', /r)
    
    first_frames_arr = (temporary(first_nod_frames)).toArray(/TRANSPOSE, /NO_COPY)
    
    ; Write with header preservation
    writefits, first_frames_outbase+'.fits', first_frames_arr, head
    
    ; Only write CR counts if not in CDS mode
    if not is_cds_mode then begin
        first_cr_arr = (temporary(first_nod_cr_counts)).toArray(/TRANSPOSE, /NO_COPY)
        writefits, first_frames_outbase+'_cr_count.fits', first_cr_arr
    endif
    
    if read_header eq 1 then begin
        first_angles_arr = (temporary(first_nod_angles)).toArray(/NO_COPY)
        first_dits_arr = (temporary(first_nod_dits)).toArray(/NO_COPY)
        first_flags_arr = (temporary(first_nod_flags)).toArray(/NO_COPY)
        save, filename=first_frames_outbase+'_parang.sav', first_angles_arr, first_dits_arr, first_flags_arr
    endif
    
    print, 'First nod frames cube saved as: ', first_frames_outbase, '.fits'
endif else if del_nodframe1 eq 1 then begin
    print, 'No first nod frames collected (del_nodframe1 was set but no frames found).'
endif

print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end
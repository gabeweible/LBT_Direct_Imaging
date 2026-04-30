pro bad_pixels_fast, output_folder, obj_name, $
    run=run, do_destripe=do_destripe, $
    just_destripe=just_destripe, debug=debug, $
    create_master_masks=create_master_masks, create_badpix_mask=create_badpix_mask, $
    use_flat_corrected=use_flat_corrected, darks_filename=darks_filename, $
    coadd=coadd, vapp=vapp, destripe_skysub=destripe_skysub, $
    pca_skysub=pca_skysub, post_pca_crop=post_pca_crop, start_nod=start_nod, $
    skip_badpix_mask=skip_badpix_mask, cube_indices=cube_indices, ncomp=ncomp

; NOTE: previously-positional `stripe` and the keywords `pre_sky_sub`,
; `bad_px_arr`, `boxhsize`, `fwhm`, `type` were unused and have been removed.
; Callers passing the positional `stripe` must drop it; callers passing any
; of the removed keywords will get an "unknown keyword" error.

compile_opt idl2, logical_predicate, strictarrsubs
newline = string(10B)

; Set default for start_nod if not provided. NOTE: must use `~` (logical not),
; not `not` (bitwise not). Under compile_opt logical_predicate, `if not 1 then`
; evaluates as `if -2 then` which is non-zero/truthy, so `if not keyword_set(...)`
; would ALWAYS fire and stomp the passed value back to the default.
if ~keyword_set(start_nod) then start_nod = 0

; Default PCA component count for the post-PCA cube filename glob.
; Accept either integer (ncomp=50) or string (ncomp='50') from caller.
if n_elements(ncomp) eq 0 then ncomp = 30
pca_suffix = '_' + string(long(ncomp[0]), format='(I0)') + 'comp.fits'

; cube_indices: optional list of nod NUMBERS (the XX in 'nodXX'). If set,
; only files whose nod number appears in cube_indices are processed. Missing
; nod numbers produce a warning and are skipped. start_nod is applied first
; so the effective set is (nods >= start_nod) AND (in cube_indices).

print, 'skip_badpix_mask: ', skip_badpix_mask, newline

; ---------- per-object parameter blocks (unchanged) ----------
if vapp eq 1 then begin
    boxehs = 60.
    g_boxehs = boxehs/5.
    owa_boxehs_in = 216
    owa_boxehs_out = 277
endif

if obj_name eq 'Alcor' then begin
    hc1_x1 = 385.-1.
    hc1_y1 = 250.-1.
    hc1_y2 = 278.-1.
    hc1_f1 = 0 * (20/coadd)
    hc1_f2 = 1545 * (20/coadd)

    hc2_x1 = 1281.-1.
    hc2_y1 = 210.-1.
    hc2_y2 = 235.-1.
    hc2_f1 = 1550 * (20/coadd)
    hc2_f2 = 2849 * (20/coadd)

    frame_max = 6550
    frame_min = 2000
    frame_max_post_destripe = 4000
    frame_min_post_destripe = -130
    hot_sigma = 2.0
    flat_thresh = [0.48, 1.15]
    max_sources = 2
    min_separation = 100.0
    min_peak_significance = 3.0
    max_cluster_size = 500000
    nantr = 15.
endif

if obj_name eq 'tyc5709' then begin
    boxehs = 80.
    hot_sigma = 2.0
    flat_thresh = [0.73, 1.14]
    frame_max = 1400
    frame_min = 665
    frame_max_post_destripe = 575
    frame_min_post_destripe = -8.25
    max_sources = 4
    min_separation = 100.0
    min_peak_significance = 3.0
    max_cluster_size = 500000
    nantr = 15.
endif

if obj_name eq 'HIP17034' then begin
    boxehs = 80.
    hot_sigma = 1.75
    flat_thresh = [0.73, 1.22]
    frame_max = 3100
    frame_min = 665
    frame_max_post_destripe = 2350
    frame_min_post_destripe = -14.
    max_sources = 8
    min_separation = 100.0
    min_peak_significance = 1.5
    max_cluster_size = 500000
    nantr = 10.
endif

if obj_name eq 'HIP17900' then begin
    boxehs = 85.
    hot_sigma = 1.75
    flat_thresh = [0.73, 1.16]
    frame_max = 6100.
    frame_min = -6100.
    frame_max_post_destripe = 6100.
    frame_min_post_destripe = -6100.
    max_sources = 4
    min_separation = 1000.0
    min_peak_significance = 4.0
    max_cluster_size = 500000
    nantr = 4.5
endif

if obj_name eq 'HIP16635' then begin
    boxehs = 85.
    hot_sigma = 1.5
    flat_thresh = [0.735, 1.16]
    frame_max = 1120
    frame_min = 590
    frame_max_post_destripe = 430.
    frame_min_post_destripe = -9.
    max_sources = 6
    min_separation = 500.0
    min_peak_significance = 3.0
    max_cluster_size = 500000
    nantr = 10.0
endif

new_filter = 27.

; ---------- debug folder (created lazily, only when debug eq 1) ----------
if debug eq 1 then begin
    debug_folder = output_folder + 'debug/'
    if ~file_test(debug_folder, /directory) then file_mkdir, debug_folder
endif

; ---------- helper: classify a file as NOD_A or NOD_B from filename ----------
; For non-PCA files the substring '_NOD_A_' / '_NOD_B_' is authoritative.
; For PCA cubes (filename: test_pca_skysub_cube_nodNN_<ncomp>comp.fits) there is
; no A/B substring, so we use the parity of the actual nod number from the
; filename. Even nod number => A, odd => B. This is read from the filename
; via extract_nod_number, NOT from the file-list position, so missing cubes
; cannot mislabel subsequent ones.
;
; (Inlined below where needed — kept as comment for clarity.)


if just_destripe ne 1 then begin

    ; ===========================================================================
    ; MAIN PATH: bad-pixel correction (and optionally destripe)
    ; ===========================================================================

    print, 'Reading in cubes for bad-px correction'
    print, 'Starting from nod number: ', start_nod

    if keyword_set(use_flat_corrected) then begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_grp??_corrected_cube.fits'
    endif else if keyword_set(destripe_skysub) and (pca_skysub eq 0) then begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_grp??_skysub_cube.fits'
    endif else if keyword_set(destripe_skysub) and (pca_skysub eq 1) then begin
        search_pattern = output_folder + 'test_pca_skysub_cube_nod??' + pca_suffix
        print, 'searching: ', search_pattern
    endif else begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_grp??_cube.fits'
    endelse

    nod_files = select_nod_cubes(search_pattern, $
                                 start_nod=start_nod, cube_indices=cube_indices, $
                                 nod_numbers=nod_numbers, count=n_nod_files, $
                                 /verbose)

    if n_nod_files eq 0 then begin
        print, 'ERROR: No NOD files found (after start_nod / cube_indices filtering).'
        return
    endif

    print, 'Found ', n_nod_files, ' NOD files (>= nod', string(start_nod, format='(I02)'), ')'
    print, 'Files sorted by nod number:'
    for i = 0, n_nod_files-1 do print, '  ', file_basename(nod_files[i])

    ; ---------- classify NOD_A / NOD_B from filename ----------
    nod_a_files = []
    nod_b_files = []
    if pca_skysub eq 0 then begin
        for i = 0, n_nod_files-1 do begin
            if strpos(nod_files[i], '_NOD_A_') ne -1 then begin
                nod_a_files = [nod_a_files, nod_files[i]]
            endif else if strpos(nod_files[i], '_NOD_B_') ne -1 then begin
                nod_b_files = [nod_b_files, nod_files[i]]
            endif
        endfor
    endif else begin
        ; PCA cubes: even nod number => NOD_A, odd => NOD_B. Use the actual
        ; nod number from the filename, not the loop index.
        for i = 0, n_nod_files-1 do begin
            if (nod_numbers[i] mod 2) eq 0 then begin
                nod_a_files = [nod_a_files, nod_files[i]]
            endif else begin
                nod_b_files = [nod_b_files, nod_files[i]]
            endelse
        endfor
    endelse

    n_nod_a = n_elements(nod_a_files)
    n_nod_b = n_elements(nod_b_files)
    print, 'Found ', n_nod_a, ' NOD_A files and ', n_nod_b, ' NOD_B files'
    if (n_nod_a eq 0) or (n_nod_b eq 0) then begin
        print, 'ERROR: Missing one or both nodding positions'
        return
    endif

    first_cube = readfits_fast(nod_a_files[0])
    cube_size = size(first_cube, /dimensions)
    x_dim = cube_size[0]
    y_dim = cube_size[1]
    frames_per_nod = cube_size[2]
    delvar, first_cube
    print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'

    ; ---------- master dark + flat construction (unchanged logic, with safer where()) ----------
    if create_master_masks eq 1 then begin
        print, 'Generating master dark...'
        darks = readfits_fast(darks_filename)
        resistant_mean, double(darks), 3.5, master_dark, dim=3, /double
        delvar, darks
        master_dark = float(master_dark)

        print, 'Saving master dark to FITS file...'
        writefits, output_folder + obj_name + '_master_dark.fits', master_dark
        print, 'Master dark FITS saved.'

        print, 'Creating outlier-resistant master flat fields for each nod position...'

        print, 'Creating outlier-resistant master flat for NOD_A...'
        nod_a_individual_flats = []
        for i = 0, n_nod_a-1 do begin
            print, 'Processing NOD_A file ', i+1, '/', n_nod_a, ' for individual flat'
            nod_cube = readfits_fast(nod_a_files[i])
            for kk = 0, frames_per_nod-1 do nod_cube[*,*,kk] -= master_dark

            normalized_frames = fltarr(x_dim, y_dim, frames_per_nod)
            valid_frame_count = 0
            for kk = 0, frames_per_nod-1 do begin
                frame = nod_cube[*,*,kk]
                frame_median = float(median(double(frame), /even, /double))
                if frame_median gt 0 then begin
                    normalized_frames[*,*,valid_frame_count] = frame / frame_median
                    valid_frame_count += 1
                endif
            endfor

            if valid_frame_count gt 0 then begin
                normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
                individual_nod_flat = dblarr(x_dim, y_dim)
                resistant_mean, double(normalized_frames), 3.5, individual_nod_flat, dim=3, /double
                individual_median = float(median(individual_nod_flat, /even, /double))
                individual_nod_flat = float(individual_nod_flat)
                if individual_median gt 0 then individual_nod_flat /= individual_median

                if n_elements(nod_a_individual_flats) eq 0 then $
                    nod_a_individual_flats = fltarr(x_dim, y_dim, n_nod_a)
                nod_a_individual_flats[*,*,i] = individual_nod_flat
            endif
            delvar, nod_cube, normalized_frames, individual_nod_flat
        endfor

        print, 'Combining individual NOD_A flats with outlier resistance...'
        nod_a_flat = dblarr(x_dim, y_dim)
        resistant_mean, double(nod_a_individual_flats), 3.5, nod_a_flat, dim=3, /double
        nod_a_median = float(median(nod_a_flat, /even, /double))
        nod_a_flat = float(nod_a_flat)
        if nod_a_median gt 0 then nod_a_flat /= nod_a_median
        delvar, nod_a_individual_flats

        print, 'Creating outlier-resistant master flat for NOD_B...'
        nod_b_individual_flats = []
        for i = 0, n_nod_b-1 do begin
            print, 'Processing NOD_B file ', i+1, '/', n_nod_b, ' for individual flat'
            nod_cube = readfits_fast(nod_b_files[i])
            for kk = 0, frames_per_nod-1 do nod_cube[*,*,kk] -= master_dark

            normalized_frames = fltarr(x_dim, y_dim, frames_per_nod)
            valid_frame_count = 0
            for kk = 0, frames_per_nod-1 do begin
                frame = nod_cube[*,*,kk]
                frame_median = float(median(double(frame), /even, /double))
                if frame_median gt 0 then begin
                    normalized_frames[*,*,valid_frame_count] = frame / frame_median
                    valid_frame_count += 1
                endif
            endfor

            if valid_frame_count gt 0 then begin
                normalized_frames = normalized_frames[*,*,0:valid_frame_count-1]
                individual_nod_flat = dblarr(x_dim, y_dim)
                resistant_mean, double(normalized_frames), 3.5, individual_nod_flat, dim=3, /double
                individual_median = float(median(individual_nod_flat, /even, /double))
                individual_nod_flat = float(individual_nod_flat)
                if individual_median gt 0 then individual_nod_flat /= individual_median

                if n_elements(nod_b_individual_flats) eq 0 then $
                    nod_b_individual_flats = fltarr(x_dim, y_dim, n_nod_b)
                nod_b_individual_flats[*,*,i] = individual_nod_flat
            endif
            delvar, nod_cube, normalized_frames, individual_nod_flat
        endfor

        print, 'Combining individual NOD_B flats with outlier resistance...'
        nod_b_flat = dblarr(x_dim, y_dim)
        resistant_mean, double(nod_b_individual_flats), 3.5, nod_b_flat, dim=3, /double
        nod_b_median = float(median(nod_b_flat, /even, /double))
        nod_b_flat = float(nod_b_flat)
        if nod_b_median gt 0 then nod_b_flat /= nod_b_median
        delvar, nod_b_individual_flats

        writefits, output_folder + obj_name + '_nod_a_master_flat_resistant.fits', nod_a_flat
        writefits, output_folder + obj_name + '_nod_b_master_flat_resistant.fits', nod_b_flat

        super_flat = nod_a_flat
        if obj_name eq 'Alcor' then super_flat[0:897, *] = nod_b_flat[0:897, *]
        if obj_name eq 'tyc5709' then super_flat[0:970, *] = nod_b_flat[0:970, *]
        if (obj_name eq 'HIP17034') or (obj_name eq 'HIP17900') or (obj_name eq 'HIP16635') then $
            super_flat[0:x_dim/2-1, *] = nod_b_flat[0:x_dim/2-1, *]

        delvar, nod_a_flat, nod_b_flat
        writefits, output_folder + obj_name + '_super_master_flat_resistant_raw.fits', super_flat

        print, 'Final normalization of stitched master flat...'
        super_median = float(median(double(super_flat), /even, /double))
        if super_median gt 0 then super_flat /= super_median
        writefits, output_folder + obj_name + '_super_master_flat_resistant_normalized.fits', super_flat
        print, 'Outlier-resistant master flat correction complete.'

        if create_badpix_mask eq 1 then begin
            print, 'Creating master bad-pixel mask from master dark and outlier-resistant master flat'
            gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask, $
                       hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
            writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
        endif

    endif else begin
        print, 'Reading in pre-generated master files...'
        resistant_flat_file = output_folder + obj_name + '_super_master_flat_resistant_normalized.fits'
        original_flat_file  = output_folder + obj_name + '_super_master_flat_normalized.fits'

        if file_test(resistant_flat_file) then begin
            print, 'Reading outlier-resistant master flat...'
            super_flat = readfits_fast(resistant_flat_file)
        endif else begin
            print, 'Outlier-resistant flat not found, reading original flat...'
            super_flat = readfits_fast(original_flat_file)
        endelse
        master_dark = readfits_fast(output_folder + obj_name + '_master_dark.fits')

        if create_badpix_mask eq 1 then begin
            print, 'Creating master bad-pixel mask from master dark and master flat'
            gen_bad_pix_mask, master_dark, super_flat, mask_out=master_badpix_mask, $
                               hot_sigma=hot_sigma, flat_thresh=flat_thresh, debug=debug
            writefits, output_folder+obj_name+'_master_badpix_mask.fits', master_badpix_mask
        endif else begin
            print, 'Reading-in existing bad pixel mask...'
            master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask_complete.fits')

            if skip_badpix_mask eq 1 then begin
                print, "Replacing master badpix mask with all 1b's (i.e., assuming all pixels are good)"
                master_badpix_mask[*,*] = 1b
                print, "skip_badpix_mask: Replacing master flat with all 1's"
                super_flat[*,*] = double(1.0)
                print, "skip_badpix_mask: Replacing master dark with all 0's"
                master_dark[*,*] = double(0.0)
            endif
        endelse
    endelse

    if skip_badpix_mask eq 0 then begin
        print, '=========================================='
        print, 'Building complete bad pixel mask...'
        print, '=========================================='

        if obj_name eq 'Alcor' then begin
            master_badpix_mask[hc1_x1, hc1_y1:hc1_y2] = 0
            master_badpix_mask[hc2_x1, hc2_y1:hc2_y2] = 0
            master_badpix_mask[931, *] = 0
        endif

        ; Always mask edges (4 px) and the bad readout channel + bad column.
        master_badpix_mask[0:4, *] = 0
        master_badpix_mask[x_dim-4:x_dim-1, *] = 0
        master_badpix_mask[*, y_dim-4:y_dim-1] = 0
        master_badpix_mask[*, 0:4] = 0
        master_badpix_mask[896:959, *] = 0
        master_badpix_mask[1984, *] = 0

        if run eq 5 then begin
            x = rebin(findgen(x_dim), x_dim, y_dim)
            y = rebin(reform(findgen(y_dim), 1, y_dim), x_dim, y_dim)

            y1 = (-105.0 / 2047.0)* x + 1130.0
            y2 = ( 460.0 / 2047.0)*x +  840.0
            y4 = (305.0 / 2047.0)*x + 810.0
            master_badpix_mask[89, 860] = 0

            y5 = 0.0*x + 983.0
            y8 = 0.0*x + 1105.0
            y3 = (173.0 / 1107.0)*x - 154.0
            y7 = (272.0 / 387.0)*x - 1150.0
            y10 = (79.0 / 1396.0)*x + 14.0
            y6 = (123.0 / 176.0)*x + 1820.0
            y9 = (-155.0 / 192.0)*x + 3455.0
            y11 = 0.0*x + 1983.0

            cs_bads = where(((y gt y2) and (y lt y1)) or $
                            ((y gt y1) and (y lt y2)) or (y lt y3) or $
                            ((y gt y4) and (y lt y2)) or $
                            ((y gt y5) and (y lt y8)) or (y gt y6) or $
                            (y lt y7) or (y gt y9) or (y gt y11) or $
                            (y lt y10), n_bad)
            if n_bad gt 0 then master_badpix_mask[cs_bads] = 0
        endif

        if obj_name eq 'tyc5709' then begin
            master_badpix_mask[1408, 228-10:244+10] = 0
            master_badpix_mask[*, 434:536] = 0
        endif

        if obj_name eq 'HIP17034' then master_badpix_mask[*, 870:1100] = 0

        if obj_name eq 'HIP17900' then begin
            master_badpix_mask[1408, 1500:1570] = 0
            master_badpix_mask[1344:1345, 470:540] = 0
            master_badpix_mask[447:448, 440:510] = 0
            master_badpix_mask[448, 1520:1590] = 0
            master_badpix_mask[575, 1525:1595] = 0
        endif

        if obj_name eq 'HIP16635' then begin
            ; (most "bad columns" only appear in the first nod groups; keep
            ; only the persistent ones — see notes in repo history.)
            master_badpix_mask[1876, *] = 0
            master_badpix_mask[1912, *] = 0
            master_badpix_mask[1301:1345, 681:730] = 0   ; bottom-right dipole
            master_badpix_mask[448, 1572:1601]    = 0    ; NOD_A top crosstalk
        endif

        print, 'Saving complete bad pixel mask...'
        writefits, output_folder+obj_name+'_master_badpix_mask_complete.fits', master_badpix_mask
    endif

    n_bad_pix = long(total(master_badpix_mask eq 0))
    print, 'Complete bad pixel mask created with ', n_bad_pix, ' bad pixels (', $
           string((100.* n_bad_pix) / (x_dim * y_dim), format='(F5.2)'), '% of each frame)'
    print, '=========================================='

    ; ===========================================================================
    ; PROCESS EACH NOD FILE (bad-pixel correction; optional destripe)
    ; ===========================================================================
    print, 'Processing NOD files individually for bad pixel correction and flat fielding...'

    total_frame_count = 0L
    ; seed intentionally undefined: randomn() will initialize from the system
    ; clock on its first call and advance the sequence on subsequent calls,
    ; producing non-reproducible noise across runs.

    for i = 0, n_nod_files-1 do begin
        print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])

        nod_cube = readfits_fast(nod_files[i])
        nod_frames = (size(nod_cube, /dimensions))[2]
        working_mask = master_badpix_mask

        print, '  Applying bad pixel correction to ', nod_frames, ' frames'

        for kk = 0, nod_frames-1 do begin
            if (kk mod 100) eq 0 then print, '    Processing frame ', kk+1, '/', nod_frames

            frame = nod_cube[*,*,kk]
            frame -= master_dark
            frame /= super_flat
            frame_mask = working_mask
            original_frame = frame

            if debug eq 1 and (i mod 100) eq 0 and (kk eq 2) then $
                writefits, debug_folder + 'original_frame_test.fits', original_frame

            ; --- mask bad/zero/out-of-range pixels (count-checked) ---
            idx = where(frame_mask eq 0, n) & if n gt 0 then frame[idx] = !values.f_nan
            idx = where(frame eq 0,        n) & if n gt 0 then frame[idx] = !values.f_nan
            idx = where(frame gt frame_max, n) & if n gt 0 then frame[idx] = !values.f_nan
            idx = where(frame lt frame_min, n) & if n gt 0 then frame[idx] = !values.f_nan

            original_nan_mask = ~finite(frame)

            ; replace small NaN clumps with weighted-24-NN interpolation
            efficient_nan_correction, frame, fixed_frame, $
                max_cluster_size=max_cluster_size, npix=24, /weight, /silent

            idx = where(fixed_frame gt frame_max, n) & if n gt 0 then fixed_frame[idx] = !values.f_nan
            idx = where(fixed_frame lt frame_min, n) & if n gt 0 then fixed_frame[idx] = !values.f_nan

            interpolated_mask = (original_nan_mask and finite(fixed_frame))
            good_pixels = where(finite(original_frame) and frame_mask eq 1, n_good)

            if n_good gt 25000 then begin
                good_values = original_frame[good_pixels]
                median_val = median(double(good_values), /even, /double)
                mad_val = float(median(abs(double(good_values) - median_val), /even, /double))
                median_val = float(median_val)
                noise_sigma = 1.4826 * mad_val

                interpolated_indices = where(interpolated_mask, n_interpolated)
                if n_interpolated gt 0 then begin
                    gaussian_noise = randomn(seed, n_interpolated) * noise_sigma
                    fixed_frame[interpolated_indices] += gaussian_noise
                endif

                if debug eq 1 and kk eq 2 and (i mod 100) eq 0 then begin
                    print, '    Noise sigma estimate: ', noise_sigma
                    print, '    Number of interpolated pixels: ', n_interpolated
                endif
            endif else begin
                if debug eq 1 and kk eq 2 and (i mod 100) eq 0 then $
                    print, '    Warning: Insufficient good pixels for noise estimation'
            endelse

            if debug eq 1 and kk eq 2 and (i mod 100) eq 0 then writefits, $
                debug_folder + strcompress('nod_'+string(i)+'_post_px_fix.fits', /r), $
                float(fixed_frame)

            nod_cube[*,*,kk] = fixed_frame
        endfor

        total_frame_count += nod_frames

        ; --- optional destripe (Block 1 of the original two destriping blocks) ---
        if (do_destripe eq 1) and ((obj_name eq 'HIP17034') or (obj_name eq 'HIP17900') or (obj_name eq 'HIP16635')) then begin
            print, '  Applying destriping to NOD cube...'
            x_arr = indgen(x_dim) # replicate(1, y_dim)
            y_arr = replicate(1, x_dim) # indgen(y_dim)

            for kk = 0, nod_frames-1 do begin
                if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
                frame_ii = nod_cube[*,*,kk]
                frame_ii_med = float(median(double(frame_ii), /double, /even))

                frame_ii_cen = frame_ii
                frame_ii_cen = frame_ii_cen > frame_ii_med
                frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)

                fin_idx = where(finite(frame_ii_cen), n_fin_cen)
                if n_fin_cen gt 0 then begin
                    noise_stats = frame_ii_cen[fin_idx]
                    noise_mean = float(mean(double(noise_stats), /double))
                    noise_stddev = 1.4826 * float(median(abs(double(noise_stats) - $
                                    median(double(noise_stats), /double, /even)), /double, /even))
                endif else begin
                    noise_mean = 0.0
                    noise_stddev = 1.0
                endelse
                min_peak_value = noise_mean + min_peak_significance * noise_stddev

                if debug eq 1 and (kk mod 100) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_cen_test_kk'+string(kk)+'.fits', /r), $
                    float(frame_ii_cen)

                peak_x = []
                peak_y = []
                peak_type = []

                temp_frame_pos = frame_ii_cen
                for source_num = 0, max_sources-1 do begin
                    mx = MAX(temp_frame_pos, location, /nan)
                    if ~finite(mx) or mx lt min_peak_value then break
                    ind = ARRAY_INDICES(temp_frame_pos, location)
                    xx = ind[0] & yy = ind[1]

                    valid_peak = 1
                    if n_elements(peak_x) gt 0 then begin
                        for prev_peak = 0, n_elements(peak_x)-1 do begin
                            sep = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
                            if sep lt min_separation then begin
                                valid_peak = 0
                                break
                            endif
                        endfor
                    endif

                    if valid_peak then begin
                        peak_x = [peak_x, xx]
                        peak_y = [peak_y, yy]
                        peak_type = [peak_type, 1]
                        if debug eq 1 and kk eq 2 then $
                            print, '    Found positive peak ', n_elements(peak_x), $
                                   ' at (', xx, ',', yy, ') with value ', mx
                    endif

                    rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                    mp = where(rad_temp le boxehs, n_mp)
                    if n_mp gt 0 then temp_frame_pos[mp] = !values.f_nan
                endfor

                frame_ii_neg = -1.0 * (frame_ii - frame_ii_med)
                frame_ii_neg = frame_ii_neg > 0
                frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
                min_neg_peak_value = min_peak_significance * noise_stddev

                temp_frame_neg = frame_ii_neg
                for source_num = 0, max_sources-1 do begin
                    mx = MAX(temp_frame_neg, location, /nan)
                    if ~finite(mx) or mx lt min_neg_peak_value then break
                    ind = ARRAY_INDICES(temp_frame_neg, location)
                    xx = ind[0] & yy = ind[1]

                    valid_peak = 1
                    if n_elements(peak_x) gt 0 then begin
                        for prev_peak = 0, n_elements(peak_x)-1 do begin
                            sep = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
                            if sep lt min_separation then begin
                                valid_peak = 0
                                break
                            endif
                        endfor
                    endif

                    if valid_peak then begin
                        peak_x = [peak_x, xx]
                        peak_y = [peak_y, yy]
                        peak_type = [peak_type, -1]
                        if debug eq 1 and kk eq 2 then begin
                            actual_neg_value = frame_ii[xx, yy]
                            print, '    Found negative peak ', n_elements(peak_x), $
                                   ' at (', xx, ',', yy, ') with value ', actual_neg_value
                        endif
                    endif

                    rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                    mp = where(rad_temp le boxehs, n_mp)
                    if n_mp gt 0 then temp_frame_neg[mp] = !values.f_nan
                endfor

                n_sources = n_elements(peak_x)
                if n_sources gt 0 then begin
                    pos_indices = where(peak_type eq 1, n_pos_sources)
                    neg_indices = where(peak_type eq -1, n_neg_sources)
                endif else begin
                    n_pos_sources = 0L
                    n_neg_sources = 0L
                    pos_indices = !null
                    neg_indices = !null
                endelse

                if debug eq 1 and kk eq 2 then begin
                    print, '    Total sources found in frame ', kk+1, ': ', n_sources
                    print, '    Positive sources: ', n_pos_sources, ', Negative sources: ', n_neg_sources
                endif

                mask = bytarr(x_dim, y_dim)
                if n_sources gt 0 then begin
                    for source_num = 0, n_sources-1 do begin
                        rad_pos = sqrt((x_arr - peak_x[source_num])^2 + (y_arr - peak_y[source_num])^2)
                        mask = mask OR (rad_pos le boxehs)
                        if debug eq 1 and kk eq 2 then begin
                            stype = peak_type[source_num] eq 1 ? 'positive' : 'negative'
                            print, '    Masking ', stype, ' source at (', $
                                   peak_x[source_num], ',', peak_y[source_num], ')'
                        endif
                    endfor
                endif else if debug eq 1 then begin
                    print, '    Warning: No significant sources found in frame ', kk+1
                endif

                ; (mask_pos / mask_neg constructed below for any downstream use)
                if n_pos_sources gt 0 then begin
                    mask_pos = bytarr(x_dim, y_dim)
                    for j = 0, n_pos_sources-1 do begin
                        idx2 = pos_indices[j]
                        rad_pos = sqrt((x_arr - peak_x[idx2])^2 + (y_arr - peak_y[idx2])^2)
                        mask_pos = mask_pos OR (rad_pos le boxehs)
                    endfor
                endif
                if n_neg_sources gt 0 then begin
                    mask_neg = bytarr(x_dim, y_dim)
                    for j = 0, n_neg_sources-1 do begin
                        idx2 = neg_indices[j]
                        rad_pos = sqrt((x_arr - peak_x[idx2])^2 + (y_arr - peak_y[idx2])^2)
                        mask_neg = mask_neg OR (rad_pos le boxehs)
                    endfor
                endif

                frame_for_stripes = frame_ii
                masked_pixels = where(mask, n_masked)
                if n_masked gt 0 then frame_for_stripes[masked_pixels] = !values.f_nan

                frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                frame_ii_stddev = 1.4826 * float(median(abs(double(frame_for_stripes) - $
                                    median(double(frame_for_stripes), /double, /even)), /double, /even))

                low_nantr_thresh  = frame_ii_median - nantr*frame_ii_stddev
                high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                if debug eq 1 and (kk mod 100) eq 0 then begin
                    print, 'low_nantr_thresh:  ', low_nantr_thresh
                    print, 'high_nantr_thresh: ', high_nantr_thresh, newline
                endif

                value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                             (frame_for_stripes gt high_nantr_thresh)
                nanmask = where(value_mask, n_nm)
                if n_nm gt 0 then begin
                    frame_for_stripes[nanmask] = !values.f_nan
                    frame_ii[nanmask] = !values.f_nan
                endif

                if debug eq 1 and (kk mod 100) eq 0 then begin
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_masked_test_kk'+string(kk)+'.fits', /r), $
                        frame_for_stripes
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_frame_ii_kk'+string(kk)+'.fits', /r), $
                        frame_ii
                endif

                if (kk mod 100) eq 0 then begin
                    n_finite_fs  = long(total(finite(frame_for_stripes)))
                    n_finite_fii = long(total(finite(frame_ii)))
                    print, '  Stripe stats nod ', string(i), ', frame ', string(kk), ':'
                    print, '  frame_for_stripes finite pixels: ', n_finite_fs
                    print, '  frame_ii finite pixels:          ', n_finite_fii
                    print, '  frame_for_stripes median: ', median(frame_for_stripes, /double)
                    print, '  frame_ii median:          ', median(frame_ii, /double), newline
                endif

                ; --- Vertical destripe with quadrant split ---
                frame_for_stripes_top    = frame_for_stripes
                frame_for_stripes_bottom = frame_for_stripes
                frame_for_stripes_top[*,    0:y_dim/2-1]    = !values.f_nan
                frame_for_stripes_bottom[*, y_dim/2:y_dim-1] = !values.f_nan

                if debug eq 1 and kk eq 2 then begin
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_masked_test_top.fits',    /r), float(frame_for_stripes_top)
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_masked_test_bottom.fits', /r), float(frame_for_stripes_bottom)
                endif

                destriped_frame_top    = destripe(frame_for_stripes_top,    frame_ii, 90., fraction=0.02, /no_fit, /nodisp)
                idx = where(destriped_frame_top gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame_top[idx] = !values.f_nan
                idx = where(destriped_frame_top lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame_top[idx] = !values.f_nan

                destriped_frame_bottom = destripe(frame_for_stripes_bottom, frame_ii, 90., fraction=0.02, /no_fit, /nodisp)
                idx = where(destriped_frame_bottom gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame_bottom[idx] = !values.f_nan
                idx = where(destriped_frame_bottom lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame_bottom[idx] = !values.f_nan

                destriped_frame = frame_ii
                destriped_frame[*, 0:y_dim/2-1]      = destriped_frame_bottom[*, 0:y_dim/2-1]
                destriped_frame[*, y_dim/2:y_dim-1]  = destriped_frame_top[*,    y_dim/2:y_dim-1]

                if debug eq 1 and (kk mod 100) eq 0 then begin
                    n_finite_df = long(total(finite(destriped_frame)))
                    n_finite_fi = long(total(finite(frame_ii)))
                    print, '  After vertical destripe (90deg) with quadrant split:'
                    print, '  destriped_frame finite pixels: ', n_finite_df
                    print, '  destriped_frame median:        ', median(destriped_frame, /double)
                    print, '  Number of NaNs added:          ', $
                           (long(x_dim)*long(y_dim) - n_finite_df) - (long(x_dim)*long(y_dim) - n_finite_fi), newline
                endif

                ; --- Horizontal destripe ---
                frame_ii = destriped_frame
                frame_for_stripes = frame_ii
                if n_masked gt 0 then frame_for_stripes[masked_pixels] = !values.f_nan

                frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                frame_ii_stddev = 1.4826 * float(median(abs(double(frame_for_stripes) - $
                                    median(double(frame_for_stripes), /double, /even)), /double, /even))

                low_nantr_thresh  = frame_ii_median - nantr*frame_ii_stddev
                high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                if debug eq 1 and (kk mod 100) eq 0 then begin
                    print, 'low_nantr_thresh:  ', low_nantr_thresh
                    print, 'high_nantr_thresh: ', high_nantr_thresh, newline
                endif

                value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                             (frame_for_stripes gt high_nantr_thresh)
                nanmask = where(value_mask, n_nm)
                if n_nm gt 0 then begin
                    frame_for_stripes[nanmask] = !values.f_nan
                    frame_ii[nanmask]          = !values.f_nan
                endif

                destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, /no_fit, /nodisp)
                idx = where(destriped_frame gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame[idx] = !values.f_nan
                idx = where(destriped_frame lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame[idx] = !values.f_nan

                nod_cube[*,*,kk] = destriped_frame

                if debug eq 1 and (kk mod 100) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_destriped_test_kk'+string(kk)+'.fits', /r), $
                    destriped_frame
            endfor
        endif

        ; --- write out the corrected cube ---
        current_file_basename = file_basename(nod_files[i])
        print, 'Debug: Processing file: ', current_file_basename

        suffix_pos = strpos(current_file_basename, '_cube.fits')
        if suffix_pos ne -1 then begin
            current_base_name = strmid(current_file_basename, 0, suffix_pos)
        endif else begin
            suffix_pos = strpos(current_file_basename, '.fits')
            if suffix_pos ne -1 then begin
                current_base_name = strmid(current_file_basename, 0, suffix_pos)
            endif else begin
                current_base_name = current_file_basename
            endelse
        endelse

        output_filename = output_folder + current_base_name + '_corrected_cube.fits'
        print, '  Saving corrected cube: ', file_basename(output_filename)
        writefits, output_filename, nod_cube
        delvar, nod_cube
    endfor

    print, 'All NOD cubes processed individually and saved.'
    print, 'Total frames processed: ', total_frame_count

    summary_file = output_folder + obj_name + '_corrected_cubes_list.txt'
    openw, lun, summary_file, /get_lun
    printf, lun, '# List of corrected cube files for object: ' + obj_name
    printf, lun, '# Generated on: ' + systime()
    printf, lun, '# Total NOD files processed: ' + string(n_nod_files)
    printf, lun, '# Total frames processed: ' + string(total_frame_count)
    if keyword_set(cube_indices) then begin
        printf, lun, '# cube_indices (nod numbers) requested: ' + $
                strjoin(string(cube_indices, format='(I02)'), ', ')
        printf, lun, '# cube_indices (nod numbers) processed: ' + $
                strjoin(string(nod_numbers, format='(I02)'), ', ')
    endif
    printf, lun, '#'
    printf, lun, '# Format: cube_file'
    for i = 0, n_nod_files-1 do begin
        base_name = file_basename(nod_files[i], '_cube.fits')
        printf, lun, base_name + '_corrected_cube.fits'
    endfor
    free_lun, lun
    print, 'Summary file created: ', file_basename(summary_file)

endif else begin

    ; ===========================================================================
    ; JUST-DESTRIPE PATH (e.g. post-PCA skysub)
    ; ===========================================================================

    print, 'Just destriping mode - processing individual NOD cubes...'
    print, 'Starting from nod number: ', start_nod

    print, 'Reading-in existing bad pixel mask...'
    master_badpix_mask = readfits_fast(output_folder+obj_name+'_master_badpix_mask_complete.fits')

    if keyword_set(use_flat_corrected) then begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_grp??_corrected_cube.fits'
    endif else if keyword_set(destripe_skysub) and (pca_skysub eq 0) then begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_grp??_skysub_cube.fits'
    endif else if keyword_set(destripe_skysub) and (pca_skysub eq 1) then begin
        search_pattern = output_folder + 'test_pca_skysub_cube_nod??' + pca_suffix
        print, 'searching: ', search_pattern
    endif else begin
        search_pattern = output_folder + obj_name + '_NOD_?_nod??_no_bad_px_pre_destripe_cube.fits'
    endelse

    nod_files = select_nod_cubes(search_pattern, $
                                 start_nod=start_nod, cube_indices=cube_indices, $
                                 nod_numbers=nod_numbers, count=n_nod_files, $
                                 /verbose)

    if n_nod_files eq 0 then begin
        print, 'ERROR: No NOD files found (after start_nod / cube_indices filtering).'
        return
    endif

    print, 'Found ', n_nod_files, ' NOD files (>= nod', string(start_nod, format='(I02)'), ')'
    print, 'Files sorted by nod number:'
    for i = 0, n_nod_files-1 do print, '  ', file_basename(nod_files[i])

    nod_a_files = []
    nod_b_files = []
    if pca_skysub eq 0 then begin
        for i = 0, n_nod_files-1 do begin
            if strpos(nod_files[i], '_NOD_A_') ne -1 then begin
                nod_a_files = [nod_a_files, nod_files[i]]
            endif else if strpos(nod_files[i], '_NOD_B_') ne -1 then begin
                nod_b_files = [nod_b_files, nod_files[i]]
            endif
        endfor
    endif else begin
        for i = 0, n_nod_files-1 do begin
            if (nod_numbers[i] mod 2) eq 0 then begin
                nod_a_files = [nod_a_files, nod_files[i]]
            endif else begin
                nod_b_files = [nod_b_files, nod_files[i]]
            endelse
        endfor
    endelse

    n_nod_a = n_elements(nod_a_files)
    n_nod_b = n_elements(nod_b_files)
    print, 'Found ', n_nod_a, ' NOD_A files and ', n_nod_b, ' NOD_B files'
    if (n_nod_a eq 0) or (n_nod_b eq 0) then begin
        print, 'ERROR: Missing one or both nodding positions'
        return
    endif

    first_cube = readfits_fast(nod_a_files[0])
    cube_size = size(first_cube, /dimensions)
    x_dim = cube_size[0]
    y_dim = cube_size[1]
    frames_per_nod = cube_size[2]
    delvar, first_cube
    print, 'Cube dimensions: ', x_dim, ' x ', y_dim, ' x ', frames_per_nod, ' per nod'

    if (do_destripe eq 1) and ((obj_name eq 'HIP17034') or (obj_name eq 'HIP17900') or (obj_name eq 'HIP16635')) then begin
        total_frame_count = 0L
        for i = 0, n_nod_files-1 do begin
            print, 'Processing NOD file ', i+1, '/', n_nod_files, ': ', file_basename(nod_files[i])

            nod_cube = readfits_fast(nod_files[i])
            nod_frames = (size(nod_cube, /dimensions))[2]
            print, 'Applying destriping to NOD cube...'

            x_arr = indgen(x_dim) # replicate(1, y_dim)
            y_arr = replicate(1, x_dim) # indgen(y_dim)

            for kk = 0, nod_frames-1 do begin
                if (kk mod 100) eq 0 then print, '    Destriping frame: ', kk+1, '/', nod_frames
                frame_ii = nod_cube[*,*,kk]
                frame_for_stripes = frame_ii

                if post_pca_crop eq 1 then begin
                    ; odd nod number => NOD_B (right): mask the LEFT half
                    ; even nod number => NOD_A (left):  mask the RIGHT half
                    if (nod_numbers[i] mod 2) eq 1 then begin
                        frame_for_stripes[0:x_dim/2-1, *]        = !values.f_nan
                    endif else begin
                        frame_for_stripes[x_dim/2:x_dim-1, *]  = !values.f_nan
                    endelse
                endif

                frame_ii_med = float(median(double(frame_for_stripes), /double, /even))

                frame_ii_cen = frame_for_stripes
                frame_ii_cen = frame_ii_cen > frame_ii_med
                frame_ii_cen = filter_image(temporary(frame_ii_cen), smooth=5.)

                fin_idx = where(finite(frame_ii_cen), n_fin_cen)
                if n_fin_cen gt 0 then begin
                    noise_stats = frame_ii_cen[fin_idx]
                    noise_mean = float(mean(double(noise_stats), /double))
                    noise_stddev = 1.4826 * float(median(abs(double(noise_stats) - $
                                    median(double(noise_stats), /double, /even)), /double, /even))
                endif else begin
                    noise_mean = 0.0
                    noise_stddev = 1.0
                endelse
                min_peak_value = noise_mean + min_peak_significance * noise_stddev

                if debug eq 1 and (kk mod 100) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_cen_test_kk'+string(kk)+'.fits', /r), $
                    frame_ii_cen

                peak_x = []
                peak_y = []
                peak_type = []

                temp_frame_pos = frame_ii_cen
                for source_num = 0, max_sources-1 do begin
                    mx = MAX(temp_frame_pos, location, /nan)
                    if ~finite(mx) or mx lt min_peak_value then break
                    ind = ARRAY_INDICES(temp_frame_pos, location)
                    xx = ind[0] & yy = ind[1]

                    valid_peak = 1
                    if n_elements(peak_x) gt 0 then begin
                        for prev_peak = 0, n_elements(peak_x)-1 do begin
                            sep = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
                            if sep lt min_separation then begin
                                valid_peak = 0
                                break
                            endif
                        endfor
                    endif

                    if valid_peak then begin
                        peak_x = [peak_x, xx]
                        peak_y = [peak_y, yy]
                        peak_type = [peak_type, 1]
                        if debug eq 1 and kk eq 2 then $
                            print, '    Found positive peak ', n_elements(peak_x), $
                                   ' at (', xx, ',', yy, ') with value ', mx
                    endif

                    rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                    mp = where(rad_temp le boxehs, n_mp)
                    if n_mp gt 0 then temp_frame_pos[mp] = !values.f_nan
                endfor

                frame_ii_neg = -1.0 * (frame_for_stripes - frame_ii_med)
                frame_ii_neg = frame_ii_neg > 0
                frame_ii_neg = filter_image(temporary(frame_ii_neg), smooth=5.)
                min_neg_peak_value = min_peak_significance * noise_stddev

                temp_frame_neg = frame_ii_neg
                for source_num = 0, max_sources-1 do begin
                    mx = MAX(temp_frame_neg, location, /nan)
                    if ~finite(mx) or mx lt min_neg_peak_value then break
                    ind = ARRAY_INDICES(temp_frame_neg, location)
                    xx = ind[0] & yy = ind[1]

                    valid_peak = 1
                    if n_elements(peak_x) gt 0 then begin
                        for prev_peak = 0, n_elements(peak_x)-1 do begin
                            sep = sqrt((xx - peak_x[prev_peak])^2 + (yy - peak_y[prev_peak])^2)
                            if sep lt min_separation then begin
                                valid_peak = 0
                                break
                            endif
                        endfor
                    endif

                    if valid_peak then begin
                        peak_x = [peak_x, xx]
                        peak_y = [peak_y, yy]
                        peak_type = [peak_type, -1]
                        if debug eq 1 and kk eq 2 then begin
                            actual_neg_value = frame_ii[xx, yy]
                            print, '    Found negative peak ', n_elements(peak_x), $
                                   ' at (', xx, ',', yy, ') with value ', actual_neg_value
                        endif
                    endif

                    rad_temp = sqrt((x_arr-xx)^2 + (y_arr-yy)^2)
                    mp = where(rad_temp le boxehs, n_mp)
                    if n_mp gt 0 then temp_frame_neg[mp] = !values.f_nan
                endfor

                n_sources = n_elements(peak_x)
                if n_sources gt 0 then begin
                    pos_indices = where(peak_type eq 1, n_pos_sources)
                    neg_indices = where(peak_type eq -1, n_neg_sources)
                endif else begin
                    n_pos_sources = 0L
                    n_neg_sources = 0L
                    pos_indices = !null
                    neg_indices = !null
                endelse

                if debug eq 1 and kk eq 2 then begin
                    print, '    Total sources found in frame ', kk+1, ': ', n_sources
                    print, '    Positive sources: ', n_pos_sources, ', Negative sources: ', n_neg_sources
                endif

                mask = bytarr(x_dim, y_dim)
                if n_sources gt 0 then begin
                    for source_num = 0, n_sources-1 do begin
                        rad_pos = sqrt((x_arr - peak_x[source_num])^2 + (y_arr - peak_y[source_num])^2)
                        mask = mask OR (rad_pos le boxehs)
                        if debug eq 1 and kk eq 2 then begin
                            stype = peak_type[source_num] eq 1 ? 'positive' : 'negative'
                            print, '    Masking ', stype, ' source at (', $
                                   peak_x[source_num], ',', peak_y[source_num], ')'
                        endif
                    endfor
                endif else if debug eq 1 then begin
                    print, '    Warning: No significant sources found in frame ', kk+1
                endif

                if n_pos_sources gt 0 then begin
                    mask_pos = bytarr(x_dim, y_dim)
                    for j = 0, n_pos_sources-1 do begin
                        idx2 = pos_indices[j]
                        rad_pos = sqrt((x_arr - peak_x[idx2])^2 + (y_arr - peak_y[idx2])^2)
                        mask_pos = mask_pos OR (rad_pos le boxehs)
                    endfor
                endif
                if n_neg_sources gt 0 then begin
                    mask_neg = bytarr(x_dim, y_dim)
                    for j = 0, n_neg_sources-1 do begin
                        idx2 = neg_indices[j]
                        rad_pos = sqrt((x_arr - peak_x[idx2])^2 + (y_arr - peak_y[idx2])^2)
                        mask_neg = mask_neg OR (rad_pos le boxehs)
                    endfor
                endif

                masked_pixels = where(mask, n_masked)
                if n_masked gt 0 then frame_for_stripes[masked_pixels] = !values.f_nan

                frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                frame_ii_stddev = 1.4826 * float(median(abs(double(frame_for_stripes) - $
                                    median(double(frame_for_stripes), /double, /even)), /double, /even))

                low_nantr_thresh  = frame_ii_median - nantr*frame_ii_stddev
                high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                if debug eq 1 and (kk mod 100) eq 0 then begin
                    print, 'low_nantr_thresh:  ', low_nantr_thresh
                    print, 'high_nantr_thresh: ', high_nantr_thresh, newline
                endif

                value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                             (frame_for_stripes gt high_nantr_thresh)
                nanmask = where(value_mask, n_nm)
                if n_nm gt 0 then begin
                    frame_for_stripes[nanmask] = !values.f_nan
                    frame_ii[nanmask]          = !values.f_nan
                endif

                bp_idx = where(master_badpix_mask eq 0, n_bp)
                if n_bp gt 0 then frame_for_stripes[bp_idx] = !values.f_nan

                ; --- Vertical destripe with quadrant split ---
                frame_for_stripes_top    = frame_for_stripes
                frame_for_stripes_bottom = frame_for_stripes
                frame_for_stripes_top[*,    0:y_dim/2-1]    = !values.f_nan
                frame_for_stripes_bottom[*, y_dim/2:y_dim-1] = !values.f_nan

                if debug eq 1 and (kk mod 200) eq 0 then begin
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_masked_test_top.fits',    /r), frame_for_stripes_top
                    writefits, debug_folder + strcompress('nod_'+string(i)+'_masked_test_bottom.fits', /r), frame_for_stripes_bottom
                endif

                destriped_frame_top = destripe(frame_for_stripes_top, frame_ii, 90., fraction=0.02, /no_fit, /nodisp)
                if debug eq 1 and (kk mod 200) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_destriped_frame_top_test_kk'+string(kk)+'.fits', /r), $
                    destriped_frame_top
                idx = where(destriped_frame_top gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame_top[idx] = !values.f_nan
                idx = where(destriped_frame_top lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame_top[idx] = !values.f_nan

                destriped_frame_bottom = destripe(frame_for_stripes_bottom, frame_ii, 90., fraction=0.02, /no_fit, /nodisp)
                if debug eq 1 and (kk mod 200) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_destriped_frame_bottom_test_kk'+string(kk)+'.fits', /r), $
                    destriped_frame_bottom
                idx = where(destriped_frame_bottom gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame_bottom[idx] = !values.f_nan
                idx = where(destriped_frame_bottom lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame_bottom[idx] = !values.f_nan

                destriped_frame = frame_ii
                destriped_frame[*, 0:y_dim/2-1]     = destriped_frame_bottom[*, 0:y_dim/2-1]
                destriped_frame[*, y_dim/2:y_dim-1] = destriped_frame_top[*,    y_dim/2:y_dim-1]

                ; --- Horizontal destripe ---
                frame_ii = destriped_frame
                frame_for_stripes = frame_ii

                if n_masked gt 0 then frame_for_stripes[masked_pixels] = !values.f_nan

                frame_ii_median = float(median(double(frame_for_stripes), /double, /even))
                frame_ii_stddev = 1.4826 * float(median(abs(double(frame_for_stripes) - $
                                    median(double(frame_for_stripes), /double, /even)), /double, /even))

                low_nantr_thresh  = frame_ii_median - nantr*frame_ii_stddev
                high_nantr_thresh = frame_ii_median + nantr*frame_ii_stddev
                if debug eq 1 and (kk mod 100) eq 0 then begin
                    print, 'low_nantr_thresh:  ', low_nantr_thresh
                    print, 'high_nantr_thresh: ', high_nantr_thresh, newline
                endif

                value_mask = (frame_for_stripes lt low_nantr_thresh) OR $
                             (frame_for_stripes gt high_nantr_thresh)
                nanmask = where(value_mask, n_nm)
                if n_nm gt 0 then begin
                    frame_for_stripes[nanmask] = !values.f_nan
                    frame_ii[nanmask]          = !values.f_nan
                endif

                bp_idx = where(master_badpix_mask eq 0, n_bp)
                if n_bp gt 0 then frame_for_stripes[bp_idx] = !values.f_nan

                destriped_frame = destripe(frame_for_stripes, frame_ii, 0., fraction=0.02, /no_fit, /nodisp)
                idx = where(destriped_frame gt frame_max_post_destripe, n) & if n gt 0 then destriped_frame[idx] = !values.f_nan
                idx = where(destriped_frame lt frame_min_post_destripe, n) & if n gt 0 then destriped_frame[idx] = !values.f_nan

                nod_cube[*,*,kk] = destriped_frame

                if debug eq 1 and (kk mod 200) eq 0 then writefits, $
                    debug_folder + strcompress('nod_'+string(i)+'_destriped_test_kk'+string(kk)+'.fits', /r), $
                    destriped_frame
            endfor

            current_file_basename = file_basename(nod_files[i])
            print, 'Debug: Processing file: ', current_file_basename

            suffix_pos = strpos(current_file_basename, '_cube.fits')
            if suffix_pos ne -1 then begin
                current_base_name = strmid(current_file_basename, 0, suffix_pos)
            endif else begin
                suffix_pos = strpos(current_file_basename, '.fits')
                if suffix_pos ne -1 then begin
                    current_base_name = strmid(current_file_basename, 0, suffix_pos)
                endif else begin
                    current_base_name = current_file_basename
                endelse
            endelse

            if keyword_set(destripe_skysub) then begin
                output_filename = output_folder + current_base_name + '_skysub_destriped_cube.fits'
            endif else begin
                output_filename = output_folder + current_base_name + '_corrected_cube.fits'
            endelse

            print, '  Saving corrected cube: ', file_basename(output_filename)
            writefits, output_filename, nod_cube
            delvar, nod_cube
        endfor
    endif
endelse

end
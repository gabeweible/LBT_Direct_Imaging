function select_nod_cubes, search_pattern, $
                            start_nod=start_nod, cube_indices=cube_indices, $
                            nod_numbers=nod_numbers_out, count=count, $
                            verbose=verbose

; Helper: locate NOD cube files matching a glob, sort by nod number, and
; optionally filter by a minimum nod number (start_nod) and/or by an explicit
; list of nod numbers (cube_indices).
;
; cube_indices is matched against the ACTUAL nod number embedded in each
; filename (the XX in 'nodXX'), not against a positional index in the file
; list. Missing nod numbers are warned about and skipped — they do NOT shift
; subsequent selections.
;
; SEMANTICS:
;   start_nod undefined  -> no minimum-nod filter
;   start_nod = 0        -> no minimum-nod filter (treated as "from the start")
;   start_nod = N (N>0)  -> keep only files with nod number >= N
;
; INPUTS:
;   search_pattern - file_search glob, e.g.:
;                    output_folder + obj_name + '_NOD_?_nod??_grp??_cube.fits'
;
; KEYWORDS (input):
;   start_nod    - minimum nod number to keep (see semantics above)
;   cube_indices - keep only files whose nod number is in this array;
;                  missing nod numbers are warned and skipped
;   verbose      - print progress and warnings (default 0)
;
; KEYWORDS (output):
;   nod_numbers  - integer array of nod numbers parallel to the returned
;                  file array
;   count        - number of files returned; 0 on any failure mode
;
; RETURNS:
;   string array of file paths, sorted by nod number, after filtering. []
;   when no files match.
;
; DEPENDENCIES:
;   extract_nod_number - user-defined function returning the integer nod
;                        number from a filename

compile_opt idl2, logical_predicate, strictarrsubs

; NOTE: `~` (logical not), not `not` (bitwise) — see comment in bad_pixels_fast.
if ~keyword_set(verbose) then verbose = 0

count           = 0L
nod_numbers_out = []

; ---------- diagnostic: what was received? ----------
if verbose then begin
    print, '------ select_nod_cubes ------'
    print, '  pattern:        ', search_pattern
    if n_elements(start_nod) eq 0 then begin
        print, '  start_nod:      <undefined>'
    endif else begin
        print, '  start_nod:      ', start_nod, $
               '   (type code = ', size(start_nod, /type), ')'
    endelse
    if n_elements(cube_indices) eq 0 then begin
        print, '  cube_indices:   <undefined>'
    endif else begin
        print, '  cube_indices:   ', $
               strjoin(string(cube_indices, format='(I0)'), ', ')
    endelse
endif

; ---------- glob ----------
all_files = file_search(search_pattern, count=n_all)
if verbose then print, '  files matched:  ', n_all
if n_all eq 0 then begin
    if verbose then print, '  (no files match pattern; returning empty)'
    return, []
endif

; ---------- extract nod numbers and sort ----------
all_nods = intarr(n_all)
for i = 0, n_all-1 do all_nods[i] = extract_nod_number(all_files[i])

sort_idx  = sort(all_nods)
all_files = all_files[sort_idx]
all_nods  = all_nods[sort_idx]

if verbose then $
    print, '  nods present:   ', strjoin(string(all_nods, format='(I02)'), ', ')

; ---------- start_nod filter ----------
; Only filter if start_nod is defined and > 0. start_nod = 0 is treated as
; "no filter" — equivalent to omitting the keyword.
apply_start_nod = 0
if n_elements(start_nod) gt 0 then begin
    sn = long(start_nod[0])
    if sn gt 0 then apply_start_nod = 1
endif

if apply_start_nod then begin
    if verbose then print, '  applying start_nod filter: keep nods >= ', sn
    keep_idx = where(all_nods ge sn, n_keep)
    if n_keep eq 0 then begin
        if verbose then print, '  (no files with nod number >= ', sn, '; returning empty)'
        return, []
    endif
    all_files = all_files[keep_idx]
    all_nods  = all_nods[keep_idx]
    if verbose then $
        print, '  nods after start_nod filter: ', $
               strjoin(string(all_nods, format='(I02)'), ', ')
endif else begin
    if verbose then print, '  start_nod filter not applied (undefined or 0)'
endelse

; ---------- cube_indices filter (match by nod number, not position) ----------
if keyword_set(cube_indices) then begin
    req   = long(cube_indices)
    n_req = n_elements(req)

    matched = lonarr(n_req) - 1L
    for j = 0, n_req-1 do begin
        m = where(all_nods eq req[j], nm)
        if nm gt 0 then matched[j] = m[0]
    endfor

    found = where(matched ge 0, n_found, complement=missing, ncomplement=n_missing)

    if n_missing gt 0 and verbose then begin
        miss_nods = req[missing]
        print, '  WARNING: requested nod number(s) not found, skipped: ', $
               strjoin(string(miss_nods, format='(I02)'), ', ')
        print, '           available nod numbers: ', $
               strjoin(string(all_nods,  format='(I02)'), ', ')
    endif

    if n_found eq 0 then begin
        if verbose then print, '  ERROR: none of the requested nod numbers exist'
        return, []
    endif

    all_files = all_files[matched[found]]
    all_nods  = all_nods[matched[found]]

    if verbose then $
        print, '  nods after cube_indices filter: ', $
               strjoin(string(all_nods, format='(I02)'), ', ')
endif

count           = n_elements(all_files)
nod_numbers_out = all_nods

if verbose then begin
    print, '  returning ', count, ' file(s)'
    print, '------------------------------'
endif

return, all_files

end
function destripe, d, d_not_masked, angle, clip_level=clip_level, comment=comment, no_fit=no_fit, $
        nodisp=nodisp, fraction=fraction

    compile_opt IDL2, HIDDEN
    
    ; Default parameters
    fraction = keyword_set(fraction) ? fraction : 0.2
    
    ; Get image dimensions efficiently
    n = size(d, /DIMENSIONS)
    nx = n[0]
    ny = n[1]
    n3 = n_elements(n) eq 3 ? n[2] : 1
    
    ; Preallocate output cube
    h = make_array(nx, ny, n3, TYPE=SIZE(d, /TYPE))
    
    ; Precompute rotation matrix
    r = angle/!RADEG
    m = [[cos(r), sin(r)], [-sin(r), cos(r)]]
    
    ; Compute extended dimensions efficiently
    corner_points = [[1,1], [1,-1], [-1,-1], [-1,1]]
    rotated_corners = m # corner_points
    nxx = ceil(nx * max(abs(rotated_corners[0,*])))
    nyy = ceil(ny * max(abs(rotated_corners[1,*])))
    
    ; Create mask and temporary array
    hhh = make_array(nxx, nyy, TYPE=SIZE(d, /TYPE), VALUE=0)
    mask = bytarr(nxx, nyy)
    
    ; Efficiently set initial mask
    x_offset = (nxx - nx) / 2
    y_offset = (nyy - ny) / 2
    mask[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] = 1
    
    ; Optional clipping
    if keyword_set(clip_level) then begin
        if n3 gt 1 then begin
            clip_mask = total(d, 3) / n3 lt clip_level
            mask[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] AND= clip_mask
        endif else begin
            clip_mask = d lt clip_level
            mask[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] AND= clip_mask
        endelse
    endif
    
    ; Optional display of mask
    if ~keyword_set(nodisp) then tvscl, mask[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1]
    
    xp = findgen(nxx)
    yp = findgen(nyy)
    
    ; Vectorized processing for different angle cases
    case angle of
        0: begin
            weight = total(mask, 1)
            for i = 0, n3-1 do begin
                hhh[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] = d[*,*,i]
                x = medcol(hhh, mask, 1)
                
                index = where(weight gt fraction*nx, cnt)
                ;if cnt lt 3 then message, 'No action possible: lower clip_level/fraction!'
                
                if keyword_set(no_fit) then $
                    st = replicate(1., nxx)#(x*(weight gt fraction*nx)) $
                else begin
                    dummy = SPL_INIT(yp[index], x[index])
                    xf = SPL_INTERP(yp[index], x[index], dummy, yp)
                    st = replicate(1., nxx)#xf
                endelse
                
                h[*,*,i] = d_not_masked[*,*,i] - st[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1]
                ;if keyword_set(comment) then print, 'loop ', i, ' finished'
            endfor
        end
        
        90: begin
            weight = total(mask, 2)
            for i = 0, n3-1 do begin
                hhh[where(hhh eq 0)] = !VALUES.F_NAN
                hhh[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] = d[*,*,i]
                x = medcol(hhh, mask, 2)
                
                index = where(weight gt fraction*ny, cnt)
                ;if cnt lt 3 then message, 'No action possible: lower clip_level/fraction!'
                
                if keyword_set(no_fit) then $
                    st = x*(weight gt fraction*ny)#replicate(1., nyy) $
                else begin
                    xf = INTERPOL(x[index], xp[index], xp)
                    st = xf#replicate(1., nyy)
                endelse
                
                h[*,*,i] = d_not_masked[*,*,i] - st[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1]
                ;if keyword_set(comment) then print, 'loop ', i, ' finished'
            endfor
        end
        
        else: begin
            mask = rot(mask, angle-90, /CUBIC)
            weight = total(mask, 2)
            
            for i = 0, n3-1 do begin
                hhh[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1] = d[*,*,i]
                hhh = rot(hhh, angle-90, CUBIC=-0.5)
                x = medcol(hhh, mask, 2)
                
                index = where(weight gt fraction*ny, cnt)
                if cnt lt 3 then message, 'No action possible: lower clip_level/fraction!'
                
                if keyword_set(no_fit) then $
                    st = x*(weight gt fraction*ny)#replicate(1., nyy) $
                else begin
                    xx = poly_fit(yp[index], x[index], 4)
                    xf = x
                    xf[index] = xx
                    st = (x-xf)#replicate(1., nyy)
                endelse
                
                st = rot(st, 90-angle, CUBIC=-0.5)
                h[*,*,i] = d_not_masked[*,*,i] - st[x_offset:x_offset+nx-1, y_offset:y_offset+ny-1]
                if keyword_set(comment) then print, 'loop ', i, ' finished'
            endfor
        end
    endcase
    
    return, h
end
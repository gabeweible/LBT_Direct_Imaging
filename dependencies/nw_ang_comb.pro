;Takes a cube with different image rotation angles and derotates and combines,
;performing a noise-weighted combination of the pixel values for the final image

;making this routine a drop-in for median or mean combination...

;assuming input has already been derotated

;possible to add running variance

; Written by Kevin Wagner
; Based off of Bottom et al. 2017: https://arxiv.org/pdf/1711.09119.pdf

function nw_ang_comb, in_cube, angles_and_tn, do_smooth=do_smooth

compile_opt IDL2
newline = string(10B)

;copy cube so as to leave original input alone
adi_cube=in_cube


xhs = (size(adi_cube))[1] / 2.
yhs = (size(adi_cube))[2] / 2.
if do_smooth gt 0 then lp_PSF = psf_Gaussian(NPIX=11, FWHM=[do_smooth, do_smooth], /normalize, /double)

;first go back to pupil stabilization (note: true north should be included in angles.)
for ii=0, (size(adi_cube))[3]-1 do begin
	adi_cube[*,*,ii]=rot(adi_cube[*,*,ii],angles_and_tn[ii],$
	1.0, xhs, yhs, cubic=-1.0, missing=median(adi_cube[*,*,ii], /double, /even), /pivot)
endfor

;need to replace NaN for the rest to work correctly
adi_cube[where(finite(adi_cube) lt 1)]=0.

; default computation
;compute variance of each pixel with time
; this may be useful - I *want* to down-weight samples at locations
; of noisy pixels - they cannot be trusted so well
; the only problem here is that the variance can tend to zero for "dead"
; pixels - it may be useful to find these values and mask or downweight them
; (see below)...
; outlier-resistance would be bad here, since noisy outliers are what I am
; trying to find here relative to "normal" pixels
var = variance(adi_cube, dim=3, /double, /nan)
var[where(finite(var) lt 1)]=0.
if do_smooth gt 0 then var = convolve(var, lp_PSF)

;copy into 3d cube
var_cube=adi_cube
for nn=0,(size(adi_cube))[3]-1 do var_cube[*,*,nn]=var

;derotate cubes
for ii=0, (size(adi_cube))[3]-1 do begin
	adi_cube[*,*,ii]=rot(temporary(adi_cube[*,*,ii]),-angles_and_tn[ii],$
		1.0, xhs, yhs, cubic=-1.0, missing=median(adi_cube[*,*,ii], /double, /even), /pivot)
	var_cube[*,*,ii]=rot(temporary(var_cube[*,*,ii]),-angles_and_tn[ii],$
		1.0, xhs, yhs, cubic=-1.0, missing=median(var_cube[*,*,ii], /double, /even), /pivot)
endfor


;generate an array of ones for division purposes
one_cube=var_cube
replicate_inplace, one_cube, 1.0

; actual weights on fluxes across index m (temporal direction)
; inverse-variance weighting
w_cube = one_cube / var_cube

; Replace NaN with zeros (weights and fluxes alike)
w_cube[where(finite(w_cube) ne 1)] = 0.0
adi_cube[where(finite(adi_cube) ne 1)] = 0.0

; weighted median (done for each pixel individually)
nwadi = var
replicate_inplace, nwadi, 0.0

for xx=0,2*xhs-1 do begin
	for yy=0,2*yhs-1 do begin
		; grab weights and star-subtracted fluxes across the third dimension
		
		w_arr = w_cube[xx, yy, *]
		flux_arr = adi_cube[xx, yy, *]
		
		; weighted median
		nwadi[xx,yy] = medianw(flux_arr, w_arr)
		
		; weighted mean
		;nwadi[xx,yy] = total(w_arr*flux_arr, /double, /nan) / total(w_arr, /double, /nan)
		
	endfor
endfor

;return the combined image
return, nwadi

end

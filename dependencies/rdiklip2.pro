function RDIKLIP2, cube, refcube, k_klip, target=target, checkpoint = checkpoint, anglemask=anglemask, $
distmask=distmask, posang=posang,scaling=scaling, diam=diam,pixelscale=pixelscale, wl=wl,angsep=angsep,hyper=hyper,$
obj=obj,nrings=nrings,wr =wr,n_ang =n_ang,annmode_inout=annmode_inout,spot_radius=spot_radius,rho=rho,phi=phi,anglemax=anglemax,$
create_basis=create_basis, kl_basis_arr=kl_basis_arr,ref_angles=ref_angles,filter=filter,bin=bin,FWHM=FWHM,$
save_basis=save_basis, short_suffix=short_suffix,output_path=output_path

compile_opt IDL2
newline = string(10B)

;Function to create PSF subtracted image using KLIP algorithm
; 2nd version for stream-lined Alcor processing with rdi.pro as a wrapper
; with a reference cube ('refcube') on a single frame with index 'target' in 'cube'
;K_KLIP - Number of KL Basis Vectors to retain

; diam - telescope diameter in m 
; wl - wavelength of the observations
; pixelscale - in 1/arcsec
; =======================

; number of target images in the full cube, number of references available in refcube
; for PSF modeling of cube index target=target
nimages = (size(cube))[3]
nrefimages = (size(refcube))[3]

; for now, only do KLIP in an annulus (of course, the inner radius may be 0, and the outer
; radius can be the half-size of the frames in 'cube')

print, newline, '--[ NUKLIP    ]-- STARTING ANNULAR KLIP for target = ',target, ' / ', nimages-1
print, 'angsep: ', angsep

;Set dimensions of single frame
n_cols = (size(cube))[1]
n_rows = n_cols; make the image a square if it is not already!
final_arr = fltarr(n_cols,n_rows, /nozero)  ; final subtracted image will be this shape
replicate_inplace, final_arr, !values.f_nan; fill it with nan!

progress_counter = 0.0          ; initialize progress counter
total_elements = nrings*n_ang   ; also used for progress update

; What pos angle difference translates to 1 FWHM in this annulus/ring? (at ring center)
arcdist = asin( FWHM / (total(annmode_inout)/2.) )/!DTOR ; Gives the minimum parallactic-angle difference
								
; initialize a list to fill for the new basis, if we are making one.
if create_basis eq 1 then kl_basis_arr = list()

; all possible indices for 'cube' and 'refcube', respectively
frames = findgen(nimages)
refframes = findgen(nrefimages)

; loop over annular segment index 'jj'
for jj=0, n_ang-1 do begin
	
	; get the min, max angle range for this annulus, with the first annulus always
	; starting at an angle of 0 and ending at 360./n_ang before the next annulus 
	; (if more than one) begins.
	ann_angles = [jj*(360./n_ang) , (jj+1)*(360./n_ang)]
	if max(ann_angles) gt 360 then stop; don't wrap further than 360 degrees!

	; Returns a nz x npix 2D array with the pixels in the given zone across the target 'cube'
	zones = getzone(cube, annmode_inout, ann_angles)
	; Same thing, but in the reference cube.
	refzones = getzone(refcube, annmode_inout, ann_angles)
	
	; grab the target zone for 'cube' index target=target and unwrap it into a 1-D vector/array
	target_seg = reform(zones.data[*,target])
	;print, newline, 'size(target_seg): ', size(target_seg)
	dense_target = target_seg[where(finite(target_seg) eq 1,/null)] ; Remove NaN from target_seg
	
	; OLD: subtract mean, divide by stdev to standardize each frame? (after removing NaNs)
	; subtract mean
	dense_target = dense_target - mean(dense_target, /double, /nan)
	
	; transpose to match how dense_data will look later (column vectors from row vectors)
	dense_target = transpose(dense_target)
	  
	; print some diagnostics on the angles, arcdist
	;print, newline, 'angsep, arcdist, anglemax: ', angsep, arcdist, anglemax
	if create_basis eq 1 then begin
		; absolute angular separation in degrees between all of the reference frames in 
		; 'refcube' and the parallactic angle of the index target=target frame in 'cube'
		angdiff = abs(ref_angles[refframes]-posang[target])
		;print,newline, 'angdiff: ', angdiff, newline
		
		; the good references (i.e., that don't include the target) have anglemax > angdiff > angsep FWHM
		; this gets their indices in 'refcube' and 'ref_angles'
		nottarget = where(angdiff ge ANGSEP*abs(arcdist) and angdiff le anglemax) ; AD
		
		; print how man frames meet the angular-differentiation criteria selected above
		print, 'Frames meeting separation criteria:', (size(nottarget))[1]
	
		; impose a minimum number of frames to use as a PSF reference — default to 4?
		; I am going to change this to use no references instead of all of them.
		if N_elements(nottarget) lt 4 then begin
			nottarget =  findgen(nimages)
			print, ' -[ NUKLIP ]- Less than 4 frames avialable for PCA - Removing angular distance constraint'
		endif
			
		; this is now the target zones in the angdiff-selected reference frames, as a 2-D array
		; and transposed	
		data_seg = refzones.data[*,nottarget]
		;print, newline,'size(data_seg): ', size(data_seg),newline
		
		for k=0,(size(nottarget))[2]-1 do begin
			data_seg_k = data_seg[*,k]; ref-image row vector k in data_seg
			data_seg_k = data_seg_k[where(finite(data_seg_k) eq 1,/null)] ; Remove NaN from data_seg
		
			; subtract mean
			data_seg[*,k] = data_seg_k - mean(data_seg_k, /double, /nan)
		endfor; nottarget for
		
		data_seg = transpose(data_seg); transpose so row vectors become column vectors	
		dense_data = reduce_matrix(data_seg) ; Condense Matrix (rid NaN, which should already be gone?)
		  
		; create a new eigenbasis from dense_data, applicable to dense_target

		print, 'creating KLIP basis!'
		
		new_basis = get_klip_basis(dense_data, k_klip)
		
		kl_basis_arr.Add, new_basis; add to the list over n_ang index jj
		;print, 'added KLIP basis for ang index: ', string(jj)
		
		; save the list of kl_bases to the desktop for future reference/debugging only
		if save_basis eq 1 then save, kl_basis_arr, filename=output_path+obj+short_suffix+'_filt_'+strcompress(string(filter),/r)+'_bin_'+strcompress(string(bin),/r)+'RDI_basis.sav'
		
	endif

	;last, we will project our target image onto our KL basis:
	kklip_basis_jj = kl_basis_arr[jj]; grab the basis for this loop over n_ang index jj
	
	; truncate the basis to only the first k_klip eigenimages, or use all of them if fewer are available
	; in kklip_basis_jj than are requested by k_klip. (columns given by (size())[1])
	; index column vecs, keep all rows
	full_basis_sz = (size(kklip_basis_jj))[1]

	kklip_basis_jj = kklip_basis_jj[0:min([k_klip-1, full_basis_sz-1]),*]; dont' try more k_klip than reference frames

	synthetic_1d = dense_target # matrix_multiply(kklip_basis_jj,kklip_basis_jj,/atranspose)  ; create psf reconstruction

	synthetic_1d = dense_target - synthetic_1d ; subtract psf reconstruction  

	final_arr[zones.indices] = synthetic_1d ; Build Final Image  -- AD: Replaced loop with single command

	progress_counter = progress_counter + 1 ; update progress
	
endfor; for angles (annular segments)
   
print, '--[ NUKLIP    ]-- COMPLETE for target = ',target, ' / ', nimages-1

return, final_arr
END

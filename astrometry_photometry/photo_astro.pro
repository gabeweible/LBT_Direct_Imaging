; for rho = ...; inject negative planet to try and cancel out the secondary psf, get the stdev of the
; square around the negative injection to be as low as possible which means that the residuals there
; are as flat as possible, as if the seondary wasn't there (loop over r, theta, and contrast to 
; refine it)
; 
; aper does aperture photometry, cntrd gets the "center of light" for the sources in the apertures

; I probably need to go full grid-search and loop over negative contrasts at all positions as well.
; I think it makes sense to pass the xx and yy values straight into the planet injection as well,
; as something could go slightly wrong with either conversion to and from radii and angles, but
; I've been having problems with this so for right now the injection is passed through a radius and
; an angle.
; The centroiding should get me within a pixel of the true best negative injection center

;stdev of differences in astrometry (picked up & injected) are uncertainties in astrometry and photometry
; make sure that there isn't a big difference in one axis and not the other, but otherwise just worrying
; about the absolute differences in position. Final results in PA and sep. More uncertainty in PA than in
; sep because of the orientation of LBTI

pro photo_astro, coadd=coadd, type=type, grid_sz=grid_sz, nod=nod, write=write
; Type is 'ADI' or 'KLIP'
; Nod is 'total', '1', '2', '3', or '4'

COMPILE_OPT IDL2; Strictarr and 32-bit ints
newline = string(10B); Make printing newlines convenient
; Get the current time in a Julian date format, used to print how long the
; photometry/astrometry processing took.
start_time = systime(/JULIAN)

;------------------------------[ Start User Input ]---------------------------------

; Where to find our files and put the results
output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_'+$
	strcompress(coadd,/r)+'_'+strcompress(nod,/r)+'/'
	
obj = 'HII1348'; Observed object

; Parameters needed to read in our total_klip or total_adi file after 
; negative injection
klip=0 ; don't run KLIP in HII1348_pipeline.pro right now.
bin = 3 & bin_type = 'mean' & combine_type = 'nwadi'
k_klip = 7 & angsep= 1. & anglemax = 360. & nrings = 4.
n_ang = 2 & do_cen_filter = 1 & filter = 17. & ct = 0.994

if ~ keyword_set(grid_sz) then grid_sz = 5; Default

; Inject_planets parameters
uncert = 0; not getting uncertainties yet
use_gauss = 1; Fit a gaussian to the pupil median PSF (more symmetric)
n_planets = 1; Only one negative injection at a time
; arcsec/pixel, updated now after Trapezium astrometry solution
pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared)
pxscale_dx = 0.010700 ; +0.000042 or -0.000051 arcsec/pixel (from Steve and Jared)

min_pxscale = min([pxscale_sx, pxscale_dx]); Both are matched to this in KLIP/ADI
; starting guess for the (negative) contrast, should be pretty close.
c_guess_total = -0.00940320
c_guess_1 = -0.00959084
c_guess_2 = -0.00915081
c_guess_3 = -0.00969534
c_guess_4 = -0.00924411
c_guess_nods = [c_guess_1, c_guess_2, c_guess_3, c_guess_4]

n_contrasts = grid_sz; Number of contrasts to test at each position, ODD

; Companion centroid guess [x,y] indices (start at 0)
guess_total = [277.301, 353.548]
guess_1 = [277.133, 353.696]
guess_2 = [277.272, 353.125]
guess_3 = [277.419, 353.635]
guess_4 = [277.377, 353.685]
guess_nods = [[guess_1], [guess_2], [guess_3], [guess_4]]

fwhm = 9.6 ; px ``width'' in reduce_lbti_HII1348.pro
; nx x ny grid around centroid result
nx = grid_sz & ny = grid_sz ; 5 x 5 grid is default


initial_hc = 1 ; plus or minus 1 %
initial_hw = 1 ; plus or minus 1 px (about 1 %)

; Let's go tiny! I wanna be done with this. Should take us to 6 sigfigs for each
hc_thresh = 0.0000005 ; percentage
hw_thresh = 0.0005

;------------------------------[ End User Input ]---------------------------------

; Default to ADI for now
if not keyword_set(type) then type = 'ADI'
if not keyword_set(nod) then nod='total'; Don't do four independent measurements
; unless specified

if nod eq 'total' then begin

	print, 'Reading in total ' + type +  ' image...'
	;if type eq 'KLIP' then begin
	;  
	;   og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
	;      +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
	;      string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	;      string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
	;      '_total_klip.fits', /rem))
	;    
	;endif

	if type eq 'ADI' then begin
		og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+$
		'filt_'+string(filter)+'_neg_inj_'+string(0)+'_uncert_0_total_adi.fits',/rem))
	endif
	
	c_guess = c_guess_total
	guess = guess_total
	
endif else begin; nod eq 'total' if
	
	print, 'Reading in nod ' + type +  nod + ' image...'
	;if type eq 'KLIP' then begin
	;  
	;   og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
	;      +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
	;      string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	;      string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
	;      '_total_klip.fits', /rem))
	;    
	;endif

	if type eq 'ADI' then begin
		og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+$
		'filt_'+string(filter)+'_neg_inj_'+string(0)+'_uncert_0_adi_nod' + nod +$
		'.fits',/rem))
	endif
	
	c_guess = c_guess_nods[fix(nod)-1]
	guess = guess_nods[*, fix(nod)-1]

endelse; nod ne 'total' else

	
x_size = (size(og_image))[1] & y_size = (size(og_image))[2]
;half-sizes
x_hsize = x_size / 2.
y_hsize = y_size / 2.

print, 'Original image read, finding companion centroid...', newline

; Calculate the centroid of the companion two ways (at least to test, for now)
; Apparently, gcntrd is better, but slower?
; For now, I'll take the average of the two results
cntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
print, 'CNTRD xcen, ycen:', XCEN, YCEN
;x_avg = XCEN & y_avg = YCEN
cen_x = XCEN & cen_y = YCEN
gcntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
print, 'GCNTRD xcen, ycen:', XCEN, YCEN
;x_avg += XCEN & x_avg *= 0.5 & y_avg += YCEN & y_avg *= 0.5
cen_x += XCEN & cen_x *= 0.5 & cen_y += YCEN & cen_y *= 0.5
print, 'Mean xcen, ycen:', cen_x, cen_y, newline
; This will get us the best place to inject our planet at
print, 'Starting loop over xx, yy around xcen, ycen'

;--------------------------------------------------------------------------------

hc = initial_hc & hw = initial_hw
x_avg = cen_x & y_avg = cen_y
con = c_guess
; Loop until we're within BOTH of our thresholds
i = 1
; make a folder to put our stuff in
file_mkdir, strcompress(output_path + 'photometry/nod_' + nod, /r)
WHILE (hc gt hc_thresh) || (hw gt hw_thresh) DO BEGIN

; Initialize arrays for our results
xxs=[] & yys=[] & cons=[] & devs=[] & means=[] & rhos=[] & thetas=[]

file_mkdir,$
strcompress(output_path + 'photometry/nod_' + nod + '/' + string(i), /r)

; Define lower bounds, upper bounds, and step sizes for our nested loops (grid search)
x_i = x_avg-hw & x_f = x_avg+hw
x_step = ((x_avg+hw)-(x_avg-hw))/(nx-1)

y_i = y_avg-hw & y_f = y_avg+hw
y_step = ((y_avg+hw)-(y_avg-hw))/(ny-1)

c_i = con*(1.0-hc) & c_f = con*(1.0+hc)
c_step = (con*(1.0+hc)-con*(1.0-hc))/(n_contrasts-1)

; Create arrays to loop through (for some reason using normal for loops didn't work...)
x_loop = [x_i : x_f : x_step]
y_loop = [y_i : y_f : y_step]
c_loop = [c_i : c_f : c_step]

trial=0; Keep track of which run we're on for a given i-value
; Instead of writing each image and then making a cube, let's just
; make a cube from the beginning to write.
;cube = list()
foreach xx, x_loop do begin; Loop over x
   foreach yy, y_loop do begin; Loop over y

      ; Something isn't working with injecting planets directly, so for now I'm
      ; using r and theta

		; radius from center of the image in arcsec
      planet_r = min_pxscale * SQRT( (xx-x_hsize)^2. + (yy-y_hsize)^2. )
      ; radians
      planet_theta = ATAN( (yy-y_hsize) / (xx-x_hsize) )
 
      foreach contrast, c_loop do begin; Loop over (negative) contrasts
         
     	    ;Inject a negative planet at (xx, yy) = (planet_r, planet_theta) with the given contrast,
     	    ;then do rotate and klip
;     	   hii1348_pipeline, planet_x=xx, planet_y=yy, contrast=contrast, pre_inj_stuff=0, neg_inj=1,$
;     	      trial=trial, coadd=coadd, use_gauss=use_gauss

     	   hii1348_pipeline, rho=planet_r,theta=planet_theta, contrast=contrast,$
     	   pre_inj=0, neg_inj=1, trial=trial, coadd=coadd, use_gauss=use_gauss,$
     	   uncert=uncert, klip=klip, fs=0, extra=nod, nod=nod; Inject and run ADI
      
     	   ; Read in the total KLIP or ADI file after the negative injection
     	   print, 'Reading in neg-injected file'
;     	   if type eq 'KLIP' then begin
;     	       
;     	      image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+$
;     	         '_type_'+bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+$
;     	         '_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'$
;     	         +string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'_neg_inj_'+string(1)+$
;     	         '_total_klip.fits', /rem))
;     	         
;     	   endif
;     	   if type eq 'ADI' then begin (only real option for now)
				if nod eq 'total' then begin
				
					image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+$
						string(ct)+'filt_'+string(filter)+'_neg_inj_'+string(1)+$
						'_uncert_0_total_adi.fits',/rem))
						
				endif else begin; nod eq 'total' if
				
					image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+$
						string(ct)+'filt_'+string(filter)+'_neg_inj_'+string(1)+$
						'_uncert_0_adi_nod'+nod+'.fits',/rem))
						
				endelse; nod neq 'total' else
					
;     	   endif
         print, 'Read in image after negative injection'
      
     	   ; Calculate the standard deviation (and the mean value) for the box around where we
   	   ; injected the negative planet over the real companion
   	   print, 'Finding stdev and mean...'
	      deviation = stdev(image[x_avg-20.:x_avg+20., y_avg-20.:y_avg+20.], mean)
	      print,'Standard deviation (40px x 40px square) =', deviation
     	   average = mean
     	   print, 'Mean =', average
     	   
     	   ; Append results to arrays
     	   print, 'Appending to arrays...'
     	   xxs=[xxs,xx] & yys=[yys,yy]
     	   cons=[cons,contrast] & devs=[devs,deviation]
     	   means=[means,average]
     	   rhos=[rhos,planet_r] & thetas=[thetas,planet_theta]
     	   ;print, 'Done.'
     	   ;print, 'Adding image to cube...'
     	   
     	   ;writefits, strcompress(output_path+'photometry/'+string(i)+'/'+obj+$
     	   	;'_trial_'+string(sigfig(trial,4))+'.fits', /rem), image
     	   	
     	   ; Add image to cube instead of writing
     	   ;cube.Add, [[image]]
     	   	
     	   print, 'Done.'+newline+'Incrementing trial...'
     	   trial += 1
     	   print, 'Done.'
     	   
     	   ; Convert to scalars so that IDL doens't throw a fit
     	   yy=yy[0]
     	   xx=xx[0]
     	   
      endforeach; contrast foreach
   endforeach; yy foreach
endforeach; xx foreach

;if write eq 1 then begin

;	print, 'Converting to array...'
;	cube = cube.toArray(/TRANSPOSE, /NO_COPY)

	; Save the results
;	save,filename=strcompress(output_path+'photometry/nod_'+nod+'/'+$
;		string(i)+'/'+obj+'_negative_inj_data_while_'+string(i)+'.sav', /r),xxs,yys,$
;		cons,devs,means,hw,hc,rhos,thetas

	; Combine all of the trials into a cube and write it to the same folder
;	print, newline, 'Saving the trials into one FITS cube'
	;folder_cube, strcompress(output_path+'photometry/'+string(i)+'/', /r),$
		;save_path=strcompress(output_path+'photometry/while_'+$
		;string(i)+'_', /r)
	
;	writefits, strcompress(output_path+'photometry/nod_'+nod+'/'+string(i)+$
;		'/'+obj+'_while_'+string(i)+'_cube.fits',/rem), cube
		
;	print, 'FITS cube created!''
		
;endif
	
print, 'Starting analysis', newline

min_devs = WHERE(devs eq MIN(devs))
median_i = median(min_devs)

; Print some good results, restrict to one value (median)
x_best = xxs[median_i]
y_best = yys[median_i]

print, newline, 'Min. STDEV:', MIN(devs), 'at (x, y): ('+string(x_best)+', '+$
	string(y_best)+')'
	
left = x_hsize - x_best; x > 250 is "minus left", since it's right (px)
up = y_best - y_hsize; How much above the center we are (px)

; Convert to RA and DEC with the plate scale
RA = min_pxscale * left
DEC = min_pxscale * up

print, "This is at RA:", RA, 'and DEC:', DEC
	
	
rho_px = SQRT( (x_best-x_hsize)^2. + (y_best-y_hsize)^2.)
theta = ATAN( (y_best-y_hsize) / (x_best-x_hsize) )

; Convert to arcsec and PA
rho_arcsec = rho_px * min_pxscale
PA = 270. + !RADEG * theta; CW from +y-axis (deg)

print, 'Or, at (rho, theta): ('+string(min_pxscale*rho_px)+', '+string(theta)+')',newline
	
print,'Best trial was trial'+string(where(devs eq min(devs)))

print, 'Compare with centroiding (x, y): ('+string(cen_x)+', '+string(cen_y)+')',newline
best_con = cons[median_i]
print, 'Contrast:', string(best_con), newline

if write eq 1 then begin

	save, filename=strcompress(output_path+'photometry/nod_'+nod+'/'+$
		string(i)+'/'+obj+'_results_while_'+string(i)+'.sav', /r), x_best, y_best,$
		left, up, RA, DEC,rho_px, theta, rho_arcsec, PA, best_con

endif; write eq 1 if

; Reset stuff for the next loop iteration:
x_avg = x_best & y_avg = y_best
hw *= 0.8 & hc *= 0.8 ; Shrink 3D ``grid'' by 20% in all dimensions.
con = best_con

i += 1; increment while loop counter
ENDWHILE; thresholds loop

;-----------------------------------------------------------------------------------
; Save the FINAL results
save,filename=strcompress(output_path+'photometry/nod_'+nod+'/'+obj+$
	'_negative_inj_data_final.sav',/r),xxs,yys,cons,devs,means
	 
photo_astro_to_csv, strcompress(output_path+'photometry/nod_'+nod+'/'+obj+$
	'_negative_inj_data_final.sav',/r)

;-----------------------------------------------------------------------------------

print, 'Completed photometry in ',(systime(/JULIAN) - start_time) * 1440.,' minutes.'

end; That's all, folks!
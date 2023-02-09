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

pro photometry, coadd=coadd, type=type, grid_sz=grid_sz
; Type is 'ADI' or 'KLIP'

COMPILE_OPT IDL2; Strictarr and 32-bit ints
newline = string(10B); Make printing newlines convenient
; Get the current time in a Julian date format, used to print how long the
; photometry/astrometry processing took.
start_time = systime(/JULIAN)

;------------------------------[ Start User Input ]---------------------------------

; Where to find our files and put the results
output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_'+$
	strcompress(coadd,/r)+'/'
	
obj = 'HII1348'; Observed object

; Parameters needed to read in our total_klip or total_adi file after 
; negative injection
bin = 3 & bin_type = 'mean' & combine_type = 'nwadi'
k_klip = 7 & angsep= 1. & anglemax = 360. & nrings = 4.
n_ang = 2 & do_cen_filter = 1 & filter = 17. & ct = 0.994

if ~ keyword_set(grid_sz) then grid_sz = 5; Default

; Inject_planets parameters
uncert = 0; not get uncertainties yet
use_gauss = 1; Fit a gaussian to the pupil median PSF (more symmetric)
n_planets = 1; Only one negative injection at a time
; arcsec/pixel, will need to be updated after Trapezium calculations
pxscale = 0.0107 
; starting guess for the (negative) contrast, should be pretty close.
c_guess = -0.008914755
n_contrasts = grid_sz; Number of contrasts to test at each position, ODD

; Companion centroid guess [x,y] indices (start at 0)
guess = [277.3, 353.5]
fwhm = 8.72059 ; px ``width'' in reduce_lbti_HII1348.pro
; nx x ny grid around centroid result
nx = grid_sz & ny = grid_sz ; 5 x 5 grid is default


initial_hc = 0.02 ; plus or minus 2%
initial_hw = 0.5 ; plus or minus half a pixel (centroiding should get within 1 px)

hc_thresh = 0.001; uncert is about 1%, so lets get within 0.1%
hw_thresh = 0.005; (5/1000 of a px on either side, uncert is about 5/100 either way)

;------------------------------[ End User Input ]---------------------------------

; Initialize arrays for our results
xxs=[] & yys=[] & cons=[] & devs=[] & means=[] & rhos=[] & thetas=[]

; Default to ADI for now
if not keyword_set(type) then type = 'ADI'

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

print, 'Original image read, finding companion centroid...', newline

; Calculate the centroid of the companion two ways (at least to test, for now)
; Apparently, gcntrd is better, but slower? For now I'll take the average of the two results
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
print, 'Starting loop over xx, yy around mean xcen, ycen (nx x ny px square) or overridden x_avg, y_avg...'

;--------------------------------------------------------------------------------

hc = 0.016 & hw = 0.4
x_avg = 277.713 & y_avg = 353.174
con = -0.00909305
; Loop until we're within BOTH of our thresholds
i = 1
WHILE (hc > hc_thresh) || (hw > hw_thresh) DO BEGIN

file_mkdir, strcompress(string(i), /r)

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
foreach xx, x_loop do begin; Loop over x
   foreach yy, y_loop do begin; Loop over y

      ; Something isn't working with injecting planets directly, so for now I'm using r and theta
      planet_r = pxscale * SQRT(((xx-250.)^2.)+((yy-250.)^2.))
      planet_theta = ATAN((250.-yy)/(250.-xx)); radians (ccw from top at 0 rad)

      foreach contrast, c_loop do begin; Loop over (negative) contrasts
         
     	    ;Inject a negative planet at (xx, yy) = (planet_r, planet_theta) with the given contrast,
     	    ;then do rotate and klip
;     	   hii1348_pipeline, planet_x=xx, planet_y=yy, contrast=contrast, pre_inj_stuff=0, neg_inj=1,$
;     	      trial=trial, coadd=coadd, use_gauss=use_gauss

     	   hii1348_pipeline, rho=planet_r, theta=planet_theta, contrast=contrast,$
     	   pre_inj=0, neg_inj=1, trial=trial, coadd=coadd, use_gauss=use_gauss,$
     	   uncert=uncert, klip=klip, fs=0; Inject and run ADI
      
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
     	      image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+string(ct)+'filt_'+ $
     	         string(filter)+'_neg_inj_'+string(1)+'_uncert_0_total_adi.fits',/rem))
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
     	   xxs=[xxs,xx] & yys=[yys,yy] & cons=[cons,contrast] & devs=[devs,deviation] & means=[means,average]
     	   rhos=[rhos,planet_r] & thetas=[thetas,planet_theta]
     	   print, 'Done.'
      
     	   ; Save our results (in the loop)
     	   print, 'Saving...'
     	   save,filename=strcompress(output_path+'photometry/'+string(i)+'/'+obj+'_negative_inj_data_trial_' + string(trial)+'.sav', /r),xxs,yys,cons,devs,means,rhos,thetas,hw,hc
     	   print, 'Done.'+newline+'Writing FITS...'
     	   writefits, strcompress(output_path+'photometry/'+string(i)+'/'+obj+'_trial_'+string(sigfig(trial,4))+'.fits', /rem), image
     	   print, 'Done.'+newline+'Incrementing trial...'
     	   trial += 1
     	   print, 'Done.'
     	   yy=yy[0]
     	   xx=xx[0]; Convert to scalars so that IDL doens't throw a fit
      endforeach; contrast foreach
   endforeach; yy foreach
endforeach; xx foreach

; Save the results
save,filename=strcompress(output_path+'photometry/'+obj+'_negative_inj_data_while_'+string(i)+'.sav', /r),xxs,yys,cons,devs,means,hw,hc,rhos,thetas

; Combine all of the trials into a cube and write it to the same folder
print, newline, 'Saving the trials into one FITS cube'
folder_cube, strcompress(output_path+'photometry/'+string(i)+'/', /r),$
	save_path=strcompress(output_path+'photometry/'+string(i)+'/while_'+string(i)+'_', /r)
	
print, 'FITS cube created! Starting analysis', newline

; Print some good results
x_best = xxs[WHERE(devs eq MIN(devs))]
y_best = yys[WHERE(devs eq MIN(devs))]

print, newline, 'Min. STDEV:', MIN(devs), 'at (x, y): ('+string(x_best)+', '+$
	string(y_best)+')
	
print, 'Or, at (rho, theta): ('+string(pxscale * SQRT(((x_best-250.)^2.)+$
	((y_best-250.)^2.)))+', '+string(ATAN((250.-y_best)/(250.-x_best)))+')', newline
	
print,'Best trial was trial'+string(where(devs eq min(devs)))

print, 'Compare with centroiding (x, y): ('+string(cen_x)+', '+string(cen_y)+')', newline
best_con = cons[WHERE(devs eq MIN(devs))]
print, 'Contrast:', string(best_con), newline

; Reset stuff for the next loop iteration:
x_avg = x_best & y_avg = y_best
hw *= 0.8 & hc *= 0.8 ; Shrink 3D ``grid'' by 20% in all dimensions.
con = best_con

i += 1; increment while loop counter
ENDWHILE; thresholds loop

;-----------------------------------------------------------------------------------
; Save the FINAL results
save,filename=strcompress(output_path+'photometry/'+obj+'_negative_inj_data_final.sav', /r),xxs,yys,cons,devs,means
folder_cube, strcompress(output_path+'photometry/', /r),$
	save_path=strcompress(output_path+'photometry/'+'final_', /r)

;-----------------------------------------------------------------------------------

print, 'Completed photometry in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!
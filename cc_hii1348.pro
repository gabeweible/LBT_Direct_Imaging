pro cc_hii1348
COMPILE_OPT IDL2

; Trying this out...
SET_PLOT, 'Z'

newline = string(10B)

;USER INPUTS
; Going being 2.0 and with more increments to see if that improves our results
; at all (though this will take longer to run, to be sure)
separations=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,$
				 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]; Arcsec

; Position angles, 18 of them all separated by 20 deg
; starting at 5 deg
pas = findgen(18)*20. + 5. 
threshold = 5.; Sigma/Significant/SNR
pxscale = 0.0107; Plate scale for our observations [arcsec/px]
contrast_guess = 0.005 ;guess at the starting point for the contrast curve

; Where to find our files and put the results
output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_25/'

; Rotate/KLIP parameters (needed to read in our total_klip file)
bin = 3
bin_type = 'mean'
combine_type = 'nwadi'
k_klip = 7
angsep= 1.
anglemax = 360.
nrings = 4.
n_ang = 2

;-------------------- END USER INPUT -------------------------------

;save start time
contstarttime=systime(/JULIAN)

;set up arrays
results=separations & replicate_inplace, results, 0.
contrasts=separations & replicate_inplace, contrasts, contrast_guess
curves=[]

;start a loop over position angle
runs = 0
for m = 0,n_elements(pas)-1 do begin

PA=separations & replicate_inplace, PA, pas[m]*!DTOR ;set all same PA for now

;start a loop over separation
for n = 0,n_elements(separations)-1 do begin
	print, 'Working on test separation',n,'/',n_elements(separations)-1
	result=99E99 ;start significance off at infinity
	
	; Loop until we get within 5% of 5-sigma
	while result gt 1.05*threshold or result lt 0.95*threshold do begin
	
		rho = separations[n] & theta = PA[n] & contrast = contrasts[n]
		if n gt 0 and finite(contrast) eq 0 then contrast=contrasts[n-1]
		print, 'Testing rho, theta, contrast = ', rho, theta/!DTOR, contrast
		print, 'Run:', runs
		runs += 1

		;set boundaries to remove the line (NEAR SPECIFIC)
		;if rho lt 0.6 then begin 
		;	aa=5. & bb=5. & ab=5. & ba=5. 
		;	if pas[m] lt 200. then ab=10.
		;endif else begin 
		;	aa=10. & bb=10.  & ab=10. & ba=10. 
		;endelse
		;aa=0. & bb=0.  & ab=0. & ba=0. ;fix to zero


		;perform reduction
		;name='S37'
		;reduce_near,rho=rho,theta=theta,contrast=contrast,/block_burn,aa=aa,bb=bb,ab=ab,ba=ba
		hii1348_pipeline, rho=rho, theta=theta, contrast=contrast, pre_inj=0,$
		uncert=0, neg_inj=0
		
		;read image
		image=readfits(strcompress(output_path + 'combined/' + 'HII1348' +$
		 	'_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+$
		 	combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
		 	string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+$
		 	'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) +$
		  	'_neg_inj_0' + '_total_klip.fits', /rem))

		; remove zeros
		image[where(image eq 0)] = !values.f_nan

		;what is rho, theta on the image?
		xs = (size(image))[1] & ys = (size(image))[2]
		xhs = xs/2. & yhs = ys/2. ;half-sizes
		print, xhs, yhs
		print, rho
		xp=xhs + ( (rho/pxscale) * Cos(theta)  ) & yp=yhs + ( (rho/pxscale) * Sin(theta)  )
		print, 'Xp,Yp = ',xp,yp

		;measure signal
		signal = max(image[xp-1:xp+1,yp-1:yp+1])

		;measure noise:
		
		nimage = image ;make a copy
		boxehs = 10. ;source exclusion half-size
		nimage[xp-boxehs:xp+boxehs-1, yp-boxehs:yp+boxehs-1] = !values.f_nan

		;exclude known companion
		boxehs=20. ;source exclusion half-size
		nimage[273-boxehs:273+boxehs-1, 356-boxehs:356+boxehs-1] = !values.f_nan

		; The noise is just the stdev of the image with the known companion removed
		noise=stddev(nimage,/nan)
		
		; S/N is well Signal/Noise !
		result = signal/noise
			
		print, 'Source recovered with significance = ', result

;		if runs gt 5 then begin & result=threshold & contrasts[n] = !values.f_nan & endif

		; Changing the gt .1 to .2 to see if that fixes my nan contrast problem...
		; Spoiler! It doesn't...but I'm debugging this...
		; Maybe gt 1.0? idk
		; 2.0?
		; Putting it to 5.0 just so that I (probably) won't have to worry about it anymore...
		if runs gt 5 and abs(1.-(result/threshold) ) gt 5.0 then begin & print, 'abs:',$
			abs(1.-(result/threshold) ) & result=threshold & contrasts[n] = !values.f_nan &$
			print, "Hey dude, I'm setting the contrast to NaN now..." & nan_problem = 1 &$
			stop & endif

		if runs gt 5 and abs(1.-(result/threshold) ) le 0.1 then result=threshold 
		
		;if result is not within 5% of the threshold, revise the injected contrast
		if abs(1.-(result/threshold) ) gt 0.05 then contrasts[n:n_elements(separations)-1]=$
			contrast/((result/threshold))	
		;if result lt 0.95*threshold then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold)) 
		results[n]=result ;store results for diagnostics
		print, contrasts
		print, results
		print, separations
		;hak

		;start plotting things
		; Do I need to manually enter the file path for cgps_open here?
		; 0 is landscape 1 is portrait
		cgps_open,strcompress(output_path + 'HII1348_contrast_curve.ps',/rem),xsize=7.25,$
		 ysize=5, /inches, portrait=0

			;open a new plot if it is the first position angle
			if m eq 0 then cgplot,separations,contrasts,xtitle='Separation (arcsec)',$
				ytitle='Contrast',/ylog,yrange=[1E-4,1E-2],title='result = '+$
				string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts),3)),$
				charsize=1	
		
			;open a new plot for the first PA and then oplot for the remaining
			if m gt 0 then begin
			
				cgplot,separations,curves[*,0],xtitle='Separation (arcsec)',$
					ytitle='Contrast',/ylog,yrange=[1E-4,1E-2],title='result = '+$
					string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts),3)),$
					charsize=1,color='gray'
					
				if m gt 1 then for l=1,m-1 do cgoplot,separations,curves[*,l],color='gray'

				cgoplot, separations, contrasts,color='gray' ;now plot the current contrasts
			endif

			;cgoplot,[0,2],[threshold*0.67E-6,threshold*0.67E-6],linestyle=2

		;plot mean contrast curve from the data reductions

		if m gt 1 then begin
			medcurve=mean(curves,dim=2,/nan)
			cgoplot,separations,medcurve,thick=4
		endif

		cgps_close,/pdf

		writefits,strcompress(output_path + 'HII1348_contrastcurve_testimage.fits',/rem),image
		writefits,strcompress(output_path + 'HII1348_contrastcurve_testnoise.fits',/rem),nimage
		;output images once the threshold has been found
		if result ge 0.95*threshold and result le 1.05*threshold then $
			if n eq 0 and m eq 0 then imgs=image else imgs=[ [[imgs]], [[image]] ]
			
		if n gt 0 or m gt 0 then writefits,strcompress(output_path +$
			'HII1348_contrastcurve_images.fits',/rem),imgs
			
	endwhile
endfor ;n

curves=[ [curves], [contrasts] ]

save,filename=strcompress(output_path + 'HII1348_contrast_output.sav',/rem),$
	curves,pas,separations,threshold
; Write curves to a CSV.
WRITE_CSV, strcompress(output_path + 'HII1348_contrast_output.csv',/rem), curves

endfor ;m

print, 'Completed contrast curve generation in ',(systime(/JULIAN)-contstarttime)*$
	86400./60.,' minutes.'

end

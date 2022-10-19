pro contrastcurve_HII1348

;USER INPUTS
separations=[ 0.1,0.2,0.3,0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
pas=findgen(10)*36.
threshold=5.
pxscale=0.0107
contrast_guess=0.005 ;guess at the starting point for the contrast curve



;-------------------- END USER INPUT -------------------------------

;save start time
contstarttime=systime(/JULIAN)

;set up arrays
results=separations & results[*]=0.
contrasts=separations & contrasts[*]=contrast_guess
curves=[]

;start a loop over position angle
for m=0,n_elements(pas)-1 do begin

PA=separations & PA[*]=pas[m]*!DTOR ;set all same PA for now

;start a loop over separation
for n=0,n_elements(separations)-1 do begin
	print, 'Working on test separation',n,'/',n_elements(separations)-1
	result=99E99 ;start significance off at infinity
	runs=0
	while result gt 1.05*threshold or result lt 0.95*threshold do begin
		rho=separations[n] & theta=PA[n] & contrast=contrasts[n]
		if n ge 1 and finite(contrast) eq 0 then contrast=contrasts[n-1]
		print, 'Testing rho, theta, contrast = ',rho,theta/!DTOR,contrast
		runs=runs+1

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
		reduce_lbti_HII1348,rho=rho,theta=theta,contrast=contrast		
		
		;read image
		image=readfits(strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_total_klip.fits',/rem))

		;remove zeros
		image[where(image eq 0)]=!values.f_nan

		;what is rho, theta on the image?
		xs=(size(image))(1) & ys=(size(image))(2)
		xhs=xs/2. & yhs=ys/2. ;half-sizes
		print, xhs, yhs
		print, rho
		xp=xhs + ( (rho/pxscale) * Cos(theta)  ) & yp=yhs + ( (rho/pxscale) * Sin(theta)  )
		print, 'Xp,Yp = ',xp,yp

		;measure signal
		signal=max(image[xp-1:xp+1,yp-1:yp+1])

		;measure noise
		nimage=image ;make a copy
		boxehs=10. ;source exclusion half-size
		nimage[xp-boxehs:xp+boxehs-1,yp-boxehs:yp+boxehs-1]=!values.f_nan

		;exclude known companion
		boxehs=20. ;source exclusion half-size
		nimage[273-boxehs:273+boxehs-1,356-boxehs:356+boxehs-1]=!values.f_nan

		noise=stddev(nimage,/nan)
		result=signal/noise
			
		print, 'Source recovered with significance = ',result

;		if runs gt 5 then begin & result=threshold & contrasts[n] = !values.f_nan & endif
		if runs gt 5 and abs(1.-(result/threshold) ) gt 0.1 then begin & result=threshold & contrasts[n] = !values.f_nan & endif
		if runs gt 5 and abs(1.-(result/threshold) ) le 0.1 then result=threshold 
		
		;if result is not within 5% of the threshold, revise the injected contrast
		if abs(1.-(result/threshold) ) gt 0.05 then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold))	
		;if result lt 0.95*threshold then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold)) 
		results[n]=result ;store results for diagnostics
		print, contrasts
		print, results
		print, separations
		;hak

		;start plotting things
		cgps_open,strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_contrast_curve.ps',/rem),xsize=7.25, ysize=5, /inches, portrait=1

			;open a new plot if it is the first position angle
			if m eq 0 then cgplot,separations,contrasts,xtitle='Separation (arcsec)',ytitle='Contrast',/ylog,yrange=[1E-4,1E-2],title='result = '+string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts),3)),charsize=1	
		
			;open a new plot for the first PA and then oplot for the remaining
			if m gt 0 then begin
				cgplot,separations,curves[*,0],xtitle='Separation (arcsec)',ytitle='Contrast',/ylog,yrange=[1E-4,1E-2],title='result = '+string(sigfig(result,3))+' '+'Floor = '+string(sigfig(min(contrasts),3)),charsize=1,color='gray'
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

		writefits,strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_contrastcurve_testimage.fits',/rem),image
		writefits,strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_contrastcurve_testnoise.fits',/rem),nimage
		;output images once the threshold has been found
		if result ge 0.95*threshold and result le 1.05*threshold then $
			if n eq 0 and m eq 0 then imgs=image else imgs=[ [[imgs]], [[image]] ]
		if n gt 0 or m gt 0 then writefits,strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_contrastcurve_images.fits',/rem),imgs
	endwhile
endfor ;n

curves=[ [curves], [contrasts] ]

save,filename=strcompress('/media/kevin/Storage/LBT/HII1348/combined/HII1348_contrast_output.sav',/rem),curves,pas,separations,threshold

endfor ;m



print, 'Completed contrast curve generation in ',(systime(/JULIAN)-contstarttime)*86400./60.,' minutes.'

end

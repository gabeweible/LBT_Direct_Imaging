pro alcor_coadd_testing2

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

; Allowed coadds that, I think, actually can run on my machine
; descending so that the coadds less likely process (and longer to process)
; are run last.
coadds = reverse(indgen(3, start=30, increment=5))
corr_threshes = findgen(10, start=0.990, increment=0.001)

foreach coadd, coadds do begin
foreach corr_thresh, corr_threshes do begin

	alcor_pipeline, neg_inj=0, pre_inj=1, coadd=coadd,$
		uncert=0, klip=0, fs=1, nod='sx_only', cds='endcds', corr_thresh=corr_thresh
		
endforeach; corr_threshes
endforeach; coadds
  
end
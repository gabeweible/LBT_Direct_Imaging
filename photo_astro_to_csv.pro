pro photo_astro_to_csv, filename
; Takes a .sav from contrast curve generation and writes a .csv file with the
; curves for easier handling with other programming languages (e.g. Python)

compile_opt IDL2; strictarr, 32 bit integers

restore, filename ; restore saved data
; Now we have curves, and separations (primarily)

; Write curves to a CSV. I think that this is all that I need
WRITE_CSV, filename + '.csv', xxs,yys,cons,devs,means

end
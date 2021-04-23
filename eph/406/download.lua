#! /usr/bin/env lua
--require 'sh'
require 'ext'
string.split([[ascm0100.406
ascm0200.406
ascm0300.406
ascm0400.406
ascm0500.406
ascm0600.406
ascm0700.406
ascm0800.406
ascm0900.406
ascm1000.406
ascm1100.406
ascm1200.406
ascm1300.406
ascm1400.406
ascm1500.406
ascm1600.406
ascm1700.406
ascm1800.406
ascm1900.406
ascm2000.406
ascm2100.406
ascm2200.406
ascm2300.406
ascm2400.406
ascm2500.406
ascm2600.406
ascm2700.406
ascm2800.406
ascm2900.406
ascm3000.406
ascp0000.406
ascp0100.406
ascp0200.406
ascp0300.406
ascp0400.406
ascp0500.406
ascp0600.406
ascp0700.406
ascp0800.406
ascp0900.406
ascp1000.406
ascp1100.406
ascp1200.406
ascp1300.406
ascp1400.406
ascp1500.406
ascp1600.406
ascp1700.406
ascp1800.406
ascp1900.406
ascp2000.406
ascp2100.406
ascp2200.406
ascp2300.406
ascp2400.406
ascp2500.406
ascp2600.406
ascp2700.406
ascp2800.406
ascp2900.406
header.406
testpo.406]], '\n'):map(function(f)
	print('downloading '..f)
	--repeat until wget ('ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de406/'..f)	-- lua-sh
	--local cmd = 'wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de406/'..f
	local cmd = 'curl -O ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de406/'..f
	print('>'..cmd)
	print(os.execute(cmd))
end)

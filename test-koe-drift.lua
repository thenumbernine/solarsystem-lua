#!/usr/bin/env luajit
--[[
test how quickly KOE drifts from the NASA Ephemeris 406 data
--]]
local table = require 'ext.table'
local range = require 'ext.range'
local string = require 'ext.string'
local ffi = require 'ffi'
local tolua = require 'ext.tolua'
local vec3d = require 'vec-ffi.vec3d'	-- is a fft metatype, not a class, so there's no "isa" for the record
local Planets = require 'planets'
local julian = require 'julian'
local KOE = require 'koe'

local t = os.date('!*t')
print('current time: '..os.date(nil, os.time(t)))
local currentDate = julian.fromCalendar(t)
print('current julian day: '..currentDate)

-- [=[
-- get the planet positions over a range of dates
-- range is inclusive, which means there's always at least 2 entries, and the table uses 1-based indexes
local startDate = currentDate
--[[ looks good, shows how long it takes to reach a full 100% error
local endDate = startDate + 100000	-- days
local n = 10000
--]]
-- [[ better for estimating error from the avarage, for comparing / fixing koe method
local endDate = startDate + 10000	-- days
local n = 1000
--]]
--[[ debugging i guess
local endDate = startDate + 1-- days
local n = 10
--]]

local planets = table()
local julianDates = table()
for i=0,n do
	local f = i / n
	local julianDate = startDate * (1 - f) + endDate * f
	julianDates:insert(julianDate)
	-- produces: per-planet: {index=number, pos=vec3d, vel=vec3d},
	planets:insert(Planets.fromEphemeris(julianDate, 406, 'eph/406'))
end
--]=]


-- for the first frame, calculate the KOE
local planets1 = planets[1]
for _,planet in ipairs(planets1) do
	planet.koe = KOE.calcKOEFromPosVel(planet, planets1, julianDates[1])
	print()
	print(planet.name)
	print(tolua(planet.koe))
end
print()

local posErrorMag = range(#Planets.planetClasses):mapi(function() return table() end)
local velErrorMag = range(#Planets.planetClasses):mapi(function() return table() end)
assert(Planets.indexes.sun == 1)
for i=1,#planets do	-- skip the 1st time entry cuz thats got our KOE in it
	local planetsi = planets[i]
	-- skip the 1st planet cuz it's the sun and isn't moving
	for j=2,#planetsi do
		local planet = planetsi[j]
		-- mind you the planets[i][j] is the ith Ephemeris frame calculated
		--  it uses this for the parent positions
		--  so maybe I should be feeding it the parent koe positions instead? :shrug:
		KOE.updatePosVel(
			planet,
			planets1[j].koe,
			julianDates[i],
			julianDates[1]
		)
		-- you gotta add the parent KOE pos/vel yourself
		if planet.parentIndex then
			local parent = planetsi[planet.parentIndex]
			while parent do 
				if parent.pos_koe then
					planet.pos_koe = planet.pos_koe + parent.pos_koe
					planet.vel_koe = planet.vel_koe + parent.vel_koe
				else
					-- doesn't have KOE, doens't have parent, is the sun / main barycenter
					--planet.pos_koe = planet.pos_koe + planets1[parentIndex].pos
					--planet.vel_koe = planet.vel_koe + planets1[parentIndex].vel / (60 * 60 * 24)
					break
				end
				parent = planetsi[parent.parentIndex]
			end
		end
		--[[ insert difference magnitude
		posErrorMag[j]:insert((planet.pos_koe - planet.pos):length())
		velErrorMag[j]:insert((planet.vel_koe - planet.vel):length())
		--]]
		-- [[ insert error, i.e. deviation relative to vector mag
		posErrorMag[j]:insert((planet.pos_koe - planet.pos):length() / planet.pos:length())
		velErrorMag[j]:insert((planet.vel_koe - planet.vel / (60 * 60 * 24)):length() / (planet.vel:length() / (60 * 60 * 24)))
		--]]
	end
end
assert(0 == #posErrorMag:remove(Planets.indexes.sun))	-- doesn't have any anyways
assert(0 == #velErrorMag:remove(Planets.indexes.sun))	-- doesn't have any anyways
print(posErrorMag:mapi(function(t) return #t end):concat', ')
print(velErrorMag:mapi(function(t) return #t end):concat', ')

-- now compare planet.pos vs KOE evaluation across the timeframes
local gnuplot = require 'gnuplot'
for _,info in ipairs{
	{src=posErrorMag, outfn='koe_pos_error.svg', name='position'},
	{src=velErrorMag, outfn='koe_vel_error.svg', name='velocity'},
} do
	print()
	print(info.name..' error avgs:')
	for i,mags in pairs(info.src) do
		print(Planets.planetClasses[i+1].name, mags:sum() / #mags)
	end
	gnuplot(
		table(
			{
				terminal = 'svg size 1024,768 background rgb "white"',
				output = info.outfn,
				key = 'left top',
				style = 'data lines',
				title = info.name..' relative error of KOE vs NASA Ephemeris 406',
				xlabel = 'days',
				ylabel = 'rel err',
				log = 'xy',
				data = table()
				:append{
					julianDates:sub(2):mapi(function(d) return d - startDate end)
				}
				:append(info.src),	--velErrorMag,
			},
			range(#info.src):mapi(function(i)
				local planetClass = Planets.planetClasses[i+1]
				-- +1 cuz sun is missing
				-- +1 to skip the date column
				return {
					using = '1:'..(i+1),
					title = planetClass.name,
					linecolor = 'rgbcolor "#'
						..('%02x'):format(math.floor(planetClass.color[1]*255))
						..('%02x'):format(math.floor(planetClass.color[2]*255))
						..('%02x'):format(math.floor(planetClass.color[3]*255))
						..'"',
				}
			end)
		)
	)
end


--[=[ uncomment to do extra checks to compare DE406 to Horizons webservice data.

-- It doesn't hurt to compare here the pos/vel of the Eph406 data vs that downloaded from JPL Horizons.
--[[
ok so it looks like these are pretty severely off.
DE406 is supposed to use the ICRF frame.
Horizons is supposed to use the ICRF frame.
What could be the problem?
The relative errors of the magnitudes of the pos's is very very small.
But the overall error is large.
So it looks like the basis vectors are different.
Hmmmmmm
TODO look at what rotation goes from one to the other in each planet, compare them all,
and compare that to the transform on that one pdf of common frames/transforms.
--]]
local file = require 'ext.file'
local json = require 'dkjson'
local function jsonstrtovec3(t)
	return vec3d(table.mapi(t, function(x) return tonumber(string.trim(x)) end):unpack())
end
local dynamicVarsFileName = 'dynamic-vars.json'
local dynamicVarsJSON = assert(file(dynamicVarsFileName):read())	-- really is js, has a var = at the front
local dynamicVars = assert(json.decode((dynamicVarsJSON:match'^.-=(.*)$')))
print()
local dvPlanets = Planets.fromEphemeris(dynamicVars.julianDate, 406, 'eph/406')
local _, hsun = table.find(dynamicVars.coords, nil, function(o) return o.name:lower() == 'sun' end)
assert(hsun, "couldn't find sun in horizons dynamicVars")
local m = 3
local n = #dvPlanets
local matrix_ffi = require 'matrix.ffi'
local Ephs = matrix_ffi{m,n}:zeros()
local Horzs = matrix_ffi{m,n}:zeros()
local earthTiltInDeg = 23.4365472133
local cosTilt = math.cos(math.rad(earthTiltInDeg))
local sinTilt = math.sin(math.rad(earthTiltInDeg))
local EclipToEquat = matrix_ffi{
	{1,0,0},
	{0,cosTilt,-sinTilt},
	{0,sinTilt,cosTilt},
}
for j,planet in ipairs(dvPlanets) do
	print()
	print(planet.name)
	print('eph406: ', planet.pos, planet.vel)
	print('eph406 magn:', planet.pos:length(), planet.vel:length())
	local _, o = table.find(dynamicVars.coords, nil, function(o)
		return o.name:lower() == planet.name
	end)
	if not o then
		print('...cannot find horizon obj')
	else
		local hpos = jsonstrtovec3(o.pos)
		local hvel = jsonstrtovec3(o.vel)
--		hpos = hpos + jsonstrtovec3(hsun.pos)
--		hvel = hvel + jsonstrtovec3(hsun.vel)
		hpos = hpos * 1000	-- km to m
		hvel = hvel * 1000	-- km/s(day?) to m/s(day?)
		local hposInEquat = vec3d((EclipToEquat * matrix_ffi{hpos:unpack()}):unpack())
		local hvelInEquat = vec3d((EclipToEquat * matrix_ffi{hvel:unpack()}):unpack())
		print('horizon:', hpos, hvel)
		print('horizon magn:', hpos:length(), hvel:length())
		print('horizon vec rel err:', 
			(planet.pos - hposInEquat):length() / planet.pos:length(),
			(planet.vel - hvelInEquat):length() / planet.vel:length())
		-- magnitude relative-error is always within 1e-5
		print('horizon magn rel err:', 
			(planet.pos:length() - hposInEquat:length()) / planet.pos:length(),
			(planet.vel:length() - hvelInEquat:length()) / planet.vel:length())
		for i=1,3 do
			Ephs[{i,j}] = planet.pos.s[i-1]
			Horzs[{i,j}] = hpos.s[i-1]
		end
	end
end
--[[
ok need linear regressin to fin the transform from eph406 to horizons
[e1|e2|...|en] = A [h1|h2|...|hn]
E = A H
E H' = A H H'
A = (E H') * (H H')^-1
--]]
print('Ephs:\n'..Ephs)
print('Horzs:\n'..Horzs)
print('E HT:\n'.. (Ephs * Horzs:T()))
print('H HT:\n'..(Horzs * Horzs:T()))
print('(H HT)^-1:\n'..(Horzs * Horzs:T()):inv())
local HorzsToEphsXform = (
	(Ephs * Horzs:T())
	* (Horzs * Horzs:T()):inv()
)
print('A = (E HT) (H HT)^-1:\n'..HorzsToEphsXform )
print('E - A H:\n'..(Ephs - HorzsToEphsXform * Horzs))
print('avg(E - A H):\n'..(Ephs - HorzsToEphsXform * Horzs):sum()/Ephs:size():prod())
print('|E - A H|/|E|:\n'..(Ephs - HorzsToEphsXform * Horzs):ediv(Ephs))
print('avg(|E - A H|/|E|):\n'..(Ephs - HorzsToEphsXform * Horzs):ediv(Ephs):sum()/Ephs:size():prod())
-- ok looks like this is the xform - which is a YZ rotation - is accurate up to 1e-5 rel error
--[[1.0000022393073, -2.3832669953094e-06, 7.4167233663047e-05],
--[-5.0763288612643e-06, 0.91748539910954, -0.39799691370649],
--[-1.4680298994207e-06, 0.39777818479152, 0.91740958678304]]
-- = +23.439269141797 degrees Y->Z axis rotation
-- so one of this is earth-relative the other is sun-relative?
-- ok in my offline viewer, rotate the grid -23.439° to rotate the grid from sun-aligned (ecliptic) to earth-aligned (equatorial, with polaris up)
-- so rotate objects +23.439° to rotate objects from sun-aligned (ecliptic) to earth-aligned (equatorial)
-- so that means Horizons is ecliptic / sun-aligned, and DE406 is equatorial / earth-aligned
-- ok I think i'm seeing that DE406 is J2000 is "ECI" "earth centered inertial"
-- https://en.wikipedia.org/wiki/Earth-centered_inertial "J2000 ... The z-axis is aligned with the Earth's rotation axis (or equivalently, the celestial North Pole) as it was at that time."
-- still having a hard time tracking down any exact DE406 reference material that says "we're earth-aligned" .... usu they do say "ICRF" and then ICRF says "J2000" ... gah.
-- how come I'm reading my own self ask on stack exchange about the ICRF / J2000 frame pointing from earth->sun as the x-axis ?  wtf?  that's ecliptic. what is up with some sources saying J2000 is equatorial?

--]=]

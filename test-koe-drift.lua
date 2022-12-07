#!/usr/bin/env luajit
--[[
test how quickly KOE drifts from the NASA Ephemeris 406 data
--]]
local table = require 'ext.table'
local ffi = require 'ffi'
local tolua = require 'ext.tolua'
local vec3d = require 'vec-ffi.vec3d'	-- is a fft metatype, not a class, so there's no "isa" for the record
local Planets = require 'planets'
local julian = require 'julian'
local KOE = require 'koe'

local t = os.date('!*t')
print('current time: '..os.date(nil, os.time(t)))
local currentDate = julian.fromCalendar(t)

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
local n = 10000
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
end
local posErrorMag = range(#Planets.planetClasses):mapi(function() return table() end)
local velErrorMag = range(#Planets.planetClasses):mapi(function() return table() end)
assert(Planets.indexes.sun == 1)
for i=2,#planets do	-- skip the 1st time entry cuz thats got our KOE in it
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
-- now compare planet.pos vs KOE evaluation across the timeframes

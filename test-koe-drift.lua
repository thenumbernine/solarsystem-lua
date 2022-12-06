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
local endDate = startDate + 1000	-- days
local n = 10000
--local endDate = startDate + 1-- days
--local n = 10

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


-- TODO put this in planets.lua?
local planetClasses = Planets.planetClasses
for i,planetClass in ipairs(planetClasses) do
	if i == Planets.indexes.sun then
	elseif i == Planets.indexes.moon then
		planetClass.parentIndex = Planets.indexes.earth
	else
		planetClass.parentIndex = Planets.indexes.sun
	end
end

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
		KOE.updatePosVel(planets[1][j], planet, planetsi, julianDates[i], julianDates[1])
		--[[ insert difference magnitude
		posErrorMag[j]:insert((planet.pos_koe - planet.pos):length())
		velErrorMag[j]:insert((planet.vel_koe - planet.vel):length())
		--]]
		-- [[ insert error, i.e. deviation relative to vector mag
		posErrorMag[j]:insert((planet.pos_koe - planet.pos):length() / planet.pos:length())
		velErrorMag[j]:insert((planet.vel_koe - planet.vel / (60 * 60 * 24)):length() / (planet.vel:length() / (60 * 60 * 24)))
		--]]

		-- position appears to be in the same range
--		if j==2 then print(planet.vel_koe, planet.vel) end
--		print(
	
	end
end
assert(0 == #posErrorMag:remove(Planets.indexes.sun))	-- doesn't have any anyways
assert(0 == #velErrorMag:remove(Planets.indexes.sun))	-- doesn't have any anyways
print(posErrorMag:mapi(function(t) return #t end):concat', ')
print(velErrorMag:mapi(function(t) return #t end):concat', ')

local gnuplot = require 'gnuplot'
for _,info in ipairs{
	{src=posErrorMag, outfn='koe_pos_error.svg'},
	{src=velErrorMag, outfn='koe_vel_error.svg'},
} do
	gnuplot(
		table(
			{
				terminal = 'svg size 1024,768',
				output = info.outfn,
				key = 'left top',
				style = 'data lines',
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
				-- +1 cuz sun is missing
				-- +1 to skip the date column
				return {using='1:'..(i+1), title=Planets.planetClasses[i+1].name}
			end)
		)
	)
end
-- now compare planet.pos vs KOE evaluation across the timeframes

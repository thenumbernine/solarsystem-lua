#!/usr/bin/env luajit
local table = require 'ext.table'
local math = require 'ext.math'
local range = require 'ext.range'
local string = require 'ext.string'
local file = require 'ext.file'
local fromlua = require 'ext.fromlua'
local Planets = require 'planets'
local KOE = require 'koe'
local julian = require 'julian'
local gnuplot = require 'gnuplot'
local CSV = require 'csv'
require 'ffi.c.time'
require 'ffi.c.sys.time'
local ffi = require 'ffi'

local timeofday = ffi.new('struct timeval[1]')
local timezone = ffi.new('struct timezone[1]')
timezone[0].tz_minuteswest = 0
timezone[0].tz_dsttime = 0
function currentDate()
	--[[
	ffi.C.gettimeofday(timeofday, timezone)
	timeofday[0].tv_sec = timeofday[0].tv_sec -- - timeZoneOffset
	--]]
	local datetable = os.date('!*t')	--, tonumber(timeofday[0].tv_sec))
	--datetable.sec = datetable.sec + tonumber(timeofday[0].tv_usec) * 1e-6
	return datetable
end

local julianDate = julian.fromCalendar(currentDate())		-- now-ish ... my conversion and my currrentDate are off by a bit ...
local planets = Planets.fromEphemeris(julianDate, 406, 'eph/406')
local earth = planets[planets.indexes.earth]

local planetNames = table.map(Planets.planetClasses, function(cl) return cl.name end)

local function castnumber(x)
	return assert(tonumber(x))
end

local earthquakeEntries = CSV.file'earthquakes_6andover.csv'
earthquakeEntries:setColumnNames(earthquakeEntries.rows:remove(1))

local magns = table()	-- list of earthquake magnitudes
local depths = table()
local validDepths = table()
--local angles = {}	-- per-planet-pair angle at earth
--local prodSinSqAngles = table()
--local avgCosSqAngles = table()
--local avgAngles = table()
local minAngles = table()
local maxAngles = table()

local angleGraphsEnabledForName = {
	sun_earth_moon = true,
}

for ei,e in ipairs(earthquakeEntries.rows) do
	local year, month, day, hour, min, sec = e.time:match'^(%d%d%d%d)%-(%d%d)%-(%d%d)T(%d%d):(%d%d):(%d%d%.%d%d%d)Z$'
	local date = {
		year = castnumber(year),
		month = castnumber(month),
		day = castnumber(day),
		hour = castnumber(hour),
		min = castnumber(min),
		sec = castnumber(sec),
	}
	local time = os.time(date)
	local lat = castnumber(e.latitude)
	local lon = castnumber(e.longitude)
	local mag = castnumber(e.mag)
	local depth = tonumber(e.depth)
	local julianDate = julian.fromCalendar(date)
	local planets = Planets.fromEphemeris(julianDate, 406, 'eph/406')

	magns[ei] = mag
	depths[ei] = depth or '-'
	if depth then validDepths:insert(ei) end

	--local count = 0
	--local prodSinSqAngle = 1
	--local sumAngle = 0
	--local sumCosSqAngle = 0
	local minAngle = math.huge
	local maxAngle = -math.huge
--[[ considering all angles of all planets
	for j=1,#planetNames do
		local nj = planetNames[j]
--]]
-- [[ considering just angles of planets to earth
	do
		local j = Planets.indexes.earth
		local nj = 'earth'
--]]
		for i=1,#planetNames-1 do
			if i ~= j then
				local ni = planetNames[i]
				for k=i+1,#planetNames do
					if k ~= j then
						local nk = planetNames[k]
						local name = table{ni,nj,nk}:concat'_'
						local va = planets[Planets.indexes[ni]].pos - planets[Planets.indexes[nj]].pos
						local vb = planets[Planets.indexes[nk]].pos - planets[Planets.indexes[nj]].pos
						va = va:normalize()
						vb = vb:normalize()
						local costh = math.clamp(va:dot(vb), -1, 1)
						local theta = math.acos(costh)
						local thetadeg = math.deg(theta)
						--[[
						if angleGraphsEnabledForName[name] then
							angles[name] = angles[name] or table()
							angles[name]:insert(thetadeg)
						end
						--]]
						-- too close to zero
						--prodSinSqAngle = prodSinSqAngle * math.sin(theta) * math.sin(theta)
						--sumAngle = sumAngle + thetadeg
						minAngle = math.min(minAngle, thetadeg)
						maxAngle = math.max(maxAngle, thetadeg)
						--sumCosAngle = sumCosAngle + costh
						--sumCosSqAngle = sumCosSqAngle + costh * costh
						--count = count + 1
					end
				end
			end
		end
	end
	--prodSinSqAngles:insert(prodSinSqAngle^(1/count))
	--avgCosSqAngles:insert(sumCosSqAngle/count)
	--avgAngles:insert(sumAngle/count)
	minAngles:insert(minAngle)
	maxAngles:insert(maxAngle)
end

-- [[ plot earthquake depth vs magnitude
local magnsForDepths = validDepths:mapi(function(i) return magns[i] end)
local depthsForDepths = validDepths:mapi(function(i) return depths[i] end)
gnuplot{
	output = 'quake-magn-vs-depth.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'depth',
	ylabel = 'magn',
	data = {depthsForDepths, magnsForDepths},
	log = 'xy',
	xrange = {depthsForDepths:inf()+.001, depthsForDepths:sup()+.001},
	yrange = {magnsForDepths:inf(), (magnsForDepths:sup())},
	{using='($1+.001):2', notitle=true},
}
--]]

--[[
how does gnuplot handle pm3d?
are the x y grid values at cell center or edges?
edges it looks.
but then what is the index range of griddata?
x y = values at edges <=> griddata only needs [1..bins] or [0..bins-1] ...
... but it looks like griddata wants [0..bins]

... but for image, it is cell-centered?
--]]
local function bin(xs, ys, xbins, ybins, xext, yext)
	assert(#xs == #ys)
	local xmin = xs:inf()
	local xmax = xs:sup()
	local ymin = ys:inf()
	local ymax = ys:sup()
	if xext then
		local dx = .5 * (xmax - xmin) / xbins
		xmin = xmin - dx
		xmax = xmax + dx
	end
	if yext then
		local dy = .5 * (ymax - ymin) / ybins
		ymin = ymin - dy
		ymax = ymax + dy
	end
	local binxs = range(xbins):mapi(function(i)
		return (i-.5)/(xbins-1) * (xmax - xmin) + xmin
	end)
	local binys = range(ybins):mapi(function(i)
		return (i-.5)/(ybins-1) * (ymax - ymin) + ymin
	end)
	local griddata = {}
	for i=1,xbins do
		griddata[i] = griddata[i] or {}
		for j=1,ybins do
			griddata[i][j] = 0
		end
	end
	for i=1,#xs do
		local xi = math.floor((xs[i] - xmin) / (xmax - xmin) * xbins) + 1
		if xi == xbins+1 then xi = xbins end
		local yi = math.floor((ys[i] - ymin) / (ymax - ymin) * ybins) + 1
		if yi == ybins+1 then yi = ybins end
		griddata[xi][yi] = griddata[xi][yi] + 1
	end
	return {
		xrange = {xmin, xmax},
		yrange = {ymin, ymax},
		griddata = {
			x = binxs,
			y = binys,
			griddata,
		},
	}
end

--[[ pointwise plot 
gnuplot{
	output = 'quake-magn-vs-angle-sun-earth-moon.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'angle',
	ylabel = 'magn',
	data = {angles.sun_earth_moon, magns},
	--log = 'xy',
	--xrange = {magnsForDepths:inf(), (magnsForDepths:sup())},
	--yrange = {depthsForDepths:inf()+.001, depthsForDepths:sup()+.001},
	{using='1:2', notitle=true},
}
--]]

--[[ density plot 
gnuplot(table(
bin(angles.sun_earth_moon, magns, 50, 30, false, false),
{
	output = 'quake-magn-vs-angle-sun-earth-moon-density.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	view = 'map',
	xlabel = 'angle',
	ylabel = 'magn',
	'max(a,b) = a > b ? a : b',
	'set palette rgb 33,13,10',
	{splot=true, using='1:2:(log(max($3,.001)))', with='image pixels', notitle=true},
}))
--]]

--[[ still zero ...
gnuplot{
	output = 'quake-magn-vs-prod-sinsq-angles.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'prod(sin(x)^2)',
	ylabel = 'magn',
	data = {prodSinSqAngles, magns},
	{using='1:2', notitle=true},
}
--]]

--[[
gnuplot{
	output = 'quake-magn-vs-avg-cossq-angles.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'avg(cos(x)^2)',
	ylabel = 'magn',
	data = {avgCosSqAngles, magns},
	{using='1:2', notitle=true},
}
--]]

--[[
gnuplot{
	output = 'quake-magn-vs-avg-angle.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'avg(angle)',
	ylabel = 'magn',
	data = {avgAngles, magns},
	{using='1:2', notitle=true},
}
--]]

-- [[
print('min angle range', minAngles:inf(), (minAngles:sup()))
gnuplot{
	output = 'quake-magn-vs-min-angle.svg',
	terminal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'min(angle)',
	ylabel = 'magn',
	data = {minAngles, magns},
	{using='1:2', notitle=true},
	savedata='minangles.txt',
}
--]]

-- [[
print('max angle range', maxAngles:inf(), (maxAngles:sup()))
gnuplot{
	output = 'quake-magn-vs-max-angle.svg',
	termaxal = 'svg size 1024,768 background rgb "white"',
	style = 'data points',
	xlabel = 'max(angle)',
	ylabel = 'magn',
	data = {maxAngles, magns},
	{using='1:2', notitle=true},
	savedata='maxangles.txt',
}
--]]

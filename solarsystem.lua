#!/usr/bin/env luajit
--[[
TODO:

- single script for downloading earthquake data from usgs

interface:
- toggle wireframe
- toggle texture vs force display
- force display options:
	- linear combinations of: 
		- tide (per planet)
		- gravity (per planet)
	- decomposition
		normal
		tangential
		
orbit locations

events and extra geometry:
	- total eclipses (and penumbra?)
	- earthquakes
	- how should the user pick a space and a time ?
	
- fix the julian/gregorian conversion routine
	- maybe add SOFA support too ...
	
- add mapm support or at least qd
--]]


local Planets = require 'planets'
local KOE = require 'koe'
local julian = require 'julian'
local ImGuiApp = require 'imguiapp'
local Mouse = require 'glapp.mouse'
local gl = require 'gl'
local ig = require 'imgui'
local sdl = require 'ffi.sdl'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
local Tex2D = require 'gl.tex2d'
local HsvTex = require 'gl.hsvtex'
require 'ffi.c.time'
require 'ffi.c.sys.time'
local ffi = require 'ffi'
local table = require 'ext.table'
local range = require 'ext.range'
local class = require 'ext.class'
local tolua = require 'ext.tolua'
local fromlua = require 'ext.fromlua'
local file = require 'ext.file'
--local shader = require 'gl.shader'

local planets = Planets()
local earth = planets[planets.indexes.earth]


--[[
local now = ffi.new('time_t[1]')
ffi.C.time(now)
local ts = ffi.C.localtime(now)[0]
--print(ts.tm_year + 1900, ts.tm_mon + 1, ts.tm_mday, ts.tm_hour, ts.tm_min, ts.tm_sec, ts.tm_gmtoff)
local timeZoneOffset = tonumber(ts.tm_gmtoff)
print(timeZoneOffset)
--]]
local timeZoneOffset = 0


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


local planets
local julianDate = 0
local targetJulianDate
--[[ astro-phys
do
	local json = require 'dkjson'
	local dateTable = currentDate()
	local dateTime = os.time(dateTable)
	local fn = 'state.json'
	local d
	if file(fn):exists() then
		d = assert(file(fn):read())
		d = assert(json.decode(d))
	else
		local http = require 'socket.http'
		local url = require 'socket.url'
		local dateStr = os.date('%Y-%m-%d %H:%M:%S', dateTime)
		d = http.request('http://www.astro-phys.com/api/de406/states?date='..url.encode(dateStr)..'&bodies='..table.concat(Planets.names, ','))
		d = assert(json.decode(d))
		d.calendarDate = dateTable
		d.linuxTime = dateTime
		file(fn):write(json.encode(d, {indent=true}))
	end
	planets = Planets.fromAstroPhys(d.results)
	julianDate = d.date	-- julian date
end
--]]
-- ephemeris
-- I'm getting these numbers from http://aa.usno.navy.mil/data/docs/JulianDate.php
--julianDate = 2456049.41235	-- matches some data I had pulled from astr-phys
--julianDate = julian.fromCalendar(currentDate())
--julianDate = 2455034.608623	--julian.fromCalendar{year=2009, month=7, day=22, hour=2, min=36, sec=25}	-- eclipse
--julianDate = 2455389.315718		--julian.fromCalendar{year=2010, month=7, day=11, hour=19, min=34, sec=38}	-- eclipse
julianDate = julian.fromCalendar(currentDate())		-- now-ish ... my conversion and my currrentDate are off by a bit ...


local startDate = julianDate	-- used for reset
planets = Planets.fromEphemeris(julianDate, 406, 'eph/406')

-- [[ get earthquake data
--local earthquakeEntries = setmetatable(assert(fromlua(assert(file'earthquakes.lua':read()))), table)
--local solarEclipseEntries = table((assert(json.decode(assert(file'solar-eclipses.json':read())))))
--]]



--[[ 
local planet = planets[planets.indexes.earth]
--planet.class.inverseFlattening = nil	-- in absence of inverse flattening, this goes closer to exactly 0.3 * circ
--local pt1 = {29.9791953,31.1320557,17}		-- giza, egypt
local pt1 = {37.2231663,38.9202232,17}		-- gobekli tepe, turkey
local pt2 = {-13.5300169,-71.9742687,14067}	-- cuzco, peru
--local pt2 = {-13.1631412,-72.5471516,17}		-- machu picchu, peru
pt1 = vec3d(planet:geodeticPosition(table.unpack(pt1)))
pt2 = vec3d(planet:geodeticPosition(table.unpack(pt2)))
print('pt1', pt1)
print('pt2', pt2)
local costh = pt1:normalize():dot(pt2:normalize())
print('cos theta', costh)
print('theta', math.acos(costh))
print('theta/(2pi)', math.acos(costh) / (2 * math.pi))

os.exit()
--]]

-- event database ... TODO list this out
local events = table()
if earthquakeEntries then
	events:append(earthquakeEntries:map(function(earthquake, index, newTable)
		local res, time = pcall(os.time, earthquake)
		if res then
			return {
				lat = earthquake.lat,
				lon = earthquake.lon,
				height = -(earthquake.depth or 0),
				julianDate = julian.fromCalendar(earthquake),	-- TODO fixme plz
				date = earthquake,
				earthquake = earthquake,
				name = 'eq '..earthquake.momentMagn..'±'..earthquake.momentMagnVar,
			}, #newTable+1
		end
	end))
end
if solarEclipseEntries then
	events:append(solarEclipseEntries:map(function(eclipse, index, newTable)
		local res, time = pcall(os.time, eclipse)
		if res then
			return {
				lat = eclipse.lat,
				lon = eclipse.lon,
				height = 0,
				julianDate = julian.fromCalendar(eclipse),
				date = eclipse,
				eclipse = eclipse,
				name = 'Total Solar Eclipse',
			}, #newTable+1
		end
	end))
end

-- track ball motion variables
local orbitPlanetIndex
local orbitGeodeticLocation
local orbitDistance
local orbitTargetDistance
local orbitZoomFactor = .9	-- upon mousewheel

--orbitGeodeticLocation = {lat=45.265837, lon=-123.834667, height=0}	--my house


-- view state & transformation
local viewPos = vec3d()
local viewAngle = quatd(0,0,0,1)

local mouse = Mouse()
local leftButtonDown


-- planet barycenter orbit periods, in days
local orbitPeriods = {
	sun = 0,
	mercury = 87.92940628134535,
	venus = 224.57167537062335,
	earth = 365.256363004,
	
	moon = 365.256363004,	-- moon around sun, or moon around earth?
	
	mars = 686.6292080074874,
	jupiter = 4331.20766313597,
	saturn = 10824.232418862814,
	uranus = 30668.842130815046,
	neptune = 60608.875026272224,
	pluto = 91856.29710876015,
}

local showTrail = {
	sun = false,
	mercury = false,
	venus = false,
	earth = false,
	
	moon = false,
	
	mars = false,
	jupiter = false,
	saturn = false,
	uranus = false,
	neptune = false,
	pluto = false,
}

-- globals
--numCycles = 12	-- more looks better esp for planet-centered rendering of other planets' orbits
numCycles = 1
showOrbitWrtPlanet = false
invertRadialDistance = false
calcKOE = false
koeInfo = nil

local integrateTimeStep		-- in days, per iteration

local gravitationConstant = 6.6738480e-11	-- m^3 / (kg * s^20

local tidalMin = 0
local tidalMax = 1
local function calcTidalForce(srcPlanet, pos)
	local accel = vec3d()
	local srcPlanetToPos = pos - srcPlanet.pos
	for _,planet in ipairs(planets) do
		if planet ~= srcPlanet then
			local x = pos - planet.pos
			local xLength = x:length()
			local xToTheSecond = xLength * xLength
			local xToTheThird = xLength * xToTheSecond
			local xToTheFourth = xLength * xToTheThird
			local xToTheFifth = xLength * xToTheFourth
			
			-- a^i = -R^i_jkl dt^j dx^k dt^l = -R^i_tjt = R^i_ttj = -phi_,ij
			-- looks like dt^j = [1,0,0,0] so that we only get the t part of R (which is zero)
			-- but what if phi changes wrt time? then phi_,tt is nonzero, and how does our Riemann metric change?
			for i=1,3 do
				for j=1,3 do
					if i == j then
						phi_ij = gravitationConstant * planet.mass * (3 * x[i] * x[i] / xToTheFifth - 1 / xToTheThird)
					else
						phi_ij = gravitationConstant * planet.mass * (3 * x[i] * x[j]) / xToTheFifth
					end
					accel[i] = accel[i] - phi_ij * srcPlanetToPos[j]
				end
			end
		end
	end
	return accel
end

local function calcGravitationForce(pos)
	local accel = vec3d()
	for _,planet in ipairs(planets) do	
		local x = pos - planet.pos
		local xLength = x:length()
		local xToTheSecond = xLength * xLength
		local xToTheThird = xLength * xToTheSecond
		accel = accel - x * (gravitationConstant * planet.mass / xToTheThird)
	end
	return accel
end

--[[
local integrationMethod = 'rkf45'
local integrationArgs = {accuracy=1, iterations=100, norm=Planets.length}
--]]
local integrationMethod = 'rk4'
local integrationArgs = {accuracy=1, iterations=100, norm=Planets.length}

local integrateFunction = function(time, planets)
	local secondsPerJulianDay = 86400
	-- m^3 / (kg * julianday^2) = 6.6738480e-11 m^3 / (kg * s^2) * (second/julianday)^2
	local G = gravitationConstant * secondsPerJulianDay * secondsPerJulianDay	-- gravitationConstant is proprotional to s^-2 and our time is in julian days, so convert
	-- f'' = G m1 m2 / r^2 = rdir G m1 m2 / r^3
	-- x1'' = rhat G m2 / r^3
	local deltaPlanets = Planets()
	for i=1,#planets do
		local pi = planets[i]
		local accel = vec3d()
		for j=1,#planets do
			local pj = planets[j]
			if i ~= j then
				local delta = pj.pos - pi.pos
				local deltaLen = delta:length()
				local scalar = G * pj.mass / (deltaLen * deltaLen * deltaLen)
				accel = accel + delta * scalar
			end
		end
		deltaPlanets[i].pos = vec3d(pi.vel)
		deltaPlanets[i].vel = accel
	end
	return deltaPlanets
end

local function planetCartesianToSolarSystemBarycentric(planet, x)
	x = planet.angle:rotate(x)		-- right now planet angles aren't stored in the planet state (for adaptive integration's sake -- how to weight and combine time + space measurements)
	
	-- now rotate by axial tilt (along
	local tiltAngle = planet.tiltAngle
	if tiltAngle then
		x = tiltAngle:rotate(x)
	end
	
	x = x + planet.pos

	return x
end

local function solarSystemBarycentricToPlanetCartesian(planet, x)
	x = x - planet.pos
	if planet.tiltAngle then
		x = planet.tiltAngle:conjugate():rotate(x)
	end
	x = planet.angle:conjugate():rotate(x)
	return x
end

local function planetGeodeticToSolarSystemBarycentric(planet, lat, lon, height)
	local x = vec3d(planet:geodeticPosition(lat, lon, height))		-- position relative to the planet center
	x = planetCartesianToSolarSystemBarycentric(planet, x)
	return x
end

local function solarSystemBarycentricToPlanetGeodetic(planet, x)
	x = solarSystemBarycentricToPlanetCartesian(planet, x)
	return planet:cartesianToGeodetic(x:unpack())
end

local hsvTex
local colorBarWidth = 50	-- in menu units
local colorBarHSVRange = 2/3	-- how much of the rainbow to use
local drawWireframe = false

showTide = false
local eventText
local mouseOverEvent


local quad = {{0,0},{0,1},{1,1},{1,0}}
local latMin, latMax, latStep = -90, 90, 5
local lonMin, lonMax, lonStep = -180, 180, 5

-- for building the display list.  only use time-independent variables
local function drawPlanetPrims(planet)
	gl.glBegin(gl.GL_QUADS)
	for baselon=lonMin,lonMax-lonStep,lonStep do
		for baselat=latMin,latMax-latStep,latStep do
			for _,ofs in ipairs(quad) do
				local lat = baselat + latStep * ofs[1]
				local lon = baselon + lonStep * ofs[2]
				gl.glTexCoord2f(lon / 360 + .5, -lat / 180 + .5)
				gl.glVertex3d(planet:geodeticPosition(lat, lon, 0))
			end
		end
	end
	gl.glEnd()
end

local function drawPlanetMesh(planet, tccoords, tcbuf)
	--[[	call list based
	if not planet.class.list then planet.class.list = {} end	-- because each planet has its own class .. and the planet objects themselves are thrown away
	glCallOrRun(planet.list, drawPlanetPrims, planet)
	--]]
	
	-- [[	array based
	gl.glVertexPointer(3, gl.GL_DOUBLE, 0, planet.vertexArray)
	gl.glTexCoordPointer(tccoords or 2, gl.GL_DOUBLE, 0, tcbuf or planet.texCoordArray)

	gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
	gl.glEnableClientState(gl.GL_TEXTURE_COORD_ARRAY)

	gl.glDrawElements(gl.GL_QUADS, planet.elementCount, gl.GL_UNSIGNED_INT, planet.elementArray)

	gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
	gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
	--]]
end

local function drawPlanet(planet)

	gl.glPushMatrix()
	gl.glTranslated(planet.pos:unpack())
	local aa = planet.angle:toAngleAxis()
	gl.glRotated(aa.w, aa.x, aa.y, aa.z)
	gl.glEnable(gl.GL_BLEND)
	gl.glDisable(gl.GL_DEPTH_TEST)
	
	if showTide then
		if planet.class.lastTideCalcDate ~= julianDate then
			planet.class.lastTideCalcDate = julianDate
			-- calc values and establish range
			tidalMin = nil
			tidalMax = nil
			for i=0,planet.vertexCount-1 do
				local planetX = vec3d(planet.vertexArray[3*i+0], planet.vertexArray[3*i+1], planet.vertexArray[3*i+2])
				local x = planetCartesianToSolarSystemBarycentric(planet, planetX)
				local accel = calcTidalForce(planet, x)
				--local accel = calcGravitationForce(x)
				local norm = (x - planet.pos):normalize()
				--local toTheMoon = (planets[planets.indexes.moon].pos - x):normalize()
				--local normCrossMoon = norm:cross(toTheMoon)	--points upwards, tangent, right angle to norm and moon
				--local tangentTowardsMoon = normCrossMoon:cross(norm)
				--local tidalAccel = accel:dot(tangentTowardsMoon)
				--local tidalAccel = accel:dot(norm)	-- normal component
				--local tidalAccel = accel:dot(toTheMoon)	-- moonward component
				local tidalAccel = (accel - norm * accel:dot(norm)):length()	-- tangent component
				--local tidalAccel = accel:length()	-- magnitude
				local t = tidalAccel
				if not tidalMin or t < tidalMin then tidalMin = t end
				if not tidalMax or t > tidalMax then tidalMax = t end
				planet.tideArray[i] = t
			end
			-- remap to color spectrum
			for i=0,planet.vertexCount-1 do
				planet.tideArray[i] = (255/256 - (planet.tideArray[i] - tidalMin) / (tidalMax - tidalMin) * 254/256) * colorBarHSVRange
			end
		end
		
		gl.glColor3f(1,1,1)
		hsvTex:bind()
		hsvTex:enable()
		drawPlanetMesh(planet, 1, planet.tideArray)
		hsvTex:disable()
	elseif planet.tex then	-- draw with tex
		gl.glColor3f(1,1,1)
		planet.tex:enable()
		planet.tex:bind()
		drawPlanetMesh(planet)
		planet.tex:disable()
	else	-- draw without tex
		gl.glColor4f(planet.color[1], planet.color[2], planet.color[3], .3)
		drawPlanetMesh(planet)
	end

	if drawWireframe and planet.visRatio > .05 then
		gl.glColor4f(1,1,1, .1)
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
		gl.glDisable(gl.GL_CULL_FACE)
		drawPlanetMesh(planet)
		gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
		gl.glEnable(gl.GL_CULL_FACE)
		gl.glColor3f(table.unpack(planet.color))
	end

	gl.glDisable(gl.GL_BLEND)
	gl.glEnable(gl.GL_DEPTH_TEST)

	gl.glPopMatrix()
end

local planetHistoryIndex = 1
local planetHistoryMax = 1000
local planetHistoryDelta = 1
local lastHistoryJulianDate = 0
local planetHistory = table()	-- circular buffer

-- gl frustum info
local zNear = 1
local zFar = 1e+15
fov = 90

local dateText
local solarSystemApp

local function mouseRay(mousePos)
	local w, h = solarSystemApp:size()
	local aspectRatio = w / h
	-- ray intersect
	local fx = mousePos.x * 2 - 1
	local fy = mousePos.y * 2 - 1
	local tanFov = math.tan(math.rad(fov/2))
	local v = vec3d(fx * aspectRatio * tanFov, fy * tanFov, -1)
	v = viewAngle:rotate(v)
	v = v:normalize()
	return v
end

local function chooseNewPlanet(mouseDir,doChoose)
	local bestDot = 0
	local bestPlanet = nil
	for i=1,#planets do
		local planet = planets[i]
		local delta = (planet.pos - viewPos):normalize()
		local dot = delta:dot(mouseDir)
		if dot > bestDot then 
			bestDot = dot
			bestPlanet = planet
		end
	end
	if bestPlanet then
		if bestPlanet.index ~= orbitPlanetIndex and doChoose then
			orbitPlanetIndex = bestPlanet.index
			orbitDistance = (viewPos - bestPlanet.pos):length()
			orbitTargetDistance = 2 * bestPlanet.radius
		end
	end
end

-- routines for calculating the earth's angle ...
-- ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/astrtime.cc
local function getLeapSec(mjd, utc)
	-- oh great there is a leap-second file...
	return 0
end
local function tai2utc(taimjd, tai)
	local mjd = taimjd
	local ls = getLeapSec(taimjd, tai)
	local utc = tai - ls / 86400
	if utc < 0 then 
		utc = utc + 1
		mjd = mjd - 1
	end
	if ls > getLeapSec(mjd, utc) + .00001 then
		utc = utc + 1/86400
		if utc > 1 then
			utc = utc - 1
			mjd = mjd + 1
		end
	end
	return mjd, utc
end
local function tt2utc(ttmjd, tt)
	local tai = tt - 32.184 / 86400
	local taimjd = ttmjd
	if tai < 0 then	
		tai = tai + 1
		taimjd = taimjd - 1
	end
	return tai2utc(taimjd, tai)
end
local function getUt1Offset(mjd, utc)
	-- yet another absurd time offset file
	return 0
end
local function utc2ut1(mjd, utc)
	local ut1mjd = mjd
	local ut1 = utc + getUt1Offset(mjd, utc) / 86400
	if ut1 > 1 then
		ut1 = ut1 - 1
		ut1mjd = ut1mjd + 1
	end
	if ut1 < 0 then
		ut1 = ut1 + 1
		ut1mjd = ut1mjd - 1
	end
	return ut1mjd, ut1
end
local function iauGmst82() return 0 end -- sofa
local function utc2gmst(mjd, utc)
	local ut1mjd, ut1 = utc2ut1(mjd, utc)
	return iauGmst82(2400000.5, ut1mjd + ut1) / (2 * math.pi)
end
local function utc2tai(mjd, utc)
	local taimjd = mjd
	local tai = utc
	tai = tai + getLeapSec(mjd, utc) / 86400
	if tai > 1 then
		tai = tai - 1
		taimjd = taimjd + 1
	end
	return taimjd, tai
end
local function utc2tt(mjd, utc)
	return utc2tai(mjd, utc)
end
local function iauEqeq94() return 0 end	-- sofa
local function utc2gast(mjd, utc)
	local gmst = utc2gmst(mjd, utc)
	local ttmjd, tt = utc2tt(mjd, utc)
	return gmst + iauEqeq94(2400000.5, ttmjd + tt) / (2 * math.pi)
end
local function utc2last(mjd,utc,lon)
	local lon = 0	-- longitude of our observatory ... from space!
	local last = utc2gast(mjd, utc) + lon / (2 * math.pi)
	if last < 0 then last = last + 1 end	--modulo?
	if last > 1 then last = last - 1 end
	return last
end
local function getEarthAngle(jd)
	local jdofs = julianDate - 2400000.5
	local ttmjd = math.floor(jdofs)
	local tt = jdofs - ttmjd
	local mjd, utc = tt2utc(ttmjd, tt)
	local arg = utc2last(mjd,utc) * math.pi * 2
	return deg
end

local function average(t)
	return t:sum() / #t
end

-- TODO also save what planet it is?
local allArcs = table()
do
	local data = file'arcs.luon':read()
	if data then
		allArcs = table(fromlua(data))
	end
end
local function saveArcs()
	file'arcs.luon':write(tolua(allArcs))
end

function drawScene(viewScale, mouseDir)
	gl.glClear(gl.GL_DEPTH_BUFFER_BIT)

	local viewWidth, viewHeight = solarSystemApp:size()
	local aspectRatio = viewWidth / viewHeight
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	local tanFov = math.tan(math.rad(fov/2))
	gl.glFrustum(-zNear * aspectRatio * tanFov, zNear * aspectRatio * tanFov, -zNear * tanFov, zNear * tanFov, zNear, zFar);

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glScaled(viewScale, viewScale, viewScale)
	local aa = viewAngle:toAngleAxis()
	gl.glRotated(-aa.w, aa.x, aa.y, aa.z)
	gl.glTranslated((-viewPos):unpack())


	-- [[
	if historyCacheDate ~= julianDate 
	or historyCacheOrbitPlanetIndex ~= orbitPlanetIndex
	then
		historyCacheDate = julianDate
		historyCacheOrbitPlanetIndex = orbitPlanetIndex
		historyCache = table.mapi(planets, function(planet,i)
			local divs = 360
			local period = assert(orbitPeriods[planet.name], "failed to find orbit period for "..planet.name)
-- how many orbits?
divs = divs * numCycles
period = period * numCycles
			local poss = range(divs):mapi(function(j)
				local tailPlanets = Planets.fromEphemeris(julianDate - ((j - 1)/(divs-1)) * period)
				local pos = tailPlanets[i].pos
				
				-- for orbits centered on a specific planet ... replace planet[i][t].pos with planet[i][t].pos - targetPlanet[t].pos + targetPlanet[now].pos
				if showOrbitWrtPlanet then
					pos = pos - tailPlanets[orbitPlanetIndex].pos 
				end

				return pos
			end)
			
			if showOrbitWrtPlanet then
				local avgLen
				if invertRadialDistance then
					avgLen = average(poss:mapi(function(pos) return pos:length() end))
				end

				poss = poss:mapi(function(pos)
					if invertRadialDistance then
						pos = pos / avgLen
						pos = pos / pos:lenSq()
						pos = pos * avgLen
					end
					pos = pos + planets[orbitPlanetIndex].pos
					return pos
				end)
			end

			return poss
		end)
	end
	--]]


	-- draw state
	gl.glPointSize(4)
	for i,planet in ipairs(planets) do
		gl.glColor3f(table.unpack(planet.color))
	
		-- [[ trace history
		gl.glLineWidth(2)
		gl.glBegin(gl.GL_LINE_STRIP)
		gl.glVertex3d(planet.pos:unpack())

		local j = (planetHistoryIndex - 2) % planetHistoryMax + 1
		while j ~= planetHistoryIndex and planetHistory[j] do
			gl.glVertex3d(planetHistory[j][i].pos:unpack())
			j = (j - 2) % planetHistoryMax + 1
		end
		gl.glEnd()
		gl.glLineWidth(1)
		--]]
	
		-- [[
		for i=1,#planets do
			if showTrail[planets[i].name] then
				local c = planets[i].color
				gl.glBegin(gl.GL_LINE_STRIP)
				local n = #historyCache[i]
				for j,pos in ipairs(historyCache[i]) do
					local f = (j-1)/(n-1)
					local s = 1-f
					gl.glColor3f(f*c[1], f*c[2],f*c[3])
					gl.glVertex3d(pos:unpack())
				end
				gl.glEnd()
			end
		end
		--]]
		
		-- koe orbit?
		if calcKOE then
			for i=1,#planets do
				local c = planets[i].color
				local koe = koeInfo.koe[i]
				local history = historyCache[i]
				if koe then
					local koeFrames = koeInfo.frames[i]
					local n = #koeFrames
					gl.glBegin(gl.GL_LINE_STRIP)
					for j,frame in ipairs(koeFrames) do
						local f = (j-1)/(n-1)
						local s = 1-f
						gl.glColor3f(f*c[1], f*c[2],f*c[3])
						gl.glVertex3d(frame.pos_koe:unpack())
					end
					gl.glEnd()
					if history then
						assert(#koeFrames == #history)
						gl.glBegin(gl.GL_LINES)
						for j=1,n do
							local f = (j-1)/(n-1)
							local s = 1-f
							gl.glColor3f(f*c[1], f*c[2],f*c[3])
							gl.glVertex3d(koeFrames[j].pos_koe:unpack())
							gl.glVertex3d(history[j]:unpack())
						end
						gl.glEnd()
					end
				end
			
			end	
		end

		gl.glColor3f(table.unpack(planet.color))

		planet.visRatio = planet.radius / (planet.pos - viewPos):length()
		if planet.visRatio >= .005 then
			-- draw sphere
			drawPlanet(planet)
		end
		
		-- [[ and a point just in case
		gl.glBegin(gl.GL_POINTS)
		gl.glVertex3d(planet.pos:unpack())
		gl.glEnd()
		--]]
	end
	
	do
		local earth = planets[planets.indexes.earth]
		if earth.visRatio >= .05 then
			local bestMouseDot = .99
			local lastMouseOverEvent = mouseOverEvent
			mouseOverEvent = nil
		
			gl.glBegin(gl.GL_POINTS)
			for _,event in ipairs(events) do
				if not event.pos then
					event.pos = vec3d(earth:geodeticPosition(event.lat, event.lon, event.height))
				end
				local ssPos = planetCartesianToSolarSystemBarycentric(earth, event.pos)
				local viewDelta = (ssPos - earth.pos):normalize()
				local viewDir = (ssPos - viewPos):normalize()
				local viewDot = viewDelta:dot(viewDir)
				if viewDot <= 0 then
					local mouseDot = viewDir:dot(mouseDir)
					if not bestMouseDot or mouseDot > bestMouseDot then
						bestMouseDot = mouseDot
						mouseOverEvent = event
					end
					if event == lastMouseOverEvent then
						gl.glColor3f(1,1,0)
					else
						gl.glColor3f(1,0,1)
					end
					gl.glVertex3d(ssPos:unpack())
				end
			end
			gl.glEnd()

			gl.glEnable(gl.GL_BLEND)
			local viewFwd = -viewAngle:zAxis()
			for _,arc in ipairs(allArcs) do
				local a = vec3d(earth:geodeticPosition(table.unpack(arc[1])))
				local b = vec3d(earth:geodeticPosition(table.unpack(arc[2])))
				local axis = a:cross(b):normalize()
				local n = 60
				local r = quatd():fromAngleAxis(axis.x, axis.y, axis.z, 360/n)
				gl.glBegin(gl.GL_LINE_LOOP)
				for i=1,n do
					local l = a:normalize():dot(viewFwd)
					local alpha = .5 - l * .5
					gl.glColor4f(1,1,1, alpha * alpha * alpha)

					gl.glVertex3d(planetCartesianToSolarSystemBarycentric(earth, a):unpack())
					a = r:rotate(a)
				end
				gl.glEnd()
				gl.glBegin(gl.GL_POINTS)
				gl.glVertex3d(planetCartesianToSolarSystemBarycentric(earth, a):unpack())
				gl.glVertex3d(planetCartesianToSolarSystemBarycentric(earth, b):unpack())
				gl.glEnd()
			end
			gl.glDisable(gl.GL_BLEND)
			
			if mouseOverEvent then
				eventText = mouseOverEvent.name..' '..os.date(nil, os.time(mouseOverEvent.date))
			else
				eventText = nil
			end
		end
		
		gl.glDisable(gl.GL_DEPTH_TEST)
		gl.glEnable(gl.GL_BLEND)
		
		-- line from moon through earth in direction of the sun (useful for eclipses)
		gl.glColor3f(1, .5, 0)
		gl.glBegin(gl.GL_LINES)
		do
			local moon = planets[planets.indexes.moon]
			local sun = planets[planets.indexes.sun]
			local sunToMoon = moon.pos - sun.pos
			local moonToEarth = moon.pos - earth.pos
			gl.glVertex3d(moon.pos:unpack())
			gl.glVertex3d((moon.pos + sunToMoon * (moonToEarth:length() / sunToMoon:length()) ):unpack())
		end
		gl.glEnd()
		gl.glDisable(gl.GL_BLEND)
		gl.glEnable(gl.GL_DEPTH_TEST)
	end
	
	gl.glPointSize(1)

end



local SolarSystemApp = class(ImGuiApp)

SolarSystemApp.title = 'Solar System Simulation'

function SolarSystemApp:initGL(gl, glname, ...)
	SolarSystemApp.super.initGL(self, gl, glname, ...)
	for _,planet in ipairs(planets) do
		-- load texture
		local fn = 'textures/'..planet.name..'.png'
		if file(fn):exists() then
			pcall(function()
				planet.class.tex = Tex2D{
					filename=fn,
					minFilter=gl.GL_LINEAR,
					magFilter=gl.GL_LINEAR,
				}
			end)
		end
		
		planet.class.angle = quatd(0,0,0,1)			-- rotation ... only used for earth at the moment
		
		-- init vertex arrays
		local latdiv = math.floor((latMax-latMin)/latStep)
		local londiv = math.floor((lonMax-lonMin)/lonStep)
		planet.class.vertexCount = (latdiv + 1) * (londiv + 1)
		planet.class.elementCount = 4 * latdiv * londiv
		planet.class.elementArray = ffi.new('unsigned int[?]', planet.class.elementCount)
		planet.class.vertexArray = ffi.new('double[?]', 3 * planet.class.vertexCount)
		planet.class.texCoordArray = ffi.new('double[?]', 2 * planet.class.vertexCount)
		planet.class.tideArray = ffi.new('double[?]', planet.class.vertexCount)
		local vertexIndex = 0
		local elementIndex = 0
		for loni=0,londiv do
			local lon = lonMin + loni * lonStep 
			for lati=0,latdiv do
				local lat = latMin + lati * latStep
				
				-- vertex
				local pos = vec3d(planet:geodeticPosition(lat, lon, 0))
				for j=0,2 do
					planet.class.vertexArray[3*vertexIndex + j] = pos.s[j]
				end
				
				-- texcoord
				planet.class.texCoordArray[2*vertexIndex + 0] = lon / 360 + .5
				planet.class.texCoordArray[2*vertexIndex + 1] = -lat / 180 + .5
				
				vertexIndex = vertexIndex + 1
				
				if loni < londiv and lati < latdiv then
					for _,ofs in ipairs(quad) do
						planet.class.elementArray[elementIndex] = (lati + ofs[1]) + (latdiv + 1) * (loni + ofs[2])
						elementIndex = elementIndex + 1
					end
				end
			end
		end
		assert(vertexIndex == planet.class.vertexCount)
		assert(elementIndex == planet.class.elementCount)
		
	end

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_CULL_FACE)
	gl.glDepthFunc(gl.GL_LEQUAL)
	gl.glClearColor(0,0,0,0)
	
	hsvTex = HsvTex(256)
			
	local earth = planets[planets.indexes.earth]
	orbitPlanetIndex = earth.index
	orbitTargetDistance = 2 * earth.radius
	orbitDistance = orbitTargetDistance
end

local clickEarthSurfaceCallback

local leftShiftDown
local rightShiftDown 
function SolarSystemApp:event(event, ...)
	SolarSystemApp.super.event(self, event, ...)
	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not ig.igGetIO()[0].WantCaptureKeyboard

	if event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if canHandleMouse then
			if event.button.button == sdl.SDL_BUTTON_WHEELUP then
				orbitTargetDistance = orbitTargetDistance * orbitZoomFactor
			elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
				orbitTargetDistance = orbitTargetDistance / orbitZoomFactor
			elseif event.button.button == sdl.SDL_BUTTON_LEFT then
				if clickEarthSurfaceCallback then
					-- line intersection from view center to Earth sphere
					-- TODO just use from last mouseDir in :update() ?
					--[[
					ray sphere intersection
					ray point = {x(t) = a + b t}
					sphere point = {|x - c| = r}
					combine: 
					(a + b t - c) . (a + b t - c) = r^2
					let d = c - a
					(b t - d) . (b t - d) = r^2
					b.b t^2 - 2 d.b t + d.d - r^2 = 0
					t = (
						2 d.b +- sqrt(4 (d.b)^2 - 4 b.b (d.d - r^2))
					) / (2 b.b)
					t = ( d.b +- sqrt((d.b)^2 - b.b (d.d - r^2)) ) / (b.b)
					--]]
					local planet = planets[planets.indexes.earth]
					local mouseDir = mouseRay(mouse.pos)
					local delta = planet.pos - viewPos
					local db = delta:dot(mouseDir)
					local dd = delta:lenSq()
					local bb = mouseDir:dot(mouseDir)
					local discr = db * db - bb * (dd - planet.radius * planet.radius)
					if discr >= 0 then
						local pm = math.sqrt(discr) / bb
						local t0 = (db - pm) / bb
						local t1 = (db + pm) / bb
						local t = t0 > 0 and t0 or (t1 > 0 and t1 or nil)
						if t then
							local pt = viewPos + mouseDir * t
							local lat, lon, height = solarSystemBarycentricToPlanetGeodetic(planet, pt)	
							clickEarthSurfaceCallback(planet, lat, lon, height)
						end
					end
				end
			end
		end
	elseif event.type == sdl.SDL_KEYDOWN or event.type == sdl.SDL_KEYUP then
		if canHandleKeyboard then
			if event.key.keysym.sym == sdl.SDLK_LSHIFT then
				leftShiftDown = event.type == sdl.SDL_KEYDOWN
			elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
				rightShiftDown = event.type == sdl.SDL_KEYDOWN
			end
		end
	end
end
	
function SolarSystemApp:update(...)
	mouse:update()
	
	local mouseDir = mouseRay(mouse.pos)
	chooseNewPlanet(mouseDir, mouse.rightPress)
	
	if mouse.leftClick then
		if mouseOverEvent then
			targetJulianDate = mouseOverEvent.julianDate
		end
	elseif mouse.leftDragging then
		if leftShiftDown or rightShiftDown then
			orbitTargetDistance = orbitTargetDistance * math.exp(100 * orbitZoomFactor * mouse.deltaPos.y)
		else
			local magn = mouse.deltaPos:length() * 1000
			if magn > 0 then
				local normDelta = mouse.deltaPos / magn
				local r = quatd():fromAngleAxis(-normDelta.y, normDelta.x, 0, -magn)
				viewAngle = (viewAngle * r):normalize()
			end
		end
	end
	
	-- track ball orbit
	
	local earth = planets[planets.indexes.earth]
	
	local orbitCenter
	if orbitGeodeticLocation then
		orbitCenter = planetGeodeticToSolarSystemBarycentric(planets[orbitPlanetIndex], orbitGeodeticLocation.lat, orbitGeodeticLocation.lon, orbitGeodeticLocation.height)
	else
		orbitCenter = planets[orbitPlanetIndex].pos
	end
	viewPos = orbitCenter + viewAngle:zAxis() * orbitDistance
	
	do
		local logDist = math.log(orbitDistance)
		local logTarget = math.log(orbitTargetDistance)
		local coeff = .05
		local newLogDist = (1 - coeff) * logDist + coeff * logTarget
		orbitDistance = math.exp(newLogDist)
	end
	
	-- setup opengl
	
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)

	drawScene(1, mouseDir)

	
	-- iterate after render 

	
	-- push old state
	if julianDate > lastHistoryJulianDate + planetHistoryDelta then
		planetHistory[planetHistoryIndex] = planets
		planetHistoryIndex = planetHistoryIndex % planetHistoryMax + 1
	end

	local datetable

	--[[ realtime update
	local datetable = currentDate()
	julianDate = julian.fromCalendar(datetable)
	--]]
	
	local lastJulianDate = julianDate
	
	if targetJulianDate then
		local logDate = math.log(julianDate)
		local logTarget = math.log(targetJulianDate)
		local coeff = .5
		local deltaLog = logTarget - logDate
		if deltaLog < 0 then
			local newLogDate = logDate + coeff * deltaLog
			julianDate = math.exp(newLogDate)
		else
			julianDate = targetJulianDate
			targetJulianDate = nil
		end
		planets = Planets.fromEphemeris(julianDate)
	elseif integrateTimeStep then
		if julianDate > 625673.500000 and julianDate < 2816787.500000 then	-- within 3000 BC to 3000 AD
			julianDate = julianDate + integrateTimeStep
			planets = Planets.fromEphemeris(julianDate)
		elseif realtime then
			-- if we're integrating ...
			local integrate = require 'integrate.ivp'
			planets = integrate(julianDate, planets, integrateTimeStep, integrateFunction, integrationMethod, integrationArgs)
			julianDate = julianDate + integrateTimeStep
		end
	end
	datetable = julian.toCalendar(julianDate)
	
	if lastJulianDate ~= julianDate then
		local deltaJulianDate = julianDate - lastJulianDate
		local deltaAngle = quatd():fromAngleAxis(0,0,1, deltaJulianDate * 360)
		local deltaPos = viewPos - orbitCenter
		--deltaAngle[4] = -deltaAngle[4]
		deltaPos = deltaAngle:rotate(deltaPos)
		--deltaAngle[4] = -deltaAngle[4]
		viewPos = viewPos + orbitCenter
		viewAngle = deltaAngle * viewAngle
	end
	

	-- TODO I'm mixing up sidereal and other kind of day, hence all this nonsense
	-- set the angle for the time
	--2455389.287573 is aligned with the eclipse of 2455389.315718
	--so scale the angle by (2455389.287573 - 2455034.608623) / (2455389.315718 - 2455034.608623)
	do
		-- .331 offset for unscaled angle
		local angleOffset = 46.5 / 360
		local angleScale = 1--(2455389.287573 - 2455034.608623) / (2455389.315718 - 2455034.608623)
		planets[planets.indexes.earth].class.angle = quatd():fromAngleAxis(0,0,1, ((julianDate * angleScale + angleOffset) % 1) * 360)
	end

	dateText = ('%f / %d-%02d-%02d %02d:%02d:%02f'):format(julianDate, datetable.year, datetable.month, datetable.day, datetable.hour, datetable.min, datetable.sec)
	
	SolarSystemApp.super.update(self, ...)
end

function guiInputFloat(name, t, k, step, stepfast, format, flags)
	step = step or .1
	stepfast = stepfast or 1
	format = format or '%.3f'
	flags = flags or ig.ImGuiInputTextFlags_EnterReturnsTrue
	return ig.luatableInputFloat(name, t, k, step, stepfast, format, flags)
end

function SolarSystemApp:updateGUI()
	local orbitPlanet = planets[orbitPlanetIndex]
	if orbitPlanet then
		ig.igText(dateText)
		ig.igText(orbitPlanet.name)
	end

	if ig.igButton'<<' then
		if not integrateTimeStep or integrateTimeStep > 0 then
			integrateTimeStep = -1/(24*60*60)
		else
			integrateTimeStep = integrateTimeStep * 2
		end
	end
	ig.igSameLine()
	if ig.igButton'||' then	
		integrateTimeStep = nil
	end
	ig.igSameLine()
	if ig.igButton'R' then
		integrateTimeStep = nil
		julianDate = startDate
		planets = Planets.fromEphemeris(julianDate)
	end
	ig.igSameLine()
	if ig.igButton'>>' then
		if not integrateTimeStep or integrateTimeStep < 0 then
			integrateTimeStep = 1/(24*60*60)
		else
			integrateTimeStep = integrateTimeStep * 2
		end
	end

	ig.luatableCheckbox('Show Tide', _G, 'showTide')
	if showTide then
		ig.igText(tostring(tidalMin)..' m/s^2')
		ig.igText(tostring(tidalMax)..' m/s^2')
	end

--[[
	local bar = gui:widget{
		parent={gui.root},
		pos={3,75},
		size={3+colorBarWidth,5},
	}
	bar.backgroundTexture = hsvTex.id
	bar.backgroundColorValue = {1,1,1,1}
	bar.backgroundOffsetValue = {colorBarHSVRange,0}
	bar.backgroundScaleValue = {-colorBarHSVRange/colorBarWidth,1}
	bar.colorValue = {0,0,0,0}
--]]	

	if guiInputFloat('cycles', _G, 'numCycles') then
		historyCacheDate = nil
	end

	guiInputFloat('fov', _G, 'fov')

	if ig.luatableCheckbox('invert radial distance', _G, 'invertRadialDistance') then
		historyCacheDate = nil
	end
	if ig.luatableCheckbox('show orbit wrt planet', _G, 'showOrbitWrtPlanet') then
		historyCacheDate = nil
	end

	if ig.igCollapsingHeader'show trails' then
		if ig.igButton'all' then	-- or checkbox?
			for _,planetClass in ipairs(Planets.planetClasses) do
				showTrail[planetClass.name] = true
			end
		end
		for _,planetClass in ipairs(Planets.planetClasses) do
			if ig.luatableCheckbox(planetClass.name, showTrail, planetClass.name) then
				historyCacheDate = nil
			end
		end
	end

	if ig.luatableCheckbox('calc koe', _G, 'calcKOE') then
		if calcKOE then
			koeInfo = {}
			koeInfo.koe = {}
			for i,planet in ipairs(planets) do
				koeInfo.koe[i] = KOE.calcKOEFromPosVel(planet, planets, julianDate)
				--print('planet', i, planet.name, koeInfo.koe[i])
			end
			koeInfo.frames = {}
			for i,planet in ipairs(planets) do
				koeInfo.frames[i] = {}
				if koeInfo.koe[i] then
					local divs = 360
					local period = orbitPeriods[planet.name]
divs = divs * numCycles
period = period * numCycles
					for j=1,divs do
						local t = julianDate - ((j-1)/(divs-1)) * period
						local frameij = {}
						koeInfo.frames[i][j] = frameij
						-- TODO 
						-- this reads the .koe created by calcKOEFromPosVel
						-- and writes .pos .vel to it
						KOE.updatePosVel(
							frameij,
							koeInfo.koe[i],
							t,
							julianDate
						)
						
						local parentIndex = planet.parentIndex
						while parentIndex
						and koeInfo.frames[parentIndex]
						and koeInfo.frames[parentIndex][j]
						and koeInfo.frames[parentIndex][j].pos_koe 
						do
							frameij.pos_koe = frameij.pos_koe + koeInfo.frames[parentIndex][j].pos_koe
							frameij.vel_koe = frameij.vel_koe + koeInfo.frames[parentIndex][j].vel_koe
							parentIndex = planets[parentIndex].parentIndex
						end
--print(i, j, julianDate - t, (koeInfo.frames[i][j].pos_koe - historyCache[i][j]):length() / historyCache[i][j]:length())
					end
--print(#koeInfo.frames[i], 'vs', #historyCache[i])
					assert(#koeInfo.frames[i] == #historyCache[i])
				end
			end
		else
			koeInfo = nil
		end
	end

	if ig.igCollapsingHeader'draw arcs' then
		if ig.igButton'new arc' then
			local pt1, pt2
			clickEarthSurfaceCallback = function(planet, lat, lon, height)
				local pt = {lat, lon, height}
				if not pt1 then
					pt1 = pt
				elseif not pt2 then
					pt2 = pt
					allArcs:insert{
						type = "great circle",
						planet = planet.name,
						pt1,
						pt2,
					}
					saveArcs()
				end
			end
		end
		for i=1,#allArcs do
			if ig.igButton('del arc '..i) then
				allArcs:remove(i)
				saveArcs()
				break
			end
		end
	end

	if eventText then
		ig.igBeginTooltip()
		ig.igText(eventText)
		ig.igEndTooltip()
	end
end

solarSystemApp = SolarSystemApp()
return solarSystemApp:run()

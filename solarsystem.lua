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


local ffi = require 'ffi'
local ImGuiApp = require 'imguiapp'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local sdl = require 'ffi.sdl'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
local Planets = require 'planets'
local Mouse = require 'glapp.mouse'


local planets = Planets()
local earth = planets[planets.indexes.earth]


local julian = require 'julian'
require 'ext'
local Tex2D = require 'gl.tex2d'
local HsvTex = require 'gl.hsvtex'
require 'ffi.c.time'
--local shader = require 'gl.shader'


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
	if io.fileexists(fn) then
		d = assert(file[fn])
		d = assert(json.decode(d))
	else
		local http = require 'socket.http'
		local url = require 'socket.url'
		local dateStr = os.date('%Y-%m-%d %H:%M:%S', dateTime)
		d = http.request('http://www.astro-phys.com/api/de406/states?date='..url.encode(dateStr)..'&bodies='..table.concat(Planets.names, ','))
		d = assert(json.decode(d))
		d.calendarDate = dateTable
		d.linuxTime = dateTime
		file[fn] = json.encode(d, {indent=true})
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
--local earthquakeEntries = setmetatable(assert(assert(load('return '..assert(file['earthquakes.lua'])))()), table)
--local solarEclipseEntries = table((assert(json.decode(assert(file['solar-eclipses.json'])))))
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
local viewAngle = quatd()

local mouse = Mouse()
local leftButtonDown


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


local function planetGeodeticToSolarSystemBarycentric(planet, lat, lon, height)
	local x = vec3d(planet:geodeticPosition(lat, lon, height))		-- position relative to the planet center
	x = planetCartesianToSolarSystemBarycentric(planet, x)
	return x
end

local hsvTex
local colorBarWidth = 50	-- in menu units
local colorBarHSVRange = 2/3	-- how much of the rainbow to use
local drawWireframe = true

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
local tanFovX = 1	--.5
local tanFovY = 1	--.5

local dateText
local solarSystemApp

local function mouseRay(mousePos)
	local w, h = solarSystemApp:size()
	local aspectRatio = w / h
	-- ray intersect
	local fx = mousePos.x * 2 - 1
	local fy = mousePos.y * 2 - 1
	local v = vec3d(fx * aspectRatio * tanFovX, fy * tanFovY, -1)
	v = viewAngle:rotate(v):normalize()
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


function drawScene(viewScale, mouseDir)
	gl.glClear(gl.GL_DEPTH_BUFFER_BIT)

	local viewWidth, viewHeight = solarSystemApp:size()
	local aspectRatio = viewWidth / viewHeight
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glFrustum(-zNear * aspectRatio * tanFovX, zNear * aspectRatio * tanFovX, -zNear * tanFovY, zNear * tanFovY, zNear, zFar);

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glScaled(viewScale, viewScale, viewScale)
	local aa = viewAngle:toAngleAxis()
	gl.glRotated(-aa.w, aa.x, aa.y, aa.z)
	gl.glTranslated((-viewPos):unpack())

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
	
		--[[
		for i=1,#planets do
			gl.glColor3f(table.unpack(planet.color))
			gl.glBegin(gl.GL_LINE_STRIP)
			for j=-180,180,10 do
				local p = Planets.fromEphemeris(julianDate + j)
				gl.glVertex3d(p[i].pos:unpack())
			end
			gl.glEnd()
		end
		--]]

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
			
			if mouseOverEvent then
				eventText = mouseOverEvent.name..' '..os.date(nil, os.time(mouseOverEvent.date))
			else
				eventText = nil
			end
		end
		
		gl.glDisable(gl.GL_DEPTH_TEST)
		gl.glEnable(gl.GL_BLEND)
		
		-- line from moon through earth in direction of the sun (useful for eclipses)
		gl.glColor3f(1,.5,0)
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


local bool = ffi.new('bool[1]', false)
local function checkbox(name, object, field) 
	bool[0] = not not object[field]
	if ig.igCheckbox(name, bool) then
		object[field] = bool[0]
		return true
	end
end




local SolarSystemApp = class(ImGuiApp)

function SolarSystemApp:initGL(gl, glname, ...)
	SolarSystemApp.super.initGL(self, gl, glname, ...)
	local colors = {
		sun={1,1,0},
		mercury={.7,0,.2},
		venus={0,1,0},
		earth={0,0,1},
		moon={.6,.6,.6},
		mars={1,0,0},
		jupiter={1,.5,0},
		saturn={1,0,.5},
		uranus={0,1,1},
		neptune={1,0,1},
		pluto={0,.5,1},
	}	
	for _,planet in ipairs(planets) do
		planet.class.color = assert(colors[planet.name])
		-- load texture
		local fn = 'textures/'..planet.name..'.png'
		if io.fileexists(fn) then
			pcall(function()
				planet.class.tex = Tex2D{
					filename=fn,
					minFilter=gl.GL_LINEAR,
					magFilter=gl.GL_LINEAR,
				}
			end)
		end
		
		planet.class.angle = quatd()			-- rotation ... only used for earth at the moment
		
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
			local integrate = require 'integrate'
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

	checkbox('Show Tide', _G, 'showTide')
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


	if eventText then
		ig.igBeginTooltip()
		ig.igText(eventText)
		ig.igEndTooltip()
	end
end

solarSystemApp = SolarSystemApp()
solarSystemApp:run()

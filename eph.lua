-- this is in equatorial frame of reference, correct?
local table = require 'ext.table'
local path = require 'ext.path'
local timer = require 'ext.timer'
local assert = require 'ext.assert'
local fromlua = require 'ext.fromlua'
local math = require 'ext.math'
local ffi = require 'ffi'
local stdio = require 'ffi.req' 'c.stdio'

local hdr			-- header data
local filestr		-- ephemeral data file (big ~ 200mb)
local filebytes		-- ... cast as uint8_t*
local filedata		-- ... cast as double*

local denum			--= 406	--422
local eph = {}

function eph.hasInitialized()
	return hdr ~= nil
end

function eph.init(denum_, dir)
	denum = denum_
	--hdr = assert(json.decode(assert(path(dir..'/header.json'):read())))
	hdr = fromlua((assert(path(dir..'/header.luaconfig'):read())))
	local fn = dir..'/f64/de'..denum..'.f64.raw'
	timer('reading ephemeris data', function()
		filestr = assert(path(fn):read())							-- keep this around so it doesn't gc
	end)
	filebytes = ffi.cast('uint8_t*', ffi.cast('char*', filestr))
	filedata = ffi.cast('double*', filebytes)					-- cast as double
end



-- with a little help from
-- http://www.cv.nrao.edu/~rfisher/Ephemerides/ephem_descr.html
-- ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/jpl_eph.cc
-- http://www.astro-phys.com/js/astro/api.js

-- also in asc2fits
-- also saved in hdr / header.luaconfig
local objNames = table{
	'Mercury',
	'Venus',
	'EM_Bary',
	'Mars',
	'Jupiter',
	'Saturn',
	'Uranus',
	'Neptune',
	'Pluto',
	'GeoCMoon',
	'Sun',
	'Nutation',
	'Libration'
};
local objIndexForName = objNames:map(function(v,k) return k,v end)

local recordEpoch1 = 1e+10
local recordEpoch2 = -1e+10

--[[
looking at mercury coefficients for 2012-05-01 (julian 2456016.5)
buffer[0] is the epoch start, 625360.5, matches header.epoch1
buffer[1] is the first record(right word?) end, 625424.5 = 64 days later, i.e. end-time of this piecewise polynomial segment
buffer[85] is the beginning of mercury's data
so what goes on from buffer[2] to buffer[84] ?

Look at header.luaconfig.
numCoeffs = 728 <=> there's 728 variables per polynomial segment.

typedef union Segment {
	double coeffs[728];
	struct {
		double t0, t1;
		struct {
			constexpr int numCoeffs = 14;
			constexpr int numComponents = 3;
			constexpr int numSubIntervals = 4;
			double coeff[numComponents][numCoeffs][numSubIntervals];
		} Mercury;
		struct {
			... same for what's in heder.luaconfig under Venus
		} Venus;
		... same for the rest of the planets
		... if you sum it all up you get 726, which is 2 less than numCoeffs=728
	};
} Segment;
--]]
local buffer
local function getCoeffBuffer(timeOrigin, timeOffset)
	local numCoeffs = hdr.numCoeffs

--print('record epoch',recordEpoch1,recordEpoch2)
	local startOffset = (timeOrigin - recordEpoch1) + timeOffset
	if startOffset < 0 or (recordEpoch1 + startOffset) > recordEpoch2 then
		startOffset = (timeOrigin - hdr.epoch1) * 1.0038	-- hardcoded value?
		startOffset = startOffset + timeOffset
--print('startOffset',startOffset)
--DEBUG:assert.ge(startOffset, 0)
--DEBUG:assert.le(hdr.epoch1 + startOffset, hdr.epoch2)
		local recordNumber = math.floor(startOffset / hdr.interval)
--print('recordNumber',recordNumber)
		local dataOffset = recordNumber * numCoeffs
--print('dataOffset',dataOffset)
		-- seek to that location and open it up!
		buffer = filedata + dataOffset
		-- make sure we have the right record
		startOffset = (timeOrigin - buffer[0]) + timeOffset
--print('dataOffset',dataOffset,'time range',buffer[0],buffer[1],'startOffset',startOffset)
		local safety = 1000
		while startOffset < 0
		and dataOffset > numCoeffs
		and safety > 0
		do
			safety = safety - 1
			dataOffset = dataOffset - numCoeffs
			buffer = filedata + dataOffset
			startOffset = (timeOrigin - buffer[0]) + timeOffset
--print('dataOffset',dataOffset,'time range',buffer[0],buffer[1],'startOffset',startOffset)
		end
		-- maybe there was a gap between records
		local safety = 1000
		local endOffset = (timeOrigin - buffer[1]) + timeOffset
--print('dataOffset',dataOffset,'time range',buffer[0],buffer[1],'endOffset',endOffset)
		while endOffset > 0 and safety > 0 do
			safety = safety - 1
			dataOffset = dataOffset + numCoeffs
			buffer = filedata + dataOffset
			endOffset = (timeOrigin - buffer[1]) + timeOffset
--print('dataOffset',dataOffset,'time range',buffer[0],buffer[1],'endOffset',endOffset)
		end
		startOffset = (timeOrigin - buffer[0]) + timeOffset
		recordEpoch1 = buffer[0]
		recordEpoch2 = buffer[1]
--[[
print('dataOffset',dataOffset)
print('final record range',buffer[0],buffer[1],'startOffset',startOffset,'endOffset',endOffset)
print('record')
for i=0,numCoeffs-1 do
	io.write('\t'..buffer[i])
end
print()
--]]
--DEBUG:assert.ge(startOffset, 0)
--DEBUG:assert.le(recordEpoch1 + startOffset, recordEpoch2)
	end

	return buffer
end

local function getCoeffSubinterval(planetIndex, coeffBuffer, timeOrigin, timeOffset)
	local planetObject = hdr.objs[planetIndex]
	local subCoeff = coeffBuffer + planetObject.offset - 1		-- pointer addition
	local startOffset = (timeOrigin - coeffBuffer[0]) + timeOffset
--print('hdr.interval',hdr.interval)
--print('coeffBuffer[0]',coeffBuffer[0],'coeffBuffer[1]',coeffBuffer[1],'delta',coeffBuffer[1]-coeffBuffer[0])
	local timeFrac = startOffset / hdr.interval
--print('timeFrac',timeFrac)
	local dnumSub = planetObject.numSubIntervals
--print('dnumSub',dnumSub)
	local dt1 = math.trunc(timeFrac)
--print('dt1',dt1)	--always zero?
	local si = math.trunc(dnumSub * timeFrac - dt1)
--print('si',si)
	local numCoeffs = planetObject.numCoeffs
	subCoeff = subCoeff + numCoeffs * planetObject.numComponents * si	-- pointer addition
	local subintervalLength = hdr.interval / dnumSub
	local subintervalFrac = (startOffset - si * subintervalLength) / subintervalLength
--DEBUG:assert.ge(subintervalFrac,0)
--DEBUG:assert.le(subintervalFrac, 1)
--[[
print('subCoeff')
for i=0,numCoeffs*planetObject.numComponents-1 do
	io.write('\t'..subCoeff[i])
end
print()
--]]
	return subCoeff, numCoeffs, subintervalLength, subintervalFrac
end

local interp
do
	local pc = ffi.new('double[18]')
	local vc = ffi.new('double[18]')
	interp = function(coeff, numCoeffs, intervalLength, time)
	--[[
		local x = 2 * time - 1
		local x2 = 2 * x
		local pos, vel
		do
			local d = 0
			local dd = 0
			local tmp
			for i = numCoeffs-1,1,-1 do
				tmp = d
				d = x2 * d - dd + coeff[i]
				dd = tmp
			end
			pos = x * d - dd + coeff[0]
		end
		do
			local k = 2 / intervalLength
			local d = 0
			local dd = 0
			local tmp
			for i=numCoeffs-1,2,-1 do
				tmp = d
				d = x2 * d - dd + coeff[i] * k
				dd = tmp
			end
			vel = x2 * d - dd + coeff[1]
		end
		return pos, vel
	--]]
	-- [[
--DEBUG:assert.ge(time, 0)
--DEBUG:assert.le(time, 1)
		-- tc is the normalized chebyshev time (-1 <= tc <= 1)
		local tc = 2 * time - 1

		pc[0] = 1
		pc[1] = tc

		local twot = tc + tc
		for i=2,numCoeffs-1 do
			pc[i] = twot * pc[i-1] - pc[i-2]
		end

		local pos = 0
		for i=numCoeffs-1,0,-1 do
			pos = pos + pc[i] * coeff[i]
		end

		vc[0] = 0
		vc[1] = 1
		vc[2] = twot + twot
		for i=3,numCoeffs-1 do
			vc[i] = twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]
		end
		local vel = 0
		for i=numCoeffs-1,1,-1 do
			vel = vel + vc[i] * coeff[i]
		end
		vel = 2 * vel / intervalLength
		return pos, vel
	--]]
	end
end

-- TODO same for libration?
function eph.nutation(timeOrigin, timeOffset)
	assert.ne(hdr.objs[objIndexForName.Nutation].numCoeffs, 0, "this dataset has no nutation information")
	timeOffset = timeOffset or 0
	local coeffBuffer = getCoeffBuffer(timeOrigin, timeOffset)
	local coeff, numCoeffs, subintervalLength, subintervalFrac = getCoeffSubinterval(objIndexForName.Nutation, coeffBuffer, timeOrigin, timeOffset)
--print(coeff, numCoeffs, subintervalLength, subintervalFrac)
	--for i=0,1 do
	local posx, velx = interp(coeff, numCoeffs, subintervalLength, subintervalFrac)
	coeff = coeff + numCoeffs
	local posy, vely = interp(coeff, numCoeffs, subintervalLength, subintervalFrac)
	--coeff = coeff + numCoeffs
	--end
	return posx, posy, velx, vely
end

function eph.posVel(planetIndex, timeOrigin, timeOffset)
	timeOffset = timeOffset or 0
	--[[
	if not timeOffset then
		local ipart = math.floor(timeOrigin)
		timeOffset = timeOrigin - ipart
		timeOrigin = ipart
	end
	--]]
--DEBUG:assert.ge(planetIndex, 1)
--DEBUG:assert.le(planetIndex, #objNames)
	local coeffBuffer = getCoeffBuffer(timeOrigin, timeOffset)
	local coeff, numCoeffs, subintervalLength, subintervalFrac = getCoeffSubinterval(planetIndex, coeffBuffer, timeOrigin, timeOffset)
--print('subinterval length',subintervalLength,'frac',subintervalFrac)
	--for i=0,2 do
--[[
print('planet coeffs: ('..numCoeffs..') dim '..i)
for j=0,numCoeffs-1 do
	io.write('\t'..coeff[j])
end
print()
--]]
--print('dim',i,'pos',pos.s[i],'vel',vel.s[i])
	local posx, velx = interp(coeff, numCoeffs, subintervalLength, subintervalFrac)
	coeff = coeff + numCoeffs
	local posy, vely = interp(coeff, numCoeffs, subintervalLength, subintervalFrac)
	coeff = coeff + numCoeffs
	local posz, velz = interp(coeff, numCoeffs, subintervalLength, subintervalFrac)
	--coeff = coeff + numCoeffs
	--end
	return posx, posy, posz, velx, vely, velz
end


function eph.mercury(...) return eph.posVel(objIndexForName.Mercury, ...) end
function eph.venus(...) return eph.posVel(objIndexForName.Venus, ...) end
function eph.mars(...) return eph.posVel(objIndexForName.Mars, ...) end
function eph.jupiter(...) return eph.posVel(objIndexForName.Jupiter, ...) end
function eph.saturn(...) return eph.posVel(objIndexForName.Saturn, ...) end
function eph.uranus(...) return eph.posVel(objIndexForName.Uranus, ...) end
function eph.neptune(...) return eph.posVel(objIndexForName.Neptune, ...) end
function eph.pluto(...) return eph.posVel(objIndexForName.Pluto, ...) end
function eph.sun(...) return eph.posVel(objIndexForName.Sun, ...) end
-- what are nutation and libration used for?

function eph.EM_Bary(...) return eph.posVel(objIndexForName.EM_Bary, ...) end
--function eph.EM_Bary(...) return eph.posVel(objIndexForName.EM_Bary, ...) end

function eph.earth(...)
	local earthMoonPosX, earthMoonPosY, earthMoonPosZ, earthMoonVelX, earthMoonVelY, earthMoonVelZ = eph.posVel(objIndexForName.EM_Bary, ...)
	local geoMoonPosX, geoMoonPosY, geoMoonPosZ, geoMoonVelX, geoMoonVelY, geoMoonVelZ = eph.posVel(objIndexForName.GeoCMoon, ...)
	local scale = 1 / (1 + hdr.emrat)
	local earthPosX = earthMoonPosX - geoMoonPosX * scale
	local earthPosY = earthMoonPosY - geoMoonPosY * scale
	local earthPosZ = earthMoonPosZ - geoMoonPosZ * scale
	local earthVelX = earthMoonVelX - geoMoonVelX * scale
	local earthVelY = earthMoonVelY - geoMoonVelY * scale
	local earthVelZ = earthMoonVelZ - geoMoonVelZ * scale
	return earthPosX, earthPosY, earthPosZ, earthVelX, earthVelY, earthVelZ
end

function eph.moon(...)
	local earthMoonPosX, earthMoonPosY, earthMoonPosZ, earthMoonVelX, earthMoonVelY, earthMoonVelZ = eph.posVel(objIndexForName.EM_Bary, ...)
	local geoMoonPosX, geoMoonPosY, geoMoonPosZ, geoMoonVelX, geoMoonVelY, geoMoonVelZ = eph.posVel(objIndexForName.GeoCMoon, ...)
	local scale = 1 / (1 + hdr.emrat)
	local earthPosX = earthMoonPosX - geoMoonPosX * scale
	local earthPosY = earthMoonPosY - geoMoonPosY * scale
	local earthPosZ = earthMoonPosZ - geoMoonPosZ * scale
	local earthVelX = earthMoonVelX - geoMoonVelX * scale
	local earthVelY = earthMoonVelY - geoMoonVelY * scale
	local earthVelZ = earthMoonVelZ - geoMoonVelZ * scale
	return geoMoonPosX + earthPosX, geoMoonPosY + earthPosY, geoMoonPosZ + earthPosZ,
		geoMoonVelX + earthVelX, geoMoonVelY + earthVelY, geoMoonVelZ + earthVelZ
end

return eph

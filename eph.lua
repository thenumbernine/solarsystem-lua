-- this is in equatorial frame of reference, correct?
local table = require 'ext.table'
local path = require 'ext.path'
local assert = require 'ext.assert'
local fromlua = require 'ext.fromlua'
local math = require 'ext.math'
local ffi = require 'ffi'
local stdio = require 'ffi.req' 'c.stdio'

local hdr		-- header data
local stream	-- ephemeral data file (big ~ 600mb)
local buffer

local denum			--= 406	--422
local eph = {}

function eph.hasInitialized()
	return hdr ~= nil
end

function eph.init(denum_, dir)
	denum = denum_
	--hdr = assert(json.decode(assert(path(dir..'/header.json'):read())))
	hdr = fromlua(path(dir..'/header.luaconfig'):read())
	local fn = dir..'/f64/de'..denum..'.f64.raw'
	stream = stdio.fopen(fn, 'rb')
	assert.ne(stream, nil, "failed to open ephemeris data file "..fn)
	buffer = ffi.new('double[?]', hdr.numCoeffs)
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
buffer[0] and buffer[1] are the epoch start and end
buffer[85] is the beginning of mercury's data
so what goes on from buffer[2] to buffer[84] ?
--]]
local function getCoeffBuffer(timeOrigin, timeOffset)

--print('record epoch',recordEpoch1,recordEpoch2)
	local startOffset = (timeOrigin - recordEpoch1) + timeOffset
	if startOffset < 0 or (recordEpoch1 + startOffset) > recordEpoch2 then
		startOffset = (timeOrigin - hdr.epoch1) * 1.0038	-- hardcoded value?
		startOffset = startOffset + timeOffset
--print('startOffset',startOffset)
		assert(not (startOffset < 0 or (hdr.epoch1 + startOffset) > hdr.epoch2))
		local recordNumber = math.floor(startOffset / hdr.interval)
--print('recordNumber',recordNumber)
		local recordLength = hdr.numCoeffs * ffi.sizeof('double')
		local fileOffset = recordNumber * recordLength
--print('fileOffset',fileOffset)
		-- seek to that location and open it up!
		-- anyone know what a fits file header size is?  I might want to store this as raw doubles ...
		assert.eq(0, stdio.fseek(stream, fileOffset, stdio.SEEK_SET))
		stdio.fread(buffer, ffi.sizeof('double'), hdr.numCoeffs, stream)
		-- make sure we have the right record
		startOffset = (timeOrigin - buffer[0]) + timeOffset
--print('fileOffset',fileOffset,'time range',buffer[0],buffer[1],'startOffset',startOffset)
		local safety = 1000
		while startOffset < 0 and fileOffset > recordLength and safety > 0 do
			safety = safety - 1
			fileOffset = fileOffset - recordLength
			assert.eq(0, stdio.fseek(stream, fileOffset, stdio.SEEK_SET))
			stdio.fread(buffer, ffi.sizeof('double'), hdr.numCoeffs, stream)
			startOffset = (timeOrigin - buffer[0]) + timeOffset
--print('fileOffset',fileOffset,'time range',buffer[0],buffer[1],'startOffset',startOffset)
		end
		-- maybe there was a gap between records
		local safety = 1000
		local endOffset = (timeOrigin - buffer[1]) + timeOffset
--print('fileOffset',fileOffset,'time range',buffer[0],buffer[1],'endOffset',endOffset)
		while endOffset > 0 and safety > 0 do
			safety = safety - 1
			fileOffset = fileOffset + recordLength
			assert.eq(0, stdio.fseek(stream, fileOffset, stdio.SEEK_SET))
			stdio.fread(buffer, ffi.sizeof('double'), hdr.numCoeffs, stream)
			endOffset = (timeOrigin - buffer[1]) + timeOffset
--print('fileOffset',fileOffset,'time range',buffer[0],buffer[1],'endOffset',endOffset)
		end
		startOffset = (timeOrigin - buffer[0]) + timeOffset
		recordEpoch1 = buffer[0]
		recordEpoch2 = buffer[1]
--[[
print('fileOffset',fileOffset)
print('final record range',buffer[0],buffer[1],'startOffset',startOffset,'endOffset',endOffset)
print('record')
for i=0,hdr.numCoeffs-1 do
	io.write('\t'..buffer[i])
end
print()
--]]
		assert(not (startOffset < 0 or recordEpoch1 + startOffset > recordEpoch2))
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
	assert(not (subintervalFrac < 0 or subintervalFrac > 1))
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
		assert.ge(time, 0)
		assert.le(time, 1)
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
	assert.ge(planetIndex, 1)
	assert.le(planetIndex, #objNames)
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

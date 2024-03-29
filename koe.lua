-- now based on the first configuration, derive Vinti args
-- TODO to derive Vinti args from KOE, this is also in solarsystem/jpl-ssd-smallbody/output_points.template.lua
-- and in JS it's also in solarsystem/planet.js 'getKOEFromSourceData'
-- but to derive Vinti args from pos/vel/mass, use in solarsystem/planet.js 'calcKOEFromPosVel'
-- I should put it all in one place ...
local math = require 'ext.math'
local vec3d = require 'vec-ffi.vec3d'	-- is a fft metatype, not a class, so there's no "isa" for the record

local KOE = {}

KOE.gravitationalConstant = 6.6738480e-11	-- m^3 / (kg * s^2)

KOE.speedOfLight = 299792458	-- m/s

--provided z, calculate x and y such that x,y,z form a basis
-- this is 1 line and constexpr-optimized in my C++ diff geom library
-- but this is what it looks like in JS .. and my JS -> Lua port ... smh
function KOE.calcBasis(x,y,z)
	-- [1,0,0] cross z
	local cxx = 0
	local cxy = z.z
	local cxz = -z.y
	-- [0,1,0] cross z
	local cyx = -z.z
	local cyy = 0
	local cyz = z.x
	-- [0,0,1] cross z
	local czx = z.y
	local czy = -z.x
	local czz = 0
	local lx = math.sqrt(cxy * cxy + cxz * cxz)
	local ly = math.sqrt(cyx * cyx + cyz * cyz)
	local lz = math.sqrt(czx * czx + czy * czy)
	if lx < ly then
		if lx < lz then	--x is smallest
			x:set(cyx, cyy, cyz)
			y:set(czx, czy, czz)
		else		--z is smallest
			x:set(cxx, cxy, cxz)
			y:set(cyx, cyy, cyz)
		end
	else
		if ly < lz then	--y is smallest
			x:set(czx, czy, czz)
			y:set(cxx, cxy, cxz)
		else		--z is smallest
			x:set(cxx, cxy, cxz)
			y:set(cyx, cyy, cyz)
		end
	end
	x = x:normalize()
	y = y:normalize()
end


function KOE.calcOrbitBasis(planet, koe, planets)
	-- based on this position and velocity, find plane of orbit
	local parentBody
	local parentBodyIndex = planet.parentIndex
	if parentBodyIndex == nil then
		--console.log(planet.name+' has no orbit parent')
		koe.orbitAxis = vec3d(0,0,1)
		koe.orbitBasis = {vec3d(1,0,0),vec3d(0,1,0),vec3d(0,0,1)}
		return
	else
		parentBody = planets[parentBodyIndex]
	end

	--consider position relative to orbiting parent
	local pos = planet.pos - parentBody.pos

	--convert from m/day to m/s to coincide with the units of our gravitational constant
	local vel = (planet.vel - parentBody.vel) / (60 * 60 * 24)

	local angularMomentum = pos:cross(vel) --m^2/s
	local angularMomentumMagSq = angularMomentum:lenSq()		--m^4/s^2
	local angularMomentumMag = math.sqrt(angularMomentumMagSq)

	if angularMomentumMag < 1e-9 then
		koe.orbitAxis = vec3d(0,0,1)
	else
		koe.orbitAxis = angularMomentum / angularMomentumMag
	end

	local basisX = vec3d()
	local basisY = vec3d()
	local basisZ = koe.orbitAxis
	KOE.calcBasis(basisX, basisY, basisZ)
	--a[j][i] = a_ij, so our indexing is backwards, but our storage is column-major
	koe.orbitBasis = {basisX, basisY, basisZ}
end


--[[
args:
	pos = vec3d pos, in meters, 1-indexed
	vel = vec3d vel, in meters-per-day, 1-indexed
	mass = (optional) mass in kg
	parent = (optional) parent body, has the same info as 'args' does
	radius = (optional) for debugging only i think?
	name = just for debugging/errors
--]]
function KOE.calcKOEFromPosVel(planet, planets, initJulianDate)

	-- based on this position and velocity, find plane of orbit
	local parentBody
	local parentBodyIndex = planet.parentIndex
	if parentBodyIndex == nil then
		print('planet',planet.name,'has no parent')
		KOE.calcOrbitBasis(planet, {}, planets)
		return
	else
		parentBody = planets[parentBodyIndex]
	end

	--consider position relative to orbiting parent
	-- should I be doing the same thing with the velocity? probably...
	local pos = planet.pos - parentBody.pos

	--convert from m/day to m/s to coincide with the units of our gravitational constant
	local vel = (planet.vel - parentBody.vel) / (60 * 60 * 24)

	local posDotVel = pos:dot(vel)	--m^2/s

	local angularMomentum = pos:cross(vel) --m^2/s
	local angularMomentumMagSq = angularMomentum:lenSq()		--m^4/s^2
	local angularMomentumMag = math.sqrt(angularMomentumMagSq)
	
	--now decompose the relative position in the coordinates of the orbit basis
	--i've eliminated all but one of the rotation degrees of freedom ...

	--http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m
	--http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

	local velSq = vel:lenSq()		--(m/s)^2
	local distanceToParent = pos:length()	--m
	local gravitationalParameter = KOE.gravitationalConstant * ((planet.mass or 0) + parentBody.mass)	--m^3 / (kg s^2) * kg = m^3 / s^2
	-- https://en.wikipedia.org/wiki/Specific_orbital_energy
	-- ε = εk + εp
	-- ε = 1/2 v^2 - μ / r
	-- ε = -1/2 (μ^2 / h^2) (1 - e^2) ... where h = L / m
	-- ε = -μ / (2 a)
	local specificOrbitalEnergy = .5 * velSq - gravitationalParameter / distanceToParent		--m^2 / s^2 - m^3 / s^2 / m = m^2/s^2, supposed to be negative for elliptical orbits
	-- a = -μ / (2 ε)
	local semiMajorAxis = -.5 * gravitationalParameter / specificOrbitalEnergy		--m^3/s^2 / (m^2/s^2) = m
	-- l = h^2 / μ
	local semiLatusRectum = angularMomentumMagSq / gravitationalParameter			--m^4/s^2 / (m^3/s^2) = m
	-- https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes
	-- h = √(a μ (1 - e^2))
	-- e = √(1 - h^2 / (a μ))
	-- e = √(1 - l / a)
	local eccentricity = math.sqrt(1 - semiLatusRectum / semiMajorAxis)		--e, unitless (assuming elliptical orbit)
	--local semiMinorAxis = semiMajorAxis * math.sqrt(1 - eccentricity * eccentricity)
	local semiMinorAxis = math.sqrt(semiMajorAxis * semiLatusRectum)

	local orbitType = nil
	local parabolicEccentricityEpsilon = 1e-7
--		if math.abs(eccentricity - 1) < parabolicEccentricityEpsilon then
--print("danger! got a parabolic orbit for ",this," based on eccentricity epsilon",math.abs(eccentricity-1))
--			orbitType = 'parabolic'
--		end else 
	if eccentricity > 1 then
		orbitType = 'hyperbolic'
	else
		orbitType = 'elliptic'
	end

	-- https://en.wikipedia.org/wiki/Orbital_elements

	local cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity						--unitless
	local sinEccentricAnomaly = posDotVel / (eccentricity * math.sqrt(gravitationalParameter * semiMajorAxis))	--m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
	local eccentricAnomaly = math.atan2(sinEccentricAnomaly, cosEccentricAnomaly)	--E, in radians (unitless)

	local cosInclination = angularMomentum.z / angularMomentumMag	--unitless
	local sinInclination = math.sqrt(angularMomentum.x * angularMomentum.x + angularMomentum.y * angularMomentum.y) / angularMomentumMag	--unitless
	local inclination = math.atan2(sinInclination, cosInclination)	--i

	local sinPericenter = ((vel.x * angularMomentum.y - vel.y * angularMomentum.x) / gravitationalParameter - pos.z / distanceToParent) / (eccentricity * sinInclination)
	local cosPericenter = (angularMomentumMag * vel.z / gravitationalParameter - (angularMomentum.x * pos.y - angularMomentum.y * pos.x) / (angularMomentumMag * distanceToParent)) / (eccentricity * sinInclination)
	local argumentOfPeriapsis = math.atan2(sinPericenter, cosPericenter)	--omega

	local cosAscending = -angularMomentum.y / (angularMomentumMag * sinInclination)
	local sinAscending = angularMomentum.x / (angularMomentumMag * sinInclination)
	local longitudeOfAscendingNode = math.atan2(sinAscending, cosAscending)	--Omega

	local semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis	--m^3
	local orbitalPeriod
	if orbitType == 'elliptic' then
		orbitalPeriod = 2 * math.pi * math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24)	--julian day
		--override for known planets? don't forget to factor in the barycenter 
		if planet.orbitalPeriod then
			print('for planet',planet.name,'orbitalPeriod calculated',orbitalPeriod,'provided',planet.orbitalPeriod)
			orbitalPeriod = planet.orbitalPeriod
		end
	end

	local longitudeOfPeriapsis = longitudeOfAscendingNode + argumentOfPeriapsis	--omega-bar
	
	local meanAnomaly
	if orbitType == 'parabolic' then
		meanAnomaly = eccentricAnomaly + eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3
	elseif orbitType == 'hyperbolic' then
		meanAnomaly = eccentricity * math.sinh(eccentricAnomaly) - eccentricAnomaly
	elseif orbitType == 'elliptic' then
		meanAnomaly = eccentricAnomaly - eccentricity * sinEccentricAnomaly
	end

	local meanLongitude = meanAnomaly + longitudeOfPeriapsis
	
	--can I do this?
	local epoch = initJulianDate
	local meanAnomalyAtEpoch = meanAnomaly

	local timeOfPeriapsisCrossing = -meanAnomaly / math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (60*60*24)	--julian day

	local A = vec3d(
		semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
		semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
		semiMajorAxis * sinPericenter * sinInclination
	)
	local B = vec3d(
		-semiMinorAxis * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
		 semiMinorAxis * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
		 semiMinorAxis * cosPericenter * sinInclination
	)

	--[[
	to convert back:
	pos[i] = A[i] * (cosEccentricAnomaly - eccentricity) + B[i] * sinEccentricAnomaly
	rDot[i] = (-A[i] * sinEccentricAnomaly + B[i] * cosEccentricAnomaly) * math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (1 - eccentricity * cosEccentricAnomaly)
	--]]
	local checkPos = A * (cosEccentricAnomaly - eccentricity) + B * sinEccentricAnomaly

	--used by planets to offset reconstructed orbit coordinates to exact position of this
	local checkPosToPos = checkPos - pos
	local checkPosToPosDist = checkPosToPos:length()
	local checkPosError = checkPosToPosDist / distanceToParent
	if not math.isfinite(checkPosError) then	--NaN? debug!
		print(planet.name..' has no orbit info. mass: '..planet.mass..' radius: '..planet.radius)
	else
		if checkPosError > 1e-5 then	--only report significant error
			print(planet.name..' error of reconstructed position '.. checkPosError)
		end
	end

	local koe = {
		relVelSq = velSq,
		gravitationalParameter = gravitationalParameter,
		specificOrbitalEnergy = specificOrbitalEnergy,
		distanceToParent = distanceToParent,
		semiMajorAxis = semiMajorAxis,
		--semiMinorAxis = semiMinorAxis,
		semiLatusRectum = semiLatusRectum,
		inclination = inclination,
		argumentOfPeriapsis = argumentOfPeriapsis,
		longitudeOfAscendingNode = longitudeOfAscendingNode,
		timeOfPeriapsisCrossing = timeOfPeriapsisCrossing,
		meanAnomaly = meanAnomaly,
		orbitalPeriod = orbitalPeriod,
		--should I do this?
		epoch = epoch,
		meanAnomalyAtEpoch = meanAnomalyAtEpoch,
		--the following are used for the orbit path shader=
		orbitType = orbitType,
		eccentricity = eccentricity,
		eccentricAnomaly = eccentricAnomaly,
		A = A,
		B = B,
		--this is used for drawing but not in the shader
		fractionOffset = 0,
	}

	--not NaN, we successfully reconstructed the position
	if not math.isfinite(checkPosError) then
		print('check-position error was NaN for planet ',this)
		KOE.calcOrbitBasis(planet, koe, planets)
		return
	end
	
	--for elliptic orbits,
	-- we can accumulate & store the largest semi major axis of all children
	--but for parabolic/hyperbolic ...
	-- there is no bound ... so maybe those objects should be stored/rendered separately?
--[[	
	if parentBody and orbitType == 'elliptic' then
		parentBody.maxSemiMajorAxisOfEllipticalSatellites = math.max(semiMajorAxis, parentBody.maxSemiMajorAxisOfEllipticalSatellites or 0)
	end
--]]

	KOE.calcOrbitBasis(planet, koe, planets)

	return koe 
end

--[[ 
reads from koe:
	orbitType
	orbitalPeriod (elliptic)
	meanAnomalyAtEpoch (elliptic)
	epoch (elliptic)
	meanAnomaly (hyperbolic)
	timeOfPeriapsisCrossing (hyperbolic)
	eccentricity
	A
	B
	semiMajorAxis
	gravitationalParameter
writes to koe:
	meanAnomaly
	eccentricAnomaly
	fractionalOffset
writes to out:
	pos_koe
	vel_koe
--]]
function KOE.updatePosVel(out, koe, julianDate, initJulianDate)
	assert(koe, "no koe")
	local timeAdvanced = julianDate - initJulianDate
	local orbitType = koe.orbitType

	-- https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
	local meanAnomaly, meanMotion
	if orbitType == 'elliptic' then
		assert(koe, 'orbitalPeriod')
		assert(koe, 'meanAnomalyAtEpoch')	--not always there.  esp not in planets.
		meanMotion = 2 * math.pi / koe.orbitalPeriod
		meanAnomaly = koe.meanAnomalyAtEpoch + meanMotion * (julianDate - koe.epoch)
	elseif orbitType == 'hyperbolic' then
		assert(koe, 'timeOfPeriapsisCrossing')
		meanAnomaly = koe.meanAnomaly
		meanMotion = koe.meanAnomaly / (julianDate - koe.timeOfPeriapsisCrossing)
	elseif orbitType == 'parabolic' then
		error('parabolic orbit')
	else
		error'here'
	end

	local eccentricity = koe.eccentricity

	--solve eccentricAnomaly from meanAnomaly via Newton Rhapson
	--for elliptical orbits:
	--	f(E) = M - E + e sin E = 0
	--	f'(E) = -1 + e cos E
	--for parabolic oribts:
	--	f(E) = M - E - E^3 / 3
	--	f'(E) = -1 - E^2
	--for hyperbolic orbits:
	--	f(E) = M - e sinh(E) - E
	--	f'(E) = -1 - e cosh(E)
	local eccentricAnomaly = meanAnomaly
	for i = 1,10 do
		local func, deriv
		if orbitType == 'parabolic' then	--parabolic
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3
			deriv = -1 - eccentricAnomaly * eccentricAnomaly
		elseif orbitType == 'elliptic' then 	--elliptical
			func = meanAnomaly - eccentricAnomaly + eccentricity * math.sin(eccentricAnomaly)
			deriv = -1 + eccentricity * math.cos(eccentricAnomaly)	--has zeroes ...
		elseif orbitType == 'hyperbolic' then	--hyperbolic
			func = meanAnomaly + eccentricAnomaly - eccentricity  * math.sinh(eccentricAnomaly)
			deriv = 1 - eccentricity * math.cosh(eccentricAnomaly)
		else
			error'here'
		end

		local delta = func / deriv
		if math.abs(delta) < 1e-15 then break end
		eccentricAnomaly = eccentricAnomaly - delta
	end

	--TODO don't use meanMotion for hyperbolic orbits
	local fractionOffset = timeAdvanced * meanMotion / (2 * math.pi) 
	local theta = timeAdvanced * meanMotion
	--local eccentricAnomalyOrig = eccentricAnomaly - theta

	local argumentOfPeriapsis = koe.argumentOfPeriapsis

	--[[
	TODO HERE - add in schwarzschild relativistic precession
	https://en.wikipedia.org//wiki/Schwarzschild_geodesics#Precession_of_elliptical_orbits
	the amount that the major/minor axii swing around the orbital plane
	change-in-angle-per-revolution:
	δφ = 6πG(M + m) / (c^2 a (1 - e^2))
	M = mass1
	m = mass2
	G = gravitational constant
	c = speed-of-light
	a = semi-major axis
	e = eccentricity
	so if the time period for one revolution is P
	then δφ = 1/P ∂φ/∂t

	μ = G (M + m) = gravitational parameter
	l = a (1 - e^2) = semi-latus rectum
	δφ = 6 π μ / (c^2 a (1 - e^2))
	δφ = 6 π μ / (c^2 l)
	
	so in terms of KOE Vinti params, what is rotating?
	if the ellipse in the orbital plane is rotating then the vectors of the semi-major and minor axii will rotate.
	this means the 'argument of periapsis' will rotate
	and the 'longitude of ascending node' will rotate
	(both proportional to the precession angle rotation?)
	how about ... meanAnomaly, meanAnomalyAtEpoch, timeOfPeriapsisCrossing?
	also wouldn't the inclination change in some way?
	seems I should be solving some differentials instead of just +='ing the angles ...
	
	sources:
	https://astronomy.stackexchange.com/questions/632/determining-effect-of-small-variable-force-on-planetary-perihelion-precession
	--]]
--[[ attempt at GR precession
	--argumentOfPeriapsis = argumentOfPeriapsis + timeAdvanced / koe.orbitalPeriod * 3 / (1 - eccentricity * eccentricity)
	--argumentOfPeriapsis = argumentOfPeriapsis + 6 * math.pi * koe.gravitationalParameter * timeAdvanced * koe.orbitalPeriod / (KOE.speedOfLight * KOE.speedOfLight * koe.semiMajorAxis * (1 - eccentricity * eccentricity))
	-- https://astronomy.stackexchange.com/questions/632/determining-effect-of-small-variable-force-on-planetary-perihelion-precession 's update #4:
	argumentOfPeriapsis = argumentOfPeriapsis + fractionOffset * 3 * (2 * math.pi / koe.orbitalPeriod)^3 * (koe.semiMajorAxis / KOE.speedOfLight)^2 / (1 - eccentricity^2)
	-- lowest error but still no significant improvement:
	--argumentOfPeriapsis = argumentOfPeriapsis + 3 * koe.gravitationalParameter * theta / (KOE.speedOfLight * KOE.speedOfLight * koe.semiLatusRectum)
--]]

	--[[ Option #1
	-- Use the stored A & B
	-- Matches close enough to NASA Ephemeris
	local A = koe.A
	local B = koe.B
	--]]
	-- [[ Option #2
	-- Recalculate them from the (arg-inferred) KOE angles
	-- NOTICE cos/sinEccentricAnomaly is based on eccentricAnomaly at KOE parameter calcuation
	-- NOT the current eccentricAnomaly which I'm constantly updating by this method.
	--  I bet there's another name for this variable ... "path eccentric anomaly" after all?
	-- but I'm recalculating here cuz I'm thinking how to do relativistic precession, which means adjusting all these vars.
	--local cosEccentricAnomaly = math.cos(eccentricAnomalyOrig)
	--local sinEccentricAnomaly = math.sin(eccentricAnomalyOrig)
	local cosInclination = math.cos(koe.inclination)
	local sinInclination = math.sin(koe.inclination)
	local cosPericenter = math.cos(argumentOfPeriapsis)
	local sinPericenter = math.sin(argumentOfPeriapsis)
	local cosAscending = math.cos(koe.longitudeOfAscendingNode)
	local sinAscending = math.sin(koe.longitudeOfAscendingNode)
	local semiMinorAxis = math.sqrt(koe.semiMajorAxis * koe.semiLatusRectum)
	--local semiMinorAxis = koe.semiMajorAxis * math.sqrt(1 - koe.eccentricity * koe.eccentricity)
	local A = vec3d(
		koe.semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
		koe.semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
		koe.semiMajorAxis * 											   sinPericenter * sinInclination
	)
	local B = vec3d(
		-semiMinorAxis * ( cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
		 semiMinorAxis * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
		 semiMinorAxis * 												 cosPericenter * sinInclination
	)
	--]]


	--matches above
	local dt_dE
	local semiMajorAxisCubed = koe.semiMajorAxis * koe.semiMajorAxis * koe.semiMajorAxis
	if orbitType == 'parabolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / koe.gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly)
	elseif orbitType == 'elliptic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / koe.gravitationalParameter) * (1 - koe.eccentricity * math.cos(eccentricAnomaly))
	elseif orbitType == 'hyperbolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / koe.gravitationalParameter) * (koe.eccentricity * math.cosh(eccentricAnomaly) - 1)
	end
	local dE_dt = 1/dt_dE
	local coeffA, coeffB
	local coeffDerivA, coeffDerivB
	if orbitType == 'parabolic' then
		--...?
	elseif orbitType == 'elliptic' then 
		coeffA = math.cos(eccentricAnomaly) - koe.eccentricity
		coeffB = math.sin(eccentricAnomaly)
		coeffDerivA = -math.sin(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cos(eccentricAnomaly) * dE_dt
	elseif orbitType == 'hyperbolic' then
		coeffA = koe.eccentricity - math.cosh(eccentricAnomaly)
		coeffB = math.sinh(eccentricAnomaly)
		coeffDerivA = -math.sinh(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cosh(eccentricAnomaly) * dE_dt
	end
	local pos = A * coeffA + B * coeffB
	local vel = A * coeffDerivA + B * coeffDerivB	--m/day

	-- don't add parent's position.  leave that for the caller to do.
	out.pos_koe = pos
	out.vel_koe = vel

	koe.meanAnomaly = meanAnomaly
	koe.eccentricAnomaly = eccentricAnomaly
	koe.fractionOffset = fractionOffset
end

return KOE

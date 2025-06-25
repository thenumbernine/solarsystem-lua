#!/usr/bin/env luajit
-- Visualize all JPL SSD small-body orbits simultaneously, allows for timestepping and Ephemeris 406 data of planets.
-- This is almost GLES3 compat except that it uses a dvec3 / GL_DOUBLE attribute in one place.

local cmdline = require 'ext.cmdline'(...)
local sdl = require 'sdl'
local gl = require 'gl.setup'(cmdline.gl)
local ffi = require 'ffi'
local vec3d = require 'vec-ffi.vec3d'
local table = require 'ext.table'
local path = require 'ext.path'
local glreport = require 'gl.report'
local GLSceneObject = require 'gl.sceneobject'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLAttribute = require 'gl.attribute'
local GLProgram = require 'gl.program'
local HSVTex2D = require 'gl.hsvtex2d'
local clnumber = require 'cl.obj.number'
local CLEnv = require 'cl.obj.env'
local template = require 'template'
local ig = require 'imgui'			-- must go after require 'imugiapp' on windows
local Julian = require 'julian'
local Planets = require 'planets'

-- [=[ this matches parse.lua
-- These *need to* match the fields in body_t_desc.lua up until bodyType etc, which are down below, and should match also.
local numberFields = table{
	'epoch',
	'perihelionDistance',		--comets
	'semiMajorAxis',		--asteroids
	'eccentricity',
	'inclination',
	'argumentOfPeriapsis',
	'longitudeOfAscendingNode',
	'meanAnomalyAtEpoch',		--asteroids
	'absoluteMagnitude',		--asteroids
	'magnitudeSlopeParameter',		--asteroids
	'timeOfPerihelionPassage',		--comets
}

local real = 'double'
-- also doubled in clenv
ffi.cdef('typedef '..real..' real;')
ffi.cdef('typedef '..real..'4 real4;')

local bodyTypeCode = template([[

enum {
	ORBIT_ELLIPTIC = 0,
	ORBIT_HYPERBOLIC = 1,
	ORBIT_PARABOLIC = 2,
};

typedef struct body_t {
<? for _,field in ipairs(numberFields) do
?>	real <?=field?>;
<? end
?>
	int bodyType;	// 0=comet, 1=numbered, 2=unnum
	int horizonID;	//for numbered asteroids only
	char name[43+1];
	char orbitSolutionReference[12+1];

	//computed parameters:
	long index;
	real pos[3], vel[3], A[3], B[3];
	real eccentricAnomaly;
	real timeOfPeriapsisCrossing;
	real meanAnomaly;
	int orbitType;	// 0=elliptic, 1=hyperbolic, 2=parabolic ... this can be derived from eccentricity
	real orbitalPeriod;
} body_t;
]], {
	real = real,
	numberFields = numberFields,
})
ffi.cdef(bodyTypeCode)
--]=] end of the code that matches solarsystem/jpl-ssd-smallbody/parse.lua

local App = require 'imgui.appwithorbit'()
App.title = 'JPL SSD Smallbody Visualizer'
App.viewDist = 2


--local calcNearLineMethod = 'shader'
local calcNearLineMethod = 'fillbuffer'

--local updateBodyMethod = 'cpu'
local updateBodyMethod = 'gpu'

-- scale for the smallbody data, which is in meters
local AU_in_m = 149597870700
local scale = 1 / AU_in_m
--local earthMoonDist_in_m = 380000000
--local scale = 1 / earthMoonDist_in_m
--earthMoonDist_in_AU = 0.002540143106462

local sunMass_kg = 1.9891e+30
local gravitationalConstant = 6.6738480e-11							-- m^3 / (kg * s^2)
local gravitationalParameter = gravitationalConstant * sunMass_kg	--assuming the comet mass is negligible, since the comet mass is not provided


local hsvTex

local distThreshold = .2	-- in AU
-- TODO this is supposed to be the spaceweather.com PHA threshold but I am not seeing any within this distance
--local distThreshold = .00256	-- AU = 1 LD (Lunar Distance)

local dateFormat = '%4d/%02d/%02d %02d:%02d:%02d'

function App:initGL(...)
	App.super.initGL(self, ...)

	self.view.znear = .001
	self.view.zfar = 100

	hsvTex = HSVTex2D(256)

	-- do a one-time snapshot, so no ffwd/rewind just yet
	local t = os.date('!*t')
	self.dateStr = dateFormat:format(t.year, t.month, t.day, t.hour, t.min, t.sec)
	self.julianDay = Julian.fromCalendar(t)
	--self.resetDay = self.julianDay
	self.resetDay = dofile'smallbody_t_desc.lua'.julianDay
	self.planets = Planets.fromEphemeris(self.julianDay, 406, 'eph/406')
assert(glreport'here')

	local earth = self.planets[self.planets.indexes.earth]
	self.view.orbit:set((earth.pos * scale):unpack())
	self.view.pos = self.view.angle:zAxis() * self.viewDist + self.view.orbit


	-- this is produced in my webgl solarsystem project,
	-- generation is done inside jpl-ssd-smallbody/output-points.template.lua
	-- the resulting file is jpl-ssd-smallbody/alldata.raw
	local data = assert(path'smallbodies.raw':read(), "failed to load smallbodies.raw")

	local bodies = ffi.cast('body_t*', data)
	self.numBodies = #data / ffi.sizeof'body_t'
--self.numBodies = math.min(self.numBodies, 1000)
	if self.numBodies ~= math.floor(self.numBodies) then
		-- structures don't align with integers, the padding on the read probably differed from the write
		error(tolua{
			['self.numBodies'] = self.numBodies,
			['math.floor(self.numBodies)'] = math.floor(self.numBodies),
		})
	end
	print('numBodies', self.numBodies)
	self.bodies = bodies

--[=[ trim the number of bodies
	local newBodies = table()
	for i=0,self.numBodies-1 do
		local body = self.bodies[i]
		local dx = (body.pos[0] - earth.pos.x) * scale	-- scale is 1/(AU in m)
		local dy = (body.pos[1] - earth.pos.y) * scale
		local dz = (body.pos[2] - earth.pos.z) * scale
		local lenSq = dx*dx + dy*dy + dz*dz
		if lenSq < 1 then
			newBodies:insert(ffi.new('body_t', self.bodies[i]))
		end
	end
print('resizing from '..self.numBodies..' to '..#newBodies)
	self.numBodies = #newBodies
	self.bodies = ffi.new('body_t[?]', self.numBodies)
	for i=0,self.numBodies-1 do
		self.bodies[i] = newBodies[i+1]
	end
--]=]

	self.env = CLEnv{size=self.numBodies, useGLSharing=false}	--, real=real}
	self.bodiesCLBuf = self.env:buffer{name='bodies', type='body_t', data=self.bodies}

	self.env.code = table{
		self.env.code,
		bodyTypeCode,
	}:concat'\n'

	self.updateCLProgram = self.env:program{
		code = template([[
kernel void update(
	global body_t* bodies,
	real julianDay,
	real resetDay,
	real4 sunPos
) {
	initKernel();

	global body_t* ke = bodies + index;

	real timeAdvanced = julianDay - resetDay;
	int orbitType = ke->orbitType;

	//https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

	real meanAnomaly, meanMotion;
	if (orbitType == ORBIT_ELLIPTIC) {
		meanMotion = 2. * M_PI / ke->orbitalPeriod;
		meanAnomaly = ke->meanAnomalyAtEpoch + meanMotion * (julianDay - ke->epoch);
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		meanAnomaly = ke->meanAnomaly;
		meanMotion = ke->meanAnomaly / (julianDay - ke->timeOfPeriapsisCrossing);
	} else if (orbitType == ORBIT_PARABOLIC) {
		//error'got a parabolic orbit'
	} else {
		//error'here'
	}

	real eccentricity = ke->eccentricity;

	//solve eccentricAnomaly from meanAnomaly via Newton Rhapson
	//for elliptical orbits:
	//	f(E) = M - E + e sin E = 0
	//	f'(E) = -1 + e cos E
	//for parabolic oribts:
	//	f(E) = M - E - E^3 / 3
	//	f'(E) = -1 - E^2
	//for hyperbolic orbits:
	//	f(E) = M - e sinh(E) - E
	//	f'(E) = -1 - e cosh(E)
	real eccentricAnomaly = meanAnomaly;
	for (int i = 0; i < 10; ++i) {
		real func, deriv;
		if (orbitType == ORBIT_PARABOLIC) {
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3.;
			deriv = -1. - eccentricAnomaly * eccentricAnomaly;
		} else if (orbitType == ORBIT_ELLIPTIC) {
			func = meanAnomaly - eccentricAnomaly + eccentricity * sin(eccentricAnomaly);
			deriv = -1. + eccentricity * cos(eccentricAnomaly);	//has zeroes ...
		} else if (orbitType == ORBIT_HYPERBOLIC) {
			func = meanAnomaly + eccentricAnomaly - eccentricity  * sinh(eccentricAnomaly);
			deriv = 1. - eccentricity * cosh(eccentricAnomaly);
		} else {
			//error'here'
		}

		real delta = func / deriv;
		if (fabs(delta) < 1e-15) break;
		eccentricAnomaly -= delta;
	}

	//TODO don't use meanMotion for hyperbolic orbits
	//real fractionOffset = timeAdvanced * meanMotion / (2 * M_PI);
	real theta = timeAdvanced * meanMotion;
	real pathEccentricAnomaly = eccentricAnomaly + theta;
	global real* A = ke->A;
	global real* B = ke->B;

	//matches above
	real dt_dE;
	real semiMajorAxisCubed = ke->semiMajorAxis * ke->semiMajorAxis * ke->semiMajorAxis;
	const real gravitationalParameter = <?=clnumber(gravitationalParameter)?>;
	if (orbitType == ORBIT_PARABOLIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (1. + pathEccentricAnomaly * pathEccentricAnomaly);
	} else if (orbitType == ORBIT_ELLIPTIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (1. - ke->eccentricity * cos(pathEccentricAnomaly));
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (ke->eccentricity * cosh(pathEccentricAnomaly) - 1.);
	}
	real dE_dt = 1. / dt_dE;
	real coeffA, coeffB;
	//real coeffDerivA, coeffDerivB;
	if (orbitType == ORBIT_PARABOLIC) {
		//...?
	} else if (orbitType == ORBIT_ELLIPTIC) {
		coeffA = cos(pathEccentricAnomaly) - ke->eccentricity;
		coeffB = sin(pathEccentricAnomaly);
		//coeffDerivA = -sin(pathEccentricAnomaly) * dE_dt;
		//coeffDerivB = cos(pathEccentricAnomaly) * dE_dt;
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		coeffA = ke->eccentricity - cosh(pathEccentricAnomaly);
		coeffB = sinh(pathEccentricAnomaly);
		//coeffDerivA = -sinh(pathEccentricAnomaly) * dE_dt;
		//coeffDerivB = cosh(pathEccentricAnomaly) * dE_dt;
	}
	real posX = A[0] * coeffA + B[0] * coeffB;
	real posY = A[1] * coeffA + B[1] * coeffB;
	real posZ = A[2] * coeffA + B[2] * coeffB;
	//real velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
	//real velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
	//real velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;

	ke->pos[0] = posX + sunPos.x;
	ke->pos[1] = posY + sunPos.y;
	ke->pos[2] = posZ + sunPos.z;
	// I'm not using this
	//ke->vel[0] = velX + parent.vel.s[0];
	//ke->vel[1] = velY + parent.vel.s[1];
	//ke->vel[2] = velZ + parent.vel.s[2];

//you only need this for calcNearLineMethod == 'shader', but not for 'fillbuffer'
//TODO CL/GL interop if you do use this
#if 0
	// now update buffers
	self.bodyToEarthArray[0+2*index].x = ke->pos[0];
	self.bodyToEarthArray[0+2*index].y = ke->pos[1];
	self.bodyToEarthArray[0+2*index].z = ke->pos[2];
	self.bodyToEarthArray[1+2*index].x = ke->pos[0];
	self.bodyToEarthArray[1+2*index].y = ke->pos[1];
	self.bodyToEarthArray[1+2*index].z = ke->pos[2];
#endif

	//ke == this.keplerianOrbitalElements;
	//ke->meanAnomaly = meanAnomaly;	// this doesn't change, does it?
	ke->eccentricAnomaly = eccentricAnomaly;	// this does change ... eccentricAnomaly and pos
	//ke->fractionOffset = fractionOffset;
}
]],		{
			gravitationalParameter = gravitationalParameter,
			clnumber = clnumber,
		}),
	}
	self.updateCLProgram:compile()
	self.updateCLKernel = self.updateCLProgram:kernel'update'

	-- TODO keep only those near the point

	self.bodyBuf = GLArrayBuffer{
		size = self.numBodies * ffi.sizeof'body_t',
		data = self.bodies,
		usage = gl.GL_STATIC_DRAW,
	}:unbind()
assert(glreport'here')

	self.bodyPosAttr = GLAttribute{
		buffer = self.bodyBuf,
		dim = 3,
		type = real == 'double' and gl.GL_DOUBLE or gl.GL_FLOAT,
		stride = ffi.sizeof'body_t',
		offset = ffi.offsetof('body_t', 'pos'),
	}
assert(glreport'here')

	-- until I get a driver with geometry shader support...
	self.bodyToEarthArray = ffi.new('vec3d_t[?]', self.numBodies*2)
	for i=0,self.numBodies-1 do
		for j=0,2 do
			self.bodyToEarthArray[0+2*i].s[j] = self.bodies[i].pos[j]
			self.bodyToEarthArray[1+2*i].s[j] = self.bodies[i].pos[j]
		end
	end
assert(glreport'here')
	self.numBodyToEarthLines = 0

	self.bodyToEarthBuf = GLArrayBuffer{
		size = ffi.sizeof(self.bodyToEarthArray),
		data = self.bodyToEarthArray,
		usage = gl.GL_STATIC_DRAW,
	}:unbind()
assert(glreport'here')

	self.bodyToEarthPosAttr = GLAttribute{
		buffer = self.bodyToEarthBuf,
		dim = 3,
		type = gl.GL_DOUBLE,
	}
assert(glreport'here')



	-- [=[ in absence of geometry buffers, I was trying to double up the body pos buffer, but that seems to be running slow
	if calcNearLineMethod == 'shader' then
		self.drawLineToEarthShader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = [[
layout(location=0) in vec4 bodyPos;
uniform vec3 earthPos;
out float lum;
uniform mat4 mvProjMat;
void main() {
	vec4 pos = bodyPos;
	float dist = length(earthPos - pos.xyz);
	lum = 1. - ]]..clnumber(scale)..[[ * dist;
	if (gl_VertexID % 2 == 0) {
		pos = vec4(earthPos, 1.);
	}
	gl_Position = mvProjMat * pos;
}
]],
			fragmentCode = [[
in float lum;
layout(location=0) out vec4 fragColor;
void main() {
	fragColor = vec4(lum, 0., 0., 1.);
}
]],
		}:useNone()
assert(glreport'here')
	elseif calcNearLineMethod == 'fillbuffer' then
		self.drawLineToEarthShader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = template([[
layout(location=0) in vec4 bodyPos;
uniform vec3 earthPos;
out vec3 color;
uniform sampler2D hsvTex;
uniform mat4 mvProjMat;
void main() {
	vec4 pos = bodyPos;
	float dist = length(earthPos - pos.xyz);
	float tc = .5 * dist / <?=clnumber(distThreshold / scale)?>;
	color = texture(hsvTex, vec2(tc, .5)).rgb;
	if (gl_VertexID % 2 == 0) {
		pos = vec4(earthPos, 1.);
	}
	gl_Position = mvProjMat * pos;
}
]], 		{
				distThreshold = distThreshold,
				scale = scale,
				clnumber = clnumber,
			}),
			fragmentCode = [[
in vec3 color;
layout(location=0) out vec4 fragColor;
void main() {
	fragColor = vec4(color, 1.);
	//fragColor = vec4(1., 0., 0., 1.);
}
]],

			uniforms = {
				hsvTex = 0,
			}
		}:useNone()
assert(glreport'here')

		-- trying another trick: just fill the buffer ... after loading planets
		-- if we aren't using one buffer for the lines then update it here after planets are loaded
		self:updateBodyToEarthLineBuf()
assert(glreport'here')
	end


	if calcNearLineMethod == 'shader' then
		self.drawLineToEarthShader:setAttrs{
			bodyPos = self.bodyPosAttr,
		}
assert(glreport'here')
	elseif calcNearLineMethod == 'fillbuffer' then
		self.drawLineToEarthShader:setAttrs{
			bodyPos = self.bodyToEarthPosAttr,
		}
assert(glreport'here')
	end


--[[
local rmin, rmax
	for i=0,self.numBodies-1 do
local x = self.bodies[i].pos[0]
local y = self.bodies[i].pos[1]
local z = self.bodies[i].pos[2]
local r = math.sqrt(x*x + y*y + z*z)
if not rmin or r < rmin then rmin = r end
if not rmax or r > rmax then rmax = r end
	end
print('r range', rmin, rmax)
--]]


	self.lineBasisObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
layout(location=0) in vec3 vertex;
in vec3 color;
out vec3 colorv;
uniform mat4 mvProjMat;
void main() {
	colorv = color;
	gl_Position = mvProjMat * vec4(vertex, 1.);
}
]],
			fragmentCode = [[
in vec3 colorv;
layout(location=0) out vec4 fragColor;
void main() {
	fragColor = vec4(colorv, 1.);
}
]],
		},
		geometry = {
			mode = gl.GL_LINES,
		},
		vertexes = {
			dim = 3,
			data = {0,0,0, 1,0,0, 0,0,0, 0,1,0, 0,0,0, 0,0,1},
		},
		attrs = {
			color = {
				buffer = {
					dim = 3,
					data = {1,0,0, 1,0,0, 0,1,0, 0,1,0, 0,0,1, 0,0,1},
				},
			},
		},
	}

	self.drawBodiesShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = template([[
layout(location=0) in <?=real == 'double' and 'dvec3' or 'vec3'?> bodyPos;
uniform mat4 mvProjMat;
void main() {
	gl_Position = mvProjMat * vec4(bodyPos, 1.);
}
]], 	{
			real = real,
		}),
		fragmentCode = [[
layout(location=0) out vec4 fragColor;
uniform vec3 color;
void main() {
	fragColor = vec4(color, 1.);
}
]],
	}:useNone()

	self.drawBodiesSceneObj = GLSceneObject{
		program = self.drawBodiesShader,
		geometry = {
			mode = gl.GL_POINTS,
			count = self.numBodies,
		},
		attrs = {
			bodyPos = {
				buffer = self.bodyBuf,
				dim = 3,
				type = real == 'double' and gl.GL_DOUBLE or gl.GL_FLOAT,
				stride = ffi.sizeof'body_t',
				offset = ffi.offsetof('body_t', 'pos'),
			},
		},
	}

	self.drawPlanetsObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
layout(location=0) in vec3 vertex;
in vec3 color;
out vec3 colorv;
uniform mat4 mvProjMat;
void main() {
	colorv = color;
	gl_Position = mvProjMat * vec4(vertex, 1.);
}
]],
			fragmentCode = [[
in vec3 colorv;
layout(location=0) out vec4 fragColor;
void main() {
	fragColor = vec4(colorv, 1.);
}
]],
		},
		geometry = {
			mode = gl.GL_POINTS,
		},
		vertexes = {
			dim = 3,
			useVec = true,
		},
		attrs = {
			color = {
				buffer = {
					dim = 3,
					useVec = true,
				},
			},
		},
	}

	gl.glEnable(gl.GL_POINT_SMOOTH)
	gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE)
	gl.glDisable(gl.GL_DEPTH_TEST)
	assert(glreport'here')
end

App.timeStep = .1

function App:update()
	self:draw()

	if self.running then
		if self.running == true then
			self.julianDay = self.julianDay + self.timeStep
		end
		local t = Julian.toCalendar(self.julianDay)
		self.dateStr = dateFormat:format(t.year, t.month, t.day, t.hour, t.min, t.sec)
		self.planets = Planets.fromEphemeris(self.julianDay)
		self:recalculateSmallBodies()
		if self.running == 'update' then
			self.running = nil
		end

		local earth = self.planets[self.planets.indexes.earth]
		self.viewDist = (self.view.pos - self.view.orbit):length()
		self.view.orbit:set((earth.pos * scale):unpack())
		self.view.pos = self.view.angle:zAxis() * self.viewDist + self.view.orbit
	end
end


App.alpha = .05
App.useBlend = true

function App:draw()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glViewport(0, 0, self.width, self.height)

	self.view:setup(self.width / self.height)
	assert(glreport'here')

	if self.useBlend then
		gl.glEnable(gl.GL_BLEND)
	end

	self.lineBasisObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	self.lineBasisObj:draw()

	local pushMat = self.view.mvMat:clone()
	self.view.mvMat:applyScale(scale, scale, scale)
	self.view.mvProjMat:mul4x4(self.view.projMat, self.view.mvMat)

	-- TODO use the original binary blob, and just pass it as a gl buffer, then use pos as a strided vertex array
	gl.glPointSize(3)
	
	if self.useBlend then
		self.drawBodiesSceneObj.uniforms.color = {self.alpha, self.alpha, self.alpha}
	else
		self.drawBodiesSceneObj.uniforms.color = {1, 1, 1}
	end
	self.drawBodiesSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	self.drawBodiesSceneObj:draw()
	assert(glreport'here')

	local earth = self.planets[self.planets.indexes.earth]
	if calcNearLineMethod == 'shader' then
		-- draw all lines and use a shader to lighten/darken them
		-- draw line to earth
		-- a geometry shader would be nice ... so I can generate one point at the earth's position
		-- until then, i'll make a second buffer
		self.drawLineToEarthShader:use()

		gl.glUniformMatrix4fv(self.drawLineToEarthShader.uniforms.mvProjMat.loc, 1, false, self.view.mvProjMat.ptr)

		if self.drawLineToEarthShader.uniforms.earthPos then
			--gl.glUniform3dv(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos.s)
			gl.glUniform3f(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos:unpack())
		end

		gl.glEnableVertexAttribArray(0)
		gl.glDrawArrays(gl.GL_LINES, 0, 2 * self.numBodies)
		gl.glDisableVertexAttribArray(0)

		self.drawLineToEarthShader:useNone()
		assert(glreport'here')
	elseif calcNearLineMethod == 'fillbuffer' then
		-- draw only those that we have filled in advance
		self.drawLineToEarthShader:use()

		gl.glUniformMatrix4fv(self.drawLineToEarthShader.uniforms.mvProjMat.loc, 1, false, self.view.mvProjMat.ptr)

		if self.drawLineToEarthShader.uniforms.earthPos then
			--gl.glUniform3dv(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos.s)
			gl.glUniform3f(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos:unpack())
		end
		hsvTex:bind()

		gl.glEnableVertexAttribArray(0)
		gl.glDrawArrays(gl.GL_LINES, 0, 2 * self.numBodyToEarthLines)
		gl.glDisableVertexAttribArray(0)

		hsvTex:unbind()
		self.drawLineToEarthShader:useNone()
	end

	gl.glPointSize(1)

	gl.glPointSize(5)
	self.drawPlanetsObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	local vertexGPU = self.drawPlanetsObj.attrs.vertex.buffer
	local colorGPU = self.drawPlanetsObj.attrs.color.buffer
	local vertexCPU = vertexGPU:beginUpdate()
	local colorCPU = colorGPU:beginUpdate()
	for _,planet in ipairs(self.planets) do
		if planet.color then
			colorCPU:emplace_back():set(table.unpack(planet.color))
		else
			colorCPU:emplace_back():set(1,1,1)
		end
		vertexCPU:emplace_back():set(planet.pos:unpack())
	end
	self.drawPlanetsObj:endUpdate()
	gl.glPointSize(1)

	if self.useBlend then
		gl.glDisable(gl.GL_BLEND)
	end

	assert(glreport'here')
	App.super.update(self)
end

local rp = ffi.new('real[1]')
local function realptr(x)
	rp[0] = assert(tonumber(x))
	return rp
end

local r4p = ffi.new('real4[1]')
local function real4ptr(x)
	r4p[0].x = x.x
	r4p[0].y = x.y
	r4p[0].z = x.z
	return r4p
end

function App:recalculateSmallBodies()
	if updateBodyMethod == 'cpu' then
		for i=0,self.numBodies-1 do
			self:recalculateSmallBody(i)
		end
	elseif updateBodyMethod == 'gpu' then
		local sun = self.planets[self.planets.indexes.sun]
		self.updateCLKernel.obj:setArg(0, self.bodiesCLBuf)
		self.updateCLKernel.obj:setArg(1, realptr(self.julianDay))
		self.updateCLKernel.obj:setArg(2, realptr(self.resetDay))
		self.updateCLKernel.obj:setArg(3, real4ptr(sun.pos))
		self.updateCLKernel()
		self.bodiesCLBuf:toCPU(self.bodies)
	end
	self.bodyBuf
		:bind()
		:updateData()
		:unbind()

	if calcNearLineMethod == 'shader' then
		self.bodyToEarthBuf
			:bind()
			:updateData()
			:unbind()
	elseif calcNearLineMethod == 'fillbuffer' then
		self:updateBodyToEarthLineBuf()
	end
end

function App:recalculateSmallBody(index)
	local ke = self.bodies[index]

	local timeAdvanced = self.julianDay - self.resetDay
	local parent = self.planets[self.planets.indexes.sun]
	local orbitType = ke.orbitType

	--https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

	local meanAnomaly, meanMotion
	if orbitType == ffi.C.ORBIT_ELLIPTIC then
		meanMotion = 2 * math.pi / ke.orbitalPeriod
		meanAnomaly = ke.meanAnomalyAtEpoch + meanMotion * (self.julianDay - ke.epoch)
	elseif orbitType == ffi.C.ORBIT_HYPERBOLIC then
		meanAnomaly = ke.meanAnomaly
		meanMotion = ke.meanAnomaly / (self.julianDay - ke.timeOfPeriapsisCrossing)
	elseif orbitType == ffi.C.ORBIT_PARABOLIC then
		error'got a parabolic orbit'
	else
		error'here'
	end

	local eccentricity = ke.eccentricity

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
	for i=0,9 do
		local func, deriv
		if orbitType == ffi.C.ORBIT_PARABOLIC then
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3
			deriv = -1 - eccentricAnomaly * eccentricAnomaly
		elseif orbitType == ffi.C.ORBIT_ELLIPTIC then
			func = meanAnomaly - eccentricAnomaly + eccentricity * math.sin(eccentricAnomaly)
			deriv = -1 + eccentricity * math.cos(eccentricAnomaly)	--has zeroes ...
		elseif orbitType == ffi.C.ORBIT_HYPERBOLIC then
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
	--local fractionOffset = timeAdvanced * meanMotion / (2 * math.pi)
	local theta = timeAdvanced * meanMotion
	local pathEccentricAnomaly = eccentricAnomaly + theta
	local A = ke.A
	local B = ke.B

	--matches above
	local dt_dE
	local semiMajorAxisCubed = ke.semiMajorAxis * ke.semiMajorAxis * ke.semiMajorAxis
	if orbitType == ffi.C.ORBIT_PARABOLIC then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + pathEccentricAnomaly * pathEccentricAnomaly)
	elseif orbitType == ffi.C.ORBIT_ELLIPTIC then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - ke.eccentricity * math.cos(pathEccentricAnomaly))
	elseif orbitType == ffi.C.ORBIT_HYPERBOLIC then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (ke.eccentricity * math.cosh(pathEccentricAnomaly) - 1)
	end
	local dE_dt = 1/dt_dE
	local coeffA, coeffB
	local coeffDerivA, coeffDerivB
	if orbitType == ffi.C.ORBIT_PARABOLIC then
		--...?
	elseif orbitType == ffi.C.ORBIT_ELLIPTIC then
		coeffA = math.cos(pathEccentricAnomaly) - ke.eccentricity
		coeffB = math.sin(pathEccentricAnomaly)
		coeffDerivA = -math.sin(pathEccentricAnomaly) * dE_dt
		coeffDerivB = math.cos(pathEccentricAnomaly) * dE_dt
	elseif orbitType == ffi.C.ORBIT_HYPERBOLIC then
		coeffA = ke.eccentricity - math.cosh(pathEccentricAnomaly)
		coeffB = math.sinh(pathEccentricAnomaly)
		coeffDerivA = -math.sinh(pathEccentricAnomaly) * dE_dt
		coeffDerivB = math.cosh(pathEccentricAnomaly) * dE_dt
	end
	local posX = A[0] * coeffA + B[0] * coeffB
	local posY = A[1] * coeffA + B[1] * coeffB
	local posZ = A[2] * coeffA + B[2] * coeffB
	local velX = A[0] * coeffDerivA + B[0] * coeffDerivB	--m/day
	local velY = A[1] * coeffDerivA + B[1] * coeffDerivB
	local velZ = A[2] * coeffDerivA + B[2] * coeffDerivB

	ke.pos[0] = posX + parent.pos.s[0]
	ke.pos[1] = posY + parent.pos.s[1]
	ke.pos[2] = posZ + parent.pos.s[2]
	-- I'm not using this
	ke.vel[0] = velX + parent.vel.s[0]
	ke.vel[1] = velY + parent.vel.s[1]
	ke.vel[2] = velZ + parent.vel.s[2]

	-- now update buffers
	self.bodyToEarthArray[0+2*index].x = ke.pos[0]
	self.bodyToEarthArray[0+2*index].y = ke.pos[1]
	self.bodyToEarthArray[0+2*index].z = ke.pos[2]
	self.bodyToEarthArray[1+2*index].x = ke.pos[0]
	self.bodyToEarthArray[1+2*index].y = ke.pos[1]
	self.bodyToEarthArray[1+2*index].z = ke.pos[2]

	--ke == this.keplerianOrbitalElements
	--ke.meanAnomaly = meanAnomaly	-- this doesn't change, does it?
	ke.eccentricAnomaly = eccentricAnomaly	-- this does change ... eccentricAnomaly and pos
	--ke.fractionOffset = fractionOffset
end

function App:updateBodyToEarthLineBuf()
	local earth = self.planets[self.planets.indexes.earth]
	local e = 0
	local minLenSq = math.huge
	for i=0,self.numBodies-1 do
		local body = self.bodies[i]
		local dx = (body.pos[0] - earth.pos.x) * scale	-- scale is 1/(AU in m)
		local dy = (body.pos[1] - earth.pos.y) * scale
		local dz = (body.pos[2] - earth.pos.z) * scale
		local lenSq = dx*dx + dy*dy + dz*dz
		if lenSq < distThreshold * distThreshold then	-- in AU
print('near body', ffi.string(body.name), 'dist', math.sqrt(lenSq))
			minLenSq = math.min(minLenSq, lenSq)
			--[[
			self.bodyToEarthArray[0+2*e].x = earth.pos.x
			self.bodyToEarthArray[0+2*e].y = earth.pos.y
			self.bodyToEarthArray[0+2*e].z = earth.pos.z
			--]]
			-- [[
			-- both 1st and 2nd line vertices are used for color determination
			-- but when rendering the 1st is overridden with the earth position
			self.bodyToEarthArray[0+2*e].x = body.pos[0]
			self.bodyToEarthArray[0+2*e].y = body.pos[1]
			self.bodyToEarthArray[0+2*e].z = body.pos[2]
			--]]

			self.bodyToEarthArray[1+2*e].x = body.pos[0]
			self.bodyToEarthArray[1+2*e].y = body.pos[1]
			self.bodyToEarthArray[1+2*e].z = body.pos[2]
			e = e + 1
		end
	end
	self.numBodyToEarthLines = e
	self.bodyToEarthBuf
		:bind()
		:updateData()	--(nil, ffi.sizeof'vec3d_t' * e)
		:unbind()
	print('dt', self.julianDay, 'minLen', math.sqrt(minLenSq))
end

local function guiInputFloat(name, t, k, flags)
	flags = flags or ig.ImGuiInputTextFlags_EnterReturnsTrue
	return ig.luatableInputFloatAsText(name, t, k, flags)
end

function App:updateGUI()
	if guiInputFloat('Julian Date', self, 'julianDay') then
		self.running = 'update'
	end
	ig.igText(self.dateStr)

	guiInputFloat('dt', self, 'timeStep')

	--ig.luatableSliderFloat('alpha', self, 'alpha', 0, 1, '%.3f', 3)
	guiInputFloat('alpha', self, 'alpha')

	if ig.igButton(self.running and 'Stop' or 'Start') then
		self.running = not self.running
	end
	if ig.igButton'Step' then
		self.running = 'update'
		self.julianDay = self.julianDay + self.timeStep
	end

	ig.luatableCheckbox('useBlend', self, 'useBlend')

	if ig.igButton'Reset' then
		self.julianDay = self.resetDay
		self.running = 'update'
	end
end

local function canHandleMouse()
	if not mouse then return false end
	if rawget(ig, 'disabled') then return false end
	return not ig.igGetIO()[0].WantCaptureMouse
end

local function canHandleKeyboard()
	if rawget(ig, 'disabled') then return false end
	return not ig.igGetIO()[0].WantCaptureKeyboard
end

function App:event(event, ...)
	if App.super.event then
		App.super.event(self, event, ...)
	end
	local shiftDown = leftShiftDown or rightShiftDown
	local guiDown = leftGuiDown or rightGuiDown
	if event[0].type == sdl.SDL_EVENT_KEY_UP then
		if canHandleKeyboard() then
			if event[0].key.key == sdl.SDLK_SPACE then
				self.running = not self.running
			elseif event[0].key.key == ('r'):byte() then
				self.julianDay = self.resetDay
				self.running = 'update'
			end
		end
	end
end

return App():run()

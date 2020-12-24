#!/usr/bin/env luajit

-- visualize all orbits simultaneously

local ffi = require 'ffi'
local vec3d = require 'vec-ffi.vec3d'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local sdl = require 'ffi.sdl'
local glreport = require 'gl.report'
local glCallOrRun = require 'gl.call'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLAttribute = require 'gl.attribute'
local GLProgram = require 'gl.program'
local HSVTex = require 'gl.hsvtex'
local clnumber = require 'cl.obj.number'
local CLEnv = require 'cl.obj.env'
local template = require 'template'
local ImGuiApp = require 'imguiapp'
local Julian = require 'julian'
local Planets = require 'planets'

local matrix_ffi = require 'matrix.ffi'
matrix_ffi.real = 'float'	-- default matrix_ffi type


-- [=[ this matches parse.lua
-- these match the fields in coldesc.lua and probably row-desc.json
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

typedef struct {
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

-- [=[ This is also in solarsystem.lua
local planetColors = {
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
--]=]

local App = class(require 'glapp.orbit'(ImGuiApp))

App.title = 'JPL SSD Smallbody Visualizer'
App.viewDist = 2


--local calcNearLineMethod = 'shader'
local calcNearLineMethod = 'fillbuffer'

--local updateBodyMethod = 'cpu'
local updateBodyMethod = 'gpu'

-- scale for the smallbody data, which is in meters
local AU_in_m = 149597870700
local earthMoonDist_in_m = 380000000
local scale = 1 / AU_in_m
--earthMoonDist_in_AU = 0.002540143106462

local sunMass_kg = 1.9891e+30 
local gravitationalConstant = 6.6738480e-11							-- m^3 / (kg * s^2)
local gravitationalParameter = gravitationalConstant * sunMass_kg	--assuming the comet mass is negligible, since the comet mass is not provided
	

local modelViewMatrix = matrix_ffi.zeros(4,4)
local projectionMatrix = matrix_ffi.zeros(4,4)
local modelViewProjectionMatrix = matrix_ffi.zeros(4,4)


local dateFormat = '%4d/%02d/%02d %02d:%02d:%02d'

function App:initGL(...)
	App.super.initGL(self, ...)
	
	self.view.znear = .001
	self.view.zfar = 100


	self.hsvTex = HSVTex(256)

	-- do a one-time snapshot, so no ffwd/rewind just yet
	local t = os.date('!*t')
	self.dateStr = dateFormat:format(t.year, t.month, t.day, t.hour, t.min, t.sec)
	self.julianDate = Julian.fromCalendar(t)
	self.resetDate = self.julianDate
	self.planets = Planets.fromEphemeris(self.julianDate, 406, 'eph/406')
	for _,planet in ipairs(self.planets) do
		planet.class.color = planetColors[planet.name]
	end
	assert(glreport'here')

	local earth = self.planets[self.planets.indexes.earth]
	self.view.orbit:set((earth.pos * scale):unpack())
	self.view.pos = self.view.angle:zAxis() * self.viewDist + self.view.orbit


	local data = file['smallbodies.raw']	-- this is produced in my webgl solarsystem project, inside jpl-ssd-smallbody/output-points.template.lua 
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

	self.env = CLEnv{size=self.numBodies, real=real}
	self.bodiesCLBuf = self.env:buffer{name='bodies', type='body_t', data=self.bodies}

	self.env.code = table{
		self.env.code,
		bodyTypeCode,
	}:concat'\n'

	self.updateCLProgram = self.env:program{
		code = template([[
kernel void update(
	global body_t* bodies,
	real julianDate,
	real resetDate,
	real4 sunPos
) {
	initKernel();

	global body_t* ke = bodies + index;

	real timeAdvanced = julianDate - resetDate;
	int orbitType = ke->orbitType;

	//https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

	real meanAnomaly, meanMotion;
	if (orbitType == ORBIT_ELLIPTIC) {
		meanMotion = 2. * M_PI / ke->orbitalPeriod;
		meanAnomaly = ke->meanAnomalyAtEpoch + meanMotion * (julianDate - ke->epoch);
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		meanAnomaly = ke->meanAnomaly;
		meanMotion = ke->meanAnomaly / (julianDate - ke->timeOfPeriapsisCrossing);
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
	}
	assert(glreport'here')
	
	self.bodyPosAttr = GLAttribute{
		buffer = self.bodyBuf,
		size = 3,
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
	}
	assert(glreport'here')
	
	self.bodyToEarthPosAttr = GLAttribute{
		buffer = self.bodyToEarthBuf,
		size = 3,
		type = gl.GL_DOUBLE,
	}
	assert(glreport'here')



	-- [=[ in absence of geometry buffers, I was trying to double up the body pos buffer, but that seems to be running slow
	if calcNearLineMethod == 'shader' then
		self.drawLineToEarthShader = GLProgram{
			vertexCode = [[
#version 460
attribute vec4 bodyPos;
uniform vec3 earthPos;
varying float lum;
uniform mat4 modelViewProjectionMatrix;
void main() {
	vec4 pos = bodyPos;
	float dist = length(earthPos - pos.xyz);
	lum = 1. - ]]..clnumber(scale)..[[ * dist;
	if (gl_VertexID % 2 == 0) {
		pos = vec4(earthPos, 1.);
	}
	gl_Position = modelViewProjectionMatrix * pos;
}
]],
			fragmentCode = [[
varying float lum;
void main() {
	gl_FragColor = vec4(lum, 0., 0., 1.);
}
]],
		}
		self.drawLineToEarthShader:useNone()
		assert(glreport'here')
	elseif calcNearLineMethod == 'fillbuffer' then
		self.drawLineToEarthShader = GLProgram{
			vertexCode = [[
#version 460
attribute vec4 bodyPos;
uniform vec3 earthPos;
varying vec3 color;
uniform sampler2D hsvTex;
uniform mat4 modelViewProjectionMatrix;
void main() {
	vec4 pos = bodyPos;
	float dist = length(earthPos - pos.xyz);
	float lum = 1. - ]]..clnumber(scale)..[[ * dist;
	color = texture2D(hsvTex, vec2(lum, .5)).rgb;
	if (gl_VertexID % 2 == 0) {
		pos = vec4(earthPos, 1.);
	}
	gl_Position = modelViewProjectionMatrix * pos;
}
]],
			fragmentCode = [[
varying vec3 color;
void main() {
	//gl_FragColor = vec4(color, 1.);
	gl_FragColor = vec4(1., 0., 0., 1.);
}
]],
		
			uniforms = {
				hsvTex = 0,
			}
		}
		self.drawLineToEarthShader:useNone()
	
		
		-- trying another trick: just fill the buffer ... after loading planets
		-- if we aren't using one buffer for the lines then update it here after planets are loaded
		self:updateBodyToEarthLineBuf()
	end


	if calcNearLineMethod == 'shader' then
		self.drawLineToEarthShader:setAttrs{
			bodyPos = self.bodyPosAttr,
		}
	elseif calcNearLineMethod == 'fillbuffer' then
		self.drawLineToEarthShader:setAttrs{
			bodyPos = self.bodyToEarthPosAttr,
		}
	end
	

	assert(glreport'here')


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


	gl.glEnable(gl.GL_POINT_SMOOTH)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE)
	gl.glDisable(gl.GL_DEPTH_TEST)
	assert(glreport'here')
end

App.timeStep = .1

function App:update()
	self:draw()

	if self.running then
		if self.running == true then
			self.julianDate = self.julianDate + self.timeStep
		end
		local t = Julian.toCalendar(self.julianDate)
		self.dateStr = dateFormat:format(t.year, t.month, t.day, t.hour, t.min, t.sec)
		self.planets = Planets.fromEphemeris(self.julianDate)
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

function App:draw()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glViewport(0, 0, self.width, self.height)

	self.view:setup(self.width / self.height)
	assert(glreport'here')


	gl.glBegin(gl.GL_LINES)
	gl.glColor3f(1,0,0) gl.glVertex3f(0,0,0) gl.glVertex3f(1,0,0)
	gl.glColor3f(0,1,0) gl.glVertex3f(0,0,0) gl.glVertex3f(0,1,0)
	gl.glColor3f(0,0,1) gl.glVertex3f(0,0,0) gl.glVertex3f(0,0,1)
	gl.glEnd()


	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glPushMatrix()
	gl.glScalef(scale, scale, scale)

	self.drawlist = self.drawlist or {}

	-- TODO use the original binary blob, and just pass it as a gl buffer, then use pos as a strided vertex array
	gl.glPointSize(3)
	gl.glColor3f(self.alpha, self.alpha, self.alpha)
	
	--[[ raw glVertex calls / with call lists
	do --glCallOrRun(self.drawlist, function()
	gl.glBegin(gl.GL_POINTS)
	for i=0,self.numBodies-1 do
		local body = self.bodies[i]
		gl.glVertex3dv(body.pos)
	end
	gl.glEnd()
	end --)
	--]]
	--[[ client state vertex pointer without array buffers
	gl.glVertexPointer(
		3, --self.numBodies * ffi.sizeof'body_t',
		real == 'double' and gl.GL_DOUBLE or gl.GL_FLOAT,
		ffi.sizeof'body_t',
		self.bodies[0].pos)
	gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
	gl.glDrawArrays(gl.GL_POINTS, 0, self.numBodies)
	gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
	--]]
	--[[ gl array buffers + client state vertex pointer
	self.bodyBuf:bind()
	gl.glVertexPointer(
		3, --self.numBodies * ffi.sizeof'body_t',
		real == 'double' and gl.GL_DOUBLE or gl.GL_FLOAT,
		ffi.sizeof'body_t',
		ffi.cast('uint8_t*', ffi.offsetof('body_t', 'pos')))
	gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
	gl.glDrawArrays(gl.GL_POINTS, 0, self.numBodies)
	gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
	self.bodyBuf:unbind()
	--]]
	-- [[ gl array buffers + vertex attrib arrays
	self.bodyBuf:bind()
	gl.glEnableVertexAttribArray(0)
	gl.glDrawArrays(gl.GL_POINTS, 0, self.numBodies)
	gl.glDisableVertexAttribArray(0)
	self.bodyBuf:unbind()
	--]]
	assert(glreport'here')
	
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, projectionMatrix.ptr)
	modelViewProjectionMatrix:mul(projectionMatrix, modelViewMatrix)

	local earth = self.planets[self.planets.indexes.earth]
	if calcNearLineMethod == 'shader' then
		-- draw all lines and use a shader to lighten/darken them
		-- draw line to earth
		-- a geometry shader would be nice ... so I can generate one point at the earth's position
		-- until then, i'll make a second buffer
		self.drawLineToEarthShader:use()

		gl.glUniformMatrix4fv(self.drawLineToEarthShader.uniforms.modelViewProjectionMatrix.loc, 1, false, modelViewProjectionMatrix.ptr)
		
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
		
		gl.glUniformMatrix4fv(self.drawLineToEarthShader.uniforms.modelViewProjectionMatrix.loc, 1, false, modelViewProjectionMatrix.ptr)
		
		if self.drawLineToEarthShader.uniforms.earthPos then
			--gl.glUniform3dv(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos.s)
			gl.glUniform3f(self.drawLineToEarthShader.uniforms.earthPos.loc, earth.pos:unpack())
		end
		self.hsvTex:bind()
		
		gl.glEnableVertexAttribArray(0)
		gl.glDrawArrays(gl.GL_LINES, 0, 2 * self.numBodyToEarthLines)
		gl.glDisableVertexAttribArray(0)
		
		self.hsvTex:unbind()
		self.drawLineToEarthShader:useNone()
	end

	gl.glPointSize(1)

	gl.glPointSize(5)
	gl.glBegin(gl.GL_POINTS)
	for _,planet in ipairs(self.planets) do
		gl.glColor3f(table.unpack(planet.color))
		gl.glVertex3dv(planet.pos.s)
	end
	gl.glEnd()
	gl.glPointSize(1)
	
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glPopMatrix()

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
		self.updateCLKernel.obj:setArg(1, realptr(self.julianDate))
		self.updateCLKernel.obj:setArg(2, realptr(self.resetDate))
		self.updateCLKernel.obj:setArg(3, real4ptr(sun.pos))
		self.updateCLKernel()
		self.bodiesCLBuf:toCPU(self.bodies)
	end
	self.bodyBuf:updateData()

	if calcNearLineMethod == 'shader' then
		self.bodyToEarthBuf:updateData()
	elseif calcNearLineMethod == 'fillbuffer' then
		self:updateBodyToEarthLineBuf()
	end
end

function App:recalculateSmallBody(index)
	local ke = self.bodies[index]
	
	local timeAdvanced = self.julianDate - self.resetDate
	local parent = self.planets[self.planets.indexes.sun]
	local orbitType = ke.orbitType

	--https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

	local meanAnomaly, meanMotion
	if orbitType == ffi.C.ORBIT_ELLIPTIC then
		meanMotion = 2 * math.pi / ke.orbitalPeriod
		meanAnomaly = ke.meanAnomalyAtEpoch + meanMotion * (self.julianDate - ke.epoch)
	elseif orbitType == ffi.C.ORBIT_HYPERBOLIC then
		meanAnomaly = ke.meanAnomaly
		meanMotion = ke.meanAnomaly / (self.julianDate - ke.timeOfPeriapsisCrossing)
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
		if lenSq < .04 then	-- 1 AU
			minLenSq = math.min(minLenSq, lenSq)
			self.bodyToEarthArray[0+2*e].x = earth.pos.x
			self.bodyToEarthArray[0+2*e].y = earth.pos.y
			self.bodyToEarthArray[0+2*e].z = earth.pos.z
			self.bodyToEarthArray[1+2*e].x = body.pos[0]
			self.bodyToEarthArray[1+2*e].y = body.pos[1]
			self.bodyToEarthArray[1+2*e].z = body.pos[2]
			e = e + 1
		end
	end
	self.numBodyToEarthLines = e
	self.bodyToEarthBuf:updateData()	--(nil, ffi.sizeof'vec3d_t' * e)
	print('dt', self.julianDate, 'minLen', math.sqrt(minLenSq))
end

-- TODO consolidate this, I use this trick in so many other projects
local fp = ffi.new('float[1]')

local function guiSliderFloat(name, t, k, vmin, vmax, format, power)
	fp[0] = assert(tonumber(t[k]))
	if ig.igSliderFloat(name, fp, vmin, vmax, format, power) then
		t[k] = fp[0]
		return true
	end
end

function guiInputFloat(name, t, k, step, stepfast, format, flags)
	step = step or .1
	stepfast = stepfast or 1
	format = format or '%.3f'
	flags = flags or ig.ImGuiInputTextFlags_EnterReturnsTrue
	fp[0] = assert(tonumber(t[k]))
	if ig.igInputFloat(name, fp, step, stepfast, format, flags) then
		t[k] = fp[0]
		return true
	end
end

function App:updateGUI()
	if guiInputFloat('Julian Date', self, 'julianDate', 1) then
		self.running = 'update'
	end
	ig.igText(self.dateStr)

	guiInputFloat('dt', self, 'timeStep')

	--guiSliderFloat('alpha', self, 'alpha', 0, 1, '%.3f', 3)
	guiInputFloat('alpha', self, 'alpha')

	if ig.igButton(self.running and 'Stop' or 'Start') then
		self.running = not self.running
	end
	if ig.igButton'Step' then
		self.running = 'update'
		self.julianDate = self.julianDate + self.timeStep
	end	

	if ig.igButton'Reset' then
		self.julianDate = self.resetDate 
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
	if event.type == sdl.SDL_KEYUP then
		if canHandleKeyboard() then
			if event.key.keysym.sym == sdl.SDLK_SPACE then
				self.running = not self.running
			elseif event.key.keysym.sym == ('r'):byte() then
				self.julianDate = self.resetDate
				self.running = 'update'
			end
		end
	end
end

return App():run()

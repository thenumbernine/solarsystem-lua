#! /usr/bin/env luajit
require 'ext'

local julian = require 'julian'
local Planets = require 'planets'

local currentDate = julian.fromCalendar(os.date('!*t'))

local function refreshGraphs(startDate, endDate)
	local planets = table()
	local julianDates = table()

	local divs = 2000	-- x resolution
	for i=1,divs do
		local dayOffset = (i-1)/(divs-1) * (endDate - startDate) + startDate
		julianDates[i] = dayOffset
		planets[i] = Planets.fromEphemeris(currentDate + dayOffset, 406, 'eph/406')
	end

	--[=[ graph angles between planets
	-- for earth, considering all other planets, find angles between all other planets
	local names = Planets.indexes:keys()
	names:removeObject('earth')
	names:removeObject('moon')

	local graphs = table()
	for i=1,#names do
		local na = names[i]
		local nb = 'moon' 
		graphs[na..' -> earth -> '..nb] = {
			enabled = true,
			julianDates,
			table.map(planets, function(planets)
				local va = planets[Planets.indexes[na]].pos - planets[Planets.indexes.earth].pos
				local vb = planets[Planets.indexes[nb]].pos - planets[Planets.indexes.earth].pos
				va = va:normalize()
				vb = vb:normalize()
				return math.deg(math.acos(math.clamp(va:dot(vb), -1, 1)))
			end),
		}
	end
	for _,k in ipairs{0,90,120,150,180} do
		graphs[tostring(k)..' degrees'] = {
			enabled = true,
			color = {1,1,1},
			{startDate, endDate},
			{k,k},
		}
	end
	--]=]
	-- [=[ graph planet distances
	local names = Planets.indexes:keys()
	names:removeObject('earth')
	local graphs = table()
	for i,name in ipairs(names) do
		graphs[name..' <-> earth'] = {
			enabled = true,
			julianDates,
			table.map(planets, function(planets)
				return (planets[Planets.indexes[name]].pos - planets[Planets.indexes.earth].pos):length()
			end),
		}
	end
	--]=]

	return graphs
end

local Plot2DApp =  require 'plot2d.app'
local box2 = require 'vec.box2'
local App = class(Plot2DApp)

App.guiScale = {8,8}
function App:init(...)
	App.super.init(self, ...)
	self.viewbboxLast = box2()
end

function App:update()
	App.super.update(self)

	if self.viewbbox.min[1] ~= self.viewbboxLast.min[1]
	or self.viewbbox.max[1] ~= self.viewbboxLast.max[1]
	then
		print('refreshing to ',self.viewbbox.min[1], self.viewbbox.max[1])
		self.dontResetView = true
		local newGraphs = refreshGraphs(self.viewbbox.min[1], self.viewbbox.max[1])
		for name,graph in pairs(self.graphs) do 
			graph[1] = newGraphs[name][1]
			graph[2] = newGraphs[name][2]
		end
		self.dontResetView = false
		self.viewbboxLast.min[1] = self.viewbbox.min[1]
		self.viewbboxLast.max[1] = self.viewbbox.max[1]
	end
end

function App:resetView()
	if self.dontResetView then return end
	App.super.resetView(self)
end

local graphs = refreshGraphs(-100,100)
App(graphs):run()

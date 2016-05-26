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

	-- for earth, considering all other planets, find angles between all other planets
	local names = Planets.indexes:keys()
	names:removeObject('earth')

	local graphs = table()
	for i=1,#names-1 do
		for j=i+1,#names do
			local na = names[i]
			local nb = names[j]
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
	end
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
		self:setGraphInfo(refreshGraphs(self.viewbbox.min[1], self.viewbbox.max[1]))
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

#! /usr/bin/env luajit
require 'ext'

local julian = require 'julian'
local Planets = require 'planets'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'

local planetNames = table.map(Planets.planetClasses, function(cl) return cl.name end)
print(planetNames:concat', ')

local function gregstr(t)
	return ('%4d/%02d/%02d %02d:%02d:%02d'):format(
		t.year,
		t.month,
		t.day,
		t.hour,
		t.min,
		t.sec)
end

local t = os.date('!*t')
print('current time: '..gregstr(t))
local currentDate = julian.fromCalendar(t)


local Plot2DApp =  require 'plot2d.app'
local box2 = require 'vec.box2'
local App = class(Plot2DApp)


function App:refreshGraphs(startDate, endDate)
	print('refreshing '
		..gregstr(julian.toCalendar(currentDate + startDate))
		..' -> '
		..gregstr(julian.toCalendar(currentDate + endDate)))
	local planets = table()
	local julianDates = table()

	local divs = 2000	-- x resolution
	for i=1,divs do
		local dayOffset = (i-1)/(divs-1) * (endDate - startDate) + startDate
		julianDates[i] = dayOffset
		planets[i] = Planets.fromEphemeris(currentDate + dayOffset, 406, 'eph/406')
	end

	-- [=[ graph angles between planets
	local graphs = table()
	for i,ni in ipairs(planetNames) do
		for j,nj in ipairs(planetNames) do
			for k,nk in ipairs(planetNames) do
				local name = table{ni,nj,nk}:concat' -> '
				if self.graphsEnabledForName[name] then
					graphs[name] = {
						enabled = true,
						julianDates,
						table.map(planets, function(planets,i)
							local va = planets[Planets.indexes[ni]].pos - planets[Planets.indexes[nj]].pos
							local vb = planets[Planets.indexes[nk]].pos - planets[Planets.indexes[nj]].pos
							va = va:normalize()
							vb = vb:normalize()
							local theta = math.deg(math.acos(math.clamp(va:dot(vb), -1, 1)))
							return theta
						end),
					}
				end
			end
		end
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
	--[=[ graph planet distances
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

function App:init()
	self.graphsEnabledForName = table()
	--local graphs = self:refreshGraphs(-100, 100)	-- before refresh ...
	App.super.init(self, {})
	self.viewsize = vec2(100, 100)
	self.viewbboxLast = box2()
end

local mouseAngle = 0
local mouseTime = 0
function App:update()
	App.super.update(self)

	if self.viewbbox.min[1] ~= self.viewbboxLast.min[1]
	or self.viewbbox.max[1] ~= self.viewbboxLast.max[1]
	or self.changedCheckbox
	then
		self.changedCheckbox = false
		self.dontResetView = true
		--[[
		local newGraphs = self:refreshGraphs(self.viewbbox.min[1], self.viewbbox.max[1])
		for name,graph in pairs(self.graphs) do 
			graph[1] = newGraphs[name][1]
			graph[2] = newGraphs[name][2]
		end
		--]]
		self:setGraphInfo(self:refreshGraphs(self.viewbbox.min[1], self.viewbbox.max[1]))
		self.dontResetView = false
		self.viewbboxLast.min[1] = self.viewbbox.min[1]
		self.viewbboxLast.max[1] = self.viewbbox.max[1]
	end

	local fx, fy = table.unpack(self.mousepos)
	mouseTime = self.viewbbox.min[1] * (1 - fx) + self.viewbbox.max[1] * fx
	mouseAngle = self.viewbbox.min[2] * (1 - fy) + self.viewbbox.max[2] * fy
end

function App:resetView()
	if self.dontResetView then return end
	App.super.resetView(self)
end

function App:getCoordText()
	local j = julian.toCalendar(mouseTime + currentDate)
	return 'mouse angle: '..tolua(mouseAngle)..'\n'
		..'mouse time (jul): '..tolua(mouseTime)..'\n'
		..'mouse time (greg): '..gregstr(j)
end

local bool = ffi.new'bool[1]'
function App:updateGUI()
	ig.igPushIDStr'planetCheckboxes'
	for j,nj in ipairs(planetNames) do
		ig.igPushIDStr(''..j)
		ig.igText(nj..':')
		ig.igSameLine()
		if ig.igCollapsingHeader'' then
			local all = true	
			for i=1,#planetNames-1 do
				if i ~= j then
					for k=i+1,#planetNames do
						if k ~= j 
						and i ~= k 
						then
							local ni = planetNames[i]
							local nk = planetNames[k]
							local name = table{ni,nj,nk}:concat' -> '
							local name2 = table{nk,nj,ni}:concat' -> '
							if not self.graphsEnabledForName[name] 
							and not self.graphsEnabledForName[name2] 
							then
								all = false
								break
							end
						end
						if not all then break end
					end
					if not all then break end
				end
			end
			bool[0] = all
			if ig.igCheckbox('all', bool) then
				all = bool[0]
				for i=1,#planetNames-1 do
					if i ~= j then
						for k=i+1,#planetNames do
							if k ~= j 
							and i ~= k
							then
								local ni = planetNames[i]
								local nk = planetNames[k]
								local name = table{ni,nj,nk}:concat' -> '
								local name2 = table{nk,nj,ni}:concat' -> '
								self.graphsEnabledForName[name] = all
								self.graphsEnabledForName[name2] = nil
								self.changedCheckbox = true
							end
						end
					end
				end
			end
			for i,ni in ipairs(planetNames) do
				if i ~= j then
					ig.igPushIDStr(''..i)
					for k,nk in ipairs(planetNames) do
						if k ~= j then
							ig.igPushIDStr(''..k)
							local name = table{ni,nj,nk}:concat' -> '
							local name2 = table{nk,nj,ni}:concat' -> '
							bool[0] = not not (self.graphsEnabledForName[name] or self.graphsEnabledForName[name2])
							if ig.igCheckbox('', bool) then
								print('toggling '..name..' to '..bool[0])
								self.graphsEnabledForName[name] = bool[0] or nil
								self.changedCheckbox = true
							end
							if ig.igIsItemHovered(ig.ImGuiHoveredFlags_Default) then
								ig.igBeginTooltip()
								ig.igText(name)
								ig.igEndTooltip()
							end
							if k ~= #planetNames then
								ig.igSameLine()
							end
							ig.igPopID()
						end
					end
					ig.igPopID()
				end
			end
		end
		ig.igPopID()
	end
	ig.igPopID()

	App.super.updateGUI(self)
end

App():run()

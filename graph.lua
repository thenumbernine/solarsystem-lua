#! /usr/bin/env luajit
local table = require 'ext.table'
local class = require 'ext.class'
local tolua = require 'ext.tolua'
local julian = require 'julian'
local Planets = require 'planets'
local ffi = require 'ffi'

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

local useLog = false 

local Plot2DApp =  require 'plot2d.app'
local box2 = require 'vec.box2'
local App = class(Plot2DApp)
App.title = 'ephemeris data graph' 

-- new problem says ffi.imgui has to be require'd after imguiapp
local ig = require 'ffi.imgui'

function App:refreshGraphs(startDate, endDate)
	print('refreshing '
		..gregstr(julian.toCalendar(currentDate + startDate))
		..' -> '
		..gregstr(julian.toCalendar(currentDate + endDate)))
	local planets = table()
	local julianDates = table()

	local divs = self.width	-- x resolution
	for i=1,divs do
		local dayOffset = (i-1)/(divs-1) * (endDate - startDate) + startDate
		julianDates[i] = dayOffset
		planets[i] = Planets.fromEphemeris(currentDate + dayOffset, 406, 'eph/406')
	end

	-- save for later, for getCoordText
	-- TODO save the angle calculations as well
	self.planets = planets

	-- [=[ graph angles between planets
	local graphs = self.graphs or table()
	for j=1,#planetNames do
		local nj = planetNames[j]
		for i=1,#planetNames-1 do
			if i ~= j then
				local ni = planetNames[i]
				for k=i+1,#planetNames do
					if k ~= j then
						local nk = planetNames[k]
						local name = table{ni,nj,nk}:concat' -> '
						local name2 = table{nk,nj,ni}:concat' -> '
						if self.angleGraphsEnabledForName[name] 
						or self.angleGraphsEnabledForName[name2]
						then
							local srcGraph = graphs[name]
							local srcEnabled = true
							if srcGraph then srcEnabled = srcGraph.enabled end
							graphs[name] = {
								enabled = srcEnabled,
								color = srcGraph and srcGraph.color or nil,
								julianDates,
								table.map(planets, function(planets,i)
									local va = planets[Planets.indexes[ni]].pos - planets[Planets.indexes[nj]].pos
									local vb = planets[Planets.indexes[nk]].pos - planets[Planets.indexes[nj]].pos
									va = va:normalize()
									vb = vb:normalize()
									local theta = math.deg(math.acos(math.clamp(va:dot(vb), -1, 1)))
									if useLog then theta = math.log(theta) end
									return theta
								end),
							}
						else
							if graphs[name] then
								graphs[name] = nil
							end
						end
					end
				end
			end
		end
	end
	
	for i=1,#planetNames-1 do
		if i ~= j then
			local ni = planetNames[i]
			for k=i+1,#planetNames do
				if k ~= j then
					local nk = planetNames[k]
					local name = table{ni,nk}:concat' <-> '
					local name2 = table{nk,ni}:concat' <-> '
					if self.distanceGraphsEnabledForName[name] 
					or self.distanceGraphsEnabledForName[name2]
					then
						local srcGraph = graphs[name]
						local srcEnabled = true
						if srcGraph then srcEnabled = srcGraph.enabled end
						graphs[name] = {
							enabled = srcEnabled,
							color = srcGraph and srcGraph.color or nil,
							julianDates,
							table.map(planets, function(planets)
								local dist = (planets[Planets.indexes[ni]].pos - planets[Planets.indexes[nk]].pos):length()
								if useLog then dist = math.log(dist) end
								return dist
							end),
						}
					else
						if graphs[name] then
							graphs[name] = nil
						end
					end
				end
			end
		end
	end
	
	for _,k in ipairs{0,90,120,150,180} do
		local name = tostring(k)..' degrees'
		if useLog then k = math.log(k) end
		local srcGraph = graphs[name]
		local srcEnabled = true
		if srcGraph then srcEnabled = srcGraph.enabled end
		graphs[name] = {
			enabled = srcEnabled,
			color = srcGraph and srcGraph.color or {1,1,1},
			{startDate, endDate},
			{k,k},
		}
	end

	return graphs
end

function App:init()
	self.angleGraphsEnabledForName = table()
	self.distanceGraphsEnabledForName = table()
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
	or self.width ~= self.lastWidth
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
		self.lastWidth = self.width
	end

	local fx, fy = table.unpack(self.mousepos)
	mouseTime = self.viewbbox.min[1] * (1 - fx) + self.viewbbox.max[1] * fx
	mouseAngle = self.viewbbox.min[2] * fy + self.viewbbox.max[2] * (1 - fy)
end

function App:resetView()
	if self.dontResetView then return end
	App.super.resetView(self)
end

function App:getCoordText()
	local j = julian.toCalendar(mouseTime + currentDate)
	local s = 
	--'graph y: '..tolua(mouseAngle)..'\n'
		'mouse time (+jul): '..tolua(mouseTime)..'\n'
		..'mouse time (greg): '..gregstr(j)
	
	local angles = table()
	if self.showAnglesInCoordText 
	-- make sure our cached data is available
	and self.planets
	-- make sure we're not mouseover the gui
	and not ig.igGetIO()[0].WantCaptureMouse
	then
		local graphs = self.graphs or table()
		for j=1,#planetNames do
			local nj = planetNames[j]
			for i=1,#planetNames-1 do
				if i ~= j then
					local ni = planetNames[i]
					for k=i+1,#planetNames do
						if k ~= j then
							local nk = planetNames[k]
							local name = table{ni,nj,nk}:concat' -> '
							local name2 = table{nk,nj,ni}:concat' -> '
							if self.angleGraphsEnabledForName[name] 
							or self.angleGraphsEnabledForName[name2]
							then
								local x = math.floor(self.mousepos[1] * self.width)
								local planetData = self.planets[x]
								if planetData then
									local va = planetData[Planets.indexes[ni]].pos - planetData[Planets.indexes[nj]].pos
									local vb = planetData[Planets.indexes[nk]].pos - planetData[Planets.indexes[nj]].pos
									va = va:normalize()
									vb = vb:normalize()
									local theta = math.deg(math.acos(math.clamp(va:dot(vb), -1, 1)))
									--if useLog then theta = math.log(theta) end
									angles:insert{name, theta}
								end
							end
						end
					end
				end
			end
		end
	
		for i=1,#planetNames-1 do
			local ni = planetNames[i]
			for k=i+1,#planetNames do
				local nk = planetNames[k]
				local name = table{ni,nk}:concat' <-> '
				if self.distanceGraphsEnabledForName[name] then
					local x = math.floor(self.mousepos[1] * self.width)
					local planetData = self.planets[x]
					if planetData then
						local dist = (planetData[Planets.indexes[ni]].pos - planetData[Planets.indexes[nk]].pos):length()
						--if useLog then dist = math.log(dist) end
						angles:insert{name, dist}
					end
				end
			end
		end
	end
	
	s = s .. '\n' .. angles
		:sort(function(a,b) return a[2] > b[2] end)
		:mapi(function(p) return p[1] .. ' = ' .. p[2] end)
		:concat'\n'
	
	return s
end

local bool = ffi.new'bool[1]'
function App:updateGUI()
	if ig.igCollapsingHeader'conjunctions:' then
		ig.igPushID_Str'angleGraphs'
		for j,nj in ipairs(planetNames) do
			ig.igPushID_Str(''..j)
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
								if not self.angleGraphsEnabledForName[name] 
								and not self.angleGraphsEnabledForName[name2] 
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
									self.angleGraphsEnabledForName[name] = all
									self.angleGraphsEnabledForName[name2] = nil
									self.changedCheckbox = true
								end
							end
						end
					end
				end
				for i,ni in ipairs(planetNames) do
					if i ~= j then
						ig.igPushID_Str(''..i)
						for k,nk in ipairs(planetNames) do
							if k ~= j then
								ig.igPushID_Str(''..k)
								local name = table{ni,nj,nk}:concat' -> '
								local name2 = table{nk,nj,ni}:concat' -> '
								bool[0] = not not (self.angleGraphsEnabledForName[name] or self.angleGraphsEnabledForName[name2])
								if ig.igCheckbox('', bool) then
									self.angleGraphsEnabledForName[name] = bool[0] or nil
									self.angleGraphsEnabledForName[name2] = bool[0] or nil
									self.changedCheckbox = true
								end
								if ig.igIsItemHovered(ig.ImGuiHoveredFlags_None) then
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
	end
	if ig.igCollapsingHeader'distances:' then
		ig.igPushID_Str'distances'
		
		local all = true	
		for i=1,#planetNames-1 do
			if i ~= j then
				for k=i+1,#planetNames do
					if k ~= j 
					and i ~= k 
					then
						local ni = planetNames[i]
						local nk = planetNames[k]
						local name = table{ni,nk}:concat' <-> '
						local name2 = table{nk,ni}:concat' <-> '
						if not self.distanceGraphsEnabledForName[name] 
						and not self.distanceGraphsEnabledForName[name2] 
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
							local name = table{ni,nk}:concat' <-> '
							local name2 = table{nk,ni}:concat' <-> '
							self.distanceGraphsEnabledForName[name] = all
							self.distanceGraphsEnabledForName[name2] = nil
							self.changedCheckbox = true
						end
					end
				end
			end
		end
		
		for i,ni in ipairs(planetNames) do
			ig.igPushID_Str(''..i)
			for k,nk in ipairs(planetNames) do
				ig.igPushID_Str(''..k)
				local name = table{ni,nk}:concat' <-> '
				local name2 = table{nk,ni}:concat' <-> '
				bool[0] = not not (self.distanceGraphsEnabledForName[name] or self.distanceGraphsEnabledForName[name2])
				if ig.igCheckbox('', bool) then
					self.distanceGraphsEnabledForName[name] = bool[0] or nil
					self.distanceGraphsEnabledForName[name2] = bool[0] or nil
					self.changedCheckbox = true
				end
				if ig.igIsItemHovered(ig.ImGuiHoveredFlags_None) then
					ig.igBeginTooltip()
					ig.igText(name)
					ig.igEndTooltip()
				end
				if k ~= #planetNames then
					ig.igSameLine()
				end
				ig.igPopID()
			end
			ig.igPopID()
		end
		ig.igPopID()
	end

	bool[0] = not not self.showAnglesInCoordText
	if ig.igCheckbox('show angles/distances in coord text', bool) then
		self.showAnglesInCoordText = bool[0] or nil
	end

	bool[0] = not not useLog
	if ig.igCheckbox('use log', bool) then
		useLog = bool[0] or nil
		self.changedCheckbox = true
	end

	App.super.updateGUI(self)
end

App():run()

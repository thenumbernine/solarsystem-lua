#!/usr/bin/env luajit

-- https://earthquake.usgs.gov/earthquakes/search/

require 'ext'

function math.fpart(x)
	return x - math.floor(x)
end

local fn = 'earthquake-catalog.txt'
local ls = assert(path(fn):read()):trim():split('[\r\n]')
local rows = table()
for _,l in ipairs(ls) do
	if not l:match('^%s') then	-- new entry
		rows:insert(l:trim():split('%s+'))
	else	-- append to last entry
		rows[#rows]:append(l:trim():split('%s+'))
	end
end

-- TODO FIXME every so often there's a single line with words like "min." ...
local entries = table()
for _,row in ipairs(rows) do
	print(row:concat(' | '))
	local year, month, day
	if row[1] then
		year, month, day = row[1]:match('^(%d%d%d%d?)(%d%d)(%d%d)$')
	else
		print("bad row: "..row:concat('\t'))
	end
	local hour, min, sec
	if row[5] then
		hour, min, sec = row[5]:match('^(%d%d?)(%d%d)(%d%d%.?%d*)$')
	end
	
	print(year,month,day,hour,min,sec)
	
	local entry = {
		year=tonumber(year),
		month=tonumber(month),
		day=tonumber(day),
		hour=tonumber(hour),
		min=tonumber(min),
		sec=tonumber(sec),
		lat=tonumber(row[3]),
		lon=tonumber(row[4]),
		depth=tonumber(row[7]),	-- km
		momentMagn=tonumber(row[8]),
		momentMagnVar=tonumber(row[9]),
	}

	-- TODO just use julian date?  since os.time() doesn't work with anything before 1970...
	if entry.year and entry.month and entry.day then
		entry.time = os.time(entry)
		if entry.time and entry.sec then entry.time = entry.time + math.fpart(entry.sec) end
	end
	
	entry.area = row[6]
	entries:insert(entry)
end

local result, json = pcall(require, 'dkjson')
if result then
	path'earthquakes.json':write(json.encode(entries, {indent=true}))
end
path'earthquakes.lua':write('{\n'..entries:map(function(entry)
	return '\t'..tolua(entry)..',\n'
end):concat()..'}')

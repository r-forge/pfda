#!/usr/bin/env  ruby

system "taskkill /F /IM rgui.exe /T"

Dfile =  File.join("pkg","DESCRIPTION")
description = IO.read(Dfile)
if description =~ /Package:\s+(\w+)/
	pkgname = $1
else
	raise "Bad description file"
end
if description =~ /Version:\s+([\d\.]+)/
	pkgversion = $1
else
	raise "Bad description file"
end

p pkgname + '_' + pkgversion + '.tar.gz'

system "Rcmd build pkg"
system "Rcmd INSTALL --build #{pkgname + '_' + pkgversion + '.tar.gz'}"

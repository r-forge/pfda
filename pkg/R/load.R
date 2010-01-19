.First.lib <- function(lib, pkg)  library.dynam("pfda", pkg, lib)
.Last.lib<-function(lib,pkg) library.dynam.unload("pfda", pkg, lib)
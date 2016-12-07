.onLoad = function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
}

.onAttach = function(libname, pkgname){
  num_gigs_ram_available = .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory") / 1e9
  packageStartupMessage(paste("Welcome to GreedyExperimentalDesign v", utils::packageVersion("GreedyExperimentalDesign"), ".\n", round(num_gigs_ram_available, 2), "GB memory available for the Java searcher.\n", sep = ""))
}
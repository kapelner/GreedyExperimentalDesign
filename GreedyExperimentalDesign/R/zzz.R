.onUnload = function(libpath) {
    try(ged_gpu_release_context_cpp(), silent = TRUE)
    try(library.dynam.unload("GreedyExperimentalDesign", libpath), silent = TRUE)
}

.onLoad = function(libname, pkgname) {
  # On Linux and macOS, try to pre-load wgpu-native so the dynamic linker can
  # find it before the package .so/.dylib is loaded (RPATH alone is unreliable
  # in R's load path on some systems).
  sysname = Sys.info()["sysname"]
  wgpu_libname = switch(sysname, Linux = "libwgpu_native.so", Darwin = "libwgpu_native.dylib", Windows = "wgpu_native.dll", NULL)
  if (!is.null(wgpu_libname)) {
      wgpu_lib = system.file(file.path("lib", wgpu_libname), package = pkgname)
      if (wgpu_lib == "") {
          proj_lib = file.path(getwd(), "..", "lib", "wgpu", "lib", wgpu_libname)
          if (file.exists(proj_lib)) {
              try(dyn.load(proj_lib, local = FALSE, now = TRUE), silent = TRUE)
          }
      } else {
          try(dyn.load(wgpu_lib, local = FALSE, now = TRUE), silent = TRUE)
      }
  }

  dll_loaded = tryCatch({
    library.dynam("GreedyExperimentalDesign", pkgname, libname)
    TRUE
  }, error = function(e) FALSE)

  .jpackage(pkgname, lib.loc = libname)

  if (dll_loaded && is.null(getOption("GreedyExperimentalDesign.use_gpu"))) {
	  options(GreedyExperimentalDesign.use_gpu = isTRUE(ged_gpu_available()))
  }

  #need to check if proper Java is installed by special request of Prof Brian Ripley
  java_spec_version = .jcall("java/lang/System", "S", "getProperty", "java.specification.version")
  java_major = suppressWarnings(as.integer(sub("^1\\.", "", sub("\\..*$", "", java_spec_version))))
  if (is.na(java_major) || java_major < 21L) {
      stop("GreedyExperimentalDesign requires Java 21 or higher; detected Java ", java_spec_version, ".")
  }

}

.onAttach = function(libname, pkgname) {
  runtime = .jcall("java/lang/Runtime", "Ljava/lang/Runtime;", "getRuntime")
  num_gigs_ram_available = .jcall(runtime, "J", "maxMemory") / 1e9
  packageStartupMessage(paste("Welcome to GreedyExperimentalDesign v", utils::packageVersion("GreedyExperimentalDesign"), ".\n", round(num_gigs_ram_available, 2), "GB memory available for the Java searcher.\n", sep = ""))
}

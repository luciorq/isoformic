.onLoad <- function(...) {
  S7::S7_on_load()
}

.onUnload <- function(libpath) {
  stop_workers()
  S7::S7_on_unload()
}

# Removes the placeholders that `S7::method<-` leaves in the namespace for
# generics of other packages (print, dim, dimnames, `[`). Must stay at the top
# level of the last collated file, after every method registration.
S7::S7_on_build()

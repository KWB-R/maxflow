

vonNeumann_krit <- function (transmiss = 20*2E-5,
                                  s = 1E-4,
                                  dt = 24*3600, 
                                  dx = 10, 
                                  dy = dx) {
  

  transmiss/s * ((dt/dx^2) + (dt/dy^2)) 
}
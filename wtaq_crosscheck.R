devtools::install_github(repo = "kwb-r/kwb.wtaq",ref = "v.0.2.1")

library(kwb.wtaq)

kf <- 2E-5*24*3600


generalConfiguration <- wtConfigureGeneral(
  ### title of the project (max. length 70 characters)
  title="Example well field, long-term pumping test of well 5"
)

kf <- 2E-5*24*3600 #m/day

aquiferConfiguration <- wtConfigureAquifer(    
  aqtype = "CONFINED", # aquifer type
  bb = 80,                # saturated aquifer thickness
  hkr = kf,            # horizontal hydraulic conductivity
  hkz = kf,          # vertical hydraulic conductivity
  ss = 1E-04,             # specific storage
  sy = 0.123              # specific yield
)

drainageConfiguration <- wtConfigureDrainage(
  idra = 0 # = instantaneous drainage in unsaturated zone
)


pumpwellConfiguration <- wtConfigurePumpwell(
  ### partially penetrating pumped well
  ipws = 0,
  ### finite diameter well
  ipwd = 1, 
  ### pumping rate of production well in (here: m3/day)
  qq = 3500, 
  ### radius of pumped well-screen (here: meter) 
  rw = 0.5, 
  ### top of filter screen below initial water table (here: meter)
  zpd = 60, 
  ### bottom of filter screen below initial water table (here: meter)
  zpl = 80, 
  ### well-bore skin parameter (dimensionless)
  sw = 0
)

observationWell1 <- wtConfigureObservationWell(
  ### name of observation well
  obname = "OW1", 
  ### distance from pumping well (here: meters)
  r = 25, 
  ### partially penetrating observation well
  iows = 0, 
  ### delayed response
  idpr = 1, 
  ### top of filter screen below initial water table (here: meters)
  z1 = 60, 
  ### bottom of filter screen below initial water table (here: meters)
  z2 = 80, 
  ### inside radius of the observation well (here: meters)
  rp = 0.2
)


timesConfiguration <- wtConfigureTimes(its = 0,tlast = 100000,nlc = 6,nox = 10)



wtaqConfiguration <- wtConfigure(
  general = generalConfiguration,
  aquifer = aquiferConfiguration, 
  drainage = drainageConfiguration, 
  times = timesConfiguration, 
  solution = wtConfigureSolution(),
  pumpwell = pumpwellConfiguration,
  obswells = list(observationWell1)
)

res <- wtRunConfiguration(configuration = wtaqConfiguration)
wtPlotResult(res,plottype = "s")

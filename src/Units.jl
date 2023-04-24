module Units


const GeV = 1.0
const eV = 1e-9 * GeV
const keV = 1e-6 * GeV
const MeV = 1e-3 * GeV
const TeV = 1e3 * GeV
const PeV = 1e6 * GeV

const GV = GeV
const V = 1e-9 * GV
const kV = 1e-6 * GV
const MV = 1e-3 * GV
const TV = 1e3 * GV
const PV = 1e6 * GV

const fm = 1.0 / (197.3269804 * MeV)
const mm = 1e12 * fm
const cm = 10.0 * mm
const cm2 = cm * cm
const cm3 = cm2 * cm
const m = 1e3 * mm
const km = 1e3 * m
const pc = 3.08567758149e16 * m
const kpc = 1e3 * pc
const Mpc = 1e6 * pc

const barn = 1e-24 * cm2
const pb = 1e-36 * cm2

const sec = 299792458 * m
const minute = 60 * sec
const hour = 60 * minute

const day = 24 * hour
# mean sidereal day
# const day = 23 * hour + 56 * minute + 4.09053 * sec
# const year = 31558149.8 * sec  # sidereal year
const year = 31556925.1 * sec  # tropical year
const yr = year


# Pl constant, reduced (h bar) = 1.0545718e-34 J sec
const kg = sec / (1.0545718e-34 * m * m)
const tonne = 1000 * kg
const gram = 1e-3 * kg
const g_cm3 = gram / cm3

"newton"
const N = kg * m / sec / sec
"joule"
const J = kg * m * m / sec / sec
"farad"
const F = m / 8.8541878128e-12
"ampere"
const A = sqrt(1.00000000055 * 4π * 1e-7 * N)
"coulomb"
const C = A * sec


end  # module Units

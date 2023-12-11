"""
Physical constants and unit conversion factors.

"""

#: Speed of light [cm/s].
const SPEED_OF_LIGHT = 2.99792458e10

#: Most probable velocity of halo DM
const V0 = 220e5 / SPEED_OF_LIGHT
# V0 = 300e5 / SPEED_OF_LIGHT
# V0 = 235e5 / SPEED_OF_LIGHT

#: Escape velocity of the Milky Way.
const V_ESC = 544e5 / SPEED_OF_LIGHT
# V_ESC = 550e5 / SPEED_OF_LIGHT

#: Velocity of the earth
const V_EARTH = 240e5 / SPEED_OF_LIGHT

const ELECTRON_MASS = 0.51099895e-3 # GeV
const ATOMIC_MASS = 0.93149410242  # GeV
const PROTON_MASS = 0.93827208816  # GeV
const HELUIM_4_MASS = 3.7284021964542418  # GeV

const NUCLEUS_NAME = Dict(
    "H" => "Hydrogen", "He" => "Helium", "Li" => "Lithium", "Be" => "Beryllium",
    "B" => "Boron", "C" => "Carbon", "N" => "Nitrogen", "O" => "Oxygen", "F" => "Fluorine",
    "Ne" => "Neon", "Na" => "Sodium", "Mg" => "Magnesium", "Al" => "Aluminium",
    "Si" => "Silicon", "P" => "Phosphorus", "S" => "Sulphur", "Cl" => "Chlorine",
    "Ar" => "Argon", "K" => "Potassium", "Ca" => "Calcium", "Sc" => "Scandium",
    "Ti" => "Titanium", "V" => "Vanadium", "Cr" => "Chromium", "Mn" => "Manganese",
    "Fe" => "Iron", "Co" => "Cobalt", "Ni" => "Nickel",
)

const NUCLEUS_GROUP = Dict(
    "L" => ("Li", "Be", "B"),
    "M" => ("C", "N", "O", "F"),
    "H" => ("Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"),
    "VH" => ("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"),
)

const NUCLEUS_NAME_REVERSE = Dict(name => brev for (brev, name) in NUCLEUS_NAME)

const NUCLEUS_Z = Dict(
    "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8,
    "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15,
    "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22,
    "V" => 23, "Cr" => 24, "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29,
    "Zn" => 30, "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36,
    "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42, "Tc" => 43,
    "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50,
    "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56, "La" => 57,
    "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64,
    "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70, "Lu" => 71,
    "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78,
    "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83, "Po" => 84, "At" => 85,
    "Rn" => 86, "Fr" => 87, "Ra" => 88, "Ac" => 89, "Th" => 90, "Pa" => 91, "U" => 92,
)

#: Nucleus number.
const NUCLEUS_A = Dict(
    "Electron" => 1, "Proton" => 1, "He3" => 3, "He4" => 4, "H" => 1, "He" => 4, "Li" => 7,
    "Be" => 9, "B" => 11, "C" => 12, "N" => 14, "O" => 16, "F" => 19, "Ne" => 20,
    "Na" => 23, "Mg" => 24, "Al" => 27, "Si" => 28, "P" => 31, "S" => 32, "Cl" => 35,
    "Ar" => 40, "K" => 39, "Ca" => 40, "Sc" => 45, "Ti" => 47, "V" => 51, "Cr" => 52,
    "Mn" => 55, "Fe" => 56, "Co" => 59, "Ni" => 59, "Cu" => 64, "Zn" => 65, "Ga" => 70,
    "Ge" => 72, "As" => 75, "Se" => 79, "Br" => 80, "Kr" => 84, "Rb" => 85, "Sr" => 88,
    "Y" => 89, "Zr" => 91, "Nb" => 93, "Mo" => 96, "Tc" => 98, "Ru" => 101, "Rh" => 103,
    "Pd" => 106, "Ag" => 108, "Cd" => 112, "In" => 115, "Sn" => 119, "Sb" => 122,
    "Te" => 128, "I" => 127, "Xe" => 131, "Cs" => 133, "Ba" => 137, "La" => 139,
    "Ce" => 140, "Pr" => 141, "Nd" => 144, "Pm" => 145, "Sm" => 150, "Eu" => 152,
    "Gd" => 157, "Tb" => 159, "Dy" => 163, "Ho" => 165, "Er" => 167, "Tm" => 169,
    "Yb" => 173, "Lu" => 175, "Hf" => 178, "Ta" => 181, "W" => 184, "Re" => 186,
    "Os" => 190, "Ir" => 192, "Pt" => 195, "Au" => 197, "Hg" => 201, "Tl" => 204,
    "Pb" => 207, "Bi" => 209, "Po" => 209, "At" => 210, "Rn" => 222, "Fr" => 223,
    "Ra" => 226, "Ac" => 227, "Th" => 232, "Pa" => 231, "U" => 238, "Xe131" => 131,
    "Xe129" => 129,
)

const Z_NUCLEUS = Dict(Z => nuc for (nuc, Z) in NUCLEUS_Z)

#: Nucleus mass [GeV].
const NUCLEUS_MASS = Dict(
    "Proton" => PROTON_MASS,
    "He3" => 3.0160293 * ATOMIC_MASS,
    "He4" => HELUIM_4_MASS,
    "Electron" => ELECTRON_MASS,
    "Xe131" => 130.9050824 * ATOMIC_MASS,
    "Xe129" => 128.9047794 * ATOMIC_MASS,
)

for nuc in keys(NUCLEUS_Z)
    NUCLEUS_MASS[nuc] = NUCLEUS_A[nuc] * ATOMIC_MASS
end

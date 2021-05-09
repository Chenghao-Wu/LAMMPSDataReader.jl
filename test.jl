include("/Users/bruce/Dropbox/code/research/LAMMPSDataReader/LAMMPSDataReader/LAMMPSDataReader.jl")
using .LAMMPSDataReader

uni=Universe("input2.data")
bonds=LAMMPSDataReader.velocities(uni)
println(bonds)
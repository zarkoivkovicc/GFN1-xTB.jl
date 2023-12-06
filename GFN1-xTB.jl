#include("src/io.jl")
include("src/main.jl")
using .GFN1_xTB
using ArgParse

function parse_commandline()
    s = ArgParseSettings(prog="GFN1-xTB Method",
    description="""This program performs a single point 
    tight binding calculation using GFN1-xTB method developed in Grimme's group.""",
    epilog="""Author: Zarko Ivkovic \n \n
    Made for the remote computational project of EMTCCM master.""")

    @add_arg_table s begin
        "--xyz"
            help = "Specify the xyz file to read geometry from."
        "--verbose"
            help = "Show more intermediate quantities in the output file."
            action = :store_true
        "name"
            help = """Name of the calculation."""
            required = true
        "--parameters", "-p"
            help = "Which parameters should be used"
            default = "parameters_best.dat"
        "--maxiter"
            help= "Maximum number of iterations in SCF cycle. Default: 50"
            default=50
            arg_type = Int64
        "--charge", "-c"
            help= "Specify charge of the molecule. Default: 0"
            default = 0
            arg_type = Int64
        "--damping"
            help= "Damping factor. Defualt: 0.4"
            default = 0.4
            arg_type = Float64
        "--tolerance"
            help= "Convergence criteria of energy for SCF cycles. Default: 1e-7
            The threshold for change in the norm of the density matrix is 1e3 times the threshold for change in energy."
            default = 1e-7
            arg_type = Float64
    end

    return parse_args(s)
end
cmd_args = parse_commandline()
if isnothing(cmd_args["xyz"])
    name = cmd_args["name"]
    cmd_args["xyz"] = "$name.xyz"
elseif !endswith(cmd_args["xyz"],".xyz")
    cmd_args["xyz"] = cmd_args["xyz"] * ".xyz"
end
if !endswith(cmd_args["parameters"],".dat")
    cmd_args["parameters"] = cmd_args["parameters"] * ".dat"
end
name = cmd_args["name"]; xyz_file = cmd_args["xyz"]; verbosity = cmd_args["verbose"];
charge = cmd_args["charge"]; maxiter = cmd_args["maxiter"]; ΔE_min = cmd_args["tolerance"]
λ_damp = cmd_args["damping"]; params = cmd_args["parameters"]
out::String = "$name.out"
rm(out, force=true)
touch(out)
main(out,xyz_file,verbosity,params,maxiter,λ_damp,charge,ΔE_min)

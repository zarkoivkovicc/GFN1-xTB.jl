include("../src/main.jl")
using .GFN1_xTB
using Test
function testscf(scf_block1::Array{String},scf_block2::Array{String})
    for (a,b) in zip(split.(scf_block1,"    "),split.(scf_block2,"    "))
        if a[begin:end-2] != b[begin:end-2]
            return false
        end
    end
    return true
end

@testset "water" begin
    rm("water.out",force=true)
    main("water.out","$binpath/test/xyz/water.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("water.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/water.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:23],out2[16:23]) #SCF
    @test out1[50:end] == out2[50:end] #Energies in charges
    rm("water.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","$binpath/test/xyz/methanol.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:24],out2[16:24]) #SCF
    @test out1[75:end] == out2[75:end] #Energies in charges
    rm("methanol.out")
end
@testset "glucose" begin
    rm("glucose.out",force=true)
    main("glucose.out","$binpath/test/xyz/glucose.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("glucose.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/glucose.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:25],out2[16:25]) #SCF
    @test out1[244:end] == out2[244:end] #Energies in charges
    rm("glucose.out")
end
@testset "co" begin
    rm("co.out",force=true)
    main("co.out","$binpath/test/xyz/co.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("co.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/co.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:26],out2[16:26]) #SCF
    @test out1[53:end] == out2[53:end] #Energies in charges
    rm("co.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","$binpath/test/xyz/methanol.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:24],out2[16:24]) #SCF
    @test out1[75:end] == out2[75:end] #Energies in charges
    rm("methanol.out")
end
@testset "methane" begin
    rm("methane.out",force=true)
    main("methane.out","$binpath/test/xyz/methane.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("methane.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methane.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:19],out2[16:19]) #SCF
    @test out1[58:end] == out2[58:end] #Energies in charges
    rm("methane.out")
end
@testset "secbutylamine" begin
    rm("secbutylamine.out",force=true)
    main("secbutylamine.out","$binpath/test/xyz/secbutylamine.xyz",false,"parameters_best.dat",50,0.4,0,1e-7,1e-3)
    open("secbutylamine.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/secbutylamine.out","r") do f2
        global out2 = readlines(f2)
    end
    @test testscf(out1[16:23],out2[16:23]) #SCF
    @test out1[152:end] == out2[152:end] #Energies in charges
    rm("secbutylamine.out")
end

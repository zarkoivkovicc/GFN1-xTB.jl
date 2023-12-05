include("../src/main.jl")
using .GFN1_xTB
using Test
using Mmap
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
    main("water.out","xyz/water.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("water.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/water.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:52] == out2[17:52] #Initial qunatities
    @test testscf(out1[56:63],out2[56:63]) #SCF
    @test out1[64:end] == out2[64:end] #Final quantities
    rm("water.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","xyz/methanol.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:85] == out2[17:85] #Initial qunatities
    @test testscf(out1[86:94],out2[86:94]) #SCF
    @test out1[95:end] == out2[95:end] #Final quantities
    rm("methanol.out")
end
@testset "glucose" begin
    rm("glucose.out",force=true)
    main("glucose.out","xyz/glucose.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("glucose.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/glucose.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:289] == out2[17:289] #Initial qunatities
    @test testscf(out1[290:299],out2[290:299]) #SCF
    @test out1[300:end] == out2[300:end] #Final quantities
    rm("glucose.out")
end
@testset "co" begin
    rm("co.out",force=true)
    main("co.out","xyz/co.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("co.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/co.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:53] == out2[17:53] #Initial qunatities
    @test testscf(out1[54:65],out2[54:65]) #SCF
    @test out1[66:end] == out2[66:end] #Final quantities
    rm("co.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","xyz/methanol.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:85] == out2[17:85] #Initial qunatities
    @test testscf(out1[86:94],out2[86:94]) #SCF
    @test out1[95:end] == out2[95:end] #Final quantities
    rm("methanol.out")
end
@testset "methane" begin
    rm("methane.out",force=true)
    main("methane.out","xyz/methane.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methane.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/methane.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:71] == out2[17:71] #Initial qunatities
    @test testscf(out1[72:75],out2[72:75]) #SCF
    @test out1[76:end] == out2[76:end] #Final quantities
    rm("methane.out")
end
@testset "secbutylamine" begin
    rm("secbutylamine.out",force=true)
    main("secbutylamine.out","xyz/secbutylamine.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("secbutylamine.out","r") do f1
        global out1 = readlines(f1)
    end
    open("reference/secbutylamine.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[17:183] == out2[17:183] #Initial qunatities
    @test testscf(out1[184:191],out2[184:191]) #SCF
    @test out1[192:end] == out2[192:end] #Final quantities
    rm("secbutylamine.out")
end

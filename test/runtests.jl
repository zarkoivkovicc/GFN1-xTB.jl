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
    @test out1[17:101] == out2[17:101] #Initial qunatities
    @test testscf(out1[102:109],out2[102:109]) #SCF
    @test out1[110:end] == out2[110:end] #Final quantities
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
    @test out1[17:163] == out2[17:163] #Initial qunatities
    @test testscf(out1[164:172],out2[164:172]) #SCF
    @test out1[173:end] == out2[173:end] #Final quantities
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
    @test out1[17:591] == out2[17:591] #Initial qunatities
    @test testscf(out1[592:601],out2[592:601]) #SCF
    @test out1[602:end] == out2[602:end] #Final quantities
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
    @test out1[17:99] == out2[17:99] #Initial qunatities
    @test testscf(out1[100:110],out2[100:110]) #SCF
    @test out1[111:end] == out2[111:end] #Final quantities
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
    @test out1[17:163] == out2[17:163] #Initial qunatities
    @test testscf(out1[164:172],out2[164:172]) #SCF
    @test out1[173:end] == out2[173:end] #Final quantities
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
    @test out1[17:133] == out2[17:133] #Initial qunatities
    @test testscf(out1[134:137],out2[134:137]) #SCF
    @test out1[138:end] == out2[138:end] #Final quantities
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
    @test out1[17:365] == out2[17:365] #Initial qunatities
    @test testscf(out1[366:373],out2[366:373]) #SCF
    @test out1[374:end] == out2[374:end] #Final quantities
    rm("secbutylamine.out")
end

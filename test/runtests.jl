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
function testorbitals(orbitals_block1::Array{String}, orbitals_block2::Array{String})
    for (c,d) in zip(orbitals_block1[begin+2:end:3],orbitals_block2[begin+2:end:3])
        a = abs.(tryparse.(Float64,split(strip(replace(c, r"\t{1,}|\s{2,}" => ' ')))))
        b = abs.(tryparse.(Float64,split(strip(replace(d, r"\t{1,}|\s{2,}" => ' ')))))
        if a!=b
            return false
        end
    end
    return true
end
@testset "water" begin
    rm("water.out",force=true)
    main("water.out","$binpath/test/xyz/water.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("water.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/water.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:57] == out2[13:57] #Initial quantities
    @test testorbitals(out1[58:81],out2[58:81]) # Initial orbitals
    @test out1[82:101] == out2[82:101] # Initial density matrix and charges
    @test testscf(out1[102:109],out2[102:109]) #SCF
    @test testorbitals(out1[111:134],out2[111:134]) # Final orbitals
    @test out1[134:end] == out2[134:end] #Final quantities
    rm("water.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","$binpath/test/xyz/methanol.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:87] == out2[13:87] #Initial quantities
    @test testorbitals(out1[88:135],out2[88:135]) # Initial orbitals
    @test out1[136:163] == out2[136:163] # Initial density matrix and charges
    @test testscf(out1[164:172],out2[164:172]) #SCF
    @test testorbitals(out1[174:221],out2[174:221]) # Final orbitals
    @test out1[222:end] == out2[222:end] #Final quantities
    rm("methanol.out")
end
@testset "glucose" begin
    rm("glucose.out",force=true)
    main("glucose.out","$binpath/test/xyz/glucose.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("glucose.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/glucose.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:291] == out2[13:291] #Initial quantities
    @test testorbitals(out1[292:507],out2[292:507]) # Initial orbitals
    @test out1[508:591] == out2[508:591] # Initial density matrix and charges
    @test testscf(out1[592:601],out2[592:601]) #SCF
    @test testorbitals(out1[603:818],out2[603:818]) # Final orbitals
    @test out1[819:end] == out2[819:end] #Final quantities
    rm("glucose.out")
end
@testset "co" begin
    rm("co.out",force=true)
    main("co.out","$binpath/test/xyz/co.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("co.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/co.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:55] == out2[13:55] #Initial quantities
    @test testorbitals(out1[56:79],out2[56:79]) # Initial orbitals
    @test out1[80:99] == out2[80:99] # Initial density matrix and charges
    @test testscf(out1[100:110],out2[100:110]) #SCF
    @test testorbitals(out1[112:135],out2[112:135]) # Final orbitals
    @test out1[136:end] == out2[136:end] #Final quantities
    rm("co.out")
end
@testset "methanol" begin
    rm("methanol.out",force=true)
    main("methanol.out","$binpath/test/xyz/methanol.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methanol.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methanol.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:87] == out2[13:87] #Initial quantities
    @test testorbitals(out1[88:135],out2[88:135]) # Initial orbitals
    @test out1[136:163] == out2[136:163] # Initial density matrix and charges
    @test testscf(out1[164:172],out2[164:172]) #SCF
    @test testorbitals(out1[174:221],out2[174:221]) # Final orbitals
    @test out1[222:end] == out2[222:end] #Final quantities
    rm("methanol.out")
end
@testset "methane" begin
    rm("methane.out",force=true)
    main("methane.out","$binpath/test/xyz/methane.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("methane.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/methane.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:73] == out2[13:73] #Initial quantities
    @test testorbitals(out1[74:109],out2[74:109]) # Initial orbitals
    @test out1[110:133] == out2[110:133] # Initial density matrix and charges
    @test testscf(out1[134:137],out2[134:137]) #SCF
    @test testorbitals(out1[139:174],out2[139:174]) # Final orbitals
    @test out1[175:end] == out2[175:end] #Final quantities
    rm("methane.out")
end
@testset "secbutylamine" begin
    rm("secbutylamine.out",force=true)
    main("secbutylamine.out","$binpath/test/xyz/secbutylamine.xyz",true,"parameters_best.dat",50,0.4,0,1e-7)
    open("secbutylamine.out","r") do f1
        global out1 = readlines(f1)
    end
    open("$binpath/test/reference/secbutylamine.out","r") do f2
        global out2 = readlines(f2)
    end
    @test out1[10] == out2[10] #E_rep
    @test out1[12] == out2[12] #E_disp
    @test out1[13:185] == out2[13:185] #Initial quantities
    @test testorbitals(out1[186:311],out2[186:311]) # Initial orbitals
    @test out1[312:365] == out2[312:365] # Initial density matrix and charges
    @test testscf(out1[366:373],out2[366:373]) #SCF
    @test testorbitals(out1[375:500],out2[375:500]) # Final orbitals
    @test out1[501:end] == out2[501:end] #Final quantities
    rm("secbutylamine.out")
end

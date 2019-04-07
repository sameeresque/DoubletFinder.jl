
using Test

@testset "linefinder" begin
test_array=[2.1,3.0,0.8,0.7,0.5,0.4,0.2,1.01,1.33,3.88,0.99,0.96,0.98,1.23,4.55]
@testset "Testing linear algebra results" begin
   @test 1 == 1
end;

@testset "findV" begin
    @test findV(3.12,3.12) == 0.0
    @test findV(3.0,3.12) ≈ 8858.92604 atol=1e-4
end;
    
@testset "findZAbs" begin
    @test findZAbs(0.0, 3.12) == 3.12
    @test findZAbs(1.0, 3.12) ≈ 3.119986 atol=1e-4
end;

@testset "getblocks1" begin
    @test getblocks1(test_array)[1] == (3, 7)
    @test getblocks1(test_array)[2] == (11, 13)
end;        
@testset "iscontained" begin
    @test iscontained(test_array,66.0) == false
    @test iscontained(test_array,3.25) == true            
end;
@testset "Regression_test" begin
    @test possible_doublets(df_clean,zEm,blocks_3s,species_dict,Dict("SiIV" => (-70000,10000)))==possible_doublets_parallel(df_clean,zEm,blocks_3s,species_dict,
        Dict("SiIV" => (-70000,10000))) 
end;
    
end; # linefinder

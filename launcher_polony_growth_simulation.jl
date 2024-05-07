###################################################
# Set Up
###################################################

using Bio, BioAlignments, BioSequences
using Distributions, StatsBase, Random
using Plots
using StaticArrays
using Profile
using FLoops
using SparseArrays
using Base.Threads
using Base.Iterators: partition
using NearestNeighbors
using Distributed
using FileIO
using Images

ncpu = Sys.CPU_THREADS
ENV["GKS_ENCODING"]="utf8"
ENV["JULIA_NUM_THREADS"] = ncpu*2
ENV["GKSwstype"]="nul"

global lp_ssDNA = 4.0 
global lbp = 0.34 
global lkuhn = 2.0*lp_ssDNA 
global scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);


Random.seed!(1234);

###################################################
# Structs
###################################################

struct OutParameter
    out_path::String
    grow_main::Bool
    grow_cross::Bool
    FASTQ_file::Bool
    Polon::Bool
    General::Bool
    graph_table::Bool
    
    # Constructor with default values
    function OutParameter(;out_path="./PCR_Simulation_results",
                           grow_main= true,
                           grow_cross=true,
                           FASTQ_file= true,
                           Polon = true,
                           General= true,
                           graph_table=true)
        new(out_path,grow_main, grow_cross, FASTQ_file, Polon, General, graph_table)
    end
end


struct param
    Output_file_name::String
    site_density::Float64
    spot_diameter::Float64
    site_area::Float64
    length_thresh::Int64
    l_prime::Int64
    n_cycles::Int64
    n_cycles_2::Int64
    l_prime_2::Int64
    Temp::Int64
    Temp_2::Int64
    rejection_threshold::Int64
    proportions::Vector{Float64}
    restriction_seq:: String
    sequences::Vector{BioSequence{DNAAlphabet{4}}}
    primers::Vector{Tuple{Any,String,String}}
    cutoff::Int64

    function param(Output_file_name::String; 
                    output="./Result_PCR_Simulation", 
                    site_density=0.07, 
                    spot_diameter=1000.0, 
                    n_cycles=30, 
                    n_cycles_2=1, 
                    length_thresh=40, 
                    cutoff=10, 
                    l_prime=15, 
                    l_prime_2=13, 
                    Temp=50, 
                    Temp_2=40, 
                    rejection_thresh=10000, 
                    proportions=[1,1,1,1,0.034,0.034], 
                    restriction_seq = "GAATTC",
                    primers = [("YOUR PRIMER SEQUENCE", "#d43737", "Alpha"),
                                ("YOUR PRIMER SEQUENCE", "#cf8c8c", "Alpha*"),
                                ("YOUR PRIMER SEQUENCE", "#26c766", "Beta"),
                                ("YOUR PRIMER SEQUENCEG", "#92d1ab", "Beta*"),
                                (400, "#000000", "Crossover")],
                    strand_sequence=[BioSequence{DNAAlphabet{4}}("YOUR PRIMER SEQUENCE"),
                                    BioSequence{DNAAlphabet{4}}("YOUR PRIMER SEQUENCE"),
                                    BioSequence{DNAAlphabet{4}}("YOUR PRIMER SEQUENCE"),
                                    BioSequence{DNAAlphabet{4}}("YOUR PRIMER SEQUENCE"),
                                    BioSequence{DNAAlphabet{4}}("YOUR STRAND SEQUENCE"),
                                    BioSequence{DNAAlphabet{4}}("YOUR STRAND SEQUENCE")])
        area = π * (spot_diameter / 2)^2
        if !isdir(output)
            mkdir(output)
        end
        Output_place = joinpath(output, Output_file_name)
        new(Output_place, site_density, spot_diameter, area, length_thresh, l_prime, n_cycles, n_cycles_2, l_prime_2, Temp, Temp_2, rejection_thresh, proportions, restriction_seq, strand_sequence, primers,cutoff)
    end
end

###################################################
# Main
###################################################


include("fun_bridge_PCR_launcher.jl")

out_param = OutParameter(out_path="./", grow_cross= true, FASTQ_file=true, Polon=true, General=true)


# Define your parameters
primer_densities =  ["PRIMER SEQUENCE AS INTEGER"] 
length_connections = ["CUTOFF AS INTEGER"]
proportion_templates = [[1,1,1,1,"TEMPLATE PROPORTION AS INTEGER","TEMPLATE PROPORTION AS INTEGER"] ]
crossover_scales = ["CUTOFF AS INTEGER"]
cycles = ["CYCLE NUMBER AS INTEGER"]

# Run Simulation
for (i, primer_density) in enumerate(primer_densities)
    for (k, cycle) in enumerate(cycles)
        for (j, length_connection) in enumerate(length_connections)
            for (m, proportion_template) in enumerate(proportion_templates)
                for (n, crossover_scale) in enumerate(crossover_scales)
                    output_filename = "PCR_SIM_new_$(primer_density)_$(proportion_template[5])_$(length_connection)_$(cycle)_$(crossover_scale)_"
                    parameter = param(output_filename, output = out_param.out_path,site_density = primer_density, cutoff=length_connection, n_cycles= cycle,n_cycles_2 = crossover_scale, proportions=proportion_template)
                try
                    result = @time PCR_via_bridge(out_param, parameter)
                catch e
                    println("Error in:", output_filename)
                    println("Error messgae:", e)
                    continue
                end
                end
            end 
        end
    end
end
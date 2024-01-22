using Bio, BioAlignments, BioSequences
using Distributions, StatsBase, Random
using Plots
using StaticArrays
using FLoops

# cd("H:/Projects_shared/Bridge_PCR_sim")
# include("bridge_PCR_launcher.jl")#
# Pkg.add.(["Bio","BioAlignments","BioSequences","Distributions","StatsBase","Random","Plots","StaticArrays","FLoops"])


###################################################
# low level function library
###################################################
struct Random_poisson_parameters # to hold parameters related to a Poisson random process.
    density::Float64 # points per square nm
    simulation_area::Float64 # a 2D disk
end

mutable struct Site
    interact_params::Vector{Float64} #parameters for the interaction probability distribution
    free_neighbors::Vector{Tuple{Tuple{Float64, Float64}, Float64,Float64}} #Vector{Neighbor_set}
    seq::BioSequence{DNAAlphabet{4}}
    Site() = new() #create an instance
end

mutable struct Uniform_bipartite_seed #represents a seed for a uniform bipartite system.
    proportions::Vector{Float64}
    sequences::Vector{BioSequence{DNAAlphabet{4}}}
    Uniform_bipartite_seed() = new() #create an instance
end

mutable struct Strand_tracker
    sequence_tag::String
    xs::Vector{Any}
    ys::Vector{Any}
    color::String
    length::Int64
    series_name::String
    Strand_tracker() = new() #create an instance
end










#using two strand_tracker functions to be flexible with the nature of the inputdata

function Strand_tracker(sequence_tag::String, color::String,series_name::String)
    out = Strand_tracker()
    out.sequence_tag = sequence_tag
    out.xs = []
    out.ys = []
    out.color = color
    out.series_name = series_name
    return out
end

function Strand_tracker(length::Int64, color::String,series_name::String)
    out = Strand_tracker()
    out.length = length
    out.xs = []
    out.ys = []
    out.color = color
    out.series_name = series_name
    return out
end

function update_strand_tracker!(tracker::Strand_tracker, coordinates, query_sequence::String)
    if occursin(tracker.sequence_tag, query_sequence)
        push!(tracker.xs, coordinates[1])
        push!(tracker.ys, coordinates[2])
    end
    return tracker
end

function update_strand_tracker!(tracker::Strand_tracker, coordinates, query_length::Int64)
    if query_length > tracker.length
        push!(tracker.xs, coordinates[1])
        push!(tracker.ys, coordinates[2])
    end
    return tracker
end

function output_tracking_plot(trackers::Vector{Strand_tracker}, site_positions, plot_name::String)
    p = scatter([i[1] for i in site_positions],[i[2] for i in site_positions], color="#E4E4E4",alpha=0.5,markerstrokewidth=0,markersize=1,label = "primer sites")
    for tracker in trackers
        scatter!(tracker.xs, tracker.ys, color=tracker.color,alpha=1,markerstrokewidth=0,markersize=2,label = tracker.series_name)
    end
    plot!(size=(400,400))
    savefig(plot_name)
end

function renormalize_site_probabilities(sites::Dict{Tuple{Float64, Float64}, Site})
    @floop for site in sites
        if isdefined(site[2],2)
            normalization_factor = 0
            for i in 1:length(site[2].free_neighbors)
                normalization_factor += site[2].free_neighbors[i][3]
            end
            for i in 1:length(site[2].free_neighbors)
                site[2].free_neighbors[i] = (site[2].free_neighbors[i][1],site[2].free_neighbors[i][2],site[2].free_neighbors[i][3]/normalization_factor)
            end
        end
    end
    return sites
end

function update_site(site::Site)
    #println("old site",site)
    normalization_factor = 0
    L_MS = length(site.seq)
    for i in 1:length(site.free_neighbors)
        distance = site.free_neighbors[i][2]
        p_interact_MS = p_end_to_end(L_MS, distance)
        site.free_neighbors[i] = (site.free_neighbors[i][1],site.free_neighbors[i][2],p_interact_MS)
        normalization_factor += p_interact_MS
    end
    for i in 1:length(site.free_neighbors)
        site.free_neighbors[i] = (site.free_neighbors[i][1],site.free_neighbors[i][2],site.free_neighbors[i][3]/normalization_factor)
    end
    #println("new site",site)
    return site
end

function bc_gen(seq)
    """
    Given a sequence, substitute "N" for random nucleotides
    Write the unique barcode "NNNN..."
    """
    barcode = BioSequence{DNAAlphabet{4}}()            # Full sequence
    unique_barcode = BioSequence{DNAAlphabet{4}}()     # Random domain NNNN "unique barcode"
    lookup = Dict{BioSequence{DNAAlphabet{4}},Array{BioSequence{DNAAlphabet{4}},1}}()
    lookup[dna"N"] = @SVector [dna"T",
                        dna"A",
                        dna"G",
                        dna"C"]

    lookup[dna"W"] = @SVector [dna"T",dna"A"]
    lookup[dna"S"] = @SVector [dna"C",dna"G"]
    lookup[dna"T"] = @SVector [dna"T"]
    lookup[dna"A"] = @SVector [dna"A"]
    lookup[dna"C"] = @SVector [dna"C"]
    lookup[dna"G"] = @SVector [dna"G"]

    for i in 1:length(seq)
        # Store random domain in "unique barcode"
        if seq[i:i] == dna"N"
            unique_barcode_letter = rand(lookup[seq[i:i]])
            unique_barcode = unique_barcode * unique_barcode_letter
        else 
            unique_barcode_letter = rand(lookup[seq[i:i]])
        end
        barcode = barcode*unique_barcode_letter
    end
    return barcode
end

function get_p_annealing(Temp, T_m; 𝜅 = 40.) # 𝜅 = 40.  the "melting constant" the higher the number, the sharper the transition
    p_anneal = 1 - 1/(1 + exp(-1*𝜅*(Temp*0.01 - 0.01*T_m)))
    return p_anneal
end

function get_thermodynamic_parameters()
    global ΔH° = Dict(
        BioSequence{DNAAlphabet{4}}("AA") => -33.1,BioSequence{DNAAlphabet{4}}("TT") => -33.1,
        BioSequence{DNAAlphabet{4}}("AT") => -30.1,BioSequence{DNAAlphabet{4}}("TA") => -30.1,
        BioSequence{DNAAlphabet{4}}("TA") => -30.1,BioSequence{DNAAlphabet{4}}("AT") => -30.1,
        BioSequence{DNAAlphabet{4}}("CA") => -35.6,BioSequence{DNAAlphabet{4}}("GT") => -35.6,
        BioSequence{DNAAlphabet{4}}("AC") => -35.1,BioSequence{DNAAlphabet{4}}("TG") => -35.1,
        BioSequence{DNAAlphabet{4}}("CT") => -32.6,BioSequence{DNAAlphabet{4}}("GA") => -32.6,
        BioSequence{DNAAlphabet{4}}("TC") => -34.3,BioSequence{DNAAlphabet{4}}("AG") => -34.3,
        BioSequence{DNAAlphabet{4}}("CG") => -44.4,BioSequence{DNAAlphabet{4}}("GC") => -44.4,
        BioSequence{DNAAlphabet{4}}("GC") => -41.0,BioSequence{DNAAlphabet{4}}("CG") => -41.0,
        BioSequence{DNAAlphabet{4}}("GG") => -33.5,BioSequence{DNAAlphabet{4}}("CC") => -33.5,
        BioSequence{DNAAlphabet{4}}("A") => 9.6,BioSequence{DNAAlphabet{4}}("T") => 9.6,
        BioSequence{DNAAlphabet{4}}("G") => 0.4,BioSequence{DNAAlphabet{4}}("C") => 0.4,
        )

    global ΔS° = Dict(
        BioSequence{DNAAlphabet{4}}("AA") => -92.9,BioSequence{DNAAlphabet{4}}("TT") => -92.9,
        BioSequence{DNAAlphabet{4}}("AT") => -85.4,BioSequence{DNAAlphabet{4}}("TA") => -85.4,
        BioSequence{DNAAlphabet{4}}("TA") => -89.1,BioSequence{DNAAlphabet{4}}("AT") => -89.1,
        BioSequence{DNAAlphabet{4}}("CA") => -95.0,BioSequence{DNAAlphabet{4}}("GT") => -95.0,
        BioSequence{DNAAlphabet{4}}("AC") => -93.7,BioSequence{DNAAlphabet{4}}("TG") => -93.7,
        BioSequence{DNAAlphabet{4}}("CT") => -87.9,BioSequence{DNAAlphabet{4}}("GA") => -87.9,
        BioSequence{DNAAlphabet{4}}("TC") => -92.9,BioSequence{DNAAlphabet{4}}("AG") => -92.9,
        BioSequence{DNAAlphabet{4}}("CG") => -113.8,BioSequence{DNAAlphabet{4}}("GC") => -113.8,
        BioSequence{DNAAlphabet{4}}("GG") => -83.3,BioSequence{DNAAlphabet{4}}("CC") => -83.3,
        BioSequence{DNAAlphabet{4}}("GC") => -102.1,BioSequence{DNAAlphabet{4}}("CG") => -102.1,
        BioSequence{DNAAlphabet{4}}("A")  => 17.2,BioSequence{DNAAlphabet{4}}("T") => 17.2,
        BioSequence{DNAAlphabet{4}}("G")  => -11.7,BioSequence{DNAAlphabet{4}}("C") => -11.7,
        )
    return ΔH°, ΔS°
end

function get_T_m(nearest_neighbor_sequence::BioSequence{DNAAlphabet{4}}; CNa = 0.5, ΔH°::Dict{DNASequence, Float64} = ΔH°, ΔS°::Dict{DNASequence, Float64} = ΔS°, cutoff_len = 18)
    ℛ = 8.31446261 # kJ⋅K−1⋅mol−1
    if length(nearest_neighbor_sequence)>cutoff_len #i.e. we only consider annealing of the last bases
        nearest_neighbor_sequence= nearest_neighbor_sequence[end-cutoff_len:end]
    end
    ΣΔH = 0
    ΣΔS = 0
    countr = 0
    for i in 1:length(nearest_neighbor_sequence)-1
        countr += 1
        ΣΔH += ΔH°[nearest_neighbor_sequence[i:i+1]]*1000
        ΣΔS += ΔS°[nearest_neighbor_sequence[i:i+1]]
        # print(nearest_neighbor_sequence[i:i+1], " ", ΔH°[nearest_neighbor_sequence[i:i+1]], " ", ΣΔH, "\n")
    end
    ΣΔH += ΔH°[nearest_neighbor_sequence[1:1]]*1000
    ΣΔH += ΔH°[nearest_neighbor_sequence[end:end]]*1000
    ΣΔS += ΔS°[nearest_neighbor_sequence[1:1]]
    ΣΔS += ΔS°[nearest_neighbor_sequence[end:end]]

    T_m = ΣΔH/(ΣΔS - ℛ*log(1000000)) + 16.6*log10(CNa) - 273
    # print(ΣΔH, " ", ΣΔS, "\n")
    # print("T_M = ", T_m)
    return T_m
end



"""
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.158105
"""
function p_end_to_end(L::Int64, r::Float64)
    s0 = L*lbp
    #p_end_to_end = (3/(2*π*L*lkuhn^2)^(3.0/2.0))*exp((-3.0*r)/(2*L*lkuhn^2))
    p_end_to_end = (r*exp(-1*s0/(8*lkuhn*(1-(r/s0)^2))))/((1-(r/s0)^2)*(2-(r/s0)^2)^2)
    #p_cumulative = 0
    #for r_i in 1:r
    #    p_cumulative += (r_i*exp(-1*s0/(8*lkuhn*(1-(r_i/s0)^2))))/((1-(r_i/s0)^2)*(2-(r_i/s0)^2)^2)
    #end
    if p_end_to_end < 0
        return 0
    else
        return p_end_to_end
    end
end


###################################################
# high level function library
###################################################
"""
Basic random Poisson distributed site generation
"""
function generate_site_positions(dist::Random_poisson_parameters)
    expected_N_points = dist.density*dist.simulation_area
    actual_N_points = rand(Poisson(expected_N_points))
    site_positions::Array{Tuple{Float64, Float64}} = []
    radius = sqrt(dist.simulation_area/π)
    println(expected_N_points)
    for i in 1:actual_N_points
        random_radius = radius*sqrt(rand())
        random_angle = rand()*2*π
        y = random_radius*sin(random_angle)
        x = random_radius*cos(random_angle)
        push!(site_positions, (x,y))
    end
    return site_positions
end

"""
takes site positions and creates Site objects and a dictionary of them
keys are the coordinates of each site
free_neighbors attribute is a list of coordinates(keys) and distance
threshold distance is used to reduce the number of neighbors
site information recorded in this step are the coordinates (also used as dict keys)
and the distance from site to neighbors for each neighbor
"""

function initialize_sites(site_positions::Array{Tuple{Float64, Float64}}, inputs::Uniform_bipartite_seed, length_thresh::Int64=40)
    total_input_normalization_factor = sum(inputs.proportions)
    weights = inputs.proportions/total_input_normalization_factor
    sequences = inputs.sequences
    N = length(site_positions)
    long_sites =  []
  
    initialized_sites = [Site() for i in 1:N]
    sites = Dict(zip(site_positions, initialized_sites))

    ### assign sequences to all sites in this loop###
    println("N is ", N)
    for i in 1:N   
        if i%1000 == 0
            print(i, " ")
        end
        main_site = sites[site_positions[i]]
        main_sequence = bc_gen(sample(sequences, Weights(weights)))
        main_site.seq = main_sequence
        
        L_MS = length(main_sequence)# length in bases 
        ls_MS = L_MS*lbp # countor length is bases*lbp
        if L_MS > length_thresh
            #println(main_site.seq)
            push!(long_sites, site_positions[i])
        end 
    end
    #distance_matrix = Array{Float64}(undef, N, N)
    println("\n assigning neighbors based on cutoff distance")
    ### check all pairwise distances and assign neighbors based on adaptive cutoff distance###
    for i in 1:N   
        if i%1000 == 0
            print(i)
        end
        x1 = site_positions[i][1]
        y1 = site_positions[i][2]
        main_site = sites[site_positions[i]]
        L_MS = length(main_site.seq)
        MS_cutoff =  136#L_MS*lbp # contour length is max - should be shorter though?
        
        @floop for j in i+1:N 
            candidate_neighbor_site = sites[site_positions[j]]
            L_CNS = length(candidate_neighbor_site.seq)
            CNS_cutoff = 136# L_CNS*lbp
            x2 = site_positions[j][1]
            y2 = site_positions[j][2]
            distance = sqrt((x1-x2)^2+(y1-y2)^2) 
            p_interact_MS = p_end_to_end(L_MS, distance)
            p_interact_CNS = p_end_to_end(L_CNS, distance)
            if MS_cutoff > distance
                # then this candidate and the main site are technically neighbors - some more likely than others though
                # since we are performing and update, we need to check first what the status of the free_neighbors lists are in both main and candidate sites
                if isdefined(sites[(x1,y1)], 2)
                    # main site already has a free_neighbors entry, but it cannot have the candidate in it yet
                    # so in this case we will push the candidate to the existing entry
                    push!(main_site.free_neighbors, ((x2,y2), distance, p_interact_MS))    
                else
                    # then this is the first iteration on the main site and it also has never been added as a neighbor by a previous main site
                    # in this case we create a new starting entry with the now-accepted neighbor candidate as the first of free_neighbors
                    main_site.free_neighbors = [((x2,y2), distance, p_interact_MS)]  
                end
            else
                # then the main site will not have this candidate as neighbor, nor will this candidate have main site as neighbor
            end
            if CNS_cutoff > distance
                # then the candidate is also eligible to have main site as one of its neighbors
                #now time to update the candidate as well
                if isdefined(sites[(x2,y2)], 2)
                    # then the candidate has an already-defined set of free_neighbors
                    # in this case we push the main site onto the existing list
                    push!(candidate_neighbor_site.free_neighbors, ((x1,y1), distance, p_interact_CNS))    
                else
                    # if not definied, then there is no free_neighbors entry for this candidate
                    # we make a new one now, assigning the main site as the first entry
                    candidate_neighbor_site.free_neighbors = [((x1,y1), distance, p_interact_CNS)] 
                end
            else
                # then the candidate is not eligible to have main site as its neighbor, regardless of whether IT may be a neighbor in the main site
            end
        end
        
    end
    return sites, long_sites#, distance_matrix
end






















###########################################
# run parameters
###########################################
site_density = 0.02546473*.01 #oligos per nm2
spot_diameter = 1000.0*10 #nm
site_area = π*(spot_diameter/2)^2
inputs = Uniform_bipartite_seed()
inputs.proportions = [10,10,10,10,.75,.75]
inputs.sequences = [
    BioSequence{DNAAlphabet{4}}("GGCCAACGGTGTCTCAATCAAC"),# strand a
    BioSequence{DNAAlphabet{4}}("GCTGCTAAGCCGGACTGAATTC"), #strand b
    BioSequence{DNAAlphabet{4}}("CGCACTGCCGCAGAATGAATTC"),#strand c
    BioSequence{DNAAlphabet{4}}("GCGGTTCCTGAACACGTTCGAA"),# strand d``
    BioSequence{DNAAlphabet{4}}("GGCCAACGGTGTCTCAATCAACCGATGGCCGAGCTCACTATGTAACTTTTGAACTAGGGCCTCTCACTGCCCTATGTACTGACATCCCTGCAAGTATTTCGATATATTCAGTGTAGCAGTGTCGGGGAGTTGCTCTACCCCGAGTTGTCGGCATTATAACGATCACGTGTGTCACGCAGTTTCATAGTTACGTGTGTGATCGCCACACATGCGTAGCTCCCCATGCGCGGACCTAACGCTGAAATATCGTATCTCAGCTTCAACAGACTGGATTCATAAGCAAATTGGCTAAACAGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNGTTGACAATACGCGAATTCAGTCCGGCTTAGCAGC"), #alpha strand #GGCCAACGGTGTCTCAATCAACCGATGGCCGAGCTCACTATGTAACTTTTGAACTAGGGCCTCTCACTGCCCTATGTACTGACATCCCTGCAAGTATTTCGATATATTCAGTGTAGCAGTGTCGGGGAGTTGCTCTACCCCGAGTTGTCGGCATTATAACGATCACGTGTGTCACGCAGTTTCATAGTTACGTGTGTGATCGCCACACATGCGTAGCTCCCCATGCGCGGACCTAACGCTGAAATATCGTATCTCAGCTTCAACAGACTGGATTCATAAGCAAATTGGCTAAACAGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNGTTGACAATACGCGAATTCAGTCCGGCTTAGCAGC
    BioSequence{DNAAlphabet{4}}("GCGGTTCCTGAACACGTTCGAAAAAAAAAAAAAAAAAAAAAATNNNNNNNNNNNNNNNNNNNNNNNNCGCGTATTGTCAACGGTCAATTGAGAGCCCCTGCACCGATCCGACGCATTTCGTTCAAGACATCCAGTTCGTAGTACCTCCCTGGAAATTCGGTGGCAGTTAAATCTAGATATCACATGGTTACAGGTCTCGTGAAGACCTGTACCCTGCCATCCGGAAAAGTTCATATGCCGTGGAAAATCTGGCACTGTTGGCAATTGTTCGACTTGCTATTCATAGTGACGTGCTGCATTTCAGATCAACTGCCAGGCAATGGGATACCCTTGGGATCAGGTAAGCAAGAATGCTATATCCGCAGGAGTAACCCTAGAGAATTCATTCTGCGGCAGTGCG") #beta strand #GCGGTTCCTGAACACGTTCGAAAAAAAAAAAAAAAAAAAAAATNNNNNNNNNNNNNNNNNNNNNNNNCGCGTATTGTCAACGGTCAATTGAGAGCCCCTGCACCGATCCGACGCATTTCGTTCAAGACATCCAGTTCGTAGTACCTCCCTGGAAATTCGGTGGCAGTTAAATCTAGATATCACATGGTTACAGGTCTCGTGAAGACCTGTACCCTGCCATCCGGAAAAGTTCATATGCCGTGGAAAATCTGGCACTGTTGGCAATTGTTCGACTTGCTATTCATAGTGACGTGCTGCATTTCAGATCAACTGCCAGGCAATGGGATACCCTTGGGATCAGGTAAGCAAGAATGCTATATCCGCAGGAGTAACCCTAGAGAATTCATTCTGCGGCAGTGCG
]
length_thresh = 40 #nt long to be considered a long strand
global l_prime = 15
n_cycles_2 = 5
global l_prime_2 = 13
Temp = 50
Temp_2 = 40
ncpu = Sys.CPU_THREADS
ENV["GKS_ENCODING"]="utf8"
ENV["JULIA_NUM_THREADS"] = ncpu*2
ENV["GKSwstype"]="nul"
n_cycles = 50

a_pol_f = Strand_tracker("TGCCCTATGTACTGACATCCCTGCAAGTATTTCG", "#d43737", "Alpha")
a_pol_r = Strand_tracker("CGAAATACTTGCAGGGATGTCAGTACATAGGGCA", "#cf8c8c", "Alpha*")
b_pol_f = Strand_tracker("CGCGTATTGTCAACGGTCAATTGAGAGCCCCTGC", "#26c766", "Beta")
b_pol_r = Strand_tracker("GCAGGGGCTCTCAATTGACCGTTGACAATACGCG", "#92d1ab", "Beta*")
tracker_400 = Strand_tracker(400, "#000000", "Crossover")

Random.seed!(1234);

rejection_thresh = 10000


open("run_info.txt", "w") do io
    write(io, "spot diameter: "*string(spot_diameter)*"nm \n")
    write(io, "spot area: "*string(site_area)*"nm2 \n")
    write(io, "site density: "*string(site_density)*"strands per nm2 \n")
    write(io, "input sequences used: "*string(inputs.sequences)*"\n")
    write(io, "input proportions of sequences: "*string(inputs.proportions)*"\n")
    
    end;








###########################################
# main launch script
###########################################
ΔH°, ΔS° = get_thermodynamic_parameters()
global lp_ssDNA = 4.0 #nm # persistence length of ssDNA  #Bernard, Tinland (1997). "Persistence Length of Single-Stranded DNA". Macromolecules. 30 (19): 5763. Bibcode:1997MaMol..30.5763T. doi:10.1021/ma970381+.
global lbp = 0.34 #nm # contour length of ssDNA by base 
global lkuhn = 2.0*lp_ssDNA # kuhn length is 2x lp for WLC model: Rubinstein and Colby (2003). Polymer Physics
global scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
println("generating site positions")
site_positions = generate_site_positions(Random_poisson_parameters(site_density,site_area)) 
println("initializing sites")
sites, long_sites = initialize_sites(site_positions, inputs, length_thresh)
println("renormalizing site probabilities")
sites = renormalize_site_probabilities(sites)
println("beginning thermal cycling")

global extension_vector = []
#######################################
#      ANNEAL w. EXTENSION PHASE      #
#######################################
extension_events = 0
for cycle in 1:n_cycles
    output_tracking_plot([a_pol_f,a_pol_r,b_pol_f,b_pol_r],site_positions,"test"*string(cycle+1000)*".png")

    println("start of cycle: ", cycle, "number of long sites: ", length(long_sites), "extensions: ", extension_events)
    push!(extension_vector, extension_events)
    global free_long_sites = deepcopy(long_sites)
    global i_no = 0
    global annealings = 0
    global reject_count = 0
    global extension_events_this_cycle = 0
    while length(free_long_sites) > 0 && reject_count < rejection_thresh
        shuffle!(free_long_sites)
        candidate = pop!(free_long_sites)
        global i_no += 1
        
        # now attempt a random interaction from the site_picks list
        if isdefined(sites[candidate],2)
            main_site = sites[candidate]
            weights = [j[3] for j in main_site.free_neighbors]
            choices = [j[1] for j in main_site.free_neighbors]
            partner_coordinates = sample(choices, Weights(weights))
            partner_site = sites[partner_coordinates]
            ############# pairwise local alignment check here! ###############
            seq_a = main_site.seq
            seq_b = partner_site.seq
            if occursin(string(seq_a[end-l_prime:end]), string(reverse_complement(seq_b))) || occursin(string(seq_b[end-l_prime:end]), string(reverse_complement(seq_a))) #because exact match much quicker to search than performing pairwise align at this stage
                full_align = pairalign(LocalAlignment(), seq_a, reverse_complement(seq_b), scoremodel)
                p_anneal = get_p_annealing(Temp, get_T_m(full_align.aln.b; CNa = .1));
                if rand() < p_anneal
                    global annealings += 1
                    # an annealing event is approved, now compute extension and update everything
                    #try
                        ### ANNEALING BLOCK ####
                        # now we need to check what type of extension should be generated given the alignment
                        seq_a_new = BioSequence{DNAAlphabet{4}}()*seq_b
                        seq_b_new = BioSequence{DNAAlphabet{4}}()*seq_b
                        len_a = length(seq_a) # 3' end of the original sequences
                        len_b = length(seq_b) # 3' end of the original sequences
                        ### fetch 3' ends of both strands partipating in the alginment ###
                        if full_align.aln.a.aln.anchors[end].op == OP_SEQ_MATCH && full_align.aln.a.aln.anchors[1].op == OP_START
                            aln_3prime_end_seq_a = full_align.aln.a.aln.anchors[end].seqpos # the end of the top strand is the 3' end of the top strand
                            aln_5prime_end_seq_b = 1 + len_b - full_align.aln.a.aln.anchors[end].refpos # the end of the bottom strand
                                ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###
                            if len_a == aln_3prime_end_seq_a && aln_5prime_end_seq_b > 1 # i.e. if the 3' end of the sequence is the same as 3' end of alignment    
                                # i.e. if the 5' end of the bottom alignment isn't also the 5' end
                                seqa_extend = true # now compute the extended portion
                                n_bases_2_fetch_seq_b = len_b - full_align.aln.a.aln.anchors[end].refpos # number 5' bases to be copied from the bottom strand seqb
                                extension_a = reverse_complement(seq_b[1:n_bases_2_fetch_seq_b])
                                seq_a_new = seq_a*extension_a # create a new sequence and update the molecule archive
                                #### MODIFY STRAND HERE WITH THE NEW SEQUENCE ####
                                sites[candidate].seq = seq_a_new
                                sites[candidate] = update_site(sites[candidate])
                                global extension_events += 1
                                global extension_events_this_cycle += 1
                                update_strand_tracker!(a_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(a_pol_r,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_r,candidate,string(seq_a_new))
                            else
                                #println("anneal with no a extension")
                                seqa_extend = false
                            end
                            aln_3prime_end_seq_b = len_b - full_align.aln.a.aln.anchors[1].refpos # the start base of the revcom is the 3' end, to convert t
                            aln_5prime_end_seq_a = 1 + full_align.aln.a.aln.anchors[1].seqpos
                            ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###
                            if len_b == aln_3prime_end_seq_b && aln_5prime_end_seq_a > 1
                                # i.e. if the 3' end of the sequence is the same as 3' end of alignment
                                seqb_extend = true # now compute the extended portion
                                n_bases_2_fetch_seq_a = full_align.aln.a.aln.anchors[1].seqpos # number 5' bases to be copied from the top strand
                                extension_b = reverse_complement(seq_a[1:n_bases_2_fetch_seq_a])
                                seq_b_new = seq_b*extension_b
                                sites[partner_coordinates].seq = seq_b_new                               
                                sites[partner_coordinates] = update_site(sites[partner_coordinates])                                
                                update_strand_tracker!(a_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(a_pol_r,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_r,partner_coordinates,string(seq_b_new))
                                global extension_events += 1
                                global extension_events_this_cycle += 1
                                if length(seq_b_new) > length_thresh
                                    push!(long_sites, partner_coordinates)
                                end   
                            else
                                #println("anneal with no b extension")
                                seqb_extend = false
                            end
                            #println("error finding 3' end of one of the sequences in alignment")
                        end
                    #catch
                        #println("warning - failure to assign 3' end to seq a alignment")
                    #end    
                else
                    global reject_count += 1
                end # annealing probability satisfied check - no annealing otherwise
            else
                global reject_count += 1
            end #priming end occurs in strand check
        else
            global reject_count += 1
        end #isdefined site check
    end #interaction loop
    println("end of cycle: ", cycle, "number of long sites: ", length(long_sites), "cumulative extensions: ", extension_events)
    println("    number of anealling rejections this cycle: ", reject_count)
    println("    number of attempted interactions this cycle: ", i_no)
    println("    number of annealing events this cycle: ", annealings)
    println("    number of extension events this cycle: ", extension_events_this_cycle)
    #println(i_no - (annealings + reject_count))
end #cycle loop

scatter(range(1,n_cycles), extension_vector, color="black",alpha=1,markerstrokewidth=.1,markersize=5,label = "total extensions", ylabel = "number of extensions", xlabel = "cycle")
plot!(size=(400,400))
savefig("extensions_per_cycle.png")


#######################################
#      quick and dirty ecori cut      #
#######################################
# go through every long_site sequence and search for the ecori recognition site
# perform a strand modification to update the new sequence
# do a search for concatemers and plot incidence of crosslinking
restriction_seq = "GAATTC"
cut_index = 1
new_long_sites = []
for i in 1:length(long_sites)
    inspected_seq = String(sites[long_sites[i]].seq)
    if occursin(restriction_seq, inspected_seq)
        #### MODIFY STRAND HERE WITH THE NEW SEQUENCE ####
        new_strand_index = findfirst(restriction_seq,inspected_seq)[1]+cut_index
        sites[long_sites[i]].seq = inspected_seq[1:new_strand_index-1]
        sites[long_sites[i]] = update_site(sites[long_sites[i]])
        println(sites[long_sites[i]].seq)
    end
    if length(sites[long_sites[i]].seq) > length_thresh
        push!(new_long_sites, long_sites[i])
    end
end
long_sites = deepcopy(new_long_sites)



######################################
# second round of annealing after cut
######################################

l_prime = l_prime_2
n_cycles = n_cycles_2
Temp = Temp_2
extension_events = 0
for cycle in 1:n_cycles
    output_tracking_plot([a_pol_f,a_pol_r,b_pol_f,b_pol_r,tracker_400],site_positions,"test"*string(cycle+2000)*".png")

    println("start of cycle: ", cycle, "number of long sites: ", length(long_sites), "extensions: ", extension_events)
    
    global free_long_sites = deepcopy(long_sites)
    global i_no = 0
    global annealings = 0
    global reject_count = 0
    global extension_events_this_cycle = 0
    while length(free_long_sites) > 0 && reject_count < rejection_thresh
        shuffle!(free_long_sites)
        candidate = pop!(free_long_sites)
        global i_no += 1
        # now attempt a random interaction from the site_picks list
        if isdefined(sites[candidate],2)
            main_site = sites[candidate]
            weights = [j[3] for j in main_site.free_neighbors]
            choices = [j[1] for j in main_site.free_neighbors]
            partner_coordinates = sample(choices, Weights(weights))
            partner_site = sites[partner_coordinates]
            ############# pairwise local alignment check here! ###############
            seq_a = main_site.seq
            seq_b = partner_site.seq
            if occursin(string(seq_a[end-l_prime:end]), string(reverse_complement(seq_b))) || occursin(string(seq_b[end-l_prime:end]), string(reverse_complement(seq_a))) #because exact match much quicker to search than performing pairwise align at this stage
                full_align = pairalign(LocalAlignment(), seq_a, reverse_complement(seq_b), scoremodel)
                p_anneal = get_p_annealing(Temp, get_T_m(full_align.aln.b; CNa = .1));
                if rand() < p_anneal
                    global annealings += 1
                    # an annealing event is approved, now compute extension and update everything
                    #try
                        ### ANNEALING BLOCK ####
                        # now we need to check what type of extension should be generated given the alignment
                        seq_a_new = BioSequence{DNAAlphabet{4}}()*seq_b
                        seq_b_new = BioSequence{DNAAlphabet{4}}()*seq_b
                        len_a = length(seq_a) # 3' end of the original sequences
                        len_b = length(seq_b) # 3' end of the original sequences
                        ### fetch 3' ends of both strands partipating in the alginment ###
                        if full_align.aln.a.aln.anchors[end].op == OP_SEQ_MATCH && full_align.aln.a.aln.anchors[1].op == OP_START
                            aln_3prime_end_seq_a = full_align.aln.a.aln.anchors[end].seqpos # the end of the top strand is the 3' end of the top strand
                            aln_5prime_end_seq_b = 1 + len_b - full_align.aln.a.aln.anchors[end].refpos # the end of the bottom strand
                                ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###
                            if len_a == aln_3prime_end_seq_a && aln_5prime_end_seq_b > 1 # i.e. if the 3' end of the sequence is the same as 3' end of alignment    
                                # i.e. if the 5' end of the bottom alignment isn't also the 5' end
                                seqa_extend = true # now compute the extended portion
                                n_bases_2_fetch_seq_b = len_b - full_align.aln.a.aln.anchors[end].refpos # number 5' bases to be copied from the bottom strand seqb
                                extension_a = reverse_complement(seq_b[1:n_bases_2_fetch_seq_b])
                                seq_a_new = seq_a*extension_a # create a new sequence and update the molecule archive
                                #### MODIFY STRAND HERE WITH THE NEW SEQUENCE ####
                                sites[candidate].seq = seq_a_new
                                sites[candidate] = update_site(sites[candidate])
                                global extension_events += 1
                                global extension_events_this_cycle += 1
                                update_strand_tracker!(a_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(a_pol_r,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_r,candidate,string(seq_a_new))
                                update_strand_tracker!(tracker_400,candidate,length(sites[candidate].seq))
                            else
                                #println("anneal with no a extension")
                                seqa_extend = false
                            end
                            aln_3prime_end_seq_b = len_b - full_align.aln.a.aln.anchors[1].refpos # the start base of the revcom is the 3' end, to convert t
                            aln_5prime_end_seq_a = 1 + full_align.aln.a.aln.anchors[1].seqpos
                            ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###
                            if len_b == aln_3prime_end_seq_b && aln_5prime_end_seq_a > 1
                                # i.e. if the 3' end of the sequence is the same as 3' end of alignment
                                seqb_extend = true # now compute the extended portion
                                n_bases_2_fetch_seq_a = full_align.aln.a.aln.anchors[1].seqpos # number 5' bases to be copied from the top strand
                                extension_b = reverse_complement(seq_a[1:n_bases_2_fetch_seq_a])
                                seq_b_new = seq_b*extension_b
                                sites[partner_coordinates].seq = seq_b_new                               
                                sites[partner_coordinates] = update_site(sites[partner_coordinates])                                
                                update_strand_tracker!(a_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(a_pol_r,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_r,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(tracker_400,partner_coordinates,length(sites[candidate].seq))
                                global extension_events += 1
                                global extension_events_this_cycle += 1
                                println(seq_b_new)
                                if length(seq_b_new) > length_thresh
                                    push!(long_sites, partner_coordinates)
                                end   
                            else
                                #println("anneal with no b extension")
                                seqb_extend = false
                            end
                            #println("error finding 3' end of one of the sequences in alignment")
                        end
                    #catch
                        #println("warning - failure to assign 3' end to seq a alignment")
                    #end    
                else
                    global reject_count += 1
                end # annealing probability satisfied check - no annealing otherwise
            else
                global reject_count += 1
            end #priming end occurs in strand check
        else
            global reject_count += 1
        end #isdefined site check
    end #interaction loop
    println("end of cycle: ", cycle, "number of long sites: ", length(long_sites), "cumulative extensions: ", extension_events)
    println("    number of anealling rejections this cycle: ", reject_count)
    println("    number of attempted interactions this cycle: ", i_no)
    println("    number of annealing events this cycle: ", annealings)
    println("    number of extension events this cycle: ", extension_events_this_cycle)
    #println(i_no - (annealings + reject_count))
end #cycle loop

search_seq = BioSequence{DNAAlphabet{4}}("AAAAAAAAAAAAAAAAAAAAAATNNNNNNNNNNNNNNNNNNNNNNNNCGCGTATTGTCAACNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGTCCCTATAGTG")
writer_r1 = FASTQ.Writer(open("edge_strands_direct.fastq", "w"))
for i in 1:length(long_sites)
    seq_i = sites[long_sites[i]].seq
    
    if length(String(seq_i)) > 400 && occursin("TTTTTTTTTTTTTTTTTTTTTT",string(seq_i))
        seq_i_string = string(seq_i)[end-200:end]
        #println(seq_i_string)
        record_object_r1 = BioSequences.FASTQ.Record(string(i), seq_i_string, repeat([1],length(seq_i_string)))
        write(writer_r1, record_object_r1)
    end
end
close(writer_r1)
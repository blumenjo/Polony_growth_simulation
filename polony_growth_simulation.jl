###################################################
# Mutable Structs
###################################################

mutable struct Site
    """
    This struct holds the important information for each site.
    For each Site the interaction parameters, neighbors, the sequence, and UMI is stored.
    """
    interact_params::Vector{Float64}
    free_neighbors::Vector{Tuple{Tuple{Float64, Float64}, Float64,Float64}} 
    seq::BioSequence{DNAAlphabet{4}}
    UMI::BioSequence{DNAAlphabet{4}}
    Site() = new()
end


mutable struct Strand_tracker
    """
    This struct keeps track of each sites origin and main charactes.
    """
    sequence_tag::String
    xs::Vector{Any}
    ys::Vector{Any}
    color::String
    length::Int64
    series_name::String
    neighbor_positions::Dict{Any, Any}
    Strand_tracker() = new()
end


###################################################
# low level function library
###################################################

# functions managing the sites 
function Strand_tracker(sequence_tag::String, color::String,series_name::String)
    """
    Create a Strand_tracker object to track a sequence with a given sequence tag, color, and series name.
    The function is mainly used to create the output plot

    # Arguments
    - `sequence_tag`: A string representing the sequence tag.
    - `color`: A string representing the color assigned to the sequence.
    - `series_name`: A string representing the name of the series.

    # Returns
    - `out`: A Strand_tracker object with the specified attributes.

    """
    out = Strand_tracker()
    out.sequence_tag = sequence_tag
    out.xs = []
    out.ys = []
    out.color = color
    out.series_name = series_name
    return out
end

function Strand_tracker(length::Int64, color::String,series_name::String)
    """
    Create a Strand_tracker object to track a sequence with a given sequence tag, color, and series name.
    The function is mainly used to create the output plot.

    # Arguments
    - `length`: A interger representing the length of the site.
    - `color`: A string representing the color assigned to the sequence.
    - `series_name`: A string representing the name of the series.
    - `neighbor_positions`:  A dictionary holding the neighbor sites

    # Returns
    - `out`: A Strand_tracker object with the specified attributes.

    """
    out = Strand_tracker()
    out.length = length
    out.xs = []
    out.ys = []
    out.color = color
    out.series_name = series_name
    out.neighbor_positions = Dict{}()
    return out
end

function update_strand_tracker!(tracker::Strand_tracker, coordinates, query_sequence::String)
    """
    This function checks if the tracker's sequence tag is present in the query sequence. If it is, it updates the tracker's coordinates with the given coordinates.

    # Arguments
    - `tracker`: A Strand_tracker object to be updated.
    - `coordinates`: A tuple representing the coordinates to be added to the tracker.
    - `query_sequence`: A string representing the site sequence to be searched for in the tracker's sequence tag.

    # Returns
    - `tracker`: The updated Strand_tracker object.
    """

    if occursin(tracker.sequence_tag, query_sequence)
        push!(tracker.xs, coordinates[1])
        push!(tracker.ys, coordinates[2])
    end
    return tracker
end

function update_strand_tracker!(tracker::Strand_tracker, coordinates, query_length::Int64, neighbor_position::Tuple{Float64, Float64})
    """
    This function checks if the length of the query is greater than the length of the tracker. 
    If it is, it updates the tracker's coordinates with the given coordinates and adds the neighbor position to the tracker's neighbor positions dictionary.

    # Arguments
    - `tracker::Strand_tracker`: A Strand_tracker object to be updated.
    - `coordinates`: A tuple representing the coordinates to be added to the tracker.
    - `query_length`: An integer representing the length of the query.
    - `neighbor_position`: A tuple representing the neighbor position to be added to the tracker.

    # Returns
    - `tracker`: The updated Strand_tracker object.
    """
    if query_length > tracker.length
        push!(tracker.xs, coordinates[1])
        push!(tracker.ys, coordinates[2])

        if haskey(tracker.neighbor_positions, coordinates)
            push!(tracker.neighbor_positions[coordinates], neighbor_position)
        else
            tracker.neighbor_positions[coordinates] = [neighbor_position]
        end

    end
    return tracker
end

# functions for seeding and priming 
function renormalize_site_probabilities(sites::Dict{Tuple{Float64, Float64}, Site})
    """
    This function iterates over each site in the dictionary and renormalizes the probabilities of free neighbors.
    It calculates the normalization factor as the sum of probabilities of all free neighbors and then divides each probability by the normalization factor.

    # Arguments
    - `sites`: A dictionary where keys are tuples representing coordinates and values are Site objects.

    # Returns
    - `sites`: The dictionary of sites with renormalized probabilities.
    """
    for site in sites
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
    """
    This function recalculates the probabilities of interaction for each free neighbor of the site based on its sequence length and distance to the neighbor. 
    It first calculates the interaction probability between the site and each free neighbor using the `p_end_to_end` function. Then, it normalizes the probabilities of interaction to ensure that they sum up to 1.

    # Arguments
    - `site`: A Site object to be updated.

    # Returns
    - `site`: The updated Site object.

    """

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
 
    return site
end

function bc_gen(seq)
    """
    This function takes a sequence as input and substitutes "N" with random nucleotides to generate a full sequence barcode.
    It also generates a unique barcode by randomly selecting nucleotides for the "N" positions.
    The function uses a lookup table to determine possible substitutions for "N" and other nucleotides.


    # Arguments
    - `seq`: A string representing the input sequence.

    # Returns
    - `barcode`: A BioSequence object representing the full sequence with "N" substituted.
    - `unique_barcode`: A BioSequence object representing the unique barcode with "N" substituted.

    """
    barcode = BioSequence{DNAAlphabet{4}}()            # Full sequence
    unique_barcode = BioSequence{DNAAlphabet{4}}()     # Random domain NNNN "unique barcode"
    lookup = Dict{BioSequence{DNAAlphabet{4}},Array{BioSequence{DNAAlphabet{4}},1}}() # an empty dictionary, created to look up possible substitution for N
    # the following defines possible substitution for the nucleotides on the left side
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

    #the following loop goes over the seqeunce and replaces the N 
    for i in 1:length(seq)
        # Store random domain in "unique barcode"
        # if sequence is N a random substitution is selected from lookup and saved to unique_barcode adn in the end to barcode 
        if seq[i:i] == dna"N"
            unique_barcode_letter = rand(lookup[seq[i:i]])
            unique_barcode = unique_barcode * unique_barcode_letter
        # if the character is not N a random substitution is chosen but not added to the unique barcode, only to barcode
        else 
            unique_barcode_letter = rand(lookup[seq[i:i]])
        end
        barcode = barcode*unique_barcode_letter
    end
    return barcode, unique_barcode
end

# functions for Annealing and extension 
function get_p_annealing(Temp, T_m; 𝜅 = 40.)
    """
    Calculate the probability of annealing based on temperature and melting temperature.
    This function calculates the probability of annealing based on the current temperature (`Temp`) and the melting temperature (`T_m`). 
    The "melting constant" (`𝜅`) controls the sharpness of the transition. A higher value of `𝜅` results in a sharper transition.
    The probability of annealing is calculated using the logistic function.

    # Arguments
    - `Temp`: The current temperature.
    - `T_m`: The melting temperature.
    - `𝜅`: The "melting constant". Default is 40.

    # Returns
    - `p_anneal`: The probability of annealing.
    """

    p_anneal = 1 - 1/(1 + exp(-1*𝜅*(Temp*0.01 - 0.01*T_m)))
    return p_anneal
end

function get_thermodynamic_parameters()
    """
    This function returns the thermodynamic parameters for DNA sequences, including the enthalpy change (ΔH°) and entropy change (ΔS°). 
    The parameters are stored in dictionaries, where the keys are BioSequence objects representing DNA sequences, and the values are the corresponding thermodynamic values.

    # Returns
    - `ΔH°`: A dictionary containing the enthalpy change values for DNA sequences.
    - `ΔS°`: A dictionary containing the entropy change values for DNA sequences.

    """
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
    """
    This function calculates the melting temperature (T_m) for a given nearest neighbor sequence based on the nearest-neighbor model. 
    It considers the enthalpy change (ΔH°) and entropy change (ΔS°) for each nearest neighbor pair in the sequence, as well as the sodium ion concentration (CNa). 
    The T_m calculation is based on the thermodynamic parameters and formula derived from the nearest-neighbor model.

    # Arguments
    - `nearest_neighbor_sequence::BioSequence{DNAAlphabet{4}}`: The nearest neighbor sequence for which T_m is calculated.
    - `CNa`: The sodium ion concentration (in mol/L). Default is 0.5.
    - `ΔH°`: A dictionary containing the enthalpy change values for DNA sequences. Default is global ΔH°.
    - `ΔS°`: A dictionary containing the entropy change values for DNA sequences. Default is global ΔS°.
    - `cutoff_len`: The cutoff length for considering annealing of the last bases. Default is 18.

    # Returns
    - `T_m`: The melting temperature.

    """

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
        
    end
    ΣΔH += ΔH°[nearest_neighbor_sequence[1:1]]*1000
    ΣΔH += ΔH°[nearest_neighbor_sequence[end:end]]*1000
    ΣΔS += ΔS°[nearest_neighbor_sequence[1:1]]
    ΣΔS += ΔS°[nearest_neighbor_sequence[end:end]]

    T_m = ΣΔH/(ΣΔS - ℛ*log(1000000)) + 16.6*log10(CNa) - 273

    return T_m
end

function p_end_to_end(L::Int64, r::Float64)
    """
    This function calculates the probability of the end-to-end distance for a sequence of length `L` given the distance `r` between the onsidered primer site. 
    It uses a model based on the worm-like chain (WLC) model to estimate the probability distribution of polymer configurations. 
    The calculation considers the contour length of the polymer chain, the Kuhn length, and the physical extension.

    # Arguments
    - `L`: The length of the sequence.
    - `r`: The distance between the sequence and primer site.

    # Returns
    - `p_end_to_end`: The probability of end-to-end distance.

    """

    s0 = L*lbp # calculating the countor length/ length at max physicall extension
    p_end_to_end = (r*exp(-1*s0/(8*lkuhn*(1-(r/s0)^2))))/((1-(r/s0)^2)*(2-(r/s0)^2)^2)
    
    if p_end_to_end < 0
        return 0
    else
        return p_end_to_end
    end
end


# functions for analysis
function polony_size(list_with_original, dic_with_neighboors)
    """
    This function calculates the size and coordinates of polonies based on the original points and their neighboring points.
    For each original point, it finds the furthest neighbor and considers it as one end of the polony. 
    Then, it finds the furthest point from this end and considers it as the other end of the polony. 
    The distance between these two points represents the size of the polony. 
    The function returns a dictionary containing the size and coordinates of both ends of polonies.

    # Arguments
    - `list_with_original`: A dictionary containing original points as keys and their coordinates as values.
    - `dic_with_neighbors`: A dictionary containing original points as keys and lists of neighboring points as values.

    # Returns
    - `polony_data`: A dictionary containing the size and coordinates of both ends of polonies, where the keys are original points and the values are tuples containing the size and coordinates.
    """

    polony_data = Dict{Int, Tuple{Float64, Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}}()
   
    for (key,value) in list_with_original
        template_distance, longest_distance = 0, 0
        X, Y, x, y = 0, 0, 0, 0
        XY, xy = (), ()

        x1 = value[1]
        y1 = value[2]

        for neighboor in dic_with_neighboors[key]
       
            x2 = neighboor[1]
            y2 = neighboor[2]
            
            distance = sqrt((x1-x2)^2+(y1-y2)^2) 

            if distance > template_distance || (distance == 0 && template_distance == 0)
                template_distance = distance 
                X = x2 
                Y = y2
            end 
        end
      
        longest_distance = template_distance
        
        for neighboor in dic_with_neighboors[key]
            x3 = neighboor[1]
            y3 = neighboor[2]
            
            distance = sqrt((X-x3)^2+(Y-y3)^2) 

            if distance > longest_distance
                longest_distance = distance 
                x = x3
                y = y3
        
            end 
        end
        xy = (x, y)
        XY = (X, Y)
        polony_data[key] = (longest_distance, (xy, XY))
    end
    return polony_data
end


function UMI_count(long_sites, sites)
    """
    This function counts UMIs and establishes polony connections based on the provided long sites and site information. 
    It iterates through the long sites, extracts UMIs, and establishes connections between UMIs found within the same polony. 
    The function returns an array of unique UMIs (`UMI_counter`) and a dictionary of polony connections (`polony_connections`).

    # Arguments
    - `long_sites`: An array of indices representing long sites.
    - `sites`: A dictionary containing site information.

    # Returns
    - `UMI_counter`: An array containing unique UMIs.
    - `polony_connections`: A dictionary containing polony connections.

    """
    
    UMI_counter = []
    polony_connections=Dict{String, Vector{String}}()


    for i in 1:length(long_sites)
    
        seq_i = sites[long_sites[i]].seq
       
        if length(sites[long_sites[i]].UMI) == 0  && occursin("ATTTTTTTTTTTTTTTTTTTTTT",string(seq_i))
            index_a = findfirst("ATTTTTTTTTTTTTTTTTTTTTT", string(seq_i))[1]
            UMI_a = seq_i[index_a - 62 : index_a - 38]
            UMI_b = seq_i[index_a - 24 : index_a]
           
    
        elseif occursin("ATTTTTTTTTTTTTTTTTTTTTT",string(seq_i))
            index_a = findfirst("ATTTTTTTTTTTTTTTTTTTTTT", string(seq_i))[1]
            UMI_a = sites[long_sites[i]].UMI
            UMI_b = seq_i[index_a - 24 : index_a]

        end
        
        if length(String(seq_i)) > 400 && occursin("ATTTTTTTTTTTTTTTTTTTTTT",string(seq_i))  
           
            if occursin(string(UMI_a), string(seq_i))
                both_UMI = string(seq_i)[findfirst(UMI_a, seq_i)[1]: findfirst("ATTTTTTTTTTTTTTTTTTTTTT", string(seq_i))[1] - 1]

                if both_UMI in UMI_counter
                    #do nothinh
                else
                    push!(UMI_counter, both_UMI)

                    if !haskey(polony_connections, string(UMI_a))
                        polony_connections[string(UMI_a)] = []
                    end

                    push!(polony_connections[string(UMI_a)], string(UMI_b))

                    if !haskey(polony_connections, string(UMI_b))
                        polony_connections[string(UMI_b)] = []
                    end

                    push!(polony_connections[string(UMI_b)], string(UMI_a))
                    
                end

            else
                continue 
            end 
        
        end
    end
    
    return UMI_counter, polony_connections
end


function distance_func(X1, Y1, X2, Y2)
    """
    This function calculates the Euclidean distance between two points specified by their x and y coordinates. 
    It applies the formula for Euclidean distance, which is the square root of the sum of the squared differences between the coordinates.

    # Arguments
    - `X1`: x-coordinate of the first point.
    - `Y1`: y-coordinate of the first point.
    - `X2`: x-coordinate of the second point.
    - `Y2`: y-coordinate of the second point.

    # Returns
    - `distance`: Euclidean distance between the two points.

    """
    distance = sqrt((X1-X2)^2+(Y1-Y2)^2)
    return distance
end


function output_tracking_plot(trackers::Vector{Strand_tracker}, site_positions, plot_name::String)
    """
    This function generates a tracking plot based on the provided trackers and site positions. 
    It creates a scatter plot with primer sites marked in gray and each tracked strand represented by colored markers. 
    The function saves the generated plot with the specified name.

    # Arguments
    - `trackers`: An array of `Strand_tracker` objects representing tracked strands.
    - `site_positions`: A collection of site positions.
    - `plot_name`: The name of the output plot file.

    # Returns
    - Nothing. The function saves the plot with the specified name.

    """
    # Create scatter plot of primer sites
    p = scatter([i[1] for i in site_positions], [i[2] for i in site_positions], color="#E4E4E4", alpha=0.5, markerstrokewidth=0, markersize=1, label="primer sites")
    
    # Add markers for each tracked strand
    for tracker in trackers
        scatter!(tracker.xs, tracker.ys, color=tracker.color, alpha=1, markerstrokewidth=0, markersize=2, label=tracker.series_name)
    end
    
    # Set plot size and save the plot
    plot!(size=(400, 400))
    savefig(plot_name)
end


function output_tracking_plot_edges(trackers::Vector{Strand_tracker}, site_positions, all_distance, plot_name::String)
    """
    This function generates a tracking plot with edges and histograms based on the provided trackers, site positions, and distances.
    It creates a scatter plot with primer sites marked in gray and plots edges between connected points for strands with specific colors.
    It also generates histograms to visualize the distribution of connection lengths for all distances and crossover connections. 
    The function saves the generated plots with the specified names.

    # Arguments
    - `trackers`: An array of `Strand_tracker` objects representing tracked strands.
    - `site_positions`: A collection of site positions.
    - `all_distance`: All distances.
    - `plot_name`: The name of the output plot file.

    # Returns
    - Nothing. The function saves the plots with the specified names.

    """
    # Create scatter plot of primer sites
    p = scatter([i[1] for i in site_positions], [i[2] for i in site_positions], color="#E4E4E4", alpha=0.5, markerstrokewidth=0, markersize=1, label="primer sites")
    
    
    connection_distance = []
    for tracker in trackers
        if tracker.color == "#000000"  # Check if the tracker color is black (indicating a specific type of strand)

            x_coords = []
            y_coords = []

            # Iterate through the neighbor positions of the tracker
            for (key, value) in tracker.neighbor_positions
                candidate_x_pos = key[1]
                candidate_y_pos = key[2]
                
                # Iterate through the partner coordinates of the neighbor positions
                for partner_coords in value
                    partner_x = partner_coords[1]
                    partner_y = partner_coords[2]

                    # Calculate the distance between candidate and partner coordinates
                    distance = distance_func(candidate_x_pos, candidate_y_pos, partner_x, partner_y)

                    # Append candidate and partner coordinates to arrays for plotting
                    push!(x_coords, [candidate_x_pos, partner_x, NaN])  # NaN creates breaks in the line
                    push!(y_coords, [candidate_y_pos, partner_y, NaN])
                    push!(connection_distance, distance)
                end
            end

            # Flatten the arrays
            x_coords = vcat(x_coords...)
            y_coords = vcat(y_coords...)

            # Plot the lines representing edges
            plot!(x_coords, y_coords, color=:black, linewidth=0.5, linestyle=:solid, label="", legend=false)
        end
    end
    
    # Save the plot with edges
    plot!(size=(400, 400))
    savefig(plot_name * ".png")
    
    # Generate and save histograms for connection lengths
    histogram(all_distance, bins=100, xlabel="Connection length", ylabel="Frequency", title="Distribution of connections")
    savefig(plot_name * "_strand_connection.png")
    histogram(connection_distance, bins=100, xlabel="Connection length", ylabel="Frequency", title="Distribution of crossover connections")
    savefig(plot_name * "_crossover_connection.png")
end

function plot_graph(cycles, extensions, Output_file_name, graph_table)
    """
    This function generates and saves a plot of extensions per cycle based on the provided data. 
    It creates a scatter plot with the x-axis representing cycle numbers and the y-axis representing the number of extensions. 
    The function also provides an option to save the data to a text file if `graph_table` is `true`.

    # Arguments
    - `cycles`: Total number of cycles.
    - `extensions`: Extensions data.
    - `Output_file_name`: Base name for the output plot and text file.
    - `graph_table`: A boolean indicating whether to generate a text file with extension data.

    # Returns
    - Nothing. The function saves the plot and optionally the data to files.

    """
    # Convert cycle and extension data to Float64 arrays
    x_data = Float64.(collect(2:cycles))
    y_data = Float64.(extensions[2:end])
    
    # Generate scatter plot
    scatter(x_data, y_data, color="black", alpha=1, markerstrokewidth=0.1, markersize=5, label="Total extensions", ylabel="Number of extensions", xlabel="Cycle", size=(400, 400))

    # Save the plot
    savefig(Output_file_name * "_per_cycle.png")
    
    # Optionally save data to a text file
    if graph_table == true 
        output_file = open(Output_file_name * "_growth_rate_data.txt", "w")
        for (x, y) in zip(x_data, y_data)
            println(output_file, "$x\t$y")
        end
        close(output_file)
    end
end


#functions for the Output

function writing_output(parameter::param, polony_data, extension_events, pc_seq, polonies_strand, long_sites, sites; Table=true, Polony=true, General=true, FASTQ_file=true)
    """

    This function writes various output information about polonies to text files and FASTQ files.

    # Arguments
    - `parameter`: A parameter object containing parameters of the experiment.
    - `polony_data`: A dictionary containing the size and coordinates of polonies.
    - `extension_events`: Number of extension events.
    - `pc_seq`: A dictionary containing polony sequences.
    - `polonies_strand`: A dictionary containing polony information.
    - `long_sites`: List of long sites.
    - `sites`: A dictionary containing site information.

    # Keyword Arguments
    - `Table`: Boolean indicating whether to write table information (default: true).
    - `Polony`: Boolean indicating whether to write polony information (default: true).
    - `General`: Boolean indicating whether to write general information (default: true).
    - `FASTQ_file`: Boolean indicating whether to write FASTQ files (default: true).
    """
    
    # Count UMIs and polony connections
    UMI_counter, polony_connections = UMI_count(long_sites, sites)

    # Initialize variables for maximum and mean distance
    max_distance = 0
    distance_one = 0
    
    # Write polony information if Polony is true
    if Polony == true
        writer = open(parameter.Output_file_name * "_polony_info.txt", "w")
        for (key, values) in pc_seq
            # Write polony information to file
            write(writer, "Polony "*string(key)*" \n")
            write(writer, "Start X, Y coordinate: "*string(values[1], values[2])*" \n")
            write(writer, "Sequence: "*string(values[3])*"\n")
            write(writer, "Distance: "*string(polony_data[key][1])*"nm \n")
            write(writer, "Coordinates: "*string(polony_data[key][2])*"\n")
            write(writer, "Number of Strands: "*string(length(polonies_strand[key]))*"\n")
            write(writer, "\n")
            
            # Update distance metrics
            distance_one += polony_data[key][1]
            if polony_data[key][1] > max_distance
                max_distance = polony_data[key][1]
            end
        end
        close(writer)
    end
    

    if Polony == false
        for (key,value) in pc_seq
            # Update distance metrics if Polony is false
            distance_one += polony_data[key][1]
            if polony_data[key][1] > max_distance
                max_distance = polony_data[key][1]
            end
        end
    end 

    # Calculate mean distance
    mean_distance = distance_one / length(polony_data)
   
    # Write table information if Table is true
    if Table == true 
        open(parameter.Output_file_name * "_tab_info.txt", "w") do io
            write(io, string(parameter.n_cycles) * "," * string(parameter.spot_diameter) * "," * string(parameter.site_area) * "," * string(parameter.site_density) * "," * string(parameter.proportions[5]) * "," * string(length(pc_seq)) * "," * string(extension_events) * "," * string(length(UMI_counter)) * "," * string(parameter.cutoff) * "," * string(mean_distance) * "," * string(max_distance) * "\n")
        end
    end

    # Write general information if General is true
    if General == true
        
        total_length = sum(length(value) for value in values(polonies_strand))
        total_connection = sum(length(value) for value in values(polony_connections))
        num_keys = length(keys(polonies_strand))
        num_connection = length(keys(polony_connections))
       
        average_length = total_length / num_keys

        average_connection = total_connection / num_connection
       
       
        open(parameter.Output_file_name * "_run_info.txt", "w") do s
            write(s, "Number of cycles "*string(parameter.n_cycles)*" \n")
            write(s, "spot diameter: "*string(parameter.spot_diameter)*"nm \n")
            write(s, "spot area: "*string(parameter.site_area)*"nm2 \n")
            write(s, "site density: "*string(parameter.site_density)*"strands per nm2 \n")
            write(s, "input sequences used: "*string(parameter.sequences)*"\n")
            write(s, "input 1_proportions of sequences: "*string(parameter.proportions)*"\n")
            write(s, "number of polonies: "*string(num_keys)*"\n")
            write(s, "number of extensions: "*string(extension_events)*"\n")
            write(s, "number of connected polonies: "*string(length(keys(polony_connections)))*"\n")
            write(s, "number of mean connection per polony: "*string(average_connection)*"\n")
            write(s, "number of mean number of strands per polony: "*string(average_length)*"\n")
            write(s, "max. length between sites to connect: "*string(parameter.cutoff)*"\n")
            write(s, "mean Polony size: "*string(mean_distance)*"\n")
            write(s, "max. Polonz size: "*string(max_distance)*"\n")
            
        end
    end 

    # Write FASTQ files if FASTQ_file is true
    if FASTQ_file == true

        writer_r1 = FASTQ.Writer(open(parameter.Output_file_name * "_edge_strands_direct.fastq", "w"))

        for i in 1:length(long_sites)
            seq_i = sites[long_sites[i]].seq
        
            if length(String(seq_i)) > 400 && occursin("TTTTTTTTTTTTTTTTTTTTTT", string(seq_i))
        
                seq_i_string = string(seq_i)[end-200:end]

                record_object_r1 = BioSequences.FASTQ.Record(string(i), seq_i_string, repeat([30], length(seq_i_string)))
                write(writer_r1, record_object_r1)
            
            end
        end
        close(writer_r1)
        
    end 

end


function reset!(list::Vector{Strand_tracker})
    for s in list
        s.sequence_tag=""
        s.xs=[]
        s.ys=[]
        s.color=""
        s.length = 0
        s.series_name=""
    end
end

###################################################
# high level function library
###################################################

# functions for seeding  and priming

function generate_site_positions(parameter::param)
    """

    This function generates random positions for sites according to a specified site density and area using basic possion distribution.

    # Arguments
    - `parameter`: A parameter object containing parameters of the PCR experiment.

    # Returns
    - `site_positions`: An array of tuples representing the positions of DNA sites on the substrate.

    """

    expected_N_points = parameter.site_density*parameter.site_area
    actual_N_points = rand(Poisson(expected_N_points)) #normalising with the possion distribution
    site_positions::Array{Tuple{Float64, Float64}} = []  # empty array created to store the generated positions , it will store tuples of 2 points
    radius = sqrt(parameter.site_area/π)
    println(expected_N_points)
    for i in 1:actual_N_points # loop to generate the random site points
        random_radius = radius*sqrt(rand())
        random_angle = rand()*2*π
        y = random_radius*sin(random_angle)
        x = random_radius*cos(random_angle)
        push!(site_positions, (x,y))
    end
    return site_positions
end

function matrix_compute_distance(site_positions::Array{Tuple{Float64, Float64}}, cutoff, N)
    """

    This function computes pairwise distances between sites in a range using a KDTree data structure for efficient nearest neighbor search.

    # Arguments
    - `site_positions`: An array of tuples representing the positions of sites.
    - `cutoff`: The cutoff distance for considering two sites as neighbors.
    - `N`: The number of  sites.

    # Returns
    - `all_distance`: A vector containing tuples of indices and corresponding distances between neighboring sites.

    """

    all_distance = Vector{Tuple{Int, Int, Float64}}()

    treeData = transpose(hcat(first.(site_positions), last.(site_positions)))
    tree = KDTree(treeData)

    for i in 1:N
        local_distance = Vector{Tuple{Int, Int, Float64}}()
        distance = inrange(tree, [treeData[1, i], treeData[2, i]], cutoff, true)

        for j in distance
            @inbounds begin
                if j <= i
                    continue
                end

                x1, y1 = treeData[1, i], treeData[2, i]
                x2, y2 = treeData[1, j], treeData[2, j]
                dis = sqrt((x1 - x2)^2 + (y1 - y2)^2)
                if dis <= cutoff && dis != 0
                    push!(local_distance, (i, j, dis))
                end
            end
        end

        # Append local_distance to all_distance
        append!(all_distance, local_distance)
    end

    return all_distance
end


function Matrix_initialize_sites(site_positions::Array{Tuple{Float64, Float64}}, parameter::param, length_thresh::Int64=40)
    """

    This function initializes  sites on a substrate and assigns sequences to each site based on a given set of sequences and their proportions. 
    Additionally, it calculates and assigns neighboring sites to each site based on a cutoff distance.

    # Arguments
    - `site_positions`: An array of tuples representing the positions of sites on a substrate.
    - `parameter`: A parameter object containing parameters of the PCR experiment.
    - `length_thresh`: Threshold length for considering a site as "long".

    # Returns
    - `sites`: A dictionary containing information about DNA sites, including sequences and neighboring sites.
    - `long_sites`: A list of long sites that have sequences longer than the specified length threshold.
    - `pc_seq`: A dictionary containing polony sequences, indexed by site positions which become the poloy ID.

    """


    total_input_normalization_factor = sum(parameter.proportions) # normalising the proportion
    weights = parameter.proportions/total_input_normalization_factor # normalise weigths
    sequences = parameter.sequences # extract sequence from object
    N =  length(site_positions) # number of sites
    long_sites = Tuple{Float64, Float64}[]
    pc_seq = Dict{Int,Tuple{Float64, Float64,DNASequence}}()
  
 
    initialized_sites = [Site() for _ in 1:N]
    sites = Dict{Tuple{Float64, Float64}, Site}()

    for i in 1:N
        sites[site_positions[i]] = initialized_sites[i]
    end

    ### assign sequences to all sites in this loop###
    println("N is ", N)

    # the following loop assigns each site with a sequence 
    for i in 1:N   
        if i%1000 == 0
            print(i, " ")
        end
        main_site = sites[site_positions[i]]
        main_sequence, UMI  = bc_gen(sample(sequences, Weights(weights))) # smaple is a function from julia and randomly picks a sequence, here weigths are added with the Weigths function which provides a vector 
        # the main_sequence contains now the sequence with the unique barcode attachet !
        main_site.seq = main_sequence
        main_site.UMI = UMI
        if length(main_sequence) > length_thresh
            X,Y = site_positions[i]
            
            push!(long_sites, site_positions[i])
            push!(pc_seq, i => (X,Y,UMI)) # to store the initial template coordinates in a list --> they are stored as a tuble 
              
        end 
    end

    
    println("\n assigning neighbors based on cutoff distance")
    ### check all pairwise distances and assign neighbors based on adaptive cutoff distance###
    
    distance_list = @time matrix_compute_distance(site_positions, parameter.cutoff, N)
    
    println("compute_done")
    
    println(length(distance_list))


    num_entries = length(distance_list)

    Threads.@threads for k in 1:num_entries 
        i, j, val = distance_list[k]
        if i <= j 
            main_site = sites[site_positions[i]]
            L_MS = length(main_site.seq)
            candidate_position = site_positions[j]
            site_po = site_positions[i]
            candidate_neighbor_site = sites[site_positions[j]]
            L_CNS = length(candidate_neighbor_site.seq)
            MS_cutoff = parameter.cutoff
            CNS_cutoff = parameter.cutoff 
            x2, y2 = candidate_position
            x1, y1 = site_po
            p_interact_MS = p_end_to_end(L_MS, val)
            p_interact_CNS = p_end_to_end(L_CNS, val)

            if MS_cutoff > val
 
                if isdefined(sites[(x1,y1)], 2)
              
                    push!(main_site.free_neighbors, (candidate_position, val, p_interact_MS))  
                else
                   
                    main_site.free_neighbors = [(candidate_position, val, p_interact_MS)] 
                end
            else
  
            end

            if CNS_cutoff > val
              
                if isdefined(sites[(x2,y2)], 2)
             
                    push!(candidate_neighbor_site.free_neighbors, (site_po, val, p_interact_CNS)) 
                else
            
                    candidate_neighbor_site.free_neighbors = [(site_po, val, p_interact_CNS)] 
                end
            else
              
            end
        end
        
    end

    return sites, long_sites, pc_seq
end

# functions for first PCR
function PCR_Run(parameter::param)  
    """
    This function simulates a PCR run by annealing and extending DNA strands.
    During each PCR cycle, the function tracks annealing, extension, and replication events for each site. 
    It also updates information about long sites and polonies.

    The PCR process involves several steps:
    1. Site Generation: Generate positions for DNA sites.
    2. Site Initialization: Initialize sites with sequences and properties.
    3. Site Renormalization: Renormalize probabilities of sites.
    4. Thermal Cycling: Iterate through PCR cycles, including annealing and extension steps.

    # Arguments
    - `parameter::param`: A parameter object containing parameters of the PCR experiment.

    # Returns
    - `sites`: A dictionary containing information about  sites, including sequences and properties.
    - `long_sites`: A list of long sites that have undergone successful extension.
    - `extension_vector`: A vector containing the number of extension events per cycle.
    - `polonies_strand`: A dictionary containing the coordinates of each strand in a polony, organized by polony ID.
    - `extension_events`: The total number of extension events that occurred throughout the PCR run.
    - `pc_seq`: A dictionary containing polony sequences, organized by polony ID.
    - `site_positions`: A list containing the positions of sites on a substrate.
    - `all_distance`: A list containing the distances between paired sites.

    #Detail
    - Extension Simulation: 
        During each PCR cycle, the function iterates over long sites (`long_sites`).
        For each site, it attempts annealing with a neighboring site chosen probabilistically based on interaction weights. 
        If annealing occurs, the function simulates DNA extension and crossover by appending a portion of the partner strand's sequence to the current strand, mimicking the behavior of DNA polymerases during PCR.
    
    - Annealing and Alignment: 
        Before extension, the function performs a pairwise local alignment check between candidate strands. 
        It ensures that the 3' end of one strand aligns with the complementary region of the other, facilitating efficient extension. 
        Annealing is governed by the temperature-dependent annealing probability calculated using the melting temperature of the alignment.
    
    - Strand Update and Monitoring: 
        Upon successful extension, the function updates the sequence of the current DNA site and records the extension event. 
        It also updates relevant trackers to monitor strand length and interactions during PCR amplification.
    
    - Cycle Management: 
        The function iterates over multiple PCR cycles, allowing for the accumulation of extension events over successive rounds of amplification.
        It tracks the number of annealing rejections, attempted interactions, and extension events per cycle, providing insights into the efficiency and progress of PCR amplification. 

    """

    # Initialize polony dictionary
    polonies_strand = Dict{Int, Vector{Tuple}}()

    println("generating site positions")
    # Generate site positions
    site_positions = generate_site_positions(parameter)
    println(length(site_positions))

    println("initializing sites")
    # Initialize sites
    sites, long_sites, pc_seq = @time Matrix_initialize_sites(site_positions, parameter, parameter.length_thresh)

    println("renormalizing site probabilities")
    # Renormalize site probabilities
    sites = renormalize_site_probabilities(sites)

    println("beginning thermal cycling")
    extension_events = 0
    a_ex = 0
    b_ex = 0
    extension_vector = []
    all_distance = []

    # Prepare the polony dictionary
    for (key,value) in pc_seq
        polonies_strand[key] = []
        push!(polonies_strand[key], (value[1], value[2]))
    end
   
    # Iterate over PCR cycles
    for cycle in 1:parameter.n_cycles

        if cycle % 10 == 0
            output_tracking_plot([a_pol_f, a_pol_r, b_pol_f, b_pol_r], site_positions, parameter.Output_file_name * string(cycle + 1000) * ".png")
        end
        
        println("start of cycle: ", cycle, "number of long sites: ", length(long_sites), "extensions: ", extension_events)
        push!(extension_vector, extension_events)
        free_long_sites = deepcopy(long_sites)
        i_no = 0
        annealings = 0
        reject_count = 0
        extension_events_this_cycle = 0
        
        # Iterate over annealing and extension events
        while length(free_long_sites) > 0 && reject_count < parameter.rejection_threshold
            shuffle!(free_long_sites)               # shuffle the free sequences
            candidate = pop!(free_long_sites)       # select last element from the free sites and store it in candidate and REMOVE it from the free sites list  
            i_no += 1
           
            # Attempt a random interaction from the site_picks list
            if isdefined(sites[candidate], 2)
                main_site = sites[candidate]
                
                weights = [j[3] for j in main_site.free_neighbors]
                choices = [j[1] for j in main_site.free_neighbors]
                partner_coordinates = sample(choices, Weights(weights))
                partner_site = sites[partner_coordinates]
                
                ############# pairwise local alignment check here! ###############
                seq_a = main_site.seq
                seq_b = partner_site.seq
                
                if occursin(string(seq_a[end-parameter.l_prime:end]), string(reverse_complement(seq_b))) || occursin(string(seq_b[end-parameter.l_prime:end]), string(reverse_complement(seq_a))) #because exact match much quicker to search than performing pairwise align at this stage
                    full_align = pairalign(LocalAlignment(), seq_a, reverse_complement(seq_b), scoremodel)
                    p_anneal = get_p_annealing(parameter.Temp, get_T_m(full_align.aln.b; CNa = .1));

                    ### ANNEALING BLOCK ####
                    if rand() < p_anneal
                        annealings += 1
                        # an annealing event is approved, now compute extension and update everything

                        ### EXTENSION BLOCK ####
                        # now we need to check what type of extension should be generated given the alignment
                        seq_a_new = BioSequence{DNAAlphabet{4}}() * seq_b # they are the same
                        seq_b_new = BioSequence{DNAAlphabet{4}}() * seq_b
                        len_a = length(seq_a)
                        len_b = length(seq_b)
                        ### fetch 3' ends of both strands partipating in the alginment ###
                     
                        if full_align.aln.a.aln.anchors[end].op == OP_SEQ_MATCH && full_align.aln.a.aln.anchors[1].op == OP_START
                            aln_3prime_end_seq_a = full_align.aln.a.aln.anchors[end].seqpos # the end of the top strand is the 3' end of the top strand
                            aln_5prime_end_seq_b = 1 + len_b - full_align.aln.a.aln.anchors[end].refpos # the end of the bottom strand
                            
                            ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###
                            
                            if len_a == aln_3prime_end_seq_a && aln_5prime_end_seq_b > 1 # i.e. if the 3' end of the sequence is the same as 3' end of alignment    
                                # i.e. if the 5' end of the bottom alignment isn't also the 5' end
                                a_ex += 1
                                seqa_extend = true # now compute the extended portion
                                n_bases_2_fetch_seq_b = len_b - full_align.aln.a.aln.anchors[end].refpos # number 5' bases to be copied from the bottom strand seqb
                                extension_a = reverse_complement(seq_b[1:n_bases_2_fetch_seq_b])
                                seq_a_new = seq_a*extension_a # create a new sequence and update the molecule archive
                                
                                #### MODIFY STRAND HERE WITH THE NEW SEQUENCE ####
                                sites[candidate].seq = seq_a_new
                                sites[candidate] = update_site(sites[candidate])
                                distance = distance_func(candidate[1],candidate[2],partner_coordinates[1],partner_coordinates[2])
                                push!(all_distance, distance)
                                extension_events += 1
                                extension_events_this_cycle += 1
                                update_strand_tracker!(a_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(a_pol_r,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_f,candidate,string(seq_a_new))
                                update_strand_tracker!(b_pol_r,candidate,string(seq_a_new))


                                if length(seq_a_new) > parameter.length_thresh
                                    push!(long_sites, partner_coordinates)
                                    
                                    for (key,value) in pc_seq
                                        if haskey(polonies_strand, key)      
                                            if candidate in polonies_strand[key]
                                                push!(polonies_strand[key], partner_coordinates) 
                                            end
                                        end
                                    end
                                end   
                                
                            else
                                seqa_extend = false
                            end
                            
                            aln_3prime_end_seq_b = len_b - full_align.aln.a.aln.anchors[1].refpos # the start base of the revcom is the 3' end, to convert t
                            aln_5prime_end_seq_a = 1 + full_align.aln.a.aln.anchors[1].seqpos
                            
                            ### if 3' annealed criterion is met for extension, we must also check the other strands 5' neighborhood for template material ###

                            if len_b == aln_3prime_end_seq_b && aln_5prime_end_seq_a > 1
                                # i.e. if the 3' end of the sequence is the same as 3' end of alignment
                             
                                seqb_extend = true # now compute the extended portion
                                b_ex += 1
                                n_bases_2_fetch_seq_a = full_align.aln.a.aln.anchors[1].seqpos # number 5' bases to be copied from the top strand
                              
                                extension_b = reverse_complement(seq_a[1:n_bases_2_fetch_seq_a])
                                seq_b_new = seq_b*extension_b
                                sites[partner_coordinates].seq = seq_b_new
                                sites[partner_coordinates] = update_site(sites[partner_coordinates])  
                                distance = distance_func(candidate[1],candidate[2],partner_coordinates[1],partner_coordinates[2])
                                push!(all_distance, distance)
                                                             
                                update_strand_tracker!(a_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(a_pol_r,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_f,partner_coordinates,string(seq_b_new))
                                update_strand_tracker!(b_pol_r,partner_coordinates,string(seq_b_new))
                                
                                extension_events += 1
                                extension_events_this_cycle += 1

                                if length(seq_b_new) > parameter.length_thresh
                                    push!(long_sites, partner_coordinates)
                                    
                                    for (key,value) in pc_seq
                                        if haskey(polonies_strand, key)      
                                            if candidate in polonies_strand[key]
                                                push!(polonies_strand[key], partner_coordinates) 
                                            end
                                        end
                                    end
                                end    
                            else
                                seqb_extend = false
                            end
                        end
                    else
                        reject_count += 1
                    end
                else
                    reject_count += 1
                end
            else
                reject_count += 1
            end
        end
        
        println("end of cycle: ", cycle, "number of long sites: ", length(long_sites), "cumulative extensions: ", extension_events)
        println("    number of annealing rejections this cycle: ", reject_count)
        println("    number of attempted interactions this cycle: ", i_no)
        println("    number of annealing events this cycle: ", annealings)
        println("    number of extension events this cycle: ", extension_events_this_cycle)
        
    end # cycle loop
    
    println("a extend:", a_ex)
    println("b extend:", b_ex)
    
    return sites, long_sites, extension_vector, polonies_strand, extension_events, pc_seq, site_positions, all_distance
end


# functions for second PCR / Crossover
function ecori_cut(long_sites, restriction_seq, cut_index, parameter::param, sites)
    """
    This function cuts the strands at EcoRI recognition sites and updates the sequences accordingly.
    
    # Arguments
    - `long_sites`: List of long sites.
    - `restriction_seq`: EcoRI recognition sequence to cut.
    - `cut_index`: Index to cut at in the recognition sequence.
    - `parameter`: A parameter object containing parameters of the PCR experiment.
    - `sites`: Dictionary containing site information.
    
    # Returns
    - `long_sites`: Updated list of long  sites after digestion.
    """

    new_long_sites = []

    # Iterate over each long site
    for i in eachindex(long_sites)
   
        inspected_seq = String(sites[long_sites[i]].seq)
        
        # Check if the restriction sequence occurs in the DNA sequence
        if occursin(restriction_seq, inspected_seq)
            # Cut the DNA sequence at the specified index after the recognition sequence
            new_strand_index = findfirst(restriction_seq, inspected_seq)[1] + cut_index
            sites[long_sites[i]].seq = inspected_seq[1:new_strand_index - 1]
            
            # Update the site if it's defined
            if isdefined(sites[long_sites[i]], 2)
                sites[long_sites[i]] = update_site(sites[long_sites[i]])
            end
        end
        
        # Check if the sequence length exceeds the threshold after cutting
        if length(sites[long_sites[i]].seq) > parameter.length_thresh
            push!(new_long_sites, long_sites[i])  # Add the site to new_long_sites
        end
    end

    long_sites = deepcopy(new_long_sites)  # Update long_sites with the new list
    return long_sites
end


function UMI_crossover(long_sites, parameter::param, sites, site_positions, all_distance)
    """

    This function simulates DNA extension events with UMI crossover during PCR cycles.

    # Arguments
    - `long_sites`: List of long  sites.
    - `parameter::param`: A parameter object containing parameters of the PCR experiment.
    - `sites`: Dictionary containing site information.
    - `site_positions`: List of site positions.
    - `all_distance`: List containing distances between paired sites.

    # Returns
    - `long_sites`: Updated list of long sites after extension events.
    - `ex_vec`: Vector containing the number of extension events per cycle.

    #Detail
    - Extension Simulation: 
        During each PCR cycle, the function iterates over long sites (`long_sites`).
        For each site, it attempts annealing with a neighboring site chosen probabilistically based on interaction weights. 
        If annealing occurs, the function simulates DNA extension and crossover by appending a portion of the partner strand's sequence to the current strand, mimicking the behavior of DNA polymerases during PCR.
    
    - Annealing and Alignment: 
        Before extension, the function performs a pairwise local alignment check between candidate strands. 
        It ensures that the 3' end of one strand aligns with the complementary region of the other, facilitating efficient extension. 
        Annealing is governed by the temperature-dependent annealing probability calculated using the melting temperature of the alignment.
    
    - Strand Update and Monitoring: 
        Upon successful extension, the function updates the sequence of the current DNA site and records the extension event. 
        It also updates relevant trackers to monitor strand length and interactions during PCR amplification.
    
    - Cycle Management: 
        The function iterates over multiple PCR cycles, allowing for the accumulation of extension events over successive rounds of amplification.
        It tracks the number of annealing rejections, attempted interactions, and extension events per cycle, providing insights into the efficiency and progress of PCR amplification.    

    """
    extension_events_link = 0
    ex_vec = []
    long_sites = long_sites
    a_ex = 0
    b_ex = 0

    # Iterate over PCR cycles
    for cycle in 1:parameter.n_cycles_2

        if cycle % 10 == 0
            output_tracking_plot([a_pol_f, a_pol_r, b_pol_f, b_pol_r, tracker_400], site_positions, parameter.Output_file_name * string(cycle + 2000) * ".png")
        end

        if cycle == parameter.n_cycles_2
            output_tracking_plot_edges([a_pol_f, a_pol_r, b_pol_f, b_pol_r, tracker_400], site_positions, all_distance, parameter.Output_file_name * string(cycle + 1000) * "_edges")
        end

        # Initialize cycle variables
        push!(ex_vec, extension_events_link)
        free_long_sites = deepcopy(long_sites)
        i_no = 0
        annealings = 0
        reject_count = 0
        extension_events_this_cycle = 0
        
        # Iterate over free long sites
        while length(free_long_sites) > 0 && reject_count < 10000000
            shuffle!(free_long_sites)
            candidate = pop!(free_long_sites)
            i_no += 1
            
            # Attempt a random interaction from the site_picks list
            if isdefined(sites[candidate], 2)
                main_site = sites[candidate]
                weights = [j[3] for j in main_site.free_neighbors]
                choices = [j[1] for j in main_site.free_neighbors]
                partner_coordinates = sample(choices, Weights(weights))
                partner_site = sites[partner_coordinates]

                ############# pairwise local alignment check here! ###############
                seq_a = main_site.seq
                seq_b = partner_site.seq
                if occursin(string(seq_a[end - parameter.l_prime_2:end]), string(reverse_complement(seq_b))) || occursin(string(seq_b[end - parameter.l_prime_2:end]), string(reverse_complement(seq_a)))
                    full_align = pairalign(LocalAlignment(), seq_a, reverse_complement(seq_b), scoremodel)
                    p_anneal = get_p_annealing(parameter.Temp_2, get_T_m(full_align.aln.b; CNa = 0.1))
                    
                    ### ANNEALING BLOCK ####
                    if rand() < p_anneal
                        annealings += 1

                        ### EXTENSION BLOCK ####
                        # Annealing event approved, compute extension
                        seq_a_new = BioSequence{DNAAlphabet{4}}() * seq_b
                        seq_b_new = BioSequence{DNAAlphabet{4}}() * seq_b
                        len_a = length(seq_a)
                        len_b = length(seq_b)
                        
                        # Fetch 3' ends of both strands participating in the alignment
                        if full_align.aln.a.aln.anchors[end].op == OP_SEQ_MATCH && full_align.aln.a.aln.anchors[1].op == OP_START
                            aln_3prime_end_seq_a = full_align.aln.a.aln.anchors[end].seqpos
                            aln_5prime_end_seq_b = 1 + len_b - full_align.aln.a.aln.anchors[end].refpos
                            
                            if len_a == aln_3prime_end_seq_a && aln_5prime_end_seq_b > 1
                                seqa_extend = true
                                a_ex += 1
                                n_bases_2_fetch_seq_b = len_b - full_align.aln.a.aln.anchors[end].refpos
                                extension_a = reverse_complement(seq_b[1:n_bases_2_fetch_seq_b])
                                seq_a_new = seq_a * extension_a
                                sites[candidate].seq = seq_a_new
                                sites[candidate] = update_site(sites[candidate])
                                extension_events_link += 1
                                extension_events_this_cycle += 1
                                update_strand_tracker!(a_pol_f, candidate, string(seq_a_new))
                                update_strand_tracker!(a_pol_r, candidate, string(seq_a_new))
                                update_strand_tracker!(b_pol_f, candidate, string(seq_a_new))
                                update_strand_tracker!(b_pol_r, candidate, string(seq_a_new))
                                update_strand_tracker!(tracker_400, candidate, length(sites[candidate].seq), partner_coordinates)
                                
                                if length(seq_a_new) > parameter.length_thresh
                                    push!(long_sites, candidate)
                                end
                            else
                                seqa_extend = false
                            end
                            
                            aln_3prime_end_seq_b = len_b - full_align.aln.a.aln.anchors[1].refpos
                            aln_5prime_end_seq_a = 1 + full_align.aln.a.aln.anchors[1].seqpos
                            
                            if len_b == aln_3prime_end_seq_b && aln_5prime_end_seq_a > 1
                                seqb_extend = true
                                b_ex += 1
                                n_bases_2_fetch_seq_a = full_align.aln.a.aln.anchors[1].seqpos
                                extension_b = reverse_complement(seq_a[1:n_bases_2_fetch_seq_a])
                                seq_b_new = seq_b * extension_b
                                sites[partner_coordinates].seq = seq_b_new
                                sites[partner_coordinates] = update_site(sites[partner_coordinates])
                                update_strand_tracker!(a_pol_f, partner_coordinates, string(seq_b_new))
                                update_strand_tracker!(a_pol_r, partner_coordinates, string(seq_b_new))
                                update_strand_tracker!(b_pol_f, partner_coordinates, string(seq_b_new))
                                update_strand_tracker!(b_pol_r, partner_coordinates, string(seq_b_new))
                                update_strand_tracker!(tracker_400, partner_coordinates, length(sites[partner_coordinates].seq), candidate)
                                extension_events_link += 1
                                extension_events_this_cycle += 1
                                
                                if length(seq_b_new) > parameter.length_thresh
                                    push!(long_sites, partner_coordinates)
                                end
                            else
                                seqb_extend = false
                            end
                        end
                    else
                        reject_count += 1
                    end
                else
                    reject_count += 1
                end
            else
                reject_count += 1
            end
        end
        
        # Print cycle information
        println("end of cycle: ", cycle, "number of long sites: ", length(long_sites), "cumulative extensions: ", extension_events_link)
        println("    number of annealing rejections this cycle: ", reject_count)
        println("    number of attempted interactions this cycle: ", i_no)
        println("    number of annealing events this cycle: ", annealings)
        println("    number of extension events this cycle: ", extension_events_this_cycle)
    end
    
    println("a extend:", a_ex)
    println("b extend:", b_ex)
    
    return long_sites, ex_vec
end


###################################################
# Main Function
###################################################


function PCR_via_bridge(out_param, all) 
    

    global a_pol_f = Strand_tracker(all.primers[1][1], all.primers[1][2], all.primers[1][3])
    global a_pol_r = Strand_tracker(all.primers[2][1], all.primers[2][2], all.primers[2][3])
    global b_pol_f = Strand_tracker(all.primers[3][1], all.primers[3][2], all.primers[3][3])
    global b_pol_r = Strand_tracker(all.primers[4][1], all.primers[4][2], all.primers[4][3])
    global tracker_400 = Strand_tracker(all.primers[5][1], all.primers[5][2], all.primers[5][3])

    get_thermodynamic_parameters()

    sites, long_sites, extension_vector, polonies_strand, extension_events, pc_seq,site_positions,all_distance = @time PCR_Run(all)

    polony_data = polony_size(pc_seq, polonies_strand)
    
    long_sites = ecori_cut(long_sites, "GAATTC", 1, all, sites)
    
    long_sites, ex_vec = UMI_crossover(long_sites,all,sites,site_positions,all_distance)

    if out_param.grow_main == true
       
        plot_graph(all.n_cycles, extension_vector, all.Output_file_name, out_param.graph_table)

    end

    if out_param.grow_cross == true
        plot_graph(all.n_cycles_2, ex_vec, all.Output_file_name*"_link",out_param.graph_table)
    end 

    writing_output(all, polony_data ,extension_events,pc_seq,polonies_strand,long_sites, sites, FASTQ_file=out_param.FASTQ_file, Polony=out_param.Polon, General = out_param.General)

    
    reset!([a_pol_f,a_pol_r,b_pol_f,b_pol_r,tracker_400])

end

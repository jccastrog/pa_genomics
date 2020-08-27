#!/usr/bin/env julia
#=
@name: simulate_operonEvolutionMarkov.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 27-Aug-2020
@version: 1.0
@license: GNU General Public License v3.0.
please type "./simulate_operonEvolutionMarkov.jl.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
using StatsBase;
using Distributions;
using LinearAlgebra;
using ProgressMeter;
# 1.2 Define functions ===================================================#
"""
parse_commandline()
Parse arguments in the command line and output help menu if required 
arguments are not provided.

#Examples

```
julia  ./simulate_operonEvolutionMarkov.jl
required option --genome_size was not provided
usage: simulate_operonEvolutionMarkov.jl -s GENOME_SIZE -n NUM_GENOMES
                        -g GAIN_PROBABILITY -l LOSS_PROBABILITY
                        -r REARRANGEMENT_PROBABILITY [-e NUM_STEPS]
                        [-b BURN_IN] [-o OUTPUT]
julia ./simulate_operonEvolutionMarkov.jl -h
usage: simulate_operonEvolutionMarkov.jl -s GENOME_SIZE -n NUM_GENOMES
                        -g GAIN_PROBABILITY -l LOSS_PROBABILITY
                        -r REARRANGEMENT_PROBABILITY [-e NUM_STEPS]
                        [-b BURN_IN] [-o OUTPUT] [-h]

optional arguments:
  -s, --genome_size GENOME_SIZE
                        Size of the genome in number of genes.
                        (default: 5)
  -n, --num_genomes NUM_GENOMES
                        Number of genomes to be simulated. (default:
                        10)
  -g, --gain_probability GAIN_PROBABILITY
                        Probability of occurence for a gene gain
                        event.
  -l, --loss_probability LOSS_PROBABILITY
                        Probability of occurence for a gene loss
                        event.
  -r, --rearrangement_probability REARRANGEMENT_PROBABILITY
                        Probability of occurence for a genome
                        rearrangement event.
  -e, --num_steps NUM_STEPS
                        Number of links in the Markov Chain. (default:
                        100)
  -b, --burn_in BURN_IN
                        Number of steps in the Markov chain ignored
                        before computation
  -o, --output OUTPUT   Output file with the distances and mutual
                        information values (default:
                        "mi_distances.tsv")
  -h, --help            show this help message and exit

"""
function parseCommandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--genome_size", "-s"
			help = "Size of the genome in number of genes."
			required = true
			default = 5
		"--num_genomes", "-n"
			help = "Number of genomes to be simulated."
			required = true
			default = 10
		"--gain_probability", "-g"
			help = "Probability of occurence for a gene gain event."
			required = true
		"--loss_probability", "-l"
			help = "Probability of occurence for a gene loss event."
			required = true
		"--rearrangement_probability", "-r"
			help = "Probability of occurence for a genome rearrangement event."
			required = true
		"--num_steps", "-e"
			help = "Number of links in the Markov Chain."
			required = false
			default = 1000
		"--burn_in", "-b"
			help = "Number of steps in the Markov chain ignored before computation"
			required = false
			default = 100
		"--output", "-o"
			help = "Output file with the distances and mutual information values"
			required = false
			default = "mi_distances.tsv"
	end
	return parse_args(s)
end
"""
"""
function simulateGenomes(num_genomes::Int64, genome_size::Int64)
	genome_dict = Dict();
	ogs_freqs = Dict();
	ogs_dist = Beta(0.5,0.5);
	sample_size = num_genomes*genome_size;
	ogs_sample_float = rand(ogs_dist, sample_size) * num_genomes;
	ogs_sample_int = [ceil.(x) for x in ogs_sample_float];
	strAdd(simp_str::String, simp_int::Int64) = try
		"$simp_str$simp_int";
	catch
		"";
	end;
	og_ids = [strAdd("og_", i) for i in 1:sample_size];
	genome_ids = [strAdd("genome_", i) for i in 1:num_genomes];
	for i in 1:length(og_ids)
		ogs_freqs[og_ids[i]] = Int(ogs_sample_int[i]);
	end
	for key in collect(keys(ogs_freqs))
		present_in = sample(genome_ids, ogs_freqs[key], replace = false);
		k = 0;
		for i in present_in
			k += 1;
			if haskey(genome_dict, i)
						push!(genome_dict[i], key);
					else
				genome_dict[i] = [key];
			end
		end
	end
	return(genome_dict, og_ids)
end
"""
"""
function calcAverageGenomeSize(genome_dict::Dict)
	size_arr =  Int64[];
	for key in collect(keys(genome_dict))
		genome_element = genome_dict[key];
		loc_size = length(genome_element);
		push!(size_arr, loc_size);
	end
	mean_size = mean(size_arr);
	return(mean_size);
end
"""
"""
function calcMinGenomeSize(genome_dict::Dict)
	size_arr =  Int64[];
	for key in collect(keys(genome_dict))
		genome_element = genome_dict[key];
		loc_size = length(genome_element);
		push!(size_arr, loc_size);
	end
	min_size = minimum(size_arr);
	return(min_size);
end
"""
"""
function parseDict2Matrix(genome_dict::Dict, og_ids::Array)
	num_genomes = length(genome_dict);
	num_ogs = length(og_ids);
	pan_matrix = zeros(num_genomes, num_ogs);
	for genome in collect(keys(genome_dict))
		genome_num_id = parse(Int64, split(genome, "_")[2]);
		for og in og_ids
			og_num_id = parse(Int64, split(og,"_")[2]);
			if og in genome_dict[genome]
				pan_matrix[genome_num_id, og_num_id] = 1;
			else
				pan_matrix[genome_num_id, og_num_id] = 0;
			end
		end
	end
	return(pan_matrix)
end
"""
"""
function calcArrayDist(genome_element::Array, og_ids::Array)
	genome_size = length(genome_element);
	distance = 0;
	if length(og_ids)>2
		error("More than 2 ogs used to calculate distance");
	else
		indexes = Int64[];
		for i in og_ids
			loc_index = findall(x -> x==i, genome_element)[1];
			push!(indexes, loc_index);
		end
		if (indexes[1] > indexes[2])
			distance1 = indexes[1] - indexes[2] - 1;
			distance2 = genome_size - indexes[1] + indexes[2] - 1;
			distance = min(distance1, distance2);
		else
			distance1 = indexes[2] - indexes[1] - 1;
			distance2 = genome_size - indexes[2] + indexes[1] - 1;
			distance = min(distance1, distance2);
		end
		return(distance)
	end
end
"""
"""
function calcAverageDistance(genome_dict::Dict, og_ids::Array)
	dict_distances = Int64[];
	for key in collect(keys(genome_dict))
		genome_element = genome_dict[key];
		if (og_ids[1] in genome_element) && (og_ids[2] in genome_element)
			loc_dist = calcArrayDist(genome_element, og_ids);
			push!(dict_distances, loc_dist)
		end
	end
	mean_distance = mean(dict_distances);
	return(mean_distance);
end
"""
"""
function calcMutualInfo(vecX::Array, vecY::Array)
	lenVec = length(vecX);
	px1 = sum(vecX)/lenVec;
	px0 = 1 - px1;
	py1 = sum(vecY)/lenVec;
	py0 = 1 - py1;
	if (px1 == 1 || px0 == 1 || py1 == 1 || py0 == 1)
		error("Failed to calculate mutual information, one
		or more elements being compared have no variation.")

	else
		xy11 = 0;
		xy10 = 0;
		xy00 = 0;
		xy01 = 0;
		for i in 1:lenVec
			if vecX[i] == 1
				if vecY[i] == 1
					xy11 += 1;
				elseif vecY[i] == 0
					xy10 += 1;
				end
			elseif vecX[i] == 0
				if vecY[i] == 1
					xy01 += 1;
				elseif vecY[i] == 0
					xy00 += 1;
				end
			end
		end
		pxy11 = xy11/lenVec;
		pxy10 = xy10/lenVec;
		pxy00 = xy00/lenVec;
		pxy01 = xy01/lenVec;
		t1 = pxy11*log((pxy11)/(px1*py1));
		t2 = pxy10*log((pxy10)/(px1*py0));
		t3 = pxy01*log((pxy01)/(px0*py1));
		t4 = pxy00*log((pxy00)/(px0*py0));
		if isnan(t1)
			ixy = t2 + t3 + t4;
		elseif isnan(t2)
			ixy = t1 + t3 + t4;
		elseif isnan(t3)
			ixy = t1 + t2 + t4;
		elseif isnan(t4)
			ixy = t1 + t2 + t3;
		else
			ixy = t1 + t2 + t3 + t4;
		end
		return(ixy);
	end
end
"""
Returns only a genome with an added set of OGs
MAKES CHANGE TO DICTIONARY ITSELF
"""
function addOGs(genome_id::String, og_ids::Array, genome_dict::Dict, replacement_size ::Int64)
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_length = length(genome_loc);
	if replacement_size >= round(genome_length/10)
		error("Attempt to add a fraction of the genome larger than 10%");
	else
		og_ids_new = filter(x -> x ∉ genome_loc, og_ids);
		if length(og_ids_new) > 0
			og_insert = sample(og_ids_new, replacement_size);
			insertion_site = sample(1:genome_length, 1)[1];
			for i in og_insert
				genome_loc = insert!(genome_loc, insertion_site, i);
			end
		else
			error("Failed to add gene, target genome has all possible OGS")
		end
	end
	return(genome_loc);
end
"""
"""
function addSelection(genome_id::String,  og_ids::Array, genome_dict::Dict, replacement_size::Int64, fit_mat::Array, all_kd::Array)
	strRemove(simp_str::String) = try
		parse(Int64, split(simp_str, "_")[2]);
	catch
		"";
	end
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_loc_int = [strRemove(x) for x in genome_loc];
	genome_length = length(genome_loc);
	og_ids_new = [];
	fit_submat = fit_mat[genome_loc_int, genome_loc_int];
	kd = all_kd[genome_loc_int];
	insert_weights_single = fit_submat * kd;
	joint_weights = Float64[];
	joint_weights_fragment = [];
	fragment_coordinates = [];
	genome_length = length(genome_loc);
	og_ids_new = filter(x -> x ∉ genome_loc, og_ids);
	for i in 1:genome_length
		loc_joint = [];
		if (i + 1) > genome_length
			j = i + 1 - genome_length;
			loc_joint = vcat(insert_weights_single[i:genome_length], insert_weights_single[1:j])
			loc_index = 1;
			loc_weights = sum(loc_joint);
			push!(joint_weights, loc_weights);
			push!(fragment_coordinates, loc_index);
		else
			j = i + 1;
			loc_joint = insert_weights_single[i:j];
			loc_index = i+1;
			loc_weights = sum(loc_joint);
			push!(joint_weights, loc_weights);
			push!(fragment_coordinates, loc_index);
		end
	end
	insert_weights_normal = joint_weights/sum(joint_weights);
	insert_probs = Weights(insert_weights_normal);
	if length(og_ids_new) > 0
		og_insert = sample(og_ids_new, replacement_size);
		insert_site = sample(fragment_coordinates, insert_probs);
		for i in og_insert
			genome_loc = insert!(genome_loc, insert_site, i);
		end
	else
		error("Failed to add gene, target genome has all possible OGS")
	end
	return(genome_loc);
end
"""
Returns only a genome with a removed set of OGs
"""
function removeOGs(genome_id::String, genome_dict::Dict, replacement_size ::Int64)
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_length = length(genome_loc);
	if replacement_size >= round(genome_length/10)
		error("Attempt to remove a fraction of the genome larger than 10%");
	else
		removal_site = sample(1:genome_length, 1)[1];
		try 
			removal_end = removal_site + replacement_size -1;
			filter!(x -> x ∉ genome_loc[removal_site:removal_end], genome_loc);
		catch
			removal_start = removal_site - replacement_size + 1;
			filter!(x -> x ∉ genome_loc[removal_start:removal_site], genome_loc);

		end
	end
	return(genome_loc);
end
"""
"""
function removeSelection(genome_id::String, genome_dict::Dict, replacement_size::Int64, fit_mat::Array, all_kd::Array)
	strRemove(simp_str::String) = try
		parse(Int64, split(simp_str, "_")[2]);
	catch
		"";
	end
	genome_id_int = strRemove(genome_id);
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_loc_int = [strRemove(x) for x in genome_loc];
	fit_submat = fit_mat[genome_loc_int, genome_loc_int];
	kd = all_kd[genome_loc_int];
	removal_weights_single = fit_submat * kd;
	genome_length = length(genome_loc);
	if replacement_size >= round(genome_length/10)
		error("Attempt to remove a fraction of the genome larger than 10%");
	else
		replacement_increment = replacement_size - 1;
		removal_weights_fragment = Float64[];
		fragment_coordinates = [];
		for i in 1:genome_length
			loc_frag = [];
			if (i + replacement_increment) > genome_length
				j = (i + replacement_increment) - genome_length;
				loc_frag = vcat(removal_weights_single[i:genome_length], removal_weights_single[1:j]);
				loc_index = vcat(i:genome_length, 1:j);
				loc_weights = sum(loc_frag);
				push!(removal_weights_fragment, loc_weights);
				push!(fragment_coordinates, loc_index);
			else
				j = i + replacement_increment;
				loc_frag = removal_weights_single[i:j];
				loc_index = i:j;
				loc_weights = sum(loc_frag);
				push!(removal_weights_fragment, loc_weights);
				push!(fragment_coordinates, loc_index);
			end
		end
		removal_weights_normal = removal_weights_fragment/sum(removal_weights_fragment);
		removal_probs = Weights(removal_weights_normal);
		removal_site = sample(fragment_coordinates, removal_probs);
		try 
			filter!(x -> x ∉ genome_loc[removal_site], genome_loc);
		catch
			error("Cannot remove requested fragment")
		end
	end
	return(genome_loc);
end
"""
"""
function shuffleGenome(genome_id::String, genome_dict::Dict, replacement_size ::Int64)
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_length = length(genome_loc);
	if replacement_size >= round(genome_length/10)
		error("Attempt to rearrange a fraction of the genome larger than 10%. Genome size is $genome_length,
		fragment size is $replacement_size");
	else
		rearrangement_site = sample(1:genome_length, 1)[1]
		fragment = [];
		try
			rearrangement_end = rearrangement_site + replacement_size - 1;
			fragment = genome_loc[rearrangement_site:rearrangement_end];
		catch 
			rearrangement_start = rearrangement_site - replacement_size + 1;
			fragment = genome_loc[rearrangement_start:rearrangement_site];
		end
		filter!(x -> x ∉ fragment, genome_loc);
		genome_length = length(genome_loc);
		insertion_site = sample(1:genome_length, 1)[1]
		for i in 1:replacement_size
			insert!(genome_loc, insertion_site, fragment[i]);
		end
	end
	return(genome_loc);
end
"""
"""
function shuffleSelection(genome_id::String, genome_dict::Dict, replacement_size::Int64, fit_mat::Array, all_kd::Array)
	strRemove(simp_str::String) = try
		parse(Int64, split(simp_str, "_")[2]);
	catch
		"";
	end
	genome_id_int = strRemove(genome_id);
	genome_loc = [];
	if haskey(genome_dict, genome_id)
		genome_loc = genome_dict[genome_id];
	else
		error("Genome id $genome_id does not exist in the current collection")
	end
	genome_loc_int = [strRemove(x) for x in genome_loc];
	fit_submat = fit_mat[genome_loc_int, genome_loc_int];
	kd = all_kd[genome_loc_int];
	shuffle_weights_single = fit_submat * kd;
	genome_length = length(genome_loc);
	joint_weights = Float64[];
	replacement_increment = replacement_size - 1;
	joint_weights_fragment = Float64[];
	shufle_weights_fragment = Float64[];
	fragment_coordinates = [];
	for i in 1:genome_length
		loc_joint = [];
		if (i + 1) > genome_length
			j = i + 1 - genome_length;
			loc_joint = vcat(shuffle_weights_single[i:genome_length], shuffle_weights_single[1:j])
			loc_index = i;
			loc_weights = sum(loc_joint);
			push!(joint_weights, loc_weights);
		else
			j = i + 1;
			loc_joint = shuffle_weights_single[i:j];
			loc_index = i;
			loc_weights = sum(loc_joint);
			push!(joint_weights, loc_weights);
		end
	end
	for i in 1:genome_length
		loc_frag = [];
		if (i + replacement_increment) > genome_length
				j = (i + replacement_increment) - genome_length;
				loc_frag = vcat(joint_weights[i], joint_weights[j]);
				loc_index = vcat(i:genome_length, 1:j); #HERE
				loc_weights = sum(loc_frag);
				push!(joint_weights_fragment, loc_weights);
				push!(fragment_coordinates, loc_index);
			else
				j = i + replacement_increment;
				loc_frag = vcat(joint_weights[i], joint_weights[j]);
				loc_index = i:j;
				loc_weights = sum(loc_frag);
				push!(joint_weights_fragment, loc_weights);
				push!(fragment_coordinates, loc_index);
		end
	end
	fragment_coordinates = append!(fragment_coordinates[2:end], [fragment_coordinates[1]]);
	joint_probs = Weights(joint_weights_fragment/sum(joint_weights_fragment));
	shuffle_fragment = sample(fragment_coordinates, joint_probs);
	fragment = genome_loc[shuffle_fragment];
	filter!(x -> x ∉ fragment, genome_loc);
	genome_length = length(genome_loc);
	genome_loc_int = [strRemove(x) for x in genome_loc];
	fit_submat = fit_mat[genome_loc_int, genome_loc_int];
	kd = all_kd[genome_loc_int];
	insert_weights_single = fit_submat * kd;
	insert_weights = Float64[]
	insertion_coordinates = [];
	for i in 1:genome_length
		loc_joint = [];
		if (i + 1) > genome_length
			j = i + 1 - genome_length;
			loc_joint = vcat(shuffle_weights_single[i:genome_length], shuffle_weights_single[1:j])
			loc_index = 1;
			loc_weights = sum(loc_joint);
			push!(insert_weights, loc_weights);
			push!(insertion_coordinates, loc_index);
		else
			j = i + 1;
			loc_joint = shuffle_weights_single[i:j];
			loc_index = i + 1;
			loc_weights = sum(loc_joint);
			push!(insert_weights, loc_weights);
			push!(insertion_coordinates, loc_index);
		end
	end
	insert_probs = Weights(insert_weights/sum(insert_weights))
	insertion_site = sample(insertion_coordinates, insert_probs)
	for i in 1:replacement_size
		insert!(genome_loc, insertion_site, fragment[i]);
	end
	return(genome_loc);
end
"""
"""
function noSelectionStep(og_ids::Array, genome_dict::Dict, replacement_size::Int64, probability_weights::Array)
	genome_set = collect(keys(genome_dict));
	genome_id = sample(genome_set, 1, replace = false)[1];
	event_select = sample([1,2,3], Weights(probability_weights));
	if event_select == 1
		try
			genome_loc = addOGs(genome_id, og_ids, genome_dict, replacement_size);
		catch
			"";
		end
	elseif event_select == 2
		try
			genome_loc = removeOGs(genome_id, genome_dict, replacement_size);
		catch
			"";
		end
	elseif event_select == 1
		try
			genome_loc = shuffleGenome(genome_id, genome_dict, replacement_size)
		catch
			"";
		end
	end
end
"""
"""
function selectionStep(og_ids::Array, genome_dict::Dict, replacement_size::Int64, probability_weights::Array, fit_mat::Array, all_kd::Array)
	genome_set = collect(keys(genome_dict));
	genome_id = sample(genome_set, 1, replace = false)[1];
	event_select = sample([1,2,3], Weights(probability_weights));
	if event_select == 1
		try
			genome_loc = addSelection(genome_id, og_ids, genome_dict, replacement_size, fit_mat, all_kd);
		catch
			"";
		end
	elseif event_select == 2
		try
			genome_loc = removeSelection(genome_id, genome_dict, replacement_size, fit_mat, all_kd);
		catch
			"";
		end
	elseif event_select == 3
		try
			genome_loc = shuffleSelection(genome_id, genome_dict, replacement_size, fit_mat, all_kd);
		catch
			"";
		end
	end
end
"""
"""
function simulateMarkovProcess(genome_dict::Dict, og_ids::Array, selection::Bool, probability_weights::Array, og_track_int::Array, num_steps::Int64, burn_in::Int64)
	strAdd(simp_str::String, simp_int::Int64) = try
		"$simp_str$simp_int";
	catch
		"";
	end;
	genomes_set = collect(keys(genome_dict));
	num_genomes = length(genomes_set);
	num_ogs = length(og_ids);
	og_track_id = [strAdd("og_", og_track_int[i]) for i in 1:2];
	initial_matrix = parseDict2Matrix(genome_dict, og_ids)
	vecX = initial_matrix[1:num_genomes, og_track_int[1]];
	vecY = initial_matrix[1:num_genomes, og_track_int[2]];
	MI_0 = 0;
	try
		MI_0 = calcMutualInfo(vecX, vecY);
	catch
		MI_0 = 0;
	end
	MI_vec = Float64[MI_0];
	dist_0 = calcAverageDistance(genome_dict, og_track_id);
	dist_vec = [dist_0];
	average_genome_size = calcAverageGenomeSize(genome_dict);
	min_genome_size = calcMinGenomeSize(genome_dict);
	genome_size_vec = Float64[average_genome_size];
	if selection
		all_kd = ones(num_ogs)*1;
		fit_mat = 1 ./(abs.(rand(Normal(0,0.1), num_ogs, num_ogs)));
		#fit_mat = reshape(ones(num_ogs^2), (num_ogs, num_ogs)); #DEBUG
		[fit_mat[x,x] = 0 for x in 1:num_ogs];
		[fit_mat[x,y] = fit_mat[y,x] for x in 1:num_ogs for y in 1:num_ogs];
		max_fitness_value = maximum(fit_mat);
		min_fitness_value = minimum(fit_mat);
		fit_mat[og_track_int[1], og_track_int[2]] = 0; #DEBUG
		fit_mat[og_track_int[2], og_track_int[1]] = 0; #DEGUB
		fit_mat[og_track_int[1], :] = zeros(num_ogs); #DEBUG
		fit_mat[:, og_track_int[1]] = zeros(num_ogs); #DEGUB
		fit_mat[og_track_int[2], :] = zeros(num_ogs); #DEGUB
		fit_mat[:, og_track_int[2]] = zeros(num_ogs); #DEGUB
		fitness_value = fit_mat[og_track_int[1], og_track_int[2]];
		@showprogress 1 "Simulating selection process......." for i in 1:num_steps
			avg_genome_size_i = calcAverageGenomeSize(genome_dict);
			replacement_size_float = (rand(1)[1] * (min_genome_size/10)) - (min_genome_size/50);
			replacement_size_int = Int(floor(replacement_size_float));
			if replacement_size_int <= 0
				replacement_size_int = 1;
			end
			e = exp(1);
			lambda = 1.44*150;
			rearrangement_probability = probability_weights[3];
			gain_probability = round((1-rearrangement_probability) * e^(-avg_genome_size_i/lambda), digits = 5);
			loss_probability = round(1 - (rearrangement_probability + gain_probability), digits = 5);
			probability_weights = [gain_probability, loss_probability, rearrangement_probability];
			selectionStep(og_ids, genome_dict, replacement_size_int, probability_weights, fit_mat, all_kd);
			if i >= burn_in
				step_matrix = parseDict2Matrix(genome_dict, og_ids);
				vecX = step_matrix[1:num_genomes, og_track_int[1]];
				vecY = step_matrix[1:num_genomes, og_track_int[2]];
				MI_i = 0;
				try
					MI_i = calcMutualInfo(vecX, vecY);
				catch
					MI_i = 0;
				end
				dist_i = calcAverageDistance(genome_dict, og_track_id);
				push!(genome_size_vec, avg_genome_size_i);
				push!(MI_vec, MI_i);
				push!(dist_vec, dist_i);
			end
		end
		return(MI_vec, dist_vec, fitness_value, max_fitness_value, min_fitness_value, genome_size_vec);
	else
		@showprogress 1 "Simulating non-selection process..." for i in 1:num_steps
			avg_genome_size_i = calcAverageGenomeSize(genome_dict);
			replacement_size_float = (rand(1)[1] * (min_genome_size/10)) - (min_genome_size/50);
			replacement_size_int = Int(floor(replacement_size_float));
			if replacement_size_int <= 0
				replacement_size_int = 1;
			end
			e = exp(1);
			lambda = 1.44*150;
			rearrangement_probability = probability_weights[3];
			gain_probability = round((1-rearrangement_probability) * e^(-avg_genome_size_i/lambda), digits = 5);
			loss_probability = round(1 - (rearrangement_probability + gain_probability), digits = 5);
			probability_weights = [gain_probability, loss_probability, rearrangement_probability];
			noSelectionStep(og_ids, genome_dict, replacement_size_int, probability_weights);
			if i >= burn_in
				step_matrix = parseDict2Matrix(genome_dict, og_ids);
				vecX = step_matrix[1:num_genomes, og_track_int[1]];
				vecY = step_matrix[1:num_genomes, og_track_int[2]];
				MI_i = 0;
				try
					MI_i = calcMutualInfo(vecX, vecY);
				catch
					MI_i = 0;
				end
				dist_i = calcAverageDistance(genome_dict, og_track_id);
				push!(MI_vec, MI_i)
				push!(dist_vec, dist_i);

			end
		end
		return(MI_vec, dist_vec);
	end
end
# 1.3 Initialize variables ===============================================#
# 1.3.1 Parser variables #
write(stdout, "Setting up!\n")
parsed_args = parseCommandline();
genome_size = parse(Int64, parsed_args["genome_size"]);
num_genomes = parse(Int64, parsed_args["num_genomes"]);
gain_probability = parse(Float64, parsed_args["gain_probability"]);
loss_probability = parse(Float64, parsed_args["loss_probability"]);
rearrangement_probability = parse(Float64, parsed_args["rearrangement_probability"]);
num_steps = parse(Int64, parsed_args["num_steps"]);
burn_in = parse(Int64, parsed_args["burn_in"]);
output = parsed_args["output"];
# 1.3.2 Gene collection variables #
genome_dict, og_ids = simulateGenomes(num_genomes,genome_size);
initial_mean_genome_size = calcAverageGenomeSize(genome_dict);
genome_dict_sel = deepcopy(genome_dict);
# 1.3.3 Gene tracking variables #
probability_weights = [gain_probability, loss_probability, rearrangement_probability];
genome_dict, og_ids = simulateGenomes(num_genomes,genome_size);
parse_mat = parseDict2Matrix(genome_dict, og_ids);
og_track_int = Int64[];
for i in 1:size(parse_mat)[2]
	if sum(parse_mat[1:size(parse_mat)[1], i]) != size(parse_mat)[1]
		if length(og_track_int) < 2
			push!(og_track_int, i);
		end
	end
end
###=== 2.0 Calculate mutual information based on a Markov simulation ===###
# 2.1 Mutual information and distance ====================================#
# 2.1.1 No Selection #
mi, dist = simulateMarkovProcess(genome_dict, og_ids, false, probability_weights, og_track_int, num_steps, burn_in);
# 2.1.2 Selection #
mi_sel, dist_sel, fitness_value, max_fitness_value, min_fitness_value, genome_size_vec = simulateMarkovProcess(genome_dict_sel, og_ids, true, probability_weights, og_track_int, num_steps, burn_in);
###========================== 3.0 Write results =========================###
# 3.1 Write No Selection ==================================================#
no_selection_ext = "_noSelection.tsv";
out_name = "$output$no_selection_ext";
out_stream = open(out_name, "w");
write(out_stream, "Step\tMI\tDistance\n");
for i in 1:length(mi)
	write(out_stream, "$i\t$(mi[i])\t$(dist[i])\n");
end
close(out_stream);
# 3.2 Write Selection ====================================================#
selection_ext = "_selection.tsv";
out_name = "$output$selection_ext";
out_stream = open(out_name, "w");
write(out_stream, "Step\tMI\tDistance\tGenome_size\n");
for i in 1:length(mi_sel)
	write(out_stream, "$i\t$(mi_sel[i])\t$(dist_sel[i])\t$(genome_size_vec[i])\n");
end
close(out_stream);
# 3.3 Write Miscelaneous stats ===========================================#
final_mean_genome_size = calcAverageGenomeSize(genome_dict);
final_mean_genome_size_selection = calcAverageGenomeSize(genome_dict_sel);
stats_ext = "_stats.txt";
out_name = "$output$stats_ext";
out_stream = open(out_name, "w");
write(out_stream, "Initial Genome Size:\t$initial_mean_genome_size\n");
write(out_stream, "Final Genome Size:\t$final_mean_genome_size\n");
write(out_stream, "Final Genome Size Selection:\t$final_mean_genome_size_selection\n");
write(out_stream, "Fitness Value:\t$fitness_value\n");
write(out_stream, "Maximum Fitness Value:\t$max_fitness_value\n");
write(out_stream, "Minimum Fitness Value:\t$min_fitness_value\n");
###=====================================================================###
write(stdout, "SUCCESS!!!\n");
###=====================================================================###
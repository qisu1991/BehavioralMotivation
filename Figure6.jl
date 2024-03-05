# simulate the stationary distribution 
# using Monte Carlo simulations
################################################################################

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using Random
using IterativeSolvers
using SparseArrays
using StatsBase: sample
rng = MersenneTwister(1234)

# generate random regular network of size N and node degree k
function random_regular(N::Int64, k::Int64)
    mat_seq = spzeros(Int64, N, N)
    if k == 2
        for i = 1:N
            h = ((i+1)%N > 0) ? ((i+1)%N) : N
            mat_seq[i,h] = 1
            mat_seq[h,i] = 1
        end
    end

    if k > 2
        net = 1  # network is being generated
        while net == 1
            mat_seq = spzeros(Int64, N, N)
            singular_unsolved = 0
            for i = 1:N
                deg_seq = sum(mat_seq, dims = 2)
                deg_seq = deg_seq[:, 1]
                deficient_seq = findall(!isequal(k), deg_seq)
                if length(deficient_seq) == 0
                    net = 0   # network is generated
                    break
                end
                host = deficient_seq[1]  #select a node to build edges
                deficient_seq = setdiff(deficient_seq, host)
                deficient_deg = k - sum(mat_seq[host, :])
                while deficient_deg > 0
                    # singular case, severing an existing edge and rewiring this sick node
                    if length(deficient_seq) == 0
                        nohost = findall(iszero, mat_seq[host,:])
                        nohost = setdiff(nohost, host)
                        mat_seq2 = mat_seq[nohost,nohost]
                        inds = findall(!iszero, mat_seq2)
                        aa = getindex.(inds, 1)
                        bb = getindex.(inds, 2)
                        rr = rand(1:length(inds))
                        a = nohost[aa[rr]]
                        b = nohost[bb[rr]]
                        mat_seq[a, b] = 0
                        mat_seq[b, a] = 0
                        mat_seq[host, a] = 1
                        mat_seq[a, host] = 1
                        mat_seq[host, b] = 1
                        mat_seq[b, host] = 1
                        deficient_deg = deficient_deg - 2
                        continue
                    end
                    # singular case, severing an existing edge and rewiring this sick node
                    # deficient_seq = [deficient_seq[i] for i = 1:length(deficient_seq)]
                    deficient_seq2 = setdiff(deficient_seq, findall(!iszero, mat_seq[host,:]))
                    if length(deficient_seq2) == 0
                        is_solved = 0
                        for node in deficient_seq
                            deficient_seq3 = setdiff(findall(!iszero, mat_seq[node,:]), findall(!iszero, mat_seq[host,:]))
                            if length(deficient_seq3) == 0
                                continue
                            end
                            guest = rand(deficient_seq3)
                            mat_seq[host, guest] = 1
                            mat_seq[guest, host] = 1
                            mat_seq[node, guest] = 0
                            mat_seq[guest, node] = 0
                            deficient_deg = deficient_deg - 1
                            is_solved = 1
                            break
                        end
                        if is_solved == 0
                            nohost = findall(iszero, mat_seq[host,:])
                            nohost = setdiff(nohost, host)
                            mat_seq2 = mat_seq[nohost,nohost]
                            inds = findall(!iszero, mat_seq2)
                            if length(inds) == 0
                                singular_unsolved = 1
                                break
                            end
                            aa = getindex.(inds, 1)
                            bb = getindex.(inds, 2)
                            rr = rand(1:length(inds))
                            a = nohost[aa[rr]]
                            b = nohost[bb[rr]]
                            mat_seq[a, b] = 0
                            mat_seq[b, a] = 0
                            mat_seq[host, a] = 1
                            mat_seq[a, host] = 1
                            mat_seq[rand(deficient_seq), b] = 1
                            mat_seq[b, deficient_seq(deficient_seq)] = 1
                            deficient_deg = deficient_deg - 1
                        end
                        continue
                    end
                    # normal case
                    guest = rand(deficient_seq)
                    mat_seq[host,guest] = 1
                    mat_seq[guest,host] = 1
                    deficient_deg = deficient_deg - 1
                    deficient_seq = setdiff(deficient_seq, guest)
                end
                if singular_unsolved == 1
                    break
                end
            end
        end
    end

    M = [mat_seq[i,j] for i = 1:N, j = 1:N]
    M = is_linked(M)

    inds = findall(!iszero, M)
    a = getindex.(inds, 1)
    b = getindex.(inds, 2)
    edge_seq = zeros(Int64, length(b), 3)
    edge_seq[:,1] = a
    edge_seq[:,2] = b
    edge_seq[:,3] = ones(Int64, length(a))

    return edge_seq
end

# judge if a network is connected and make adjustment if not connected
function is_linked(M::Array{Int64,2})
    N = size(M,2)
    seq_node = [i for i = 1:N]

    while true
        seed = rand(1:N)
        # node sequence to be checked
        seq_unchecked = setdiff(seq_node, seed)
        # node sequence checked
        seq_checked = [seed]
        # node sequence as seeds
        seq_seed = findall(!iszero, M[seed,:])
        # node sequence connected
        seq_connected = vcat(seq_checked, seq_seed)
        seq_unchecked = setdiff(seq_unchecked, seq_seed)

        while length(seq_seed) > 0
            host = seq_seed[1]
            seq_seed = setdiff(seq_seed, host)
            seq_checked = vcat(seq_checked, host)
            seq = findall(!iszero, M[host,:])
            seq = setdiff(seq, seq_connected)
            seq_connected = vcat(seq_connected, seq)
            seq_seed = vcat(seq_seed, seq)
            seq_unchecked = setdiff(seq_unchecked, seq)
        end

        if length(seq_checked) < N/2 continue end
        if length(seq_checked) < N
            node_unchecked1 = seq_unchecked[1]
            neg_seq = findall(!iszero, M[node_unchecked1,:])
            # node_unchecked1 is an isolated node
            if length(neg_seq) == 0
                inds = findall(!iszero, M[seq_checked, :])
                a = getindex.(inds, 1)
                b = getindex.(inds, 2)
                c = rand(1:length(a))
                node_checked1 = seq_checked[a[c]]
                node_checked2 = b[c]
                M[node_unchecked1, node_checked1] = 1
                M[node_checked1, node_unchecked1] = 1
                M[node_checked1, node_checked2] = 0
                M[node_checked2, node_checked1] = 0
                continue
            end

            # node_unchecked1 is not an isolated node
            # and not being connected to seq_connected yet
            node_unchecked2 = rand(neg_seq)
            inds = findall(!iszero, M[seq_checked, :])
            a = getindex.(inds, 1)
            b = getindex.(inds, 2)
            c = rand(1:length(a))
            node_checked1 = seq_checked[a[c]]
            node_checked2 = b[c]
            M[node_unchecked1, node_checked1] = 1
            M[node_checked1, node_unchecked1] = 1
            M[node_unchecked2, node_checked2] = 1
            M[node_checked2, node_unchecked2] = 1
            M[node_checked1, node_checked2] = 0
            M[node_checked2, node_checked1] = 0
            M[node_unchecked1, node_unchecked2] = 0
            M[node_unchecked2, node_unchecked1] = 0
        end
        if length(seq_checked) == N break end
    end

    return M
end

# get the trajectory by Monte Carlo Simulations
# output one's lambda and alpha throughtout the evolutionary process
# on networks
function lambda_alpha_trajectory(edge_seq::Array{Int64,2}, b::Float64, c::Float64, lambda_alpha::Array{Float64,2}, mutation::Float64, generation::Int64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg

    # initial action configuration
    action =  rand(0:1,N)
    generation_action = 10
    sample_interval = 10000
    lambda_alpha_output = zeros(Float64, Int64(generation/sample_interval), 3)

    for g = 1:generation
        payoff_g = zeros(Float64, N)
        action_g = 0.0

        # action updating for generation_action times 
        for ga = 1:generation_action
            payoff = zeros(Float64, N)
            action_g = action_g+sum(action)

            # calculate the payoff
            payoff = P*action*b-action*c
            payoff_g = payoff_g + payoff

            # action updating
            for i = 1:N 
                probability = 1/(1+exp(-lambda_alpha[i,1]*(payoff[i]-lambda_alpha[i,2])))
                if rand() < probability
                    action[i] = 1 
                else
                    action[i] = 0
                end
            end
        end

        payoff_g = payoff_g/generation_action

        # selecting nodes to be updated
        sample_number = 1
        node_updated = sample(1:N, sample_number, replace=false)
        lambda_alpha_updated = lambda_alpha[node_updated,:]
        for i = 1:sample_number 
            node_competitor = findall(!iszero, P[node_updated[i],:])
            payoff_competitor = P[node_updated[i],:].*payoff_g
            threshold = rand(1)*sum(payoff_competitor)
            sum1 = 0
            for j = 1:length(node_competitor)
                sum1 = sum1 + payoff_competitor[node_competitor[j]]
                if sum1 >= threshold[1]
                    if rand() >= mutation
                        lambda_alpha_updated[i,:] = lambda_alpha[node_competitor[j],:]
                    else
                        lambda_alpha_updated[i,1] = lambda_alpha[node_competitor[j],1] + 0.01*(2*rand()-1)
                        lambda_alpha_updated[i,2] = lambda_alpha[node_competitor[j],2] + 0.1*(2*rand()-1)
                    end
                    break
                end
            end
        end
        lambda_alpha[node_updated,:] = lambda_alpha_updated[:,:]

        if g%sample_interval == 0
            lambda_alpha_output[Int64(g/sample_interval),1:2] = sum(lambda_alpha, dims = 1)/N
            lambda_alpha_output[Int64(g/sample_interval),3] = action_g/(N*generation_action)
        end
    end

    return lambda_alpha_output
end


N = 100  # network size
k = 6    # node degree
edge_seq = zeros(Int64,k*N,3)   #edge list of networks
edge_seq = random_regular(N, k)
mutation = 0.01   # mutation rate of behavioral motivations 
generation = 50000000   # generations of each simulation

# b/c < (b/c)*
c = 1.0
b = 2.0
lambda_alpha_list = [-0.05 0.25;
                     -0.05 0.75;
                     0.05 0.25;
                     0.05 0.75]

lambda_alpha = zeros(Float64,N,2)
for ll = 1:300
        for i = 1:size(lambda_alpha_list,1)
                lambda_alpha[:,1] = ones(Float64,N)*lambda_alpha_list[i,1]
                lambda_alpha[:,2] = ones(Float64,N)*lambda_alpha_list[i,2]

                lambda_alpha_output = lambda_alpha_trajectory(edge_seq, b, c, lambda_alpha, mutation, generation)

                g = open("lambda_belowthreshold.txt","a")
                writedlm(g, [lambda_alpha_list[i,1] transpose(lambda_alpha_output[:,1])])
                close(g)

                g = open("alpha_belowthreshold.txt","a")
                writedlm(g, [lambda_alpha_list[i,2] transpose(lambda_alpha_output[:,2])])
                close(g)

                g = open("xC_belowthreshold.txt","a")
                writedlm(g, transpose(lambda_alpha_output[:,3]))
                close(g)
        end
end


# b/c > (b/c)*
c = 0.2
b = 2.0
lambda_alpha_list = [-0.05 0.6;
                     -0.05 1.2;
                     0.05 0.6;
                     0.05 1.2]

lambda_alpha = zeros(Float64,N,2)
for ll = 1:300
        for i = 1:size(lambda_alpha_list,1)
                lambda_alpha[:,1] = ones(Float64,N)*lambda_alpha_list[i,1]
                lambda_alpha[:,2] = ones(Float64,N)*lambda_alpha_list[i,2]

                lambda_alpha_output = lambda_alpha_trajectory(edge_seq, b, c, lambda_alpha, mutation, generation)

                g = open("lambda_abovethreshold.txt","a")
                writedlm(g, [lambda_alpha_list[i,1] transpose(lambda_alpha_output[:,1])])
                close(g)

                g = open("alpha_abovethreshold.txt","a")
                writedlm(g, [lambda_alpha_list[i,2] transpose(lambda_alpha_output[:,2])])
                close(g)

                g = open("xC_abovethreshold.txt","a")
                writedlm(g, transpose(lambda_alpha_output[:,3]))
                close(g)
        end
end

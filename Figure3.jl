# simulate the stationary distribution 
# using Monte Carlo simulations
# cooperation frequency as a function of need threshold
################################################################################

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using IterativeSolvers
using SparseArrays
using Random
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

# generate er random network of size N and average node degree k
function er_random(N::Int64, k::Int64)
    mat_seq = spzeros(Int64, N, N)
    edge_num = Int64(N*k*0.5)

    edge_all = zeros(Int64, Int64(N*(N-1)/2))
    id = 1
    for i = 1:N, j = i+1:N
        edge_all[id] = (i-1)*N+j
        id = id + 1
    end
    edge_all = shuffle(edge_all)
    ii = [Int64(ceil(edge_all[i]/N)) for i = 1:edge_num]
    jj = [edge_all[i]-(ii[i]-1)*N for i = 1:edge_num]
    mat_seq = sparse(ii, jj, ones(Int64,edge_num), N, N)+sparse(jj, ii, ones(Int64,edge_num), N, N)

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

# generate sw small-world network with rewiring probability 0.1
function sw_small_world(N::Int64, k::Int64, p::Float64)
    mat_seq = spzeros(Int64, N, N)

    for i = 1:N, j = 1:Int64(0.5*k)
        h = ((i+j)%N > 0) ? ((i+j)%N) : N
        mat_seq[i,h] = 1
        mat_seq[h,i] = 1
    end
    for i = 1:N, j = 1:Int64(0.5*k)
        h = ((i+j)%N > 0) ? ((i+j)%N) : N
        toss = rand(1)
        if p < toss[1] continue end
        choice = findall(iszero, mat_seq[i,:])
        choice = setdiff(choice, i)
        if length(choice) == 0 continue end
        hh = rand(choice)
        mat_seq[i,h] = 0
        mat_seq[h,i] = 0
        mat_seq[i,hh] = 1
        mat_seq[hh,i] = 1
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

# generate ba scale-free network of size N and average node degree k
function ba_scale_free(N::Int64, k::Int64)
    mat_seq = spzeros(Int64, N, N)
    mat = ones(Int64, k+1, k+1)
    mat[diagind(mat)] = zeros(Int64,k+1)
    mat_seq[1:k+1,1:k+1] = mat

    for i = k+2:N
        for j = 1:Int64(k/2)
            choice = findall(iszero, mat_seq[i,:])
            choice = setdiff(choice, i)
            deg_seq = sum(mat_seq, dims = 2)
            deg_seq = deg_seq[:,1]
            toss = rand(1)
            threshold = toss[1]*sum(deg_seq[choice])
            deg_acc = 0.0
            for ll in choice
                deg_acc = deg_acc + deg_seq[ll]
                if deg_acc > threshold
                    mat_seq[i,ll] = 1
                    mat_seq[ll,i] = 1
                    break
                end
            end
        end
    end

    inds = findall(!iszero, mat_seq)
    a = getindex.(inds, 1)
    b = getindex.(inds, 2)
    edge_seq = zeros(Int64, length(b), 3)
    edge_seq[:,1] = a
    edge_seq[:,2] = b
    edge_seq[:,3] = ones(Int64, length(a))

    return edge_seq
end

# generate gkk scale-free network with exponent 2.5
# universal behavior of load distribution in scale-Free networks
# the degree distribution is p(k) = k^{-gamma}
function gkk_scale_free(N::Int64, k::Int64, gamma::Float64)
    mat_seq = spzeros(Int64, N, N)
    alpha = 1/(gamma-1)
    weight = [i^(-alpha) for i = 1:N]

    edge_num = Int64(N*k/2)
    while edge_num > 0
        deg_seq = sum(mat_seq, dims = 2)
        deg_seq = deg_seq[:,1]
        start_seq = findall(!isequal(N-1), deg_seq)
        start_toss = rand(1)
        start_node = bisection(weight[start_seq], start_toss[1]*sum(weight[start_seq]))
        start_node = start_seq[start_node]

        end_seq = findall(iszero, mat_seq[start_node,:])
        end_seq = setdiff(end_seq, start_node)
        end_toss = rand(1)
        end_node = bisection(weight[end_seq], end_toss[1]*sum(weight[end_seq]))
        end_node = end_seq[end_node]

        mat_seq[start_node, end_node] = 1
        mat_seq[end_node, start_node] = 1
        edge_num = edge_num - 1
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

# generate hk scale-free network with triad formation probability p
# growing scale-free networks with tunable clustering
function hk_scale_free(N::Int64, k::Int64, p::Float64)
    mat_seq = spzeros(Int64, N, N)
    mat = ones(Int64, k+1, k+1)
    mat[diagind(mat)] = zeros(Int64,k+1)
    mat_seq[1:k+1,1:k+1] = mat

    node_chosen = 1
    for i = k+2:N
        for j = 1:Int64(k/2)
            toss = rand(1)
            neg_i = findall(!iszero, mat_seq[i,:])
            neg_chosen = findall(!iszero, mat_seq[node_chosen,:])
            neg_diff = setdiff(neg_chosen, neg_i)
            neg_diff = setdiff(neg_diff, i)

            if j == 1 || toss[1] > p || length(neg_diff) == 0
                choice = findall(iszero, mat_seq[i,:])
                choice = setdiff(choice, i)
                deg_seq = sum(mat_seq, dims = 2)
                deg_seq = deg_seq[:,1]
                toss = rand(1)
                threshold = toss[1]*sum(deg_seq[choice])
                deg_acc = 0.0
                for ll in choice
                    deg_acc = deg_acc + deg_seq[ll]
                    if deg_acc > threshold
                        mat_seq[i,ll] = 1
                        mat_seq[ll,i] = 1
                        node_chosen = ll
                        break
                    end
                end
                continue
            end

            trid_node = rand(neg_diff)
            mat_seq[i,trid_node] = 1
            mat_seq[trid_node,i] = 1
            node_chosen = trid_node
        end
    end

    inds = findall(!iszero, mat_seq)
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

# bisection
function bisection(weight::Array{Float64}, toss::Float64)
    N = length(weight)
    boundary = 0
    a = 0
    b = N
    while true
        c = Int64(ceil((a+b)/2))
        if (sum(weight[1:c]) >= toss) b = c end
        if (sum(weight[1:c]) < toss) a = c end
        if b-a == 1 || b == a
            boundary = b
            break
        end
    end
    return boundary
end

# get the frequency of a specified type of behavior by Monte Carlo Simulations
function frequency_simulation(edge_seq::Array{Int64,2}, M::Array{Float64,2}, lambda_alpha::Array{Float64,2}, generation::Int64, sample_times::Int64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    L = size(M,1)

    frequency_list = zeros(Float64, sample_times)
    sample_generation = 0.2*generation

    # simulation
    for t = 1:sample_times  
        # initial action configuration
        str = rand(1:L,N)

        for g = 1:generation
            # selecting the node to be updated
            node_updated = rand(1:N)

            # calculate the payoff
            payoff_list = [M[str[node_updated],str[j]] for j = 1:N]
            payoff = transpose(P[node_updated,:])*payoff_list

            # action updating
            probability_list = zeros(Float64,L)
            for ell = 1:L 
                lambda = lambda_alpha[node_updated,1]
                if ell < str[node_updated]
                    alpha = ell*lambda_alpha[node_updated,2]/(L-1)
                    probability_list[ell] = 2/(L+L*exp(lambda*(payoff-alpha)))
                elseif ell > str[node_updated]
                    alpha = (ell-1)*lambda_alpha[node_updated,2]/(L-1)
                    probability_list[ell] = 2/(L+L*exp(-lambda*(payoff-alpha)))
                end
            end
            probability_list[str[node_updated]] = 1-sum(probability_list)

            # action to be used
            threshold = rand(1)
            sum1 = 0
            for i = 1:L
                sum1 = sum1 + probability_list[i]
                if sum1 >= threshold[1]
                    str[node_updated] = i
                    break
                end
            end

            if g > generation-sample_generation
                frequency_list[t] = frequency_list[t] + (sum(str)-N)/(N*(L-1))
            end
        end
        frequency_list[t] = frequency_list[t]/sample_generation
    end

    return sum(frequency_list)/sample_times
end

N = 100  # network size
k = 6    # average node degree
sample_times = 100   # number of simulations
generation = 100000   # generations of each simulation

net_list = ["RR","SW","ER","baSF","gkkSF","hkSF"]

for game = 1:4
        for net = 1:6
                edge_seq = zeros(Int64,k*N,3)
                # generate random regular networks
                if net == 1
                        edge_seq = random_regular(N, k)
                end
                # generate sw small-world networks
                if net == 2
                        edge_seq = sw_small_world(N, k, 0.1)
                end
                # generate er random networks
                if net == 3
                        edge_seq = er_random(N, k)
                end
                # generate ba scale-free networks
                if net == 4
                        edge_seq = ba_scale_free(N, k)
                end
                # generate gkk scale-free networks
                if net == 5
                        edge_seq = gkk_scale_free(N, k, 2.5)
                end
                # generate hk scale-free networks
                if net == 6
                        edge_seq = hk_scale_free(N, k, 0.1)
                end

                frequency_list = zeros(Float64, 8)
                for i = 1:8
                    L = 3
                    if game == 1 
                            L = 2
                    end
                    M = spzeros(Float64, L, L)
                    c = 1.0
                    b = 6.0

                    # two-action donation game
                    if game == 1
                            M = [0 b; -c b-c]
                    end 
                    # three-action linear donation game
                    if game == 2 
                            w = 1
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end
                    # three-action donation game with w=1.5
                    if game == 3 
                            w = 0.5
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end
                    # three-action donation game with w=0.5
                    if game == 4 
                            w = 1.5
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end

                    for j = 1:sample_times
                            # set independent and random behavioral motivation for all individuals, [motivation intensity, need threshold]
                            lambda_alpha = zeros(Float64, N, 2)
                            lambda_alpha[:,1] = 0.01*(rand(N).+0.5) # philanthropic motivation
                            lambda_alpha[:,2] = rand(N).-2.5.+i
                            frequency_list[i] = frequency_list[i] + frequency_simulation(edge_seq, M, lambda_alpha, generation, 1)
                    end
                end
                
                # output the behavior frequency 
                str_output = string("xC_philanthropic","_",net_list[net],"_game",string(game),".txt")
                g = open(str_output,"a")
                writedlm(g, transpose(frequency_list/sample_times))
                close(g)

                frequency_list = zeros(Float64, 8)
                for i = 1:8
                    L = 3
                    if game == 1 
                            L = 2
                    end
                    M = spzeros(Float64, L, L)
                    c = 1.0
                    b = 6.0

                    if game == 1
                            M = [0 b; -c b-c]
                    end 
                    if game == 2 
                            w = 1
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end
                    if game == 3 
                            w = 0.5
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end
                    if game == 4 
                            w = 1.5
                            M = [0 b/(1+w) b;-c/2 -c/2+b/(1+w) -c/2+b;-c -c+b/(1+w) -c+b]
                    end

                    for j = 1:sample_times
                            lambda_alpha = zeros(Float64, N, 2)
                            lambda_alpha[:,1] = 0.01*(rand(N).-1.5)  # aspirational motivation
                            lambda_alpha[:,2] = rand(N).-2.5.+i
                            frequency_list[i] = frequency_list[i] + frequency_simulation(edge_seq, M, lambda_alpha, generation, 1)
                    end
                end
                
                str_output = string("xC_aspirational","_",net_list[net],"_game",string(game),".txt")
                g = open(str_output,"a")
                writedlm(g, transpose(frequency_list/sample_times))
                close(g)
        end
end


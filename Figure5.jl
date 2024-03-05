# simulate the stationary distribution 
# using Monte Carlo simulations
################################################################################

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using IterativeSolvers
using SparseArrays
using Random
rng = MersenneTwister(1234)

# random regular network
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

# ba scale-free network
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

# function to get reproductive value of each node
# input
# edge list, edge_seq 
function pi_DB(edge_seq::Array{Int64,2})
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    M = M./sum(M,dims=1)

    MatA = spzeros(Float64, N, N)
    MatA = M - sparse(I,N,N)
    MatB_reduced = -MatA[1:N-1,N]
    MatA = MatA - MatA[:,N]*ones(Float64,1,N)
    MatA_reduced = MatA[1:N-1,1:N-1]

    pi_solution_reduced = idrs(MatA_reduced,MatB_reduced)

    pi_solution = zeros(Float64,N)
    pi_solution[1:N-1] = pi_solution_reduced
    pi_solution[N] = 1.0-sum(pi_solution)

    return pi_solution
end

# solve a set of equations to obtain eta_ij
# for rhoA > rhoB under a uniform distribution of initial mutant
# input
# edge list, edge_seq
function eta_DB(edge_seq::Array{Int64,2})
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    M = M./sum(M,dims=1)

    MT = transpose(M)
    inds = findall(!iszero, MT)
    a = getindex.(inds, 1)
    b = getindex.(inds, 2)
    vec = [i for i = 1:N]
    vec = reshape(vec,1,N)
    X1 = N*(a.-1)*ones(Int64,1,N)+ones(Int64,length(a))*vec
    Y1 = N*(b.-1)*ones(Int64,1,N)+ones(Int64,length(b))*vec
    W1 = MT[inds]*ones(Int64,1,N)

    vec = [(i-1)*N for i = 1:N]
    vec = reshape(vec,1,N)
    X2 = a*ones(Int64,1,N)+ones(Int64,size(a,1))*vec
    Y2 = b*ones(Int64,1,N)+ones(Int64,size(b,1))*vec
    W2 = MT[inds]*ones(Int64,1,N)

    X11 = reshape(X1,:,1)
    Y11 = reshape(Y1,:,1)
    W11 = reshape(W1,:,1)
    X22 = reshape(X2,:,1)
    Y22 = reshape(Y2,:,1)
    W22 = reshape(W2,:,1)
    matA = sparse(X11[:,1],Y11[:,1],W11[:,1],N^2,N^2)+sparse(X22[:,1],Y22[:,1],W22[:,1],N^2,N^2)

    matA = matA/2 - sparse(I, N^2, N^2)
    vec = setdiff(1:N^2,[(i-1)*N+i for i = 1:N])
    matA_reduced = matA[vec, vec]
    matB_reduced = ones(Float64,N^2-N)

    eta_solution_reduced = idrs(matA_reduced,matB_reduced)
    eta_solution = spzeros(Float64, N^2)
    eta_solution[vec] = eta_solution_reduced
    eta_solution = reshape(eta_solution, N, N)

    return transpose(eta_solution)
end

# get (b/c)^* under death-birth updating
# input
# edge list, edge_seq
function bc_DB(edge_seq::Array{Int64,2})
    pi = pi_DB(edge_seq)
    eta = eta_DB(edge_seq)
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    Mp = M./sum(M,dims=2)

    n1 = sum(pi.*eta.*Mp)
    n3 = sum(pi.*eta.*Mp^3)
    n2 = sum(pi.*eta.*Mp^2)

    return n2/(n3-n1)
end

# get fixation probability of A-individuals by Monte Carlo Simulations
function rhoA_simulation_expression_DB(edge_seq::Array{Int64,2}, M::Array{Float64,2}, lambda_alpha::Array{Float64,2}, generation::Int64, sample_times::Int64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    L = size(M,1)

    fixation_times = 0
    # simulation
    for t = 1:sample_times 
        # initial strategy configuration
        str = zeros(Int64,N)
        mutant = rand(1:N)
        str[mutant] = 1

        for g = 1:generation
            # selecting the node to be updated
            node_updated = rand(1:N)
            node_competitor = findall(!iszero, P[node_updated,:])

            # check if the population enters into an aborbing state
            if g % 1000 == 0
                if sum(str) == N
                    fixation_times = fixation_times + 1
                    break
                end
                if sum(str) == 0
                    break
                end
            end
            if sum(str[node_competitor]) == length(node_competitor)*str[node_updated]
                continue
            end

            # calculate the payoff
            b = M[1,2]
            c = -M[2,1]
            A = lambda_alpha[1,1]*(b-c-2*lambda_alpha[1,2])
            B = lambda_alpha[2,1]*(b-c-2*lambda_alpha[2,2])
            payoff_list = P*str*b*(A-B) - str*c*(A-B) .+ b*B .- c*B
            payoff_list = payoff_list/8 .+ (b-c)/2

            payoff_competitor = P[node_updated,:].*payoff_list
            threshold = rand(1)*sum(payoff_competitor)
            sum1 = 0
            for i = 1:length(node_competitor)
                sum1 = sum1 + payoff_competitor[node_competitor[i]]
                if sum1 >= threshold[1]
                    str[node_updated] = str[node_competitor[i]]
                    break
                end
            end
        end
    end

    return fixation_times/sample_times
end

N = 100
k = 6
generation = 20000000
sample_times = 5000000

# random regular network
edge_seq = zeros(Int64,k*N,3)
edge_seq = random_regular(N, k)

lambda_alpha = [0.01 1; 0.01 3]
c = 1.0
b_critical = bc_DB(edge_seq)
for i = 1:17
        b = b_critical+(i-9)*0.5
        M = [0 b; -c b-c]

        rhoA = rhoA_simulation_expression_DB(edge_seq, M, lambda_alpha, generation, sample_times)

        g = open("RR_rhoA_case1.txt","a")
        writedlm(g, [b rhoA])
        close(g)
end

lambda_alpha = [0.03 2; 0.01 2]
c = 1.0
b_critical = bc_DB(edge_seq)
for i = 1:17
        b = b_critical+(i-9)*0.5
        M = [0 b; -c b-c]

        rhoA = rhoA_simulation_expression_DB(edge_seq, M, lambda_alpha, generation, sample_times)

        g = open("RR_rhoA_case2.txt","a")
        writedlm(g, [b rhoA])
        close(g)
end

# ba scale-free networks
edge_seq = zeros(Int64,k*N,3)
edge_seq = ba_scale_free(N, k)

lambda_alpha = [0.01 1; 0.01 3]
c = 1.0
b_critical = bc_DB(edge_seq)
for i = 1:17
        b = b_critical+(i-9)*0.5
        M = [0 b; -c b-c]

        rhoA = rhoA_simulation_expression_DB(edge_seq, M, lambda_alpha, generation, sample_times)

        g = open("SF_rhoA_case1.txt","a")
        writedlm(g, [b rhoA])
        close(g)
end

lambda_alpha = [0.03 2; 0.01 2]
c = 1.0
b_critical = bc_DB(edge_seq)
for i = 1:17
        b = b_critical+(i-9)*0.5
        M = [0 b; -c b-c]

        rhoA = rhoA_simulation_expression_DB(edge_seq, M, lambda_alpha, generation, sample_times)

        g = open("SF_rhoA_case2.txt","a")
        writedlm(g, [b rhoA])
        close(g)
end

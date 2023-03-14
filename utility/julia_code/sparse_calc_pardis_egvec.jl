using DelimitedFiles, LinearAlgebra, JSON
using HDF5
using ArgParse
using SparseArrays
using Pardiso, Arpack, LinearMaps
using JLD


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir", "-i"
            help = "path of rlat.dat, orbital_types.dat, site_positions.dat, hamiltonians_pred.h5, and overlaps.h5"
            arg_type = String
            default = "./"
        "--output_dir", "-o"
            help = "path of output openmx.Band"
            arg_type = String
            default = "./"
        "--config"
            help = "config file in the format of JSON"
            arg_type = String
    end
    return parse_args(s)
end
parsed_args = parse_commandline()


function _create_dict_h5(filename::String)
    fid = h5open(filename, "r")
    T = eltype(fid[keys(fid)[1]])
    d_out = Dict{Array{Int64,1}, Array{T, 2}}()
    for key in keys(fid)
        data = read(fid[key])
        nk = map(x -> parse(Int64, convert(String, x)), split(key[2 : length(key) - 1], ','))
        d_out[nk] = permutedims(data)
    end
    close(fid)
    return d_out
end


# The function construct_linear_map below is come from https://discourse.julialang.org/t/smallest-magnitude-eigenvalues-of-the-generalized-eigenvalue-equation-for-a-large-sparse-matrix/75485/11
function construct_linear_map(H, S)
    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.COMPLEX_HERM_INDEF)
    pardisoinit(ps)
    fix_iparm!(ps, :N)
    H_pardiso = get_matrix(ps, H, :N)
    b = rand(ComplexF64, size(H, 1))
    set_phase!(ps, Pardiso.ANALYSIS)
    pardiso(ps, H_pardiso, b)
    set_phase!(ps, Pardiso.NUM_FACT)
    pardiso(ps, H_pardiso, b)
    return (
        LinearMap{ComplexF64}(
            (y, x) -> begin
                set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
                pardiso(ps, y, H_pardiso, S * x)
            end,
            size(H, 1);
            ismutating=true
        ),
        ps
    )
end

const ev2Hartree = 0.036749324533634074
const Bohr2Ang = 0.529177249

function genlist(x)
    return collect(range(x[1], stop = x[2], length = Int64(x[3])))
end

function k_data2num_ks(kdata::AbstractString)
    return parse(Int64,split(kdata)[1])
end

function k_data2kpath(kdata::AbstractString)
    return map(x->parse(Float64,x), split(kdata)[2:7])
end

function std_out_array(a::AbstractArray)
    return string(map(x->string(x," "),a)...)
end

default_dtype = Complex{Float64}

println(parsed_args["config"])
config = JSON.parsefile(parsed_args["config"])
calc_job = config["calc_job"]

if isfile(joinpath(parsed_args["input_dir"],"info.json"))
    spinful = JSON.parsefile(joinpath(parsed_args["input_dir"],"info.json"))["isspinful"]
else
    spinful = false
end

site_positions = readdlm(joinpath(parsed_args["input_dir"], "site_positions.dat"))
nsites = size(site_positions, 2)

orbital_types_f = open(joinpath(parsed_args["input_dir"], "orbital_types.dat"), "r")
site_norbits = zeros(nsites)
orbital_types = Vector{Vector{Int64}}()
for index_site = 1:nsites
    orbital_type = parse.(Int64, split(readline(orbital_types_f)))
    push!(orbital_types, orbital_type)
end
site_norbits = (x->sum(x .* 2 .+ 1)).(orbital_types) * (1 + spinful)
norbits = sum(site_norbits)
site_norbits_cumsum = cumsum(site_norbits)

rlat = readdlm(joinpath(parsed_args["input_dir"], "rlat.dat"))


if isfile(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"))
    @info string("read sparse matrix from ", parsed_args["input_dir"], "/sparse_matrix.jld")
    H_R = load(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "H_R")
    S_R = load(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "S_R")
else
    if true
        @info "read h5"
        begin_time = time()
        hamiltonians_pred = _create_dict_h5(joinpath(parsed_args["input_dir"], "hamiltonians_pred.h5"))
        overlaps = _create_dict_h5(joinpath(parsed_args["input_dir"], "overlaps.h5"))
        println("Time for reading h5: ", time() - begin_time, "s")

        I_R = Dict{Vector{Int64}, Vector{Int64}}()
        J_R = Dict{Vector{Int64}, Vector{Int64}}()
        H_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()
        S_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()

        @info "construct sparse matrix in the format of COO"
        begin_time = time()
        for key in collect(keys(hamiltonians_pred))
            hamiltonian_pred = hamiltonians_pred[key]
            if (key ∈ keys(overlaps))
                overlap = overlaps[key]
            else
                # continue
                overlap = zero(hamiltonian_pred)
            end
            if spinful
                overlap = vcat(hcat(overlap,zeros(size(overlap))),hcat(zeros(size(overlap)),overlap)) # the readout overlap matrix only contains the upper-left block # TODO maybe drop the zeros?
            end
            R = key[1:3]; atom_i=key[4]; atom_j=key[5]

            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(hamiltonian_pred)
            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(overlap)
            if !(R ∈ keys(I_R))
                I_R[R] = Vector{Int64}()
                J_R[R] = Vector{Int64}()
                H_V_R[R] = Vector{default_dtype}()
                S_V_R[R] = Vector{default_dtype}()
            end
            for block_matrix_i in 1:site_norbits[atom_i]
                for block_matrix_j in 1:site_norbits[atom_j]
                    coo_i = site_norbits_cumsum[atom_i] - site_norbits[atom_i] + block_matrix_i
                    coo_j = site_norbits_cumsum[atom_j] - site_norbits[atom_j] + block_matrix_j
                    push!(I_R[R], coo_i)
                    push!(J_R[R], coo_j)
                    push!(H_V_R[R], hamiltonian_pred[block_matrix_i, block_matrix_j])
                    push!(S_V_R[R], overlap[block_matrix_i, block_matrix_j])
                end
            end
        end
        println("Time for constructing sparse matrix in the format of COO: ", time() - begin_time, "s")
    else
        @info "construct sparse matrix in the format of COO from h5"
        begin_time = time()
        I_R = Dict{Vector{Int64}, Vector{Int64}}()
        J_R = Dict{Vector{Int64}, Vector{Int64}}()
        H_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()
        S_V_R = Dict{Vector{Int64}, Vector{default_dtype}}()
        fid_H = h5open(joinpath(parsed_args["input_dir"], "hamiltonians_pred.h5"), "r")
        fid_S = h5open(joinpath(parsed_args["input_dir"], "overlaps.h5"), "r")

        for key in keys(fid_H)
            hamiltonian_pred = permutedims(read(fid_H[key]))
            overlap = permutedims(read(fid_S[key]))
            Rij = map(x -> parse(Int64, convert(String, x)), split(key[2 : length(key) - 1], ','))
            R = Rij[1:3]; atom_i=Rij[4]; atom_j=Rij[5]

            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(hamiltonian_pred)
            @assert (site_norbits[atom_i], site_norbits[atom_j]) == size(overlap)
            if !(R ∈ keys(I_R))
                I_R[R] = Vector{Int64}()
                J_R[R] = Vector{Int64}()
                H_V_R[R] = Vector{default_dtype}()
                S_V_R[R] = Vector{default_dtype}()
            end
            for block_matrix_i in 1:site_norbits[atom_i]
                for block_matrix_j in 1:site_norbits[atom_j]
                    coo_i = site_norbits_cumsum[atom_i] - site_norbits[atom_i] + block_matrix_i
                    coo_j = site_norbits_cumsum[atom_j] - site_norbits[atom_j] + block_matrix_j
                    push!(I_R[R], coo_i)
                    push!(J_R[R], coo_j)
                    push!(H_V_R[R], hamiltonian_pred[block_matrix_i, block_matrix_j])
                    push!(S_V_R[R], overlap[block_matrix_i, block_matrix_j])
                end
            end
        end
        println("Time for reading h5 in the format of COO: ", time() - begin_time, "s")
    end

    @info "convert sparse matrix to the format of CSC"
    begin_time = time()
    H_R = Dict{Vector{Int64}, SparseMatrixCSC{default_dtype, Int64}}()
    S_R = Dict{Vector{Int64}, SparseMatrixCSC{default_dtype, Int64}}()

    for R in keys(I_R)
        H_R[R] = sparse(I_R[R], J_R[R], H_V_R[R], norbits, norbits)
        S_R[R] = sparse(I_R[R], J_R[R], S_V_R[R], norbits, norbits)
    end
    println("Time for converting to the format of CSC: ", time() - begin_time, "s")

    save(joinpath(parsed_args["input_dir"], "sparse_matrix.jld"), "H_R", H_R, "S_R", S_R)
end

if calc_job == "band"
    which_k = config["which_k"] # which k point to calculate, start counting from 1, 0 for all k points
    fermi_level = config["fermi_level"]
    max_iter = config["max_iter"]
    num_band = config["num_band"]
    k_data = config["k_data"]

    @info "calculate bands"
    num_ks = k_data2num_ks.(k_data)
    kpaths = k_data2kpath.(k_data)
    println(num_ks)
    println(size(kpaths,1))
    egvals = zeros(Float64, num_band, sum(num_ks)[1])

    begin_time = time()
    idx_k = 1
    for i = 1:size(kpaths, 1)
        begin_time = time()
        kpath = kpaths[i]
        pnkpts = num_ks[i]
        kxs = LinRange(kpath[1], kpath[4], pnkpts)
        kys = LinRange(kpath[2], kpath[5], pnkpts)
        kzs = LinRange(kpath[3], kpath[6], pnkpts)
        for (kx, ky, kz) in zip(kxs, kys, kzs)
            global idx_k
            if which_k == 0 || which_k == idx_k
                H_k = spzeros(default_dtype, norbits, norbits)
                S_k = spzeros(default_dtype, norbits, norbits)
                for R in keys(H_R)
                    H_k += H_R[R] * exp(im*2π*([kx, ky, kz]⋅R))
                    S_k += S_R[R] * exp(im*2π*([kx, ky, kz]⋅R))
                end
                # H_k = (H_k + H_k') / 2
                flush(stdout)
                lm, ps = construct_linear_map(H_k - (fermi_level) * S_k, S_k)
                println("Time for No.$idx_k matrix factorization: ", time() - begin_time, "s")
                egval_inv, X = eigs(lm, nev=num_band, which=:LM, ritzvec=true, maxiter=max_iter)
                set_phase!(ps, Pardiso.RELEASE_ALL)
                pardiso(ps)
                egval = real(1 ./ egval_inv) .+ (fermi_level)
                # runzhang
                egvec = X
                # Tangzechen
                # Pband based on i equals Re(\sum_j conj(c_j)*c_i*S_{ij})
                # Ritz vectors aren't normalized, so we must normalize for each vector
                println(num_band)
                println(typeof(S_k))
                println(typeof(H_k))
                println(typeof(egval))
                println(typeof(egvec))
                for ivec = 1:num_band
                    normalization_factor = egvec[:, ivec]' * S_k * egvec[:, ivec]
                    # println(ivec)
                    # println(egval[ivec])
                    # println(egval_inv[ivec])
                    # println(normalization_factor)
                    # println(norm(H_k * egvec[:, ivec] - egval[ivec] * S_k * egvec[:, ivec]))
                    # for jvec = 1:num_band
                    #     println(ivec)
                    #     println(jvec)
                    #     println(egvec[:, jvec]' * S_k * egvec[:, ivec])
                    # end
                    egvec[:, ivec] /= sqrt(real(normalization_factor))
                    # println( egvec[:, ivec]' * S_k * egvec[:, ivec])
                end

                # begin_time = time()
                # pband = zero(egvec)
                # for iorb = 1:norbits
                #     println(iorb)
                #     flush(stdout)
                #     for jvec = 1:num_band
                #         pband[iorb, jvec] += conj(egvec[iorb, jvec]) * (S_k * egvec[:, jvec])[iorb]
                #     end
                # end
                # println("Old method: ", time() - begin_time, "s")
                
                begin_time = time()
                pband = zero(egvec)
                Segvec = S_k * egvec
                pband = conj(egvec) .* Segvec
                println("New method: ", time() - begin_time, "s")

                
                # egval = real(eigs(H_k, S_k, nev=num_band, sigma=(fermi_level + lowest_band), which=:LR, ritzvec=false, maxiter=max_iter)[1])
                if which_k == 0
                    println(egval .- fermi_level)
                else
                    open(joinpath(parsed_args["output_dir"], "kpoint.dat"), "w") do f
                        writedlm(f, [kx, ky, kz])
                    end
                    open(joinpath(parsed_args["output_dir"], "egval.dat"), "w") do f
                        writedlm(f, egval)
                    end
                    # runzhang
                    # open(joinpath(parsed_args["output_dir"], "egvec.dat"), "w") do f
                    #     writedlm(f, egvec)
                    # end
                    # Tang Zechen
                    open(joinpath(parsed_args["output_dir"], "pband.dat"), "w") do f
                        writedlm(f, real(pband))
                    end
                end
                egvals[:, idx_k] = egval
                println("Time for solving No.$idx_k eigenvalues at k = ", [kx, ky, kz], ": ", time() - begin_time, "s")
            end
            idx_k += 1
        end
    end

    # output in openmx band format
    f = open(joinpath(parsed_args["output_dir"], "openmx.Band"),"w")
    println(f, num_band, " ", 0, " ", ev2Hartree * fermi_level)
    openmx_rlat = reshape((rlat .* Bohr2Ang), 1, :)
    println(f, std_out_array(openmx_rlat))
    println(f, length(k_data))
    for line in k_data
        println(f,line)
    end
    idx_k = 1
    for i = 1:size(kpaths, 1)
        pnkpts = num_ks[i]
        kstart = kpaths[i][1:3]
        kend = kpaths[i][4:6]
        k_list = zeros(Float64,pnkpts,3)
        for alpha = 1:3
            k_list[:,alpha] = genlist([kstart[alpha],kend[alpha],pnkpts])
        end
        for j = 1:pnkpts
            global idx_k
            kvec = k_list[j,:]
            println(f, num_band, " ", std_out_array(kvec))
            println(f, std_out_array(ev2Hartree * egvals[:, idx_k]))
            idx_k += 1
        end
    end
    close(f)
end

### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 04c7fbc4-4d94-11ed-1276-f58f8d304231
using Random

# ╔═╡ 68dcf61a-3162-4f11-b9d1-335ada739e43
using Combinatorics

# ╔═╡ d74ce4a0-6f9b-4e63-9728-0208632acf7d
#This file creates simulation results in Table 1, Column (Algorithm-4.3) in the paper.

# ╔═╡ 7f963ab5-958a-4bb1-abf5-82ccfe403b23
function neyman_estimator(data)
	m = floor(Int64,length(data)/2)
	T = sum(data[1:m])/m - sum(data[m+1:2*m])/m
end

# ╔═╡ 00b56923-ce85-4a1e-8ed7-68d238d2e07f
function tau_true(table)
	return (sum(table[1,:]) - sum(table[2,:]))/length(table[1,:])
end

# ╔═╡ ea6bc7de-5112-4b5a-9175-01cbcd81ea78
function p_value(T_obs, w,  num_perm)
    num_perm = num_perm + 1
    p = 0
    m = floor(Int64, length(w[1, :]) / 2)
    tau_0 = tau_true(w)

    # Number of subjects
    n = length(w[1, :])

    # Loop over all (n choose n/2) combinations of subjects
    for comb in combinations(1:n, m)
        # Place the selected subjects in the first m positions
        w_comb = w[:, comb]  # Select the chosen subjects
        w_rest = setdiff(1:n, comb)  # Remaining subjects
        
        # Construct the new data array
        w_data = vcat(w[1, comb], w[2, w_rest])

        T_w = neyman_estimator(w_data)
        if abs(T_w - tau_0) >= abs(T_obs - tau_0)
            p += 1
        end
    end

    return p / binomial(n, m), num_perm
end

# ╔═╡ fa8efaef-18cf-4f80-adcf-d31da0bd5b01
function create_table(i,j,k,l)
	y  = zeros(Int64, 2, 0)
	y = hcat(y, repeat([1,1],1,i), repeat([1,0],1,j), repeat([0,1],1,k), repeat([0,0],1,l))
	return y
end

# ╔═╡ 29577845-c090-4edc-a953-af08d303a592
function check_possible(N,V)
	return ( max(0, N[1] - V[2], V[1] - N[3], V[1] + V[3] - N[2] - N[3]) 
		<= min(V[1], N[1], V[1] + V[3] - N[3], sum(N) - V[2] - N[3] - N[2]) )
end

# ╔═╡ c7432f87-83bf-404d-9618-db05592eb071
function median(a::Int64, b::Int64)
	return floor(Int64,(a+b)/2)
end

# ╔═╡ b4e52381-2725-4cd0-b44c-e17830c5dd52
function is_possible(N,j,t0)

	#This function checks whether the potential outcome table, indexed by j, is compatible with actual data, N, and the hypothesis to be tested, t0

	flag1 = (j >= t0 + N[3]) && (j >= N[1]) && (sum(N)>=j+N[2]) && (N[1] + t0 + N[2] + N[3] >= j)
	flag2 = max(t0, j - N[1] - N[3], N[1] + N[3] + t0 -j) <= min(j, N[1] + N[4], N[2] + N[3] + t0, sum(N) + t0 -j)

	return flag1 && flag2
end

# ╔═╡ f2cb7f02-6623-4421-b8c0-7e0830d10257
function find_v10(N,j,t0)
	return max(t0, j - N[1] - N[3], N[1] + N[3] + t0 -j,0)
end 

# ╔═╡ 8fb19f7e-5d04-4f96-9270-9ffed94147d1
function find_v10_upper(N,j,t0)
	return min(j, N[1] + N[4], N[2] + N[3] + t0, sum(N) + t0 - j)
end

# ╔═╡ 966da19d-6713-45de-aa74-4bb743c0b721
function tau_compatible(N,t0,alpha,num_perm) #here t0 is an integer
	
	###inputs
	# T_obs: estimator using the actual data
	# samples: number of Monte Carlo draws
	# t0: hypothesized value to be tested
	# alpha: coverage of the interval
	# num_perm: (global variable) to keep track of the number of permutation tests
	# multi_thread: whether the MC procedures are paralleled or not
	###
	
	compatible = false 
	T_obs_rescaled = 2*(N[1] - N[3])
	n = sum(N)
	for j in 0:n
		#println(j)
		#println(is_possible(N,j,t0))
		if is_possible(N,j,t0)
			v10 = find_v10(N,j,t0)
			v11 = j - v10
			v01 = v10 - t0
			v00 = n - j - v10 + t0
			w = create_table(v11,v10,v01,v00)
			p_val, num_perm = p_value(T_obs_rescaled/n, w,num_perm)
			
			if p_val >= alpha
				compatible = true
				break 
			end

			if v10 == 0 && v01 == 0
				v10 = 1
				v11 = j - v10
				v01 = v10 - t0
				v00 = n - j - v10 + t0
				if v11 >= 0 && v11 >= 0 && v00 >= 0 && v10 <= find_v10_upper(N,j,t0)

					#output p value and num of permutation sofar
					p_val, num_perm = p_value(T_obs_rescaled/n, w, samples,num_perm,multi_thread)

					if p_val >= alpha
					compatible = true
					break 
					end
				end
			end
		end
	end
	return compatible, num_perm
end

# ╔═╡ a2798fc9-24f1-4a76-bb0a-9ae31d09f491
function efficient_perm(N, alpha,num_perm)

	T_obs_rescaled = 2*(N[1] - N[3])
	#m = floor(Int64, sum(N)/2)
	#T_obs = T_obs_rescaled/m
	max_compatible_tau = N[1] + N[4] #rescaled by n
	min_compatible_tau = -N[2] - N[3]
	#upper = 1
	#lower = -1

	bin_lower = T_obs_rescaled
	bin_upper = max_compatible_tau
	
	counter = 0
	while counter < 50
		#println(counter)
		counter += 1 
		if bin_lower == bin_upper 
			break
		end

		#last step
		if bin_lower == bin_upper - 1
			if_compatible, num_perm =tau_compatible(N,bin_upper,alpha,num_perm)
			if if_compatible == true
				bin_lower = bin_upper
				break
			else 
				bin_upper = bin_lower
				break
			end
		end
			

		t0 = median(bin_lower, bin_upper)

		if_compatible, num_perm = tau_compatible(N,t0,alpha,num_perm)
		
		if if_compatible == true
			bin_lower = t0
		end

		if if_compatible == false
			bin_upper = t0
		end
		
	end

	n = sum(N)
	upper = bin_lower/n

	bin_lower = min_compatible_tau
	bin_upper = T_obs_rescaled

	counter = 0
	while counter < 50
		#println(counter)
		counter += 1 
		if bin_lower == bin_upper 
			break
		end

		if bin_lower == bin_upper - 1
			if_compatible, num_perm = tau_compatible(N,bin_lower,alpha,num_perm)
			if if_compatible == true
				bin_upper = bin_lower
				break
			else 
				bin_lower = bin_upper
				break
			end
		end
			

		t0 = median(bin_lower, bin_upper)

		if_compatible, num_perm = tau_compatible(N,t0,alpha,num_perm)
		
		if if_compatible == true
			bin_upper = t0
		end

		if if_compatible == false
			bin_lower = t0
		end
		
	end

	lower = bin_lower/n
	return n*[lower,upper], num_perm
	
end

# ╔═╡ 2fca2bd9-8cf3-4cb5-8379-364f81fef446
function perm_test(n11,n10,n01,n00,alpha)
	
	#printing basic information
	n=n11+n10+n01+n00
	println("Number of units in the experiment: ",n )
	println("Number of treated with outcome 1: ",n11)
	println("Number of treated with outcome 0: ",n10)
	println("Number of untreated with outcome 1: ",n01)
	println("Number of untreated with outcome 0: ",n00)

	println("Coverage probability of the two-sided interval: ", alpha)

	#nsample = 1/ (epsilon^2) * log(8 * n * log2(n)/epsilon)
	#println("Number of Monte Carlo draws for each permutation test: ", nsample)

	#this variable counts the number of permutation tests conducted
	num_perm ::Int64 = 0

	N_tuple= (n11,n10,n01,n00)

	#decrease alpha by epsilon as equation (35) of the paper
	#alpha = alpha - epsilon
	interval, num_perm = efficient_perm(N_tuple, alpha, num_perm)

	println("")
	println("95% Confidence Interval:", interval)
	println("Number of permutation tests: ", num_perm)
end

# ╔═╡ ceb20b1c-a5d3-426d-b546-e76d3c62972f
@time perm_test(2,6,8,0,0.05)

# ╔═╡ 08a698d2-c4b0-45f9-8a1a-1f8e4f8825bb
@time perm_test(6,4,4,6,0.05)

# ╔═╡ 409ef880-13fc-4704-9675-73abdb7f951e
@time perm_test(8,4,5,7,0.05)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Combinatorics = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "ccc200d4efd8889f0534966e4807c3d3e8e2fbee"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
"""

# ╔═╡ Cell order:
# ╠═d74ce4a0-6f9b-4e63-9728-0208632acf7d
# ╠═04c7fbc4-4d94-11ed-1276-f58f8d304231
# ╠═68dcf61a-3162-4f11-b9d1-335ada739e43
# ╠═7f963ab5-958a-4bb1-abf5-82ccfe403b23
# ╠═00b56923-ce85-4a1e-8ed7-68d238d2e07f
# ╠═ea6bc7de-5112-4b5a-9175-01cbcd81ea78
# ╠═fa8efaef-18cf-4f80-adcf-d31da0bd5b01
# ╠═29577845-c090-4edc-a953-af08d303a592
# ╠═c7432f87-83bf-404d-9618-db05592eb071
# ╠═b4e52381-2725-4cd0-b44c-e17830c5dd52
# ╠═f2cb7f02-6623-4421-b8c0-7e0830d10257
# ╠═8fb19f7e-5d04-4f96-9270-9ffed94147d1
# ╠═966da19d-6713-45de-aa74-4bb743c0b721
# ╠═a2798fc9-24f1-4a76-bb0a-9ae31d09f491
# ╠═2fca2bd9-8cf3-4cb5-8379-364f81fef446
# ╠═ceb20b1c-a5d3-426d-b546-e76d3c62972f
# ╠═08a698d2-c4b0-45f9-8a1a-1f8e4f8825bb
# ╠═409ef880-13fc-4704-9675-73abdb7f951e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

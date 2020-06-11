function ADMM(problem, num_iter)
	X = zeros(problem.n*(problem.N+1))
	U = zeros(problem.m*problem.N)
	Y = zeros(problem.m*problem.N)
	Λ̄ = ones(problem.m*problem.N)

	for i=1:num_iter
		(X, U) = LQR(ρ, Y - Λ̄, problem)
		# (X, U) = LQRCVX(ρ, Y,  Λ̄, problem)
		for k=1:N
			β = problem.α/problem.ρ
			v = U[SelectControl(k)] + Λ̄[SelectControl(k)]
			Y[SelectControl(k)] = L2Prox(β, v, problem)
			Λ̄[SelectControl(k)] += U[SelectControl(k)] - Y[SelectControl(k)]

			# @show norm(U[SelectControl(k)] - Y[SelectControl(k)])
		end
		println("Iter ", i, " -- Norm(U-Y): ", norm(U-Y))
	end

	return X, U, Y
end

function ADMMGroupConsensus(problem, num_iter)
	X = zeros(problem.n*(problem.N+1))
	U = zeros(problem.m*problem.N)
	Y = 1e-2*ones(problem.m*problem.N)
	Z = 1e-2*ones(problem.m*problem.N)

	W = 1/3*(U + Y + Z)

	λ = ones(problem.m*problem.N) # Lagrange multiplier for LQR
	ν = ones(problem.m*problem.N) # Lagrange multiplier for L2
	μ = ones(problem.m*problem.N) # Lagrange multiplier for control constraints

	for i=1:num_iter
		U_ref = W - (1/problem.ρ)*λ
		(X, U) = LQR(ρ, U_ref, problem)
		for k=1:N
			β = problem.α/problem.ρ
			v = W[SelectControl(k)] - 1/problem.ρ * ν[SelectControl(k)]
			Y[SelectControl(k)] = L2Prox(β, v, problem)

			Z[SelectControl(k)] = ControlProjection(W[SelectControl(k)] - 1/problem.ρ * μ[SelectControl(k)], problem)
			# @show norm(U[SelectControl(k)] - Y[SelectControl(k)])
		end
		W = 1/3*(U + Y + Z)
		λ = λ + problem.ρ*(U - W)
		ν = ν + problem.ρ*(Y - W)
		μ = μ + problem.ρ*(Z - W)

		println("Iter ", i, " -- Norm(U-Y): ", norm(U-Y))
	end

	return X, U, Y, Z
end

function objective(X, U, p)
	cost = 0
	for i=1:(p.N)
		cost += p.α*norm(U[SelectControl(i)])
		if i != 1
			cost += 0.5*X[SelectState(i)]'*p.Q_k*X[SelectState(i)]
		end
	end

	cost += 0.5*X[SelectState(p.N+1)]'*p.Q_f*X[SelectState(p.N+1)]

	return cost
end#

function LQR(ρ, U_ref, p)
	# Cost to go function over time
	V = zeros(p.n, p.n, p.N+1)
	V[:,:, p.N+1] = p.Q_f
	v = zeros(p.n, p.N+1)

	K = zeros(p.m, p.n, p.N)
	d = zeros(p.m, p.N)

	R = ρ/2*Matrix{Float64}(I, p.m, p.m)

	k = p.N
	while k > 0
		K[:,:,k] = -inv(R + p.B'*V[:,:,k+1]*p.B)*p.B'*V[:,:,k+1]*p.A
		d[:,k] = inv(R + p.B'*V[:,:,k+1]*p.B)*(R*U_ref[SelectControl(k)] - p.B'*v[:,k+1])

		V[:,:, k] = p.Q_k + K[:,:,k]'*R*K[:,:,k] + (p.A + p.B*K[:,:,k])'*V[:,:,k+1]*(p.A+p.B*K[:,:,k])
		v[:, k] = ((d[:,k] - U_ref[SelectControl(k)])'*R*K[:,:,k] +
					d[:,k]'*p.B'*V[:,:,k+1]*(p.A+p.B*K[:,:,k]) + v[:, k+1]'*(p.A+p.B*K[:,:,k]))'

		k -= 1
	end

	X = zeros(p.n*(p.N+1))
	U = zeros(p.m*p.N)
	X[SelectState(1)] = p.x0

	for i=2:(p.N+1)
		U[SelectControl(i-1)] = K[:,:, i-1]*X[SelectState(i-1)] + d[:, i-1]
		X[SelectState(i)] = p.A*X[SelectState(i-1)] + p.B*U[SelectControl(i-1)]
	end

	return (X, U)
end

function L2Prox(β, v, p)
	return max(0, 1-β/norm(v, 2))*v
end

function L1Prox(β, v, p)
	return max.(0, v .- β) - max.(0, -v .- β)
end

function ControlProjection(v, p)
	if norm(v) > p.u_max
		return v/norm(v)*p.u_max
	else
		return v
	end
end

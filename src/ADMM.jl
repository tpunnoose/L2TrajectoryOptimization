function ADMM(problem, num_iter)
	n = 6
	N = 10
	m = 3

	X = zeros(problem.n*(problem.N+1))
	U = zeros(problem.m*problem.N)
	Y = zeros(problem.m*problem.N)
	Λ̄ = ones(problem.m*problem.N)

	for i=1:num_iter
		# (X, U) = LQR(ρ, Y, Λ̄, problem)
		(X, U) = LQRCVX(ρ, Y,  Λ̄, problem)
		for k=1:N
			β = problem.α/problem.ρ
			v = U[SelectControl(k)] + Λ̄[SelectControl(k)]
			Y[SelectControl(k)] = L1Prox(β, v, p)
			Λ̄[SelectControl(k)] += U[SelectControl(k)] - Y[SelectControl(k)]

			# @show norm(U[SelectControl(k)] - Y[SelectControl(k)])
		end
		println("Iter ", i, " -- Norm(U-Y): ", norm(U-Y))
	end

	return X, U, Y
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
end

function LQR(ρ, Y, Λ̄, p)
	# Cost to go function over time
	V = zeros(p.n, p.n, p.N+1)
	V[:,:, p.N+1] = p.Q_f
	v = zeros(p.n, p.N+1)

	K = zeros(p.m, p.n, p.N)
	d = zeros(p.m, p.N)

	R = ρ/2*Matrix{Float64}(I, p.m, p.m)
	U_ref = Y - Λ̄

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
	l2_prox = max(0, 1-β/norm(v, 2))*v
	if norm(l2_prox) > p.u_max
		return l2_prox/norm(l2_prox)*p.u_max
	else
		return l2_prox
	end
end

function L1Prox(β, v, p)
	return max.(0, v .- β) - max.(0, -v .- β)
end

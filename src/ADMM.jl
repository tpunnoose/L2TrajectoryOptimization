function ADMM(problem, num_iter)
	n = 6
	N = 10
	m = 3

	X = zeros(problem.n*(problem.N+1))
	U = zeros(problem.m*problem.N)
	Y = ones(problem.m*problem.N)
	Λ̄ = ones(problem.m*problem.N)

	for i=1:num_iter
		(X, U) = LQR(ρ, Y, Λ̄, problem)
		for k=1:N
			β = problem.α/problem.ρ
			v = U[Select3(k)] + Λ̄[Select3(k)]
			Y[Select3(k)] = L2Prox(β, v)

			Λ̄[Select3(k)] += U[Select3(k)] - Y[Select3(k)]

			@show norm(U[Select3(k)] - Y[Select3(k)])
		end

	end

	return (X, U)
end

function LQR(ρ, Y, Λ̄, p)
	# Cost to go function over time
	V = zeros(p.n, p.n, p.N+1)
	V[:,:, p.N+1] = p.Q_f
	v = zeros(p.n, p.N+1)
	c = zeros(p.N+1)

	K = zeros(p.m, p.n, p.N)
	d = zeros(p.m, p.N)

	R = ρ/2*Matrix{Float64}(I, p.m, p.m)
	U_ref = Y - Λ̄

	k = p.N
	while k > 0
		K[:,:,k] = -inv(R + p.B'*V[:,:,k+1]*p.B)*p.B'*V[:,:,k+1]*p.A
		d[:,k] = inv(R + p.B'*V[:,:,k+1]*p.B)*(R*U_ref[Select3(k)] - 0.5*p.B'*v[:,k+1])

		V[:,:, k] = p.Q_k + K[:,:,k]'*R*K[:,:,k] + (p.A + p.B*K[:,:,k])'*V[:,:,k+1]*(p.A+p.B*K[:,:,k])
		v[:, k] = 2*((d[:,k] - U_ref[Select3(k)])'*R*K[:,:,k] + d[:,k]'*p.B'*V[:,:,k+1]*(p.A+p.B*K[:,:,k]))'
		c[k] = (d[:,k]-U_ref[Select3(k)])'*R*(d[:,k]-U_ref[Select3(k)]) + d[:,k]'*p.B'*V[:,:,k+1]*p.B*d[:,k]

		k -= 1
	end

	X = zeros(p.n*(p.N+1))
	U = zeros(p.m*p.N)
	X[Select6(1)] = p.x0

	for i=2:(p.N+1)
		U[Select3(i-1)] = K[:,:, i-1]*X[Select6(i-1)] + d[:, i-1]
		X[Select6(i)] = p.A*X[Select6(i-1)] + p.B*U[Select3(i-1)]
	end

	return (X, U)
end

function L2Prox(β, v)
	return max(0, 1-β/norm(v, 2))*v
end

function L1Prox(β, v)
	return max.(0, v .- β) - max.(0, -v .- β)
end
